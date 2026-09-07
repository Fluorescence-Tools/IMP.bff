/**
 * \file ModelPrecision.cpp
 * \brief How far a docked body wanders between independent runs.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/ModelPrecision.h>

#include <IMP/bff/internal/OutputView.h>

#include <IMP/bff/Base.h>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

// Named, not anonymous: IMP compiles this module as one translation unit.
namespace precision {

//! `(lines, xyz, chains)` for the ATOM/HETATM records of a PDB file.
bool read_pdb_atoms(const std::string& path, std::vector<std::string>& lines,
                    std::vector<double>& xyz,
                    std::vector<std::string>& chains) {
    std::ifstream in(path.c_str());
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.compare(0, 4, "ATOM") != 0 &&
            line.compare(0, 6, "HETATM") != 0) {
            continue;
        }
        if (line.size() < 54) continue;
        lines.push_back(line);
        xyz.push_back(std::atof(line.substr(30, 8).c_str()));
        xyz.push_back(std::atof(line.substr(38, 8).c_str()));
        xyz.push_back(std::atof(line.substr(46, 8).c_str()));
        chains.push_back(line.substr(21, 1));
    }
    return true;
}

//! The rigid transform minimising \f$\|R x + t - y\|\f$ (Kabsch).
void kabsch(const Eigen::MatrixXd& mobile, const Eigen::MatrixXd& target,
            Eigen::Matrix3d& r, Eigen::Vector3d& t) {
    const Eigen::Vector3d mc = mobile.colwise().mean();
    const Eigen::Vector3d tc = target.colwise().mean();
    const Eigen::Matrix3d h =
            (mobile.rowwise() - mc.transpose()).transpose() *
            (target.rowwise() - tc.transpose());
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(
            h, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d d = Eigen::Matrix3d::Identity();
    // The reflection guard: an SVD of a near-degenerate H can hand back a
    // determinant of -1, which superposes the structure onto its mirror image
    // and reports a plausible-looking RMSF for a chirality flip.
    d(2, 2) = (svd.matrixV() * svd.matrixU().transpose()).determinant() < 0.0
                      ? -1.0
                      : 1.0;
    r = svd.matrixV() * d * svd.matrixU().transpose();
    t = tc - r * mc;
}

Eigen::MatrixXd as_matrix(const std::vector<double>& xyz) {
    const int n = static_cast<int>(xyz.size() / 3);
    Eigen::MatrixXd out(n, 3);
    for (int i = 0; i < n; ++i) {
        out(i, 0) = xyz[3 * i + 0];
        out(i, 1) = xyz[3 * i + 1];
        out(i, 2) = xyz[3 * i + 2];
    }
    return out;
}

}  // namespace precision

double PositionUncertainty::get_rmsf_mean() const {
    if (rmsf.empty()) return std::nan("");
    double total = 0.0;
    for (std::size_t i = 0; i < rmsf.size(); ++i) total += rmsf[i];
    return total / rmsf.size();
}

double PositionUncertainty::get_rmsf_max() const {
    if (rmsf.empty()) return std::nan("");
    double best = rmsf[0];
    for (std::size_t i = 1; i < rmsf.size(); ++i) {
        if (rmsf[i] > best) best = rmsf[i];
    }
    return best;
}

double PositionUncertainty::get_mobile_rmsf_mean() const {
    if (rmsf.empty()) return std::nan("");
    double total = 0.0;
    int n = 0;
    for (std::size_t i = 0; i < rmsf.size(); ++i) {
        if (i < fixed_mask.size() && fixed_mask[i]) continue;
        total += rmsf[i];
        ++n;
    }
    // Every atom fixed means there is no mobile body; the whole-structure RMSF
    // is then what "the mobile body wandered" can mean.
    return n > 0 ? total / n : get_rmsf_mean();
}

void PositionUncertainty::get_rmsf(double** out_view, int* n_out_view) const {
    internal::copy_to_view(rmsf, out_view, n_out_view);
}

void PositionUncertainty::get_mean_coords(double** out_view,
                                          int* n_out_view) const {
    internal::copy_to_view(mean_coords, out_view, n_out_view);
}

PositionUncertainty estimate_position_uncertainty(
        const std::vector<std::string>& pdb_paths,
        const std::vector<std::string>& fixed_chains) {
    PositionUncertainty out;

    std::vector<std::string> lines0, chains0;
    std::vector<double> xyz0;
    std::vector<std::string> readable;
    for (std::size_t i = 0; i < pdb_paths.size(); ++i) {
        if (pdb_paths[i].empty()) continue;
        std::ifstream probe(pdb_paths[i].c_str());
        if (probe) readable.push_back(pdb_paths[i]);
    }
    if (readable.size() < 2) {
        out.n_models = static_cast<int>(readable.size());
        return out;
    }
    precision::read_pdb_atoms(readable[0], lines0, xyz0, out.chains);
    const int n_atoms = static_cast<int>(out.chains.size());
    if (n_atoms == 0) {
        out.n_models = 0;
        return out;
    }

    out.fixed_mask.assign(n_atoms, 0);
    int n_fixed = 0;
    for (int i = 0; i < n_atoms; ++i) {
        for (std::size_t c = 0; c < fixed_chains.size(); ++c) {
            if (out.chains[i] == fixed_chains[c]) {
                out.fixed_mask[i] = 1;
                ++n_fixed;
                break;
            }
        }
    }
    if (n_fixed == 0) out.fixed_mask.assign(n_atoms, 1);

    const Eigen::MatrixXd reference = precision::as_matrix(xyz0);
    Eigen::MatrixXd ref_fixed(n_fixed > 0 ? n_fixed : n_atoms, 3);
    std::vector<int> fixed_rows;
    for (int i = 0; i < n_atoms; ++i) {
        if (out.fixed_mask[i]) fixed_rows.push_back(i);
    }
    for (std::size_t i = 0; i < fixed_rows.size(); ++i) {
        ref_fixed.row(i) = reference.row(fixed_rows[i]);
    }

    std::vector<Eigen::MatrixXd> aligned;
    aligned.push_back(reference);
    for (std::size_t p = 1; p < readable.size(); ++p) {
        std::vector<std::string> lines, chains;
        std::vector<double> xyz;
        precision::read_pdb_atoms(readable[p], lines, xyz, chains);
        if (static_cast<int>(chains.size()) != n_atoms) continue;
        const Eigen::MatrixXd m = precision::as_matrix(xyz);
        Eigen::MatrixXd m_fixed(fixed_rows.size(), 3);
        for (std::size_t i = 0; i < fixed_rows.size(); ++i) {
            m_fixed.row(i) = m.row(fixed_rows[i]);
        }
        Eigen::Matrix3d r;
        Eigen::Vector3d t;
        precision::kabsch(m_fixed, ref_fixed, r, t);
        aligned.push_back((m * r.transpose()).rowwise() + t.transpose());
    }

    out.n_models = static_cast<int>(aligned.size());
    Eigen::MatrixXd mean = Eigen::MatrixXd::Zero(n_atoms, 3);
    for (std::size_t i = 0; i < aligned.size(); ++i) mean += aligned[i];
    mean /= static_cast<double>(aligned.size());

    out.mean_coords.assign(static_cast<std::size_t>(n_atoms) * 3, 0.0);
    for (int i = 0; i < n_atoms; ++i) {
        for (int k = 0; k < 3; ++k) out.mean_coords[3 * i + k] = mean(i, k);
    }

    out.rmsf.assign(n_atoms, 0.0);
    for (int i = 0; i < n_atoms; ++i) {
        double total = 0.0;
        for (std::size_t m = 0; m < aligned.size(); ++m) {
            const double dx = aligned[m](i, 0) - mean(i, 0);
            const double dy = aligned[m](i, 1) - mean(i, 1);
            const double dz = aligned[m](i, 2) - mean(i, 2);
            total += dx * dx + dy * dy + dz * dz;
        }
        out.rmsf[i] = std::sqrt(total / aligned.size());
    }
    return out;
}

std::vector<std::string> pdb_chain_ids(const std::string& path) {
    std::vector<std::string> lines, chains;
    std::vector<double> xyz;
    if (!precision::read_pdb_atoms(path, lines, xyz, chains)) {
        IMP_THROW("Cannot read " << path, IOException);
    }
    std::sort(chains.begin(), chains.end());
    chains.erase(std::unique(chains.begin(), chains.end()), chains.end());
    return chains;
}

void write_position_uncertainty_pdb(const PositionUncertainty& uncertainty,
                                    const std::string& template_pdb,
                                    const std::string& out_pdb) {
    std::vector<std::string> lines, chains;
    std::vector<double> xyz;
    if (!precision::read_pdb_atoms(template_pdb, lines, xyz, chains)) {
        IMP_THROW("Cannot read " << template_pdb, IOException);
    }
    std::ofstream fh(out_pdb.c_str());
    if (!fh) IMP_THROW("Cannot write " << out_pdb, IOException);
    const std::size_t n = std::min(lines.size(), uncertainty.rmsf.size());
    char buffer[128];
    for (std::size_t i = 0; i < n; ++i) {
        // 999.99 is what the six-column B-factor field holds; a wider value
        // would run into the occupancy column and shift every field after it.
        const double b = std::min(999.99, uncertainty.rmsf[i]);
        std::snprintf(buffer, sizeof(buffer), "%8.3f%8.3f%8.3f",
                      uncertainty.mean_coords[3 * i + 0],
                      uncertainty.mean_coords[3 * i + 1],
                      uncertainty.mean_coords[3 * i + 2]);
        const std::string& line = lines[i];
        fh << line.substr(0, 30) << buffer
           << (line.size() > 60 ? line.substr(54, 6) : std::string("      "));
        std::snprintf(buffer, sizeof(buffer), "%6.2f", b);
        fh << buffer << (line.size() > 66 ? line.substr(66) : std::string())
           << "\n";
    }
    fh << "END\n";
}

void write_position_uncertainty_csv(const PositionUncertainty& uncertainty,
                                    const std::string& out_csv) {
    std::ofstream fh(out_csv.c_str());
    if (!fh) IMP_THROW("Cannot write " << out_csv, IOException);
    fh << "atom_index,chain,rmsf\n";
    char buffer[64];
    for (std::size_t i = 0; i < uncertainty.rmsf.size(); ++i) {
        std::snprintf(buffer, sizeof(buffer), "%.3f", uncertainty.rmsf[i]);
        fh << i << ","
           << (i < uncertainty.chains.size() ? uncertainty.chains[i]
                                             : std::string())
           << "," << buffer << "\n";
    }
}

IMPBFF_END_NAMESPACE
