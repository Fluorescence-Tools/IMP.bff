#include <IMP/bff/FRETRotamer.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/Text.h>
#include <IMP/bff/ProbeLibrary.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>


// -------- from RotamerFret.cpp --------
/**
 * (formerly RotamerFret.cpp, now a section of this file)
 * \brief FRET over a trajectory from two screened rotamer libraries.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */




IMPBFF_BEGIN_NAMESPACE

namespace {
//! The `(n, 4)` points and the `(n, 3)` dipoles of any states.
/*! The dipoles come back empty when the states carry none *or* carry a number
    that does not match their points -- an AV cloud has none, and a mismatched
    orientation array cannot be paired with a point by index, which is the only
    way \f$\kappa^2\f$ can use it. */
void states_arrays(const States& s, std::vector<double>& points,
                   std::vector<double>& mu) {
    internal::OwnedView p, o;
    s.get_points(&p.data, &p.size);
    s.get_orientations(&o.data, &o.size);
    points = p.vector();
    mu.clear();
    const std::size_t n = points.size() / 4;
    if (o.size > 0 && static_cast<std::size_t>(o.size) == n * 3) {
        mu = o.vector();
    }
}

}  // namespace

FRETPairGeometry ProbeRotamerEnsemble::pair_geometry(const States& other,
                                                bool use_dipoles) const {
    std::vector<double> other_points, other_mu;
    states_arrays(other, other_points, other_mu);
    if (!use_dipoles) other_mu.clear();

    const int n1 = get_n_points();
    const std::size_t n2 = other_points.size() / 4;
    std::vector<double> c1(static_cast<std::size_t>(n1) * 3), w1(n1);
    for (int i = 0; i < n1; ++i) {
        c1[i * 3] = points_[i * 4];
        c1[i * 3 + 1] = points_[i * 4 + 1];
        c1[i * 3 + 2] = points_[i * 4 + 2];
        w1[i] = points_[i * 4 + 3];
    }
    std::vector<double> c2(n2 * 3), w2(n2);
    for (std::size_t i = 0; i < n2; ++i) {
        c2[i * 3] = other_points[i * 4];
        c2[i * 3 + 1] = other_points[i * 4 + 1];
        c2[i * 3 + 2] = other_points[i * 4 + 2];
        w2[i] = other_points[i * 4 + 3];
    }
    return fret_pair_geometry(
            c1, w1, c2, w2,
            use_dipoles ? orientations_ : std::vector<double>(), other_mu);
}

FRETPairEfficiencies ProbeRotamerEnsemble::pair_distribution(
        const States& other, double forster_radius, double tau0) const {
    return fret_pair_efficiencies(pair_geometry(other), forster_radius, tau0);
}

FRETPairEfficiencies ProbeRotamerEnsemble::pair_distribution_from_probes(
        const States& other, const std::string& donor,
        const std::string& acceptor, double tau0) const {
    const FRETPairGeometry geometry = pair_geometry(other);
    // Angstrom, like every other length here; the spectra are nanometres and
    // `forster_radius_from_spectra` converts once, inside.
    const double r0 = forster_radius_from_spectra(donor, acceptor,
                                                  geometry.kappa2_avg);
    return fret_pair_efficiencies(geometry, r0, tau0);
}



namespace {

//! `np.savetxt`'s default: one row per line, `%.18e` per column.
void save_txt(const std::string& path, const std::vector<double>& values,
              int columns, const std::string& header = "") {
    std::ofstream out(path.c_str());
    if (!out) IMP_THROW("cannot write " << path, IOException);
    if (!header.empty()) {
        std::istringstream lines(header);
        std::string line;
        while (std::getline(lines, line)) out << "# " << line << "\n";
    }
    out << std::scientific << std::setprecision(18);
    for (std::size_t i = 0; i < values.size(); ++i) {
        out << values[i];
        const bool row_end = (columns <= 1) ||
                             ((i + 1) % static_cast<std::size_t>(columns) == 0);
        out << (row_end ? "\n" : " ");
    }
}

//! Every number of a whitespace-separated text file, in file order.
std::vector<double> load_txt(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) IMP_THROW("cannot read " << path, IOException);
    std::vector<double> out;
    std::string line;
    while (std::getline(in, line)) {
        const std::string trimmed = internal::trimmed(line);
        if (trimmed.empty() || trimmed[0] == '#') continue;
        std::istringstream fields(trimmed);
        double v;
        while (fields >> v) out.push_back(v);
    }
    return out;
}

//! The finite entries of \p values and their weights, both renormalised.
void finite_pairs(const std::vector<double>& values,
                  const std::vector<double>& weights,
                  std::vector<double>& out_values,
                  std::vector<double>& out_weights) {
    out_values.clear();
    out_weights.clear();
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (!std::isfinite(values[i])) continue;
        out_values.push_back(values[i]);
        out_weights.push_back(i < weights.size() ? weights[i] : 0.0);
    }
}

}  // namespace

void FRETRotamer::load_libraries() {
    if (residues_.size() != 2) {
        IMP_THROW("the residue list must contain exactly 2 residue numbers, "
                  "not " << residues_.size(),
                  ValueException);
    }
    while (chains_.size() < 2) chains_.push_back(std::string());
    lib_1_ = load_probe_rotamer_library(libname_1_);
    lib_2_ = load_probe_rotamer_library(libname_2_);
}

FRETRotamer::FRETRotamer(
        const std::string& protein, const std::vector<int>& residues,
        const std::vector<std::string>& chains, const std::string& donor,
        const std::string& acceptor, const std::string& libname_1,
        const std::string& libname_2, double temperature, bool electrostatic,
        const std::string& potential, bool ign_H, double sigma_scaling,
        double epsilon_scaling, bool fixed_R0, double r0,
        const std::string& r0lib, double z_cutoff,
        const std::string& output_prefix, int max_frames,
        const std::vector<double>& user_weights)
    : frames_(load_protein_frames(protein, max_frames)), residues_(residues),
      chains_(chains), donor_(donor), acceptor_(acceptor),
      libname_1_(libname_1), libname_2_(libname_2), r0lib_(r0lib),
      output_prefix_(output_prefix),
      site_(temperature, electrostatic, potential, ign_H, sigma_scaling,
            epsilon_scaling),
      fixed_r0_(fixed_R0), r0_(r0), z_cutoff_(z_cutoff),
      user_weights_(user_weights) {
    load_libraries();
}

FRETRotamer::FRETRotamer(const std::vector<ProteinFrame>& frames,
                         const std::vector<int>& residues,
                         const std::vector<std::string>& chains,
                         const std::string& donor,
                         const std::string& acceptor,
                         const std::string& libname_1,
                         const std::string& libname_2,
                         const ProbeRotamerSiteOptions& site, bool fixed_R0,
                         double r0, const std::string& r0lib, double z_cutoff,
                         const std::string& output_prefix,
                         const std::vector<double>& user_weights)
    : frames_(frames), residues_(residues), chains_(chains), donor_(donor),
      acceptor_(acceptor), libname_1_(libname_1), libname_2_(libname_2),
      r0lib_(r0lib), output_prefix_(output_prefix), site_(site),
      fixed_r0_(fixed_R0), r0_(r0), z_cutoff_(z_cutoff),
      user_weights_(user_weights) {
    load_libraries();
}

FRETRotamer FRETRotamer::from_frames(
        const std::vector<ProteinFrame>& frames,
        const std::vector<int>& residues,
        const std::vector<std::string>& chains, const std::string& donor,
        const std::string& acceptor, const std::string& libname_1,
        const std::string& libname_2, double temperature, bool electrostatic,
        const std::string& potential, bool ign_H, double sigma_scaling,
        double epsilon_scaling, bool fixed_R0, double r0,
        const std::string& r0lib, double z_cutoff,
        const std::string& output_prefix,
        const std::vector<double>& user_weights) {
    return FRETRotamer(frames, residues, chains, donor, acceptor, libname_1,
                       libname_2,
                       ProbeRotamerSiteOptions(temperature, electrostatic,
                                          potential, ign_H, sigma_scaling,
                                          epsilon_scaling),
                       fixed_R0, r0, r0lib, z_cutoff, output_prefix,
                       user_weights);
}

FRETRotamerFrameResult FRETRotamer::frame_fret(const ProteinFrame& frame) {
    const ProbeRotamerEnsemble donor = ProbeRotamerEnsemble::from_frame(
            frame, chains_[0], residues_[0], lib_1_, site_);
    const ProbeRotamerEnsemble acceptor = ProbeRotamerEnsemble::from_frame(
            frame, chains_[1], residues_[1], lib_2_, site_);

    FRETRotamerFrameResult out;
    out.z_donor = donor.get_partition();
    out.z_acceptor = acceptor.get_partition();

    const FRETPairGeometry geometry = donor.pair_geometry(acceptor);
    out.kappa2_avg = geometry.kappa2_avg;
    if (!fixed_r0_) {
        // FRETpredict's convention: R0 is recomputed per frame at that
        // frame's <kappa2> and then applied with the isotropic formula.
        r0_ = forster_radius_from_spectra(donor_, acceptor_,
                                          geometry.kappa2_avg, r0lib_);
        if (r0_ == 0.0) {
            const double nan = std::numeric_limits<double>::quiet_NaN();
            out.kappa2_avg = nan;
            out.static_efficiency = nan;
            out.dynamic1 = nan;
            out.dynamic2 = nan;
            return out;
        }
    }
    const FRETPairEfficiencies eff = fret_pair_efficiencies(geometry, r0_);
    out.static_efficiency = eff.static_efficiency;
    out.dynamic1 = eff.dynamic1;
    out.dynamic2 = eff.dynamic2;
    return out;
}

void FRETRotamer::trajectory_analysis() {
    const std::size_t n = frames_.size();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    z_values_.assign(n * 2, 0.0);
    k2_values_.assign(n, nan);
    estatic_.assign(n, nan);
    edynamic1_.assign(n, nan);
    edynamic2_.assign(n, nan);

    for (std::size_t i = 0; i < n; ++i) {
        const FRETRotamerFrameResult result = frame_fret(frames_[i]);
        z_values_[i * 2] = result.z_donor;
        z_values_[i * 2 + 1] = result.z_acceptor;
        if (result.z_donor <= z_cutoff_ || result.z_acceptor <= z_cutoff_) {
            continue;
        }
        k2_values_[i] = result.kappa2_avg;
        estatic_[i] = result.static_efficiency;
        edynamic1_[i] = result.dynamic1;
        edynamic2_[i] = result.dynamic2;
    }
}

void FRETRotamer::write_summary(const std::string& prefix,
                                const std::vector<double>& k2,
                                const std::vector<double>& es,
                                const std::vector<double>& ed1,
                                const std::vector<double>& ed2,
                                const std::vector<double>& weights) const {
    std::vector<double> rows;
    const std::vector<double>* quantities[4] = {&k2, &es, &ed1, &ed2};
    if (k2.size() == 1) {
        // One frame has no spread: the average is the value and the SD and SE
        // are undefined rather than zero.
        const double nan = std::numeric_limits<double>::quiet_NaN();
        for (int q = 0; q < 4; ++q) {
            rows.push_back((*quantities[q])[0]);
            rows.push_back(nan);
            rows.push_back(nan);
        }
    } else {
        for (int q = 0; q < 4; ++q) {
            std::vector<double> values, w;
            finite_pairs(*quantities[q], weights, values, w);
            const std::vector<double> stats = weighted_average_sd_se(
                    values.empty() ? NULL : &values[0],
                    static_cast<int>(values.size()),
                    w.empty() ? NULL : &w[0], static_cast<int>(w.size()));
            rows.insert(rows.end(), stats.begin(), stats.end());
        }
    }
    std::ostringstream path;
    path << prefix << "-data-" << residues_[0] << "-" << residues_[1] << ".dat";
    save_txt(path.str(), rows, 3,
             "quantity Average SD SE\nk2 Estatic Edynamic1 Edynamic2");
}

void FRETRotamer::save(const std::string& output_prefix) {
    if (z_values_.empty()) trajectory_analysis();
    const std::string prefix =
            output_prefix.empty() ? output_prefix_ : output_prefix;
    std::ostringstream tail;
    tail << "-" << residues_[0] << "-" << residues_[1] << ".dat";

    internal::OwnedView ws;
    frame_weights_from_partitions(z_values_.empty() ? NULL : &z_values_[0],
                                  static_cast<int>(z_values_.size() / 2), 2,
                                  &ws.data, &ws.size);

    save_txt(prefix + "-Z" + tail.str(), z_values_, 2);
    save_txt(prefix + "-w_s" + tail.str(), ws.vector(), 1);
    save_txt(prefix + "-k2" + tail.str(), k2_values_, 1);
    save_txt(prefix + "-Es" + tail.str(), estatic_, 1);
    save_txt(prefix + "-Ed1" + tail.str(), edynamic1_, 1);
    save_txt(prefix + "-Ed2" + tail.str(), edynamic2_, 1);

    std::vector<double> weights(k2_values_.size(), 1.0);
    if (!user_weights_.empty()) {
        if (user_weights_.size() != k2_values_.size()) {
            IMP_THROW("Weights array has size "
                              << user_weights_.size()
                              << " whereas the number of frames is "
                              << k2_values_.size(),
                      ValueException);
        }
        weights = user_weights_;
    }
    double total = 0.0;
    for (std::size_t i = 0; i < weights.size(); ++i) total += weights[i];
    if (total > 0.0) {
        for (std::size_t i = 0; i < weights.size(); ++i) weights[i] /= total;
    }
    write_summary(prefix, k2_values_, estatic_, edynamic1_, edynamic2_,
                  weights);
}

void FRETRotamer::reweight(bool boltzmann_weights,
                           const std::vector<double>& user_weights,
                           const std::string& output_prefix) {
    const std::string prefix =
            output_prefix.empty() ? output_prefix_ : output_prefix;
    std::ostringstream tail;
    tail << "-" << residues_[0] << "-" << residues_[1] << ".dat";

    if (boltzmann_weights) {
        std::vector<double> z = load_txt(output_prefix_ + "-Z" + tail.str());
        internal::OwnedView ws;
        frame_weights_from_partitions(z.empty() ? NULL : &z[0],
                                      static_cast<int>(z.size() / 2), 2,
                                      &ws.data, &ws.size);
        weights_ = ws.vector();
    }
    if (!user_weights.empty()) user_weights_ = user_weights;

    const std::vector<double> k2 = load_txt(output_prefix_ + "-k2" +
                                            tail.str());
    const std::vector<double> es = load_txt(output_prefix_ + "-Es" +
                                            tail.str());
    const std::vector<double> ed1 = load_txt(output_prefix_ + "-Ed1" +
                                             tail.str());
    const std::vector<double> ed2 = load_txt(output_prefix_ + "-Ed2" +
                                             tail.str());

    std::vector<double> weights(k2.size(), 1.0);
    if (weights_.size() == k2.size()) weights = weights_;
    if (!user_weights_.empty()) {
        if (user_weights_.size() != k2.size()) {
            IMP_THROW("Weights array has size "
                              << user_weights_.size()
                              << " whereas the number of frames is "
                              << k2.size(),
                      ValueException);
        }
        for (std::size_t i = 0; i < weights.size(); ++i) {
            weights[i] *= user_weights_[i];
        }
    }
    double total = 0.0;
    for (std::size_t i = 0; i < weights.size(); ++i) total += weights[i];
    if (total > 0.0) {
        for (std::size_t i = 0; i < weights.size(); ++i) weights[i] /= total;
    }
    write_summary(prefix, k2, es, ed1, ed2, weights);
}

void FRETRotamer::run() {
    trajectory_analysis();
    save();
}

void FRETRotamer::get_z_values(double** out_view, int* n_out_view) const {
    internal::copy_to_view(z_values_, out_view, n_out_view);
}
void FRETRotamer::get_k2_values(double** out_view, int* n_out_view) const {
    internal::copy_to_view(k2_values_, out_view, n_out_view);
}
void FRETRotamer::get_estatic_values(double** out_view, int* n_out_view) const {
    internal::copy_to_view(estatic_, out_view, n_out_view);
}
void FRETRotamer::get_edynamic1_values(double** out_view,
                                       int* n_out_view) const {
    internal::copy_to_view(edynamic1_, out_view, n_out_view);
}
void FRETRotamer::get_edynamic2_values(double** out_view,
                                       int* n_out_view) const {
    internal::copy_to_view(edynamic2_, out_view, n_out_view);
}

IMPBFF_END_NAMESPACE
