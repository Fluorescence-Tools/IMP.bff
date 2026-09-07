/**
 * \file RotamerFret.cpp
 * \brief FRET over a trajectory from two screened rotamer libraries.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/RotamerFret.h>

#include <IMP/bff/ProbeLibrary.h>
#include <IMP/bff/FRETPair.h>
#include <IMP/bff/RotamerFps.h>
#include <IMP/bff/RotamerSite.h>
#include <IMP/bff/RotamerEnsemble.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/Text.h>

#include <IMP/bff/Base.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

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

void RotamerFRET::load_libraries() {
    if (residues_.size() != 2) {
        IMP_THROW("the residue list must contain exactly 2 residue numbers, "
                  "not " << residues_.size(),
                  ValueException);
    }
    while (chains_.size() < 2) chains_.push_back(std::string());
    lib_1_ = load_rotamer_library(libname_1_);
    lib_2_ = load_rotamer_library(libname_2_);
}

RotamerFRET::RotamerFRET(
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

RotamerFRET::RotamerFRET(const std::vector<ProteinFrame>& frames,
                         const std::vector<int>& residues,
                         const std::vector<std::string>& chains,
                         const std::string& donor,
                         const std::string& acceptor,
                         const std::string& libname_1,
                         const std::string& libname_2,
                         const RotamerSiteOptions& site, bool fixed_R0,
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

RotamerFRET RotamerFRET::from_frames(
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
    return RotamerFRET(frames, residues, chains, donor, acceptor, libname_1,
                       libname_2,
                       RotamerSiteOptions(temperature, electrostatic,
                                          potential, ign_H, sigma_scaling,
                                          epsilon_scaling),
                       fixed_R0, r0, r0lib, z_cutoff, output_prefix,
                       user_weights);
}

FRETFrameResult RotamerFRET::frame_fret(const ProteinFrame& frame) {
    const RotamerEnsemble donor = RotamerEnsemble::from_frame(
            frame, chains_[0], residues_[0], lib_1_, site_);
    const RotamerEnsemble acceptor = RotamerEnsemble::from_frame(
            frame, chains_[1], residues_[1], lib_2_, site_);

    FRETFrameResult out;
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

void RotamerFRET::trajectory_analysis() {
    const std::size_t n = frames_.size();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    z_values_.assign(n * 2, 0.0);
    k2_values_.assign(n, nan);
    estatic_.assign(n, nan);
    edynamic1_.assign(n, nan);
    edynamic2_.assign(n, nan);

    for (std::size_t i = 0; i < n; ++i) {
        const FRETFrameResult result = frame_fret(frames_[i]);
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

void RotamerFRET::write_summary(const std::string& prefix,
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

void RotamerFRET::save(const std::string& output_prefix) {
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

void RotamerFRET::reweight(bool boltzmann_weights,
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

void RotamerFRET::run() {
    trajectory_analysis();
    save();
}

void RotamerFRET::get_z_values(double** out_view, int* n_out_view) const {
    internal::copy_to_view(z_values_, out_view, n_out_view);
}
void RotamerFRET::get_k2_values(double** out_view, int* n_out_view) const {
    internal::copy_to_view(k2_values_, out_view, n_out_view);
}
void RotamerFRET::get_estatic_values(double** out_view, int* n_out_view) const {
    internal::copy_to_view(estatic_, out_view, n_out_view);
}
void RotamerFRET::get_edynamic1_values(double** out_view,
                                       int* n_out_view) const {
    internal::copy_to_view(edynamic1_, out_view, n_out_view);
}
void RotamerFRET::get_edynamic2_values(double** out_view,
                                       int* n_out_view) const {
    internal::copy_to_view(edynamic2_, out_view, n_out_view);
}

RotamerFRET rotamer_fret_from_fps(const std::string& fps_path,
                                  const std::string& protein,
                                  const std::string& distance_name,
                                  const std::string& output_prefix,
                                  bool fixed_R0, double r0, double temperature,
                                  bool electrostatic) {
    const RotamerFpsSelection selection = read_rotamer_fps(fps_path,
                                                           distance_name);
    std::vector<int> residues;
    residues.push_back(selection.donor.residue);
    residues.push_back(selection.acceptor.residue);
    std::vector<std::string> chains;
    chains.push_back(selection.donor.chain);
    chains.push_back(selection.acceptor.chain);

    const std::string donor = selection.distance.donor.empty()
                                      ? selection.donor.dye
                                      : selection.distance.donor;
    const std::string acceptor = selection.distance.acceptor.empty()
                                         ? selection.acceptor.dye
                                         : selection.distance.acceptor;
    const std::string libname_1 = selection.distance.libname_1.empty()
                                          ? selection.donor.library
                                          : selection.distance.libname_1;
    const std::string libname_2 = selection.distance.libname_2.empty()
                                          ? selection.acceptor.library
                                          : selection.distance.libname_2;
    return RotamerFRET(protein, residues, chains, donor, acceptor, libname_1,
                       libname_2, temperature, electrostatic, "lj", true, 0.5,
                       1.0, fixed_R0, r0, "", 0.05, output_prefix);
}

IMPBFF_END_NAMESPACE
