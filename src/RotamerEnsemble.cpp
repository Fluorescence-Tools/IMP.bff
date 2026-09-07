/**
 * \file RotamerEnsemble.cpp
 * \brief A rotamer library placed and screened at one labelling site.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/RotamerEnsemble.h>

#include <IMP/bff/ProbeLibrary.h>
#include <IMP/bff/FPSIO.h>
#include <IMP/bff/RotamerSite.h>
#include <IMP/bff/Scoring.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/Text.h>
#include <IMP/bff/internal/json.h>

#include <cmath>
#include <cstdlib>
#include <sstream>

#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

const char* const SIMULATION_TYPE_R1 = "R1";

namespace {

//! A metadata selector entry, which the registry writes as a string or a list.
std::vector<std::string> selector_list(const nlohmann::json& meta,
                                       const char* key) {
    std::vector<std::string> out;
    if (!meta.is_object() || !meta.contains(key) || meta[key].is_null()) {
        return out;
    }
    const nlohmann::json& v = meta[key];
    if (v.is_string()) {
        out.push_back(v.get<std::string>());
    } else if (v.is_array()) {
        for (nlohmann::json::const_iterator it = v.begin(); it != v.end();
             ++it) {
            if (it->is_string()) out.push_back(it->get<std::string>());
        }
    }
    return out;
}

//! `str()` of a number, spelled as Python spells it -- `298.15`, not
//! `298.150000`. The params map is provenance a reader reads back.
std::string number_text(double v) {
    std::ostringstream s;
    s << v;
    return s.str();
}

std::string bool_text(bool v) { return v ? "True" : "False"; }

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

// --------------------------------------------------------------------------
// RotamerEnsemble
// --------------------------------------------------------------------------

RotamerEnsemble::RotamerEnsemble(
        const std::vector<double>& points,
        const std::vector<double>& attachment_point,
        const std::vector<double>& orientations,
        const std::string& position_name,
        const std::map<std::string, std::string>& params,
        const std::vector<double>& atoms,
        const std::vector<std::string>& atom_names,
        const std::vector<std::string>& resnames,
        const std::vector<double>& energies, double partition,
        const std::string& library, const std::string& chain, int residue)
    : States(points, attachment_point, orientations, position_name, params),
      atoms_(atoms), atom_names_(atom_names), resnames_(resnames),
      energies_(energies), partition_(partition), library_(library),
      chain_(chain), residue_(residue) {}

void RotamerEnsemble::get_atoms(double** out_view, int* n_out_view) const {
    internal::copy_to_view(atoms_, out_view, n_out_view);
}

void RotamerEnsemble::get_energies(double** out_view, int* n_out_view) const {
    internal::copy_to_view(energies_, out_view, n_out_view);
}

void RotamerEnsemble::get_centres(double** out_view, int* n_out_view) const {
    const int n = get_n_points();
    double* buffer = internal::new_double_view(static_cast<std::size_t>(n) * 3,
                                               out_view, n_out_view);
    if (buffer == NULL) return;
    for (int i = 0; i < n; ++i) {
        buffer[i * 3] = points_[i * 4];
        buffer[i * 3 + 1] = points_[i * 4 + 1];
        buffer[i * 3 + 2] = points_[i * 4 + 2];
    }
}

void RotamerEnsemble::get_weights(double** out_view, int* n_out_view) const {
    const int n = get_n_points();
    double* buffer = internal::new_double_view(static_cast<std::size_t>(n),
                                               out_view, n_out_view);
    if (buffer == NULL) return;
    for (int i = 0; i < n; ++i) buffer[i] = points_[i * 4 + 3];
}

int RotamerEnsemble::get_n_atoms() const {
    const int n = get_n_points();
    if (n == 0) return 0;
    return static_cast<int>(atoms_.size() / (static_cast<std::size_t>(n) * 3));
}

double RotamerEnsemble::get_effective_sample_size() const {
    double sum = 0.0, sum_sq = 0.0;
    for (int i = 0; i < get_n_points(); ++i) {
        const double w = points_[i * 4 + 3];
        if (!(w > 0.0) || !std::isfinite(w)) continue;
        sum += w;
        sum_sq += w * w;
    }
    if (sum_sq <= 0.0) return 0.0;
    return sum * sum / sum_sq;
}

FRETPairGeometry RotamerEnsemble::pair_geometry(const States& other,
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

FRETPairEfficiencies RotamerEnsemble::pair_distribution(
        const States& other, double forster_radius, double tau0) const {
    return fret_pair_efficiencies(pair_geometry(other), forster_radius, tau0);
}

FRETPairEfficiencies RotamerEnsemble::pair_distribution_from_probes(
        const States& other, const std::string& donor,
        const std::string& acceptor, double tau0) const {
    const FRETPairGeometry geometry = pair_geometry(other);
    // Angstrom, like every other length here; the spectra are nanometres and
    // `forster_radius_from_spectra` converts once, inside.
    const double r0 = forster_radius_from_spectra(donor, acceptor,
                                                  geometry.kappa2_avg);
    return fret_pair_efficiencies(geometry, r0, tau0);
}

// --------------------------------------------------------------------------
// Placing a library at a site
// --------------------------------------------------------------------------

RotamerEnsemble RotamerEnsemble::from_frame(
        const ProteinFrame& frame, const std::string& chain, int residue,
        const RotamerLibrary& library, const RotamerSiteOptions& options,
        const std::string& position_name) {
    internal::OwnedView backbone;
    resolve_backbone_site(frame.coords, frame.atom_names, frame.chain_ids,
                          frame.residue_indices, chain, residue,
                          &backbone.data, &backbone.size);
    const std::vector<double> ca(backbone.data, backbone.data + 3);
    const std::vector<double> n_atom(backbone.data + 3, backbone.data + 6);
    const std::vector<double> c_atom(backbone.data + 6, backbone.data + 9);

    internal::OwnedView placed;
    transform_library_to_site(library.coords, ca, n_atom, c_atom, &placed.data,
                              &placed.size);
    const std::vector<double> rotamers = placed.vector();

    nlohmann::json meta = nlohmann::json::object();
    if (!library.metadata.empty()) {
        meta = nlohmann::json::parse(library.metadata, NULL, false);
        if (meta.is_discarded()) meta = nlohmann::json::object();
    }

    const RotamerScoreResult score = compute_rotamer_score(
            rotamers, frame.coords, frame.atom_names, frame.resnames,
            library.atom_names, selector_list(meta, "positive"),
            selector_list(meta, "negative"), library.resnames,
            frame.residue_indices, frame.chain_ids, residue, chain,
            library.weights, options.temperature, options.ignore_h,
            options.electrostatic, options.potential, options.sigma_scaling,
            options.epsilon_scaling);

    const int n_atoms = static_cast<int>(library.atom_names.size());
    const int n_rotamers =
            n_atoms > 0 ? static_cast<int>(rotamers.size() / (n_atoms * 3)) : 0;

    // Which atom is the chromophore centre, and which two span the transition
    // dipole, are the registry's (`r` and `mu`). A library outside the
    // registry has neither: its first atom stands for the centre and the
    // vector from its first atom to its second for the dipole -- a direction,
    // not the dye's, and an off-registry library should carry selectors.
    const std::vector<int> centre_idx = selector_atom_indices(
            library.atom_names, selector_list(meta, "r"), library.resnames);
    const int centre = centre_idx.empty() ? 0 : centre_idx[0];
    const std::vector<int> mu_idx = selector_atom_indices(
            library.atom_names, selector_list(meta, "mu"), library.resnames);
    const int mu_from = mu_idx.size() >= 2 ? mu_idx[0] : 0;
    const int mu_to = mu_idx.size() >= 2 ? mu_idx[1] : 1;

    std::vector<double> points(static_cast<std::size_t>(n_rotamers) * 4, 0.0);
    std::vector<double> mu(static_cast<std::size_t>(n_rotamers) * 3, 0.0);
    for (int k = 0; k < n_rotamers; ++k) {
        const std::size_t base = static_cast<std::size_t>(k) * n_atoms * 3;
        for (int d = 0; d < 3; ++d) {
            points[k * 4 + d] = rotamers[base + centre * 3 + d];
        }
        points[k * 4 + 3] =
                k < static_cast<int>(score.weights.size()) ? score.weights[k]
                                                           : 0.0;
        if (n_atoms > mu_to && n_atoms > mu_from) {
            double norm = 0.0;
            for (int d = 0; d < 3; ++d) {
                const double v = rotamers[base + mu_to * 3 + d] -
                                 rotamers[base + mu_from * 3 + d];
                mu[k * 3 + d] = v;
                norm += v * v;
            }
            norm = std::sqrt(norm);
            // Two atoms at one position have no direction between them. The
            // Python divided anyway and put three NaNs in the array, which
            // reach kappa2 and poison every pair the rotamer takes part in.
            if (norm > 0.0) {
                for (int d = 0; d < 3; ++d) mu[k * 3 + d] /= norm;
            }
        }
    }

    std::string library_name;
    if (meta.is_object() && meta.contains("library_name") &&
        meta["library_name"].is_string()) {
        library_name = meta["library_name"].get<std::string>();
    } else if (meta.is_object() && meta.contains("name") &&
               meta["name"].is_string()) {
        library_name = meta["name"].get<std::string>();
    } else {
        library_name = library.path;
    }

    std::map<std::string, std::string> params;
    params["simulation_type"] = SIMULATION_TYPE_R1;
    params["library"] = library_name;
    params["chain"] = chain;
    params["residue"] = number_text(residue);
    params["temperature"] = number_text(options.temperature);
    params["electrostatic"] = bool_text(options.electrostatic);
    params["potential"] = options.potential;
    params["ignore_h"] = bool_text(options.ignore_h);
    params["sigma_scaling"] = number_text(options.sigma_scaling);
    params["epsilon_scaling"] = number_text(options.epsilon_scaling);
    params["partition"] = number_text(score.partition);

    std::string name = position_name;
    if (name.empty()) name = chain + number_text(residue);

    return RotamerEnsemble(points, ca, mu, name, params, rotamers,
                           library.atom_names, library.resnames,
                           score.energies, score.partition, library_name,
                           chain, residue);
}

RotamerEnsemble RotamerEnsemble::from_site(const std::string& structure,
                                           const std::string& chain,
                                           int residue,
                                           const std::string& library,
                                           const RotamerSiteOptions& options,
                                           const std::string& position_name,
                                           int frame_index) {
    const std::vector<ProteinFrame> frames =
            load_protein_frames(structure, frame_index + 1);
    if (frame_index >= static_cast<int>(frames.size())) {
        IMP_THROW(structure << " has " << frames.size() << " frame(s), frame "
                            << frame_index << " requested",
                  IMP::ValueException);
    }
    return from_frame(frames[frame_index], chain, residue,
                      load_rotamer_library(library), options, position_name);
}

// --------------------------------------------------------------------------
// Weighted averaging over an ensemble
// --------------------------------------------------------------------------

void frame_weights_from_partitions(double* z_values, int n_frames, int n_pair,
                           double** out_view, int* n_out_view) {
    if (n_pair != 2) {
        throw std::invalid_argument("Z must have shape (n_frames, 2)");
    }
    double* out = internal::new_double_view(n_frames < 0 ? 0 : n_frames,
                                           out_view, n_out_view);
    if (out == nullptr) return;

    double total = 0.0;
    for (int i = 0; i < n_frames; ++i) {
        out[i] = z_values[2 * i] * z_values[2 * i + 1];
        total += out[i];
    }
    if (total == 0.0) {
        // uniform, not undefined: a frame where neither dye has an accessible
        // conformer says nothing about the others
        const double u = n_frames > 0 ? 1.0 / (double) n_frames : 0.0;
        for (int i = 0; i < n_frames; ++i) out[i] = u;
        return;
    }
    for (int i = 0; i < n_frames; ++i) out[i] /= total;
}

std::vector<double> weighted_average_sd_se(double* values, int n_values,
                                           double* weights, int n_weights) {
    if (n_values != n_weights) {
        throw std::invalid_argument("values and weights must be the same length");
    }
    std::vector<double> v, w;
    double total = 0.0;
    for (int i = 0; i < n_values; ++i) {
        // a frame where the dye could not be placed contributes nothing rather
        // than poisoning the mean
        if (!std::isfinite(values[i])) continue;
        v.push_back(values[i]);
        w.push_back(weights[i]);
        total += weights[i];
    }
    std::vector<double> out(3, std::nan(""));
    if (v.empty()) return out;

    double mean = 0.0;
    for (size_t i = 0; i < v.size(); ++i) mean += v[i] * (w[i] / total);
    double variance = 0.0;
    for (size_t i = 0; i < v.size(); ++i) {
        const double d = v[i] - mean;
        variance += d * d * (w[i] / total);
    }
    out[0] = mean;
    out[1] = std::sqrt(variance);
    // by the surviving frame count, not the effective count
    out[2] = std::sqrt(variance / (double) v.size());
    return out;
}

double effective_frame_fraction(double* weights, int n_weights) {
    std::vector<double> w;
    for (int i = 0; i < n_weights; ++i)
        if (weights[i] != 0.0) w.push_back(weights[i]);
    if (w.empty()) return 0.0;

    const double uniform = 1.0 / (double) w.size();
    double entropy = 0.0;
    for (size_t i = 0; i < w.size(); ++i)
        entropy -= w[i] * std::log(w[i] / uniform);
    return std::exp(entropy);
}

IMPBFF_END_NAMESPACE
