/**
 *  \file Scoring.cpp
 *  \brief Stage-2 scoring orchestration: CHARMM36, LJ, Boltzmann, AABB, and
 *         the end-to-end rotamer score.
 *
 * The inner kernels are in Rotamer.cpp (the RotamerEnergy section); this file is the layer that
 * builds parameters, calls them, and turns the energies into weights.
 */

#include <IMP/bff/Scoring.h>
#include <IMP/bff/Mol2IO.h>
#include <IMP/bff/Rotamer.h>
#include <IMP/bff/ZMatrix.h>

#include <IMP/algebra/vector_generators.h>
#include <IMP/constants.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <IMP/bff/internal/OutputView.h>

#include <IMP/bff/internal/Text.h>

#include <IMP/bff/Base.h>
#include <IMP/core/AngleRestraint.h>
#include <IMP/core/DihedralRestraint.h>
#include <IMP/core/DistanceRestraint.h>
#include <IMP/core/Harmonic.h>
#include <IMP/core/HarmonicLowerBound.h>
#include <IMP/container/ListPairContainer.h>
#include <IMP/container/PairsRestraint.h>
#include <IMP/core/SphereDistancePairScore.h>
#include <IMP/core/XYZ.h>
#include <IMP/core/XYZR.h>
#include <IMP/core/internal/dihedral_helpers.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>

IMPBFF_BEGIN_NAMESPACE

using internal::upper;

namespace {

std::string strip(const std::string& s) {
    const std::size_t a = s.find_first_not_of(" \t");
    if (a == std::string::npos) return std::string();
    const std::size_t b = s.find_last_not_of(" \t");
    return s.substr(a, b - a + 1);
}

}  // namespace

std::string atom_type(const std::string& atom_name) {
    const std::string name = upper(strip(atom_name));
    for (const char* atom : {"S", "N", "O", "H", "C"}) {
        if (name.rfind(atom, 0) == 0) return atom;
    }
    return "C";
}

namespace {

const std::map<std::string, std::pair<double, double>>& charmm36_table() {
    static const std::map<std::string, std::pair<double, double>> t = {
        {"C", {2.02446316, -0.06394724}},
        {"N", {1.89285714, -0.15428571}},
        {"O", {1.693,      -0.12642017}},
        {"S", {2.1,         -0.47}},
        {"H", {0.98357778,  -0.03466645}},
    };
    return t;
}

const std::pair<double, double>& charmm36_entry(const std::string& element) {
    auto& t = charmm36_table();
    auto it = t.find(element);
    return it == t.end() ? t.at("C") : it->second;
}

}  // namespace

void charmm36_lj(const std::string& element,
                 double** out_view, int* n_out_view) {
    const std::pair<double, double>& p = charmm36_entry(element);
    double* out = internal::new_double_view(2, out_view, n_out_view);
    if (out == NULL) return;
    out[0] = p.first;
    out[1] = p.second;
}

void lj_cross(const std::string& elem_i, const std::string& elem_j,
              double** out_view, int* n_out_view) {
    const std::pair<double, double>& pi = charmm36_entry(elem_i);
    const std::pair<double, double>& pj = charmm36_entry(elem_j);
    double* out = internal::new_double_view(2, out_view, n_out_view);
    if (out == NULL) return;
    out[0] = pi.first + pj.first;
    out[1] = std::sqrt(pi.second * pj.second);
}

LJArrays lj_parameter_arrays(const std::vector<std::string>& elements) {
    LJArrays out;
    out.rmin_half.resize(elements.size());
    out.epsilon.resize(elements.size());
    for (std::size_t i = 0; i < elements.size(); i++) {
        const auto& p = charmm36_entry(elements[i]);
        out.rmin_half[i] = p.first;
        out.epsilon[i] = p.second;
    }
    return out;
}

LJArrays scaled_parameters(const std::vector<std::string>& elements,
                           double sigma_scaling, double epsilon_scaling) {
    LJArrays out = lj_parameter_arrays(elements);
    for (auto& v : out.rmin_half) v *= sigma_scaling;
    for (auto& v : out.epsilon) v *= epsilon_scaling;
    return out;
}

double lj_score(double r, double rmin, double eps) {
    if (rmin == 0.0 || eps == 0.0) return 0.0;
    const double ri = std::max(r, 0.01);
    const double ratio6 = std::pow(rmin / ri, 6.0);
    const double e = eps * (ratio6 * ratio6 - 2.0 * ratio6);
    return r < rmin ? e : 0.0;
}

LJArrays cross_lj_params(const std::vector<std::string>& elements_a,
                         const std::vector<std::string>& elements_b) {
    LJArrays out;
    out.rmin_half.reserve(elements_a.size() * elements_b.size());
    out.epsilon.reserve(elements_a.size() * elements_b.size());
    for (const auto& ea : elements_a) {
        for (const auto& eb : elements_b) {
            const auto& pi = charmm36_entry(ea);
            const auto& pj = charmm36_entry(eb);
            out.rmin_half.push_back(pi.first + pj.first);
            out.epsilon.push_back(std::sqrt(pi.second * pj.second));
        }
    }
    return out;
}

double lj_pairs_sum(const std::vector<double>& coords_a,
                    const std::vector<double>& coords_b,
                    const std::vector<double>& rmin_ij,
                    const std::vector<double>& eps_ij,
                    double r_cutoff) {
    const std::size_t na = coords_a.size() / 3;
    const std::size_t nb = coords_b.size() / 3;
    double total = 0.0;
    std::size_t k = 0;
    for (std::size_t i = 0; i < na; i++) {
        for (std::size_t j = 0; j < nb; j++, k++) {
            const double dx = coords_a[i * 3] - coords_b[j * 3];
            const double dy = coords_a[i * 3 + 1] - coords_b[j * 3 + 1];
            const double dz = coords_a[i * 3 + 2] - coords_b[j * 3 + 2];
            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (r >= r_cutoff || r >= rmin_ij[k]) continue;
            const double ri = std::max(r, 0.01);
            const double ratio6 = std::pow(rmin_ij[k] / ri, 6.0);
            total += eps_ij[k] * (ratio6 * ratio6 - 2.0 * ratio6);
        }
    }
    return total;
}

void pair_energy_matrix(const std::vector<double>& coords_a,
                        const std::vector<double>& coords_b,
                        const std::vector<std::string>& elements_a,
                        const std::vector<std::string>& elements_b,
                        int n_a_conf, int n_b_conf,
                        double** out_view, int* n_out_view,
                        double r_cutoff, double aabb_pad) {
    const LJArrays cross = cross_lj_params(elements_a, elements_b);
    const int n_a_atoms = static_cast<int>(elements_a.size());
    const int n_b_atoms = static_cast<int>(elements_b.size());
    std::vector<double> ca = coords_a, cb = coords_b;
    pair_energy_matrix_kernel(ca, cb, cross.rmin_half, cross.epsilon,
                               n_a_conf, n_a_atoms, n_b_conf, n_b_atoms,
                               out_view, n_out_view, r_cutoff, aabb_pad);
}

std::vector<double> lj_energy(const std::vector<double>& r,
                                const std::vector<double>& rmin,
                                const std::vector<double>& eps,
                                bool repulsive_only, double cutoff,
                                double r_floor) {
    size_t n = r.size();
    std::vector<double> energy(n);
    bool has_cutoff = cutoff > 0.0;
    for (size_t i = 0; i < n; i++) {
        double ri = std::max(r[i], r_floor);
        double ratio6 = std::pow(rmin[i] / ri, 6.0);
        double e = eps[i] * (ratio6 * ratio6 - 2.0 * ratio6);
        if (repulsive_only && r[i] >= rmin[i]) e = 0.0;
        if (has_cutoff && r[i] >= cutoff) e = 0.0;
        energy[i] = e;
    }
    return energy;
}

std::vector<double> boltzmann_weights(const std::vector<double>& energies,
                                       double temperature) {
    const double KB = 0.0019872041;
    double kt = KB * temperature;
    double e_min = *std::min_element(energies.begin(), energies.end());
    std::vector<double> w(energies.size());
    double sum = 0.0;
    for (size_t i = 0; i < energies.size(); i++) {
        w[i] = std::exp(-(energies[i] - e_min) / kt);
        sum += w[i];
    }
    if (sum > 0.0) {
        for (auto& v : w) v /= sum;
    }
    return w;
}

std::vector<double> cluster_weights(
        const std::vector<int>& assignments,
        const std::vector<double>& frame_weights,
        int n_clusters) {
    std::vector<double> w(n_clusters, 0.0);
    for (size_t i = 0; i < assignments.size(); i++) {
        if (assignments[i] >= 0 && assignments[i] < n_clusters) {
            w[assignments[i]] += frame_weights[i];
        }
    }
    double sum = std::accumulate(w.begin(), w.end(), 0.0);
    if (sum > 0.0) {
        for (auto& v : w) v /= sum;
    }
    return w;
}

std::vector<double> aabb_build(const std::vector<double>& coords,
                                 int n_frames, int n_atoms, double pad) {
    std::vector<double> boxes(n_frames * 6);
    for (int f = 0; f < n_frames; f++) {
        const double* p = &coords[f * n_atoms * 3];
        double mn[3] = {p[0], p[1], p[2]};
        double mx[3] = {p[0], p[1], p[2]};
        for (int a = 1; a < n_atoms; a++) {
            for (int d = 0; d < 3; d++) {
                double v = p[a * 3 + d];
                if (v < mn[d]) mn[d] = v;
                if (v > mx[d]) mx[d] = v;
            }
        }
        boxes[f * 6 + 0] = mn[0] - pad;
        boxes[f * 6 + 1] = mn[1] - pad;
        boxes[f * 6 + 2] = mn[2] - pad;
        boxes[f * 6 + 3] = mx[0] + pad;
        boxes[f * 6 + 4] = mx[1] + pad;
        boxes[f * 6 + 5] = mx[2] + pad;
    }
    return boxes;
}

std::vector<double> aabb_build_single(const std::vector<double>& coords,
                                        int n_atoms, double pad) {
    double mn[3] = {coords[0], coords[1], coords[2]};
    double mx[3] = {coords[0], coords[1], coords[2]};
    for (int a = 1; a < n_atoms; a++) {
        for (int d = 0; d < 3; d++) {
            double v = coords[a * 3 + d];
            if (v < mn[d]) mn[d] = v;
            if (v > mx[d]) mx[d] = v;
        }
    }
    return {mn[0] - pad, mn[1] - pad, mn[2] - pad,
            mx[0] + pad, mx[1] + pad, mx[2] + pad};
}

std::vector<int> aabb_intersects_reference(
        const std::vector<double>& rot_boxes,
        const std::vector<double>& ref_box) {
    int n_frames = rot_boxes.size() / 6;
    std::vector<int> mask(n_frames, 0);
    for (int f = 0; f < n_frames; f++) {
        const double* b = &rot_boxes[f * 6];
        bool ox = b[3] >= ref_box[0] && b[0] <= ref_box[3];
        bool oy = b[4] >= ref_box[1] && b[1] <= ref_box[4];
        bool oz = b[5] >= ref_box[2] && b[2] <= ref_box[5];
        mask[f] = (ox && oy && oz) ? 1 : 0;
    }
    return mask;
}

// --------------------------------------------------------------------------
// The selector mini-language and the masks built from it
// --------------------------------------------------------------------------

namespace {

//! The clauses of `A [and resname R] ...`, whitespace-stripped.
std::vector<std::string> selector_parts(const std::string& selector) {
    std::vector<std::string> parts;
    std::size_t start = 0;
    while (true) {
        const std::size_t at = selector.find(" and ", start);
        parts.push_back(strip(selector.substr(
                start, at == std::string::npos ? std::string::npos : at - start)));
        if (at == std::string::npos) break;
        start = at + 5;
    }
    return parts;
}

//! The atom-name clause of a selector: its first part, sans a `name ` prefix.
std::string selector_atom_clause(const std::string& selector) {
    std::string atom = selector_parts(selector).empty()
            ? strip(selector) : selector_parts(selector)[0];
    atom = upper(atom);
    if (atom.rfind("NAME ", 0) == 0) atom = strip(atom.substr(5));
    return atom;
}

}  // namespace

bool selector_matches(const std::string& selector,
                      const std::string& atom_name,
                      const std::string& resname) {
    if (selector_atom_clause(selector) != upper(atom_name)) return false;
    std::string selector_resname;
    const std::vector<std::string> parts = selector_parts(selector);
    for (std::size_t i = 1; i < parts.size(); i++) {
        std::istringstream tokens(parts[i]);
        std::string head, value;
        tokens >> head >> value;
        if (upper(head) == "RESNAME") selector_resname = upper(value);
    }
    return selector_resname.empty() ||
            (!resname.empty() && upper(resname) == selector_resname);
}

std::vector<std::string> selector_resnames(
        const std::vector<std::string>& selectors) {
    std::vector<std::string> out;
    for (const auto& selector : selectors) {
        for (const auto& part : selector_parts(selector)) {
            std::istringstream tokens(part);
            std::string head, value;
            tokens >> head >> value;
            if (upper(head) == "RESNAME" && !value.empty()) {
                const std::string r = upper(value);
                if (std::find(out.begin(), out.end(), r) == out.end()) {
                    out.push_back(r);
                }
            }
        }
    }
    return out;
}

std::vector<std::string> selector_atom_names(
        const std::vector<std::string>& selectors,
        const std::vector<std::string>& atom_names,
        const std::vector<std::string>& resnames) {
    std::vector<std::string> out;
    if (selectors.empty()) return out;
    if (selector_resnames(selectors).size() > 1) return out;  // ambiguous

    auto add_unique = [&out](const std::string& name) {
        if (std::find(out.begin(), out.end(), name) == out.end()) {
            out.push_back(name);
        }
    };
    if (resnames.empty()) {
        for (const auto& selector : selectors) {
            add_unique(selector_atom_clause(selector));
        }
        return out;
    }
    for (std::size_t i = 0; i < atom_names.size() && i < resnames.size(); i++) {
        for (const auto& selector : selectors) {
            if (selector_matches(selector, atom_names[i], resnames[i])) {
                add_unique(upper(atom_names[i]));
                break;
            }
        }
    }
    return out;
}

std::vector<int> hydrogen_mask(const std::vector<std::string>& atom_names) {
    std::vector<int> out(atom_names.size(), 0);
    for (std::size_t i = 0; i < atom_names.size(); i++) {
        out[i] = atom_type(atom_names[i]) == "H" ? 1 : 0;
    }
    return out;
}

std::vector<int> site_mask(const std::vector<std::string>& atom_names,
                           const std::vector<int>& residue_indices,
                           int site_residue, const std::string& site_chain,
                           const std::vector<std::string>& chain_ids,
                           bool mask_backbone) {
    const std::size_t n = atom_names.size();
    std::vector<int> out(n, 0);
    for (std::size_t i = 0; i < n; i++) {
        const int residx = i < residue_indices.size() ? residue_indices[i] : -1;
        const bool same_chain = site_chain.empty() || i >= chain_ids.size() ||
                upper(chain_ids[i]) == upper(site_chain);
        if (site_residue >= 0 && same_chain && residx == site_residue) {
            out[i] = 1;
            continue;
        }
        if (mask_backbone) {
            const std::string name = upper(atom_names[i]);
            if (name == "CA" || name == "C" || name == "N" || name == "O") {
                out[i] = 1;
            }
        }
    }
    return out;
}

std::vector<double> protein_charge_mask(
        const std::vector<std::string>& atom_names,
        const std::vector<std::string>& resnames) {
    std::vector<double> q(atom_names.size(), 0.0);
    for (std::size_t i = 0; i < atom_names.size() && i < resnames.size(); i++) {
        const std::string name = upper(atom_names[i]);
        const std::string res = upper(resnames[i]);
        if ((name == "CZ" && res == "ARG") || (name == "NZ" && res == "LYS")) {
            q[i] = 1.0;
        } else if ((name == "CG" && res == "ASP") ||
                   (name == "CD" && res == "GLU")) {
            q[i] = -1.0;
        } else if ((name == "ND1" || name == "NE2") && res == "HIS") {
            q[i] = 0.25;
        }
    }
    return q;
}

std::vector<double> rotamer_charge_mask(
        const std::vector<std::string>& atom_names,
        const std::vector<std::string>& positive,
        const std::vector<std::string>& negative,
        const std::vector<std::string>& resnames) {
    std::vector<std::string> pos = positive, neg = negative;
    if (selector_resnames(pos).size() > 1) pos.clear();
    if (selector_resnames(neg).size() > 1) neg.clear();
    std::vector<double> q(atom_names.size(), 0.0);
    for (std::size_t i = 0; i < atom_names.size(); i++) {
        const std::string resname =
                i < resnames.size() ? resnames[i] : std::string();
        bool is_negative = false, is_positive = false;
        for (const auto& s : neg) {
            if (selector_matches(s, atom_names[i], resname)) {
                is_negative = true;
                break;
            }
        }
        if (!is_negative) {
            for (const auto& s : pos) {
                if (selector_matches(s, atom_names[i], resname)) {
                    is_positive = true;
                    break;
                }
            }
        }
        if (is_negative) q[i] = -1.0;
        else if (is_positive) q[i] = 0.5;
    }
    return q;
}

// --------------------------------------------------------------------------
// BoundingBoxFilter
// --------------------------------------------------------------------------

std::vector<double> BoundingBoxFilter::build(const std::vector<double>& coords,
                                             int n_frames, int n_atoms) const {
    return aabb_build(coords, n_frames, n_atoms, pad_);
}

std::vector<double> BoundingBoxFilter::build_single(
        const std::vector<double>& coords, int n_atoms) const {
    return aabb_build_single(coords, n_atoms, pad_);
}

bool BoundingBoxFilter::intersects(const std::vector<double>& a,
                                   const std::vector<double>& b) {
    return a[3] >= b[0] && a[0] <= b[3] && a[4] >= b[1] && a[1] <= b[4] &&
            a[5] >= b[2] && a[2] <= b[5];
}

std::vector<int> BoundingBoxFilter::intersects_reference(
        const std::vector<double>& rot_boxes,
        const std::vector<double>& ref_box) const {
    return aabb_intersects_reference(rot_boxes, ref_box);
}

AABBFilterResult BoundingBoxFilter::filter_frames(
        const std::vector<double>& coords, int n_frames, int n_atoms,
        const std::vector<double>& reference_coords, int n_ref_atoms) const {
    const std::vector<double> boxes = aabb_build(coords, n_frames, n_atoms, pad_);
    const std::vector<double> ref =
            aabb_build_single(reference_coords, n_ref_atoms, pad_);
    const std::vector<int> mask = aabb_intersects_reference(boxes, ref);

    AABBFilterResult out;
    out.mask = mask;
    int n_kept = 0;
    for (int f = 0; f < n_frames; f++) n_kept += mask[f];
    out.coords.resize(static_cast<std::size_t>(n_kept) * n_atoms * 3);
    int k = 0;
    for (int f = 0; f < n_frames; f++) {
        if (!mask[f]) continue;
        for (int a = 0; a < n_atoms * 3; a++) {
            out.coords[static_cast<std::size_t>(k) * n_atoms * 3 + a] =
                    coords[static_cast<std::size_t>(f) * n_atoms * 3 + a];
        }
        k++;
    }
    return out;
}

RotamerScoreResult get_rotamer_score(
        const std::vector<double>& rotamer_coords,
        const std::vector<double>& protein_coords,
        const std::vector<std::string>& protein_atom_names,
        const std::vector<std::string>& protein_resnames,
        const std::vector<std::string>& rotamer_atom_names,
        const std::vector<std::string>& positive,
        const std::vector<std::string>& negative,
        const std::vector<std::string>& rotamer_resnames,
        const std::vector<int>& protein_residue_indices,
        const std::vector<std::string>& protein_chain_ids,
        int site_residue, const std::string& site_chain,
        const std::vector<double>& rotamer_weights,
        double temperature, bool ignore_h, bool electrostatic,
        const std::string& potential, double sigma_scaling,
        double epsilon_scaling) {
    const int potential_code = potential == "lj" ? ROTAMER_POTENTIAL_LJ
            : potential == "gauss" ? ROTAMER_POTENTIAL_GAUSS
            : -1;
    if (potential_code < 0) {
        IMP_THROW("Unknown potential '" << potential << "'", ValueException);
    }

    const int n_dye = static_cast<int>(rotamer_atom_names.size());
    const int n_prot = static_cast<int>(protein_atom_names.size());
    const int n_rotamers = n_dye == 0
            ? 0
            : static_cast<int>(rotamer_coords.size() /
                               (static_cast<std::size_t>(n_dye) * 3));

    RotamerScoreResult result;
    result.partition = 0.0;
    const auto uniform = [&]() {
        result.weights.assign(n_rotamers > 0 ? n_rotamers : 1,
                              n_rotamers > 0 ? 1.0 / n_rotamers : 1.0);
    };

    // The masks: hydrogens and (on the dye side) backbone names are never
    // sterically scored; the labelled residue is masked on the protein side.
    std::vector<int> protein_keep(n_prot, 1), rotamer_keep(n_dye, 1);
    {
        const std::vector<int> p_h = hydrogen_mask(protein_atom_names);
        const std::vector<int> p_site = site_mask(
                protein_atom_names, protein_residue_indices, site_residue,
                site_chain, protein_chain_ids, false);
        for (int i = 0; i < n_prot; i++) {
            if ((ignore_h && p_h[i]) || p_site[i]) protein_keep[i] = 0;
        }
        const std::vector<int> r_backbone =
                site_mask(rotamer_atom_names, std::vector<int>(), -1, "",
                          std::vector<std::string>(), true);
        const std::vector<int> r_h = hydrogen_mask(rotamer_atom_names);
        for (int i = 0; i < n_dye; i++) {
            if (r_backbone[i] || (ignore_h && r_h[i])) rotamer_keep[i] = 0;
        }
    }
    std::vector<int> protein_idx, rotamer_idx;
    for (int i = 0; i < n_prot; i++) {
        if (protein_keep[i]) protein_idx.push_back(i);
    }
    for (int i = 0; i < n_dye; i++) {
        if (rotamer_keep[i]) rotamer_idx.push_back(i);
    }

    if (protein_idx.empty() || rotamer_idx.empty() || n_rotamers == 0) {
        uniform();
        result.partition = 1.0;
        result.energies.assign(n_rotamers, 0.0);
        return result;
    }

    // Lorentz-Berthelot parameters over the surviving atoms.
    std::vector<std::string> protein_types, rotamer_types;
    for (int i : protein_idx) protein_types.push_back(atom_type(protein_atom_names[i]));
    for (int i : rotamer_idx) rotamer_types.push_back(atom_type(rotamer_atom_names[i]));
    const LJArrays p = scaled_parameters(protein_types, sigma_scaling, epsilon_scaling);
    const LJArrays r = scaled_parameters(rotamer_types, sigma_scaling, epsilon_scaling);
    const int nd = static_cast<int>(rotamer_idx.size());
    const int np = static_cast<int>(protein_idx.size());
    std::vector<double> rmin_ij(static_cast<std::size_t>(nd) * np),
            eps_ij(static_cast<std::size_t>(nd) * np);
    for (int i = 0; i < nd; i++) {
        for (int j = 0; j < np; j++) {
            rmin_ij[static_cast<std::size_t>(i) * np + j] =
                    r.rmin_half[i] + p.rmin_half[j];
            eps_ij[static_cast<std::size_t>(i) * np + j] =
                    std::sqrt(r.epsilon[i] * p.epsilon[j]);
        }
    }

    std::vector<double> q_rotamer, q_protein;
    if (electrostatic) {
        std::vector<std::string> dye_names, dye_res;
        for (int i : rotamer_idx) {
            dye_names.push_back(rotamer_atom_names[i]);
            dye_res.push_back(i < static_cast<int>(rotamer_resnames.size())
                                      ? rotamer_resnames[i]
                                      : std::string());
        }
        std::vector<std::string> prot_names, prot_res;
        for (int i : protein_idx) {
            prot_names.push_back(protein_atom_names[i]);
            prot_res.push_back(protein_resnames[i]);
        }
        q_rotamer = rotamer_charge_mask(dye_names, positive, negative, dye_res);
        q_protein = protein_charge_mask(prot_names, prot_res);
    }

    // Gather the dye atoms that survived the mask.
    std::vector<double> selected(
            static_cast<std::size_t>(n_rotamers) * nd * 3);
    for (int f = 0; f < n_rotamers; f++) {
        for (int i = 0; i < nd; i++) {
            for (int d = 0; d < 3; d++) {
                selected[(static_cast<std::size_t>(f) * nd + i) * 3 + d] =
                        rotamer_coords[(static_cast<std::size_t>(f) * n_dye +
                                        rotamer_idx[i]) * 3 + d];
            }
        }
    }
    std::vector<double> prot_pos(static_cast<std::size_t>(np) * 3);
    for (int j = 0; j < np; j++) {
        for (int d = 0; d < 3; d++) {
            prot_pos[static_cast<std::size_t>(j) * 3 + d] =
                    protein_coords[static_cast<std::size_t>(protein_idx[j]) * 3 + d];
        }
    }

    const std::vector<double> energies = rotamer_interaction_energies(
            selected.data(), static_cast<int>(selected.size()),
            prot_pos.data(), static_cast<int>(prot_pos.size()),
            rmin_ij.data(), static_cast<int>(rmin_ij.size()),
            eps_ij.data(), static_cast<int>(eps_ij.size()),
            q_rotamer, q_protein, n_rotamers, nd, np, potential_code);

    const double GAS_CONSTANT = 1.9858775e-3;
    const double rt = GAS_CONSTANT * temperature;
    std::vector<double> boltzmann(n_rotamers);
    result.energies.assign(n_rotamers, 0.0);
    double partition = 0.0;
    for (int i = 0; i < n_rotamers; i++) {
        const double pot = energies[i * 2];
        const double dh = energies[i * 2 + 1];
        result.energies[i] = pot + dh;
        double b = std::exp(-pot / rt - dh);
        if (!rotamer_weights.empty() &&
            i < static_cast<int>(rotamer_weights.size())) {
            b *= rotamer_weights[i];
        }
        if (std::isnan(b) || std::isinf(b)) b = 0.0;
        boltzmann[i] = b;
        partition += b;
    }
    if (partition <= 0.0) {
        uniform();
        result.partition = 0.0;
        return result;
    }
    result.partition = partition;
    result.weights.resize(n_rotamers);
    for (int i = 0; i < n_rotamers; i++) result.weights[i] = boltzmann[i] / partition;
    return result;
}

// --------------------------------------------------------------------------
// Mean-field weights
// --------------------------------------------------------------------------

std::vector<double> rotamer_mean_field_weights(
        const std::vector<double>& rotamer_coords,
        const std::vector<double>& initial_weights,
        const std::vector<double>& protein_coords,
        const std::vector<std::string>& probe_elements,
        const std::vector<std::string>& protein_elements,
        double K, int n_iter, double aabb_pad, double r_cutoff) {
    const int n_clusters = static_cast<int>(initial_weights.size());
    const int n_probe_atoms = static_cast<int>(probe_elements.size());
    std::vector<double> q = initial_weights;

    // (n_clusters, 1): the protein is one "conformer"
    double* e_bb = NULL;
    int n_e_bb = 0;
    pair_energy_matrix(rotamer_coords, protein_coords, probe_elements,
                       protein_elements, n_clusters, 1, &e_bb, &n_e_bb,
                       r_cutoff, aabb_pad);
    const std::vector<double> E_bb(e_bb, e_bb + (n_e_bb > 0 ? n_e_bb : 0));
    std::free(e_bb);

    for (int it = 0; it < n_iter; it++) {
        double log_q_max = -std::numeric_limits<double>::infinity();
        std::vector<double> log_q(n_clusters);
        for (int c = 0; c < n_clusters; c++) {
            log_q[c] = std::log(std::max(q[c], 1e-300)) - K * E_bb[c];
            log_q_max = std::max(log_q_max, log_q[c]);
        }
        double denom = 0.0;
        std::vector<double> q_new(n_clusters);
        for (int c = 0; c < n_clusters; c++) {
            q_new[c] = std::exp(log_q[c] - log_q_max);
            denom += q_new[c];
        }
        if (denom > 0.0) {
            for (int c = 0; c < n_clusters; c++) q[c] = q_new[c] / denom;
        } else {
            q.assign(n_clusters, 1.0 / n_clusters);
        }
    }
    return q;
}

// --------------------------------------------------------------------------
// Typed-system walkers
// --------------------------------------------------------------------------

std::map<std::string, std::string> site_element_map(
        const ProbeForceFieldSystem& system) {
    std::map<std::string, std::string> out;
    for (const auto& s : system.get_sites()) {
        // The site's own element when it has one -- a producer that read a
        // real element column knows better than any name rule -- and the name
        // rule otherwise, from `element_from_atom_name` rather than from a
        // fourth copy of it inline here. On every structure in the tree the
        // two agree, because the MOL2 reader derives the element with exactly
        // this rule; they part company on a halogen, where the name rule
        // calls CL3 carbon and a `_atom_site.type_symbol` does not.
        out[s.id] = s.element.empty() ? element_from_atom_name(s.atom_name)
                                      : s.element;
    }
    return out;
}

std::vector<LJSitePair> get_lj_pair_sites(
        const ProbeForceFieldSystem& system) {
    const std::set<std::pair<std::string, std::string>> excluded =
            system.get_exclusions();
    const std::map<std::string, std::string> elem = site_element_map(system);
    const std::vector<FFSite>& sites = system.get_sites();

    std::vector<LJSitePair> out;
    for (std::size_t i = 0; i < sites.size(); i++) {
        for (std::size_t j = i + 1; j < sites.size(); j++) {
            const std::string& a = sites[i].id;
            const std::string& b = sites[j].id;
            if (excluded.count({a, b}) || excluded.count({b, a})) continue;
            auto ea = elem.find(a);
            auto eb = elem.find(b);
            const auto& pi = charmm36_entry(ea == elem.end() ? "C" : ea->second);
            const auto& pj = charmm36_entry(eb == elem.end() ? "C" : eb->second);
            LJSitePair pair;
            pair.site_a = a;
            pair.site_b = b;
            pair.rmin = pi.first + pj.first;
            pair.eps = std::sqrt(pi.second * pj.second);
            out.push_back(pair);
        }
    }
    return out;
}

std::map<std::string, FFLJType> get_lj_type_table(
        const std::vector<std::string>& elements) {
    std::map<std::string, FFLJType> out;
    for (const auto& elem : elements) {
        const auto& p = charmm36_entry(elem);
        FFLJType t;
        t.element = elem;
        t.rmin_half = p.first;
        t.epsilon = p.second;
        out["LJ_" + elem] = t;
    }
    return out;
}

IntramolecularEnergy::IntramolecularEnergy(
        const ProbeForceFieldSystem& system) {
    pairs_ = get_lj_pair_sites(system);
    std::map<std::string, int> id_to_idx;
    const std::vector<FFSite>& sites = system.get_sites();
    for (std::size_t i = 0; i < sites.size(); i++) id_to_idx[sites[i].id] = i;
    idx_a_.reserve(pairs_.size());
    idx_b_.reserve(pairs_.size());
    rmin_.reserve(pairs_.size());
    eps_.reserve(pairs_.size());
    for (const auto& p : pairs_) {
        idx_a_.push_back(id_to_idx[p.site_a]);
        idx_b_.push_back(id_to_idx[p.site_b]);
        rmin_.push_back(p.rmin);
        eps_.push_back(p.eps);
    }
}

double IntramolecularEnergy::evaluate(
        const std::vector<double>& coords, int n_atoms) const {
    if (pairs_.empty()) return 0.0;
    std::vector<double> c = coords;
    return lj_pair_energies(c.data(), static_cast<int>(c.size()), idx_a_, idx_b_,
                            rmin_, eps_, 1, n_atoms,
                            static_cast<int>(pairs_.size()), true)[0];
}

std::vector<double> IntramolecularEnergy::evaluate_batch(
        const std::vector<double>& coords, int n_frames, int n_atoms) const {
    if (pairs_.empty()) return std::vector<double>(n_frames, 0.0);
    std::vector<double> c = coords;
    return lj_pair_energies(c.data(), static_cast<int>(c.size()), idx_a_, idx_b_,
                            rmin_, eps_, n_frames, n_atoms,
                            static_cast<int>(pairs_.size()), true);
}

EnergyMaskResult IntramolecularEnergy::evaluate_batch_filtered(
        const std::vector<double>& coords, int n_frames, int n_atoms,
        const std::vector<double>& reference_coords, int n_ref_atoms,
        double pad) const {
    const BoundingBoxFilter bbf(pad);
    const AABBFilterResult kept =
            bbf.filter_frames(coords, n_frames, n_atoms, reference_coords,
                              n_ref_atoms);
    EnergyMaskResult out;
    out.energies.assign(n_frames, 0.0);
    out.mask = kept.mask;
    if (!kept.coords.empty()) {
        const int n_kept = static_cast<int>(kept.coords.size() / (n_atoms * 3));
        const std::vector<double> kept_energies =
                evaluate_batch(kept.coords, n_kept, n_atoms);
        int k = 0;
        for (int f = 0; f < n_frames; f++) {
            if (kept.mask[f]) out.energies[f] = kept_energies[k++];
        }
    }
    return out;
}

namespace {
//! `pair_energy_matrix` into a vector, so the caller is not handed a buffer.
std::vector<double> pair_energies(const std::vector<double>& a,
                                  const std::vector<double>& b,
                                  const std::vector<std::string>& ea,
                                  const std::vector<std::string>& eb,
                                  int n_a, int n_b, double r_cutoff,
                                  double aabb_pad) {
    double* buffer = NULL;
    int n = 0;
    pair_energy_matrix(a, b, ea, eb, n_a, n_b, &buffer, &n, r_cutoff,
                       aabb_pad);
    std::vector<double> out(buffer, buffer + n);
    std::free(buffer);
    return out;
}
}  // namespace

std::vector<std::vector<double> > rotamer_mean_field_weights_multi_probe(
        const std::vector<std::vector<double> >& rotamer_coords_list,
        const std::vector<std::vector<double> >& initial_weights_list,
        const std::vector<double>& protein_coords,
        const std::vector<std::vector<std::string> >& probe_elements_list,
        const std::vector<std::string>& protein_elements,
        double K, int n_iter, double aabb_pad, double r_cutoff) {
    const std::size_t n_dyes = rotamer_coords_list.size();
    if (initial_weights_list.size() != n_dyes ||
        probe_elements_list.size() != n_dyes) {
        IMP_THROW("rotamer_mean_field_weights_multi_probe: " << n_dyes
                  << " coordinate sets, " << initial_weights_list.size()
                  << " weight sets and " << probe_elements_list.size()
                  << " element sets", IMP::ValueException);
    }
    std::vector<std::vector<double> > q = initial_weights_list;

    // Each dye against the protein: one column per conformer.
    std::vector<std::vector<double> > e_bb(n_dyes);
    for (std::size_t d = 0; d < n_dyes; ++d) {
        e_bb[d] = pair_energies(rotamer_coords_list[d], protein_coords,
                                probe_elements_list[d], protein_elements,
                                static_cast<int>(q[d].size()), 1, r_cutoff,
                                aabb_pad);
    }
    // And each dye against every other: a matrix per ordered pair, the
    // transpose shared rather than recomputed.
    std::vector<std::vector<double> > e_sc(n_dyes * n_dyes);
    for (std::size_t d1 = 0; d1 < n_dyes; ++d1) {
        for (std::size_t d2 = d1 + 1; d2 < n_dyes; ++d2) {
            const std::vector<double> m = pair_energies(
                    rotamer_coords_list[d1], rotamer_coords_list[d2],
                    probe_elements_list[d1], probe_elements_list[d2],
                    static_cast<int>(q[d1].size()),
                    static_cast<int>(q[d2].size()), r_cutoff, aabb_pad);
            e_sc[d1 * n_dyes + d2] = m;
            std::vector<double> t(m.size());
            for (std::size_t i = 0; i < q[d1].size(); ++i) {
                for (std::size_t j = 0; j < q[d2].size(); ++j) {
                    t[j * q[d1].size() + i] = m[i * q[d2].size() + j];
                }
            }
            e_sc[d2 * n_dyes + d1] = t;
        }
    }

    for (int iteration = 0; iteration < n_iter; ++iteration) {
        std::vector<std::vector<double> > e_p = e_bb;
        for (std::size_t d1 = 0; d1 < n_dyes; ++d1) {
            for (std::size_t d2 = 0; d2 < n_dyes; ++d2) {
                if (d1 == d2) continue;
                const std::vector<double>& m = e_sc[d1 * n_dyes + d2];
                for (std::size_t i = 0; i < q[d1].size(); ++i) {
                    double acc = 0.0;
                    for (std::size_t j = 0; j < q[d2].size(); ++j) {
                        acc += m[i * q[d2].size() + j] * q[d2][j];
                    }
                    e_p[d1][i] += acc;
                }
            }
        }
        for (std::size_t d = 0; d < n_dyes; ++d) {
            std::vector<double> log_q(q[d].size());
            double top = -std::numeric_limits<double>::infinity();
            for (std::size_t i = 0; i < q[d].size(); ++i) {
                log_q[i] = std::log(std::max(q[d][i], 1e-300)) - K * e_p[d][i];
                top = std::max(top, log_q[i]);
            }
            double denom = 0.0;
            for (std::size_t i = 0; i < q[d].size(); ++i) {
                log_q[i] = std::exp(log_q[i] - top);
                denom += log_q[i];
            }
            for (std::size_t i = 0; i < q[d].size(); ++i) {
                q[d][i] = denom > 0.0 ? log_q[i] / denom
                                      : 1.0 / static_cast<double>(q[d].size());
            }
        }
    }
    return q;
}

IMP::core::Cosine* torsion_cosine(const FFTorsionType& type) {
    // CHARMM's k(1 + cos(n phi - delta)) against Cosine's k(1 - cos(...)):
    // the same curve with the phase shifted by pi.
    return new IMP::core::Cosine(type.k, type.periodicity,
                                 type.phase + IMP::algebra::PI);
}

IMP::Restraints create_probe_restraints(
        IMP::Model* model, const ProbeForceFieldSystem& system,
        const std::vector<std::string>& site_ids,
        const IMP::ParticleIndexes& particles, bool nonbonded) {
    if (site_ids.size() != particles.size()) {
        IMP_THROW("create_probe_restraints: " << site_ids.size() << " site ids "
                  << "against " << particles.size() << " particles",
                  IMP::ValueException);
    }
    std::map<std::string, IMP::ParticleIndex> site;
    for (std::size_t i = 0; i < site_ids.size(); ++i) {
        site[site_ids[i]] = particles[i];
    }
    const bool has_all = true;
    IMP::Restraints out;

    // -- bonds ------------------------------------------------------------
    const std::map<std::string, double>& bt = system.get_bond_types();
    for (std::size_t i = 0; i < system.get_bonds().size(); ++i) {
        const FFBond& b = system.get_bonds()[i];
        std::map<std::string, IMP::ParticleIndex>::const_iterator a =
                site.find(b.site_a), c = site.find(b.site_b);
        if (a == site.end() || c == site.end()) continue;
        std::map<std::string, double>::const_iterator k = bt.find(b.type_id);
        if (k == bt.end()) continue;
        // A length of zero is a length nobody set: restrain about the
        // geometry as it stands rather than pulling the sites together.
        const double length =
                b.length > 0.0
                        ? b.length
                        : IMP::core::get_distance(
                                  IMP::core::XYZ(model, a->second),
                                  IMP::core::XYZ(model, c->second));
        out.push_back(new IMP::core::DistanceRestraint(
                model, new IMP::core::Harmonic(length, k->second),
                model->get_particle(a->second),
                model->get_particle(c->second)));
    }

    // -- angles -----------------------------------------------------------
    const std::map<std::string, double>& at = system.get_angle_types();
    for (std::size_t i = 0; i < system.get_angles().size(); ++i) {
        const FFAngle& an = system.get_angles()[i];
        std::map<std::string, IMP::ParticleIndex>::const_iterator a =
                site.find(an.site_a), b = site.find(an.site_b),
                c = site.find(an.site_c);
        if (a == site.end() || b == site.end() || c == site.end()) continue;
        std::map<std::string, double>::const_iterator k = at.find(an.type_id);
        if (k == at.end()) continue;
        const double theta =
                an.theta > 0.0
                        ? an.theta
                        : bond_angle_rad(
                                  IMP::core::XYZ(model, a->second)
                                          .get_coordinates(),
                                  IMP::core::XYZ(model, b->second)
                                          .get_coordinates(),
                                  IMP::core::XYZ(model, c->second)
                                          .get_coordinates());
        out.push_back(new IMP::core::AngleRestraint(
                model, new IMP::core::Harmonic(theta, k->second),
                model->get_particle(a->second),
                model->get_particle(b->second),
                model->get_particle(c->second)));
    }

    // -- torsions ---------------------------------------------------------
    const std::map<std::string, FFTorsionType>& tt = system.get_torsion_types();
    for (std::size_t i = 0; i < system.get_dihedrals().size(); ++i) {
        const FFTorsion& t = system.get_dihedrals()[i];
        std::map<std::string, IMP::ParticleIndex>::const_iterator
                a = site.find(t.site_a), b = site.find(t.site_b),
                c = site.find(t.site_c), d = site.find(t.site_d);
        if (a == site.end() || b == site.end() || c == site.end() ||
            d == site.end()) continue;
        std::map<std::string, FFTorsionType>::const_iterator ty =
                tt.find(t.type_id);
        if (ty == tt.end()) continue;
        out.push_back(new IMP::core::DihedralRestraint(
                model, torsion_cosine(ty->second),
                model->get_particle(a->second),
                model->get_particle(b->second),
                model->get_particle(c->second),
                model->get_particle(d->second)));
    }

    // -- impropers: harmonic about the geometry as it stands now ----------
    const std::map<std::string, FFTorsionType>& it =
            system.get_improper_types();
    for (std::size_t i = 0; i < system.get_impropers().size(); ++i) {
        const FFTorsion& t = system.get_impropers()[i];
        std::map<std::string, IMP::ParticleIndex>::const_iterator
                a = site.find(t.site_a), b = site.find(t.site_b),
                c = site.find(t.site_c), d = site.find(t.site_d);
        if (a == site.end() || b == site.end() || c == site.end() ||
            d == site.end()) continue;
        std::map<std::string, FFTorsionType>::const_iterator ty =
                it.find(t.type_id);
        if (ty == it.end()) continue;
        const IMP::core::XYZ xa(model, a->second), xb(model, b->second),
                xc(model, c->second), xd(model, d->second);
        const double theta0 = IMP::core::get_dihedral(xa, xb, xc, xd);
        out.push_back(new IMP::core::DihedralRestraint(
                model, new IMP::core::Harmonic(theta0, ty->second.k),
                model->get_particle(a->second),
                model->get_particle(b->second),
                model->get_particle(c->second),
                model->get_particle(d->second)));
    }

    // -- repulsion --------------------------------------------------------
    // One soft-sphere restraint over every non-excluded pair: the same term a
    // Monte-Carlo step is scored against, and differentiable, so a dynamics
    // run uses it too. It was per-pair Lennard-Jones lower bounds here and
    // soft spheres there -- two implementations of one piece of physics, and
    // thousands of restraints where one does.
    if (nonbonded) {
        IMP::Restraint* steric =
                create_steric_restraint(model, system, site_ids, particles);
        if (steric != NULL) out.push_back(steric);
    }
    (void)has_all;
    return out;
}

namespace {

//! The mean position of a set of particles.
IMP::algebra::Vector3D centre_of(IMP::Model* model,
                                 const IMP::ParticleIndexes& ps) {
    IMP::algebra::Vector3D c(0.0, 0.0, 0.0);
    for (std::size_t i = 0; i < ps.size(); ++i) {
        c += IMP::core::XYZ(model, ps[i]).get_coordinates();
    }
    return ps.empty() ? c : c / static_cast<double>(ps.size());
}

}  // namespace

double place_guest_by_score(IMP::ScoringFunction* scoring_function,
                            IMP::Model* model,
                            const IMP::ParticleIndexes& host,
                            const IMP::ParticleIndexes& guest, double distance,
                            int n_trials, int seed) {
    if (guest.empty() || scoring_function == NULL) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const IMP::algebra::Vector3D host_centre = centre_of(model, host);
    const IMP::algebra::Vector3D guest_centre = centre_of(model, guest);
    IMP::algebra::Vector3Ds shape;
    for (std::size_t i = 0; i < guest.size(); ++i) {
        shape.push_back(IMP::core::XYZ(model, guest[i]).get_coordinates() -
                        guest_centre);
    }

    boost::mt19937 rng(static_cast<boost::uint32_t>(seed));
    boost::uniform_real<double> unit(0.0, 1.0);

    IMP::algebra::Vector3Ds best = shape;
    IMP::algebra::Vector3D best_offset = host_centre;
    double best_score = std::numeric_limits<double>::infinity();
    for (int trial = 0; trial < std::max(1, n_trials); ++trial) {
        const IMP::algebra::Rotation3D rotation =
                IMP::algebra::get_random_rotation_3d();
        // A direction drawn uniformly on the sphere: z uniform in [-1, 1] and
        // the azimuth uniform, which is the one construction that does not
        // crowd the poles.
        const double z = 2.0 * unit(rng) - 1.0;
        const double azimuth = 2.0 * IMP::PI * unit(rng);
        const double r = std::sqrt(std::max(0.0, 1.0 - z * z));
        const IMP::algebra::Vector3D direction(r * std::cos(azimuth),
                                               r * std::sin(azimuth), z);
        const IMP::algebra::Vector3D offset = host_centre + direction * distance;

        IMP::algebra::Vector3Ds pose;
        for (std::size_t i = 0; i < shape.size(); ++i) {
            pose.push_back(offset + rotation.get_rotated(shape[i]));
        }
        for (std::size_t i = 0; i < guest.size(); ++i) {
            IMP::core::XYZ(model, guest[i]).set_coordinates(pose[i]);
        }
        const double score = scoring_function->evaluate(false);
        if (score < best_score) {
            best_score = score;
            best = pose;
        }
    }
    for (std::size_t i = 0; i < guest.size(); ++i) {
        IMP::core::XYZ(model, guest[i]).set_coordinates(best[i]);
    }
    (void)best_offset;
    return best_score;
}

IMP::Restraints create_go_restraints(
        IMP::Model* model, const ProbeForceFieldSystem& system,
        const std::vector<std::string>& site_ids,
        const IMP::ParticleIndexes& particles,
        const std::map<std::string, std::string>& site_atom_names,
        const std::string& component,
        const std::vector<std::string>& only_sites, double k, double cutoff) {
    IMP::Restraints out;
    if (site_ids.size() != particles.size()) {
        IMP_THROW("create_go_restraints: " << site_ids.size() << " site ids "
                  << "against " << particles.size() << " particles",
                  IMP::ValueException);
    }
    std::map<std::string, IMP::ParticleIndex> site;
    for (std::size_t i = 0; i < site_ids.size(); ++i) {
        site[site_ids[i]] = particles[i];
    }
    const std::set<std::string> released(only_sites.begin(), only_sites.end());

    // The component's heavy sites, sorted, so the restraints come out in the
    // site ids' order and a repeated run builds the same set.
    std::vector<std::string> heavy;
    for (std::size_t i = 0; i < system.get_sites().size(); ++i) {
        const FFSite& s = system.get_sites()[i];
        if (s.component != component) continue;
        if (site.find(s.id) == site.end()) continue;
        std::map<std::string, std::string>::const_iterator name =
                site_atom_names.find(s.id);
        if (name != site_atom_names.end() && !name->second.empty() &&
            (name->second[0] == 'H' || name->second[0] == 'h')) {
            continue;
        }
        heavy.push_back(s.id);
    }
    std::sort(heavy.begin(), heavy.end());

    for (std::size_t i = 0; i < heavy.size(); ++i) {
        for (std::size_t j = i + 1; j < heavy.size(); ++j) {
            if (!released.empty() && released.count(heavy[i]) == 0 &&
                released.count(heavy[j]) == 0) {
                continue;
            }
            // Within two bonds the bonded terms already say what the distance
            // is; a contact there would be a second opinion.
            if (system.is_within_bonds(heavy[i], heavy[j], 2)) continue;
            const IMP::core::XYZ a(model, site[heavy[i]]),
                    b(model, site[heavy[j]]);
            const double d = IMP::core::get_distance(a, b);
            if (d > cutoff) continue;
            IMP::Restraint* r = new IMP::core::DistanceRestraint(
                    model, new IMP::core::Harmonic(std::max(d, 1.0), k),
                    model->get_particle(site[heavy[i]]),
                    model->get_particle(site[heavy[j]]));
            r->set_name("go_" + component + "_" + heavy[i] + "_" + heavy[j]);
            out.push_back(r);
        }
    }
    return out;
}

IMP::Restraint* create_steric_restraint(IMP::Model* model,
                                       const ProbeForceFieldSystem& system,
                                       const std::vector<std::string>& site_ids,
                                       const IMP::ParticleIndexes& particles,
                                       double k) {
    if (site_ids.size() != particles.size()) {
        IMP_THROW("create_steric_restraint: " << site_ids.size()
                  << " site ids against " << particles.size() << " particles",
                  IMP::ValueException);
    }
    std::map<std::string, IMP::ParticleIndex> site;
    for (std::size_t i = 0; i < site_ids.size(); ++i) {
        site[site_ids[i]] = particles[i];
    }
    // A system that says its non-bonded term is off has no steric restraint,
    // rather than one nobody asked for.
    if (!system.get_nonbonded().enabled) return NULL;
    const std::set<std::pair<std::string, std::string> > excluded =
            system.get_exclusions();

    // Sorted, so the pair order is the site ids' and not a hash's: a run that
    // is repeated has to build the same container.
    std::vector<std::string> ids;
    for (std::map<std::string, IMP::ParticleIndex>::const_iterator it =
                 site.begin();
         it != site.end(); ++it) {
        ids.push_back(it->first);
    }

    // A soft sphere is a sphere: a site the caller decorated without a radius
    // would contribute nothing at all, silently. The system says what radius
    // each site has, so give it that rather than score an empty term.
    for (std::map<std::string, IMP::ParticleIndex>::const_iterator it =
                 site.begin();
         it != site.end(); ++it) {
        IMP::Particle* p = model->get_particle(it->second);
        if (IMP::core::XYZR::get_is_setup(p)) continue;
        double radius = 1.7;
        for (std::size_t i = 0; i < system.get_sites().size(); ++i) {
            if (system.get_sites()[i].id != it->first) continue;
            radius = system.get_sites()[i].radius;
            break;
        }
        IMP::core::XYZR::setup_particle(p, radius);
    }

    IMP::ParticleIndexPairs pairs;
    for (std::size_t i = 0; i < ids.size(); ++i) {
        for (std::size_t j = i + 1; j < ids.size(); ++j) {
            if (excluded.count(std::make_pair(ids[i], ids[j])) > 0) continue;
            pairs.push_back(IMP::ParticleIndexPair(site[ids[i]], site[ids[j]]));
        }
    }
    if (pairs.empty()) return NULL;

    IMP_NEW(IMP::container::ListPairContainer, container, (model, pairs));
    const double strength = k >= 0.0 ? k : system.get_nonbonded().k;
    return new IMP::container::PairsRestraint(
            new IMP::core::SoftSpherePairScore(strength), container,
            "steric");
}

IMPBFF_END_NAMESPACE
