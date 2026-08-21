/**
 *  \file Scoring.cpp
 *  \brief Stage-2 scoring orchestration: CHARMM36, LJ, Boltzmann, AABB, and
 *         the end-to-end rotamer score.
 *
 * The inner kernels are in RotamerEnergy.cpp; this file is the layer that
 * builds parameters, calls them, and turns the energies into weights.
 */

#include <IMP/bff/Scoring.h>
#include <IMP/bff/RotamerEnergy.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <string>

IMPBFF_BEGIN_NAMESPACE

namespace {

const std::map<std::string, std::array<double, 2>>& charmm36_table() {
    static const std::map<std::string, std::array<double, 2>> t = {
        {"C", {2.02446316, -0.06394724}},
        {"N", {1.89285714, -0.15428571}},
        {"O", {1.693,      -0.12642017}},
        {"S", {2.1,         -0.47}},
        {"H", {0.98357778,  -0.03466645}},
    };
    return t;
}

}  // namespace

std::array<double, 2> charmm36_lj(const std::string& element) {
    auto& t = charmm36_table();
    auto it = t.find(element);
    if (it == t.end()) return t.at("C");
    return it->second;
}

std::array<double, 2> lj_cross(const std::string& elem_i,
                                const std::string& elem_j) {
    auto pi = charmm36_lj(elem_i);
    auto pj = charmm36_lj(elem_j);
    double rmin = pi[0] + pj[0];
    double eps = std::sqrt(pi[1] * pj[1]);
    return {rmin, eps};
}

void lj_parameter_arrays(const std::vector<std::string>& elements,
                          std::vector<double>& rmin_half,
                          std::vector<double>& epsilon) {
    rmin_half.resize(elements.size());
    epsilon.resize(elements.size());
    for (size_t i = 0; i < elements.size(); i++) {
        auto p = charmm36_lj(elements[i]);
        rmin_half[i] = p[0];
        epsilon[i] = p[1];
    }
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

std::vector<double> rotamer_cluster_weights(
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

RotamerScoreResult compute_rotamer_score(
        const std::vector<double>& rotamer_coords,
        const std::vector<double>& protein_coords,
        const std::vector<double>& rmin_ij,
        const std::vector<double>& eps_ij,
        const std::vector<double>& q_dye,
        const std::vector<double>& q_protein,
        const std::vector<double>& library_weights,
        int n_rotamers, int n_dye_atoms, int n_protein_atoms,
        int potential, double temperature) {
    RotamerScoreResult result;
    result.partition = 0.0;
    result.weights.assign(n_rotamers, 1.0 / n_rotamers);
    result.energies.assign(n_rotamers, 0.0);

    if (n_rotamers == 0 || n_dye_atoms == 0 || n_protein_atoms == 0) {
        return result;
    }

    // The all-pairs inner loop is the existing C++ kernel.
    std::vector<double> rotamer_coords_mut = rotamer_coords;
    std::vector<double> protein_coords_mut = protein_coords;
    std::vector<double> rmin_mut = rmin_ij;
    std::vector<double> eps_mut = eps_ij;

    auto energies = rotamer_interaction_energies(
            rotamer_coords_mut.data(), rotamer_coords_mut.size(),
            protein_coords_mut.data(), protein_coords_mut.size(),
            rmin_mut.data(), rmin_mut.size(),
            eps_mut.data(), eps_mut.size(),
            q_dye, q_protein,
            n_rotamers, n_dye_atoms, n_protein_atoms,
            potential);

    // Two values per rotamer: steric, then electrostatic
    const double GAS_CONSTANT = 1.9858775e-3;
    double rt = GAS_CONSTANT * temperature;

    std::vector<double> boltzmann(n_rotamers);
    double partition = 0.0;
    for (int i = 0; i < n_rotamers; i++) {
        double pot = energies[i * 2];
        double dh = energies[i * 2 + 1];
        result.energies[i] = pot + dh;
        double b = std::exp(-pot / rt - dh);
        if (!library_weights.empty() && i < (int) library_weights.size()) {
            b *= library_weights[i];
        }
        if (std::isnan(b) || std::isinf(b)) b = 0.0;
        boltzmann[i] = b;
        partition += b;
    }

    result.partition = partition;
    if (partition <= 0.0) {
        result.weights.assign(n_rotamers, 1.0 / n_rotamers);
        result.partition = 0.0;
        return result;
    }

    result.weights.resize(n_rotamers);
    for (int i = 0; i < n_rotamers; i++) {
        result.weights[i] = boltzmann[i] / partition;
    }
    return result;
}

IMPBFF_END_NAMESPACE
