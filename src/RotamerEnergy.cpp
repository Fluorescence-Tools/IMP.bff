/**
 * \file RotamerEnergy.cpp
 * \brief The interaction energy of each rotamer with its protein.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/RotamerEnergy.h>

#include <algorithm>
#include <cmath>

IMPBFF_BEGIN_NAMESPACE

std::vector<double> rotamer_interaction_energies(
        const std::vector<double>& rotamer_coords,
        const std::vector<double>& protein_coords,
        const std::vector<double>& rmin_ij,
        const std::vector<double>& eps_ij,
        const std::vector<double>& q_dye,
        const std::vector<double>& q_protein,
        int n_rotamers, int n_dye_atoms, int n_protein_atoms,
        int potential,
        double steric_cutoff, double coulomb_cutoff,
        double debye_length, double coulomb_prefactor) {
    std::vector<double> out(static_cast<std::size_t>(std::max(0, n_rotamers)) * 2, 0.0);
    if (n_rotamers <= 0 || n_dye_atoms <= 0 || n_protein_atoms <= 0) return out;

    const bool electrostatic =
            static_cast<int>(q_dye.size()) == n_dye_atoms &&
            static_cast<int>(q_protein.size()) == n_protein_atoms;
    const double steric_cut2 = steric_cutoff * steric_cutoff;
    const double coulomb_cut2 = coulomb_cutoff * coulomb_cutoff;

#pragma omp parallel for schedule(static)
    for (int r = 0; r < n_rotamers; ++r) {
        double steric = 0.0, coulomb = 0.0;
        const double* conf =
                &rotamer_coords[static_cast<std::size_t>(r) * n_dye_atoms * 3];
        for (int a = 0; a < n_dye_atoms; ++a) {
            const double ax = conf[3 * a + 0];
            const double ay = conf[3 * a + 1];
            const double az = conf[3 * a + 2];
            const double qa = electrostatic ? q_dye[a] : 0.0;
            const double* rmin_row =
                    &rmin_ij[static_cast<std::size_t>(a) * n_protein_atoms];
            const double* eps_row =
                    &eps_ij[static_cast<std::size_t>(a) * n_protein_atoms];
            for (int b = 0; b < n_protein_atoms; ++b) {
                const double dx = ax - protein_coords[3 * b + 0];
                const double dy = ay - protein_coords[3 * b + 1];
                const double dz = az - protein_coords[3 * b + 2];
                const double d2 = dx * dx + dy * dy + dz * dz;

                // Squared distances until a cutoff passes: the square root is
                // the expensive part and most pairs never need it.
                if (d2 < steric_cut2) {
                    const double d = std::sqrt(d2);
                    const double ratio = rmin_row[b] / d;
                    if (potential == ROTAMER_POTENTIAL_GAUSS) {
                        steric += eps_row[b] * std::exp(-0.5 / (ratio * ratio));
                    } else {
                        const double r6 = ratio * ratio * ratio;
                        const double ratio6 = r6 * r6;
                        steric += eps_row[b] * (ratio6 * ratio6 - 2.0 * ratio6);
                    }
                }
                if (electrostatic && qa != 0.0 && q_protein[b] != 0.0 &&
                    d2 < coulomb_cut2) {
                    const double d = std::sqrt(d2);
                    coulomb += qa * q_protein[b] * coulomb_prefactor / d *
                               std::exp(-d / debye_length);
                }
            }
        }
        out[2 * r + 0] = steric;
        out[2 * r + 1] = coulomb;
    }
    return out;
}

namespace {
//! Padded axis-aligned bounding box of one conformer: xmin,ymin,zmin,xmax,...
inline void conformer_box(const double* xyz, int n_atoms, double pad, double* box) {
    for (int k = 0; k < 3; ++k) { box[k] = xyz[k]; box[3 + k] = xyz[k]; }
    for (int a = 1; a < n_atoms; ++a) {
        for (int k = 0; k < 3; ++k) {
            const double v = xyz[3 * a + k];
            if (v < box[k]) box[k] = v;
            if (v > box[3 + k]) box[3 + k] = v;
        }
    }
    for (int k = 0; k < 3; ++k) { box[k] -= pad; box[3 + k] += pad; }
}

inline bool boxes_overlap(const double* a, const double* b) {
    for (int k = 0; k < 3; ++k) {
        if (a[3 + k] < b[k] || b[3 + k] < a[k]) return false;
    }
    return true;
}
}  // namespace

std::vector<double> rotamer_pair_energy_matrix(
        const std::vector<double>& coords_a, const std::vector<double>& coords_b,
        const std::vector<double>& rmin, const std::vector<double>& eps,
        int n_a_conf, int n_a_atoms, int n_b_conf, int n_b_atoms,
        double r_cutoff, double aabb_pad, double r_floor) {
    std::vector<double> out(
            static_cast<std::size_t>(std::max(0, n_a_conf)) * std::max(0, n_b_conf), 0.0);
    if (n_a_conf <= 0 || n_b_conf <= 0 || n_a_atoms <= 0 || n_b_atoms <= 0) return out;
    if (rmin.size() != static_cast<std::size_t>(n_a_atoms) * n_b_atoms ||
        eps.size() != rmin.size()) {
        return out;   // no parameters means no interaction, as in the Python
    }

    std::vector<double> box_a(static_cast<std::size_t>(n_a_conf) * 6);
    std::vector<double> box_b(static_cast<std::size_t>(n_b_conf) * 6);
    for (int i = 0; i < n_a_conf; ++i) {
        conformer_box(&coords_a[static_cast<std::size_t>(i) * n_a_atoms * 3],
                      n_a_atoms, aabb_pad, &box_a[static_cast<std::size_t>(i) * 6]);
    }
    for (int j = 0; j < n_b_conf; ++j) {
        conformer_box(&coords_b[static_cast<std::size_t>(j) * n_b_atoms * 3],
                      n_b_atoms, aabb_pad, &box_b[static_cast<std::size_t>(j) * 6]);
    }

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n_a_conf; ++i) {
        const double* ca = &coords_a[static_cast<std::size_t>(i) * n_a_atoms * 3];
        for (int j = 0; j < n_b_conf; ++j) {
            if (!boxes_overlap(&box_a[static_cast<std::size_t>(i) * 6],
                               &box_b[static_cast<std::size_t>(j) * 6])) continue;
            const double* cb = &coords_b[static_cast<std::size_t>(j) * n_b_atoms * 3];
            double e = 0.0;
            for (int a = 0; a < n_a_atoms; ++a) {
                const double ax = ca[3 * a + 0], ay = ca[3 * a + 1], az = ca[3 * a + 2];
                const std::size_t row = static_cast<std::size_t>(a) * n_b_atoms;
                for (int b = 0; b < n_b_atoms; ++b) {
                    const double dx = ax - cb[3 * b + 0];
                    const double dy = ay - cb[3 * b + 1];
                    const double dz = az - cb[3 * b + 2];
                    const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
                    const double rm = rmin[row + b];
                    // Repulsive-only, and inside the cutoff. Both tests are on
                    // the unclamped distance, as the Python has them.
                    if (!(d < rm) || !(d < r_cutoff)) continue;
                    const double safe = d > r_floor ? d : r_floor;
                    const double ratio = rm / safe;
                    const double r6 = ratio * ratio * ratio;
                    const double ratio6 = r6 * r6;
                    e += eps[row + b] * (ratio6 * ratio6 - 2.0 * ratio6);
                }
            }
            out[static_cast<std::size_t>(i) * n_b_conf + j] = e;
        }
    }
    return out;
}

std::vector<double> lj_pair_energies(
        const std::vector<double>& coords,
        const std::vector<int>& index_a, const std::vector<int>& index_b,
        const std::vector<double>& rmin, const std::vector<double>& eps,
        int n_frames, int n_atoms, int n_pairs,
        bool repulsive_only, double r_floor) {
    std::vector<double> out(static_cast<std::size_t>(std::max(0, n_frames)), 0.0);
    if (n_frames <= 0 || n_pairs <= 0 || n_atoms <= 0) return out;

#pragma omp parallel for schedule(static)
    for (int f = 0; f < n_frames; ++f) {
        const double* xyz = &coords[static_cast<std::size_t>(f) * n_atoms * 3];
        double e = 0.0;
        for (int p = 0; p < n_pairs; ++p) {
            const int ia = index_a[p], ib = index_b[p];
            const double dx = xyz[3 * ia + 0] - xyz[3 * ib + 0];
            const double dy = xyz[3 * ia + 1] - xyz[3 * ib + 1];
            const double dz = xyz[3 * ia + 2] - xyz[3 * ib + 2];
            const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
            // Both tests on the unclamped distance, as the Python has them;
            // only the ratio uses the floor.
            if (repulsive_only && !(d < rmin[p])) continue;
            const double safe = d > r_floor ? d : r_floor;
            const double ratio = rmin[p] / safe;
            const double r6 = ratio * ratio * ratio;
            const double ratio6 = r6 * r6;
            e += eps[p] * (ratio6 * ratio6 - 2.0 * ratio6);
        }
        out[f] = e;
    }
    return out;
}

IMPBFF_END_NAMESPACE
