/**
 * \file RotamerEnergy.cpp
 * \brief The interaction energy of each rotamer with its protein.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/RotamerEnergy.h>

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

IMPBFF_END_NAMESPACE
