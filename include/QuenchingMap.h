/**
 *  \file IMP/bff/QuenchingMap.h
 *  \brief Mobility, quenching-rate and FRET-rate fields on an AV grid.
 *
 * Ported from Python by PRD-113: numba is a prototyping tool in this package,
 * not a runtime dependency, so every numerical kernel is C++.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_QUENCHINGMAP_H
#define IMPBFF_QUENCHINGMAP_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Slow the mobility field wherever the dye is in contact with an atom.
/*!
    The factor is applied **once per contacting atom**, so it compounds: a voxel
    touching 45 atoms at a factor of 0.985 ends at 0.5, and at 0.9 it ends at
    9e-3. Only values within a whisker of 1.0 are meaningful, which is a
    property of the model rather than a taste.

    \param[in] d_map the unmodified mobility field, flat ng^3
    \param[in] density accessible-volume density; zero voxels stay zero
    \param[in] axis voxel offsets from the grid centre, length ng
    \param[in] r0 grid origin
    \param[in] atoms_xyz obstacle coordinates, flat, three per atom
    \param[in] min_distance_sq squared contact distance
    \param[in] factor per-contact slowing factor
*/
IMPBFFEXPORT void slow_near_atoms(
        const std::vector<double>& d_map,
        const std::vector<double>& density,
        const std::vector<double>& axis,
        const std::vector<double>& r0,
        const std::vector<double>& atoms_xyz,
        double min_distance_sq,
        double factor,
        double** out_view, int* n_out_view
);

//! Quenching-rate field: \f$k(r) = 1/\tau_0 + \sum_a k_a e^{-(|r-r_a| - r_{dye})/r_{C,a}}\f$.
/*!
    Voxels outside the accessible volume stay at **zero**, not at \f$1/\tau_0\f$:
    the dye cannot be there, so it has no decay rate there.

    \param[in] density accessible-volume density
    \param[in] axis voxel offsets, length ng
    \param[in] r0 grid origin
    \param[in] atoms_xyz obstacle coordinates, flat
    \param[in] kQ per-atom rate constant; zero means the atom does not quench
    \param[in] rC per-atom attenuation length
    \param[in] dye_radius subtracted from the centre-to-centre distance
    \param[in] inv_tau0 the radiative floor
*/
IMPBFFEXPORT void quenching_map(
        const std::vector<double>& density,
        const std::vector<double>& axis,
        const std::vector<double>& r0,
        const std::vector<double>& atoms_xyz,
        const std::vector<double>& kQ,
        const std::vector<double>& rC,
        double dye_radius,
        double inv_tau0,
        double** out_view, int* n_out_view
);

//! FRET-rate field of a donor volume against a whole acceptor volume.
/*!
    **This is the harmonic mean** -- the acceptor-weighted mean transfer *time*
    is accumulated and inverted, not the mean rate. That is the static limit, in
    which the acceptor does not move within the donor's excited-state lifetime.
    The arithmetic mean of rates is the fast-exchange limit and is what
    ``fret_rate_trace`` computes along a trajectory. The two are different
    physics and differ measurably; neither is a variant of the other.

    \param[in] density_d,density_a the two densities, flat
    \param[in] axis_d,axis_a voxel offsets of each grid
    \param[in] r0_d,r0_a the two grid origins
    \param[in] r0_6 the Forster radius to the sixth power
    \param[in] kf radiative rate
    \param[in] step stride over the acceptor grid; 1 visits every voxel
*/
IMPBFFEXPORT void fret_map(
        const std::vector<double>& density_d,
        const std::vector<double>& density_a,
        const std::vector<double>& axis_d,
        const std::vector<double>& axis_a,
        const std::vector<double>& r0_d,
        const std::vector<double>& r0_a,
        double r0_6,
        double kf,
        int step,
        double** out_view, int* n_out_view
);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_QUENCHINGMAP_H
