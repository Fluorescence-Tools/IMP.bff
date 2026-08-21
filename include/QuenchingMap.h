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

//! Voxel-centre offsets from the grid anchor, in Angstrom.
/*!
    \f$(i - (ng-1)/2)\,dg\f$ with the **integer** centre offset, the same one
    grid_center_index() uses, so a map built here indexes the way the Brownian
    walk samples it.
*/
IMPBFFEXPORT std::vector<double> grid_axis(int ng, double dg);

//! A diffusion map from a radial profile, evaluated on an integer-radius axis.
/*!
    \f$f\f$ is a profile of the distance from the grid anchor, sampled once per
    integer Angstrom: \p radial[k] is the value at \f$r = k\f$. Every voxel gets
    \f$\p radial[\operatorname{round}(|r|)]\f$, clamped to the table -- the old
    Python took a *callable* and there is no C++ spelling of one that does not
    call back into the interpreter per voxel, so the caller lands \f$f\f$ on a
    radius grid first. This is what a turnover lamp (``stretched linker'')
    mobility looks like before #slow_near_atoms adds the local crowding.

    \param[in] density binary occupancy, flat ng^3, whose grid side gives ng
    \param[in] dg voxel edge, A
    \param[in] radial \f$f\f$ at integer radii in Angstrom; length 0 is the
               identity (returns 1.0 everywhere)
    \param[out] out_view,n_out_view flat ng^3
*/
IMPBFFEXPORT void radial_diffusion_map(
        const std::vector<double>& density, double dg,
        const std::vector<double>& radial, double** out_view, int* n_out_view);

//! The mobility field: a base coefficient, slowed by nearby atoms.
/*!
    \param[in] density binary occupancy of the accessible volume, flat ng^3
    \param[in] r0 the grid anchor
    \param[in] dg voxel edge, A
    \param[in] atoms_xyz obstacle coordinates, flat
    \param[in] free_diffusion the unhindered dye diffusion coefficient, A^2/ns,
               used wherever \p base is not given
    \param[in] min_distance contact distance, roughly dye radius + 2 vdW
    \param[in] slow_factor factor applied once per contacting atom, in [0, 1]
    \param[in] base an optional per-voxel base coefficient replacing the
               constant \p free_diffusion -- a radial profile, say. Length 0
               takes the constant.
*/
IMPBFFEXPORT void diffusion_coefficient_map(
        const std::vector<double>& density,
        const std::vector<double>& r0,
        double dg,
        const std::vector<double>& atoms_xyz,
        double free_diffusion,
        double min_distance,
        double slow_factor,
        const std::vector<double>& base,
        double** out_view, int* n_out_view);

//! Total decay rate per voxel: intrinsic plus exponential PET.
/*!
    \f$k(r) = 1/\tau_0 + \sum_a k_{Q,a}\,e^{-(|r - r_a| - r_{dye})/r_{C,a}}\f$

    The distance is measured from the dye **surface**. Unlike the contact-sphere
    law in quenching_rate_grid() there is no cut-off: a distant atom contributes
    exponentially little rather than nothing. Voxels outside the accessible
    volume stay at zero -- **not** at \f$1/\tau_0\f$; the solver reads this as
    the rate field on a domain masked by the same occupancy, so an unreachable
    voxel carries no rate at all.

    \param[in] tau0 unquenched lifetime, ns; non-positive drops the floor
*/
IMPBFFEXPORT void quenching_rate_map(
        const std::vector<double>& density,
        const std::vector<double>& r0,
        double dg,
        const std::vector<double>& atoms_xyz,
        const std::vector<double>& kQ,
        const std::vector<double>& rC,
        double tau0,
        double dye_radius,
        double** out_view, int* n_out_view);

//! An effective FRET rate for every donor voxel, from the acceptor cloud.
/*!
    A fixed donor position does not have *a* FRET rate: it has a distribution of
    them, one per accessible acceptor position. This approximates that sum by a
    single exponential whose rate is the reciprocal of the **mean transfer
    time** -- the *harmonic* mean, dominated by the slow, distant acceptor
    positions. fret_rate_trace() takes the arithmetic mean instead, because it
    models a fast-exchanging acceptor where the *rates* average.

    \param[in] kf the donor's radiative rate, \f$1/\tau_0\f$
    \param[in] acceptor_step stride over the acceptor grid; the cost is the
               product of the two grids' sizes
    \throw ValueException for a non-positive \p kf
*/
IMPBFFEXPORT void fret_rate_map(
        const std::vector<double>& density_donor,
        const std::vector<double>& density_acceptor,
        const std::vector<double>& r0_donor,
        const std::vector<double>& r0_acceptor,
        double dg_donor,
        double dg_acceptor,
        double forster_radius,
        double kf,
        int acceptor_step,
        double** out_view, int* n_out_view);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_QUENCHINGMAP_H
