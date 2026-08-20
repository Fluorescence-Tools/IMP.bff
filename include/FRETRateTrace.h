/**
 *  \file IMP/bff/FRETRateTrace.h
 *  \brief FRET rate along a dye trajectory.
 *
 * Ported from Python by PRD-113: numba is a prototyping tool in this package,
 * not a runtime dependency, so every numerical kernel is C++.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_FRETRATETRACE_H
#define IMPBFF_FRETRATETRACE_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/OrientationFactor.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! FRET rate per frame, donor trajectory against a static acceptor cloud.
/*!
    **The arithmetic mean of rates**, averaged over the acceptor cloud -- the
    fast-exchange limit, in which the acceptor re-randomises within the donor's
    excited-state lifetime. ``fret_map`` accumulates the mean transfer *time*
    instead and inverts it, which is the static limit. Different physics, not
    two spellings of one thing.

    \param[in] trajectory donor positions, flat, three per frame
    \param[in] acceptor_points acceptor cloud, flat, three per point
    \param[in] R0 Forster radius
    \param[in] tau0 unquenched donor lifetime
    \param[in] r_min2 squared closest approach; pairs nearer than this are
        clamped, because \f$1/r^6\f$ diverges and two dyes cannot interpenetrate
    \param[in] kappa2_scale multiplier on \f$R_0^6\f$ from the orientation factor
    \param[out] out_view,n_out_view the rate at each frame, as a numpy view over
                the kernel's buffer -- a trajectory runs to millions of frames
                and a returned `std::vector` costs ~35-40 ns each to hand back
                (internal/OutputView.h carries the measured table).
                See internal/OutputView.h.
*/
IMPBFFEXPORT void fret_rate_trace_kernel(
        const std::vector<double>& trajectory,
        const std::vector<double>& acceptor_points,
        double R0,
        double tau0,
        double r_min2,
        double kappa2_scale,
        double** out_view, int* n_out_view
);

//! FRET rate per frame from two trajectories, paired frame by frame.
/*!
    Both dyes resolved in time. The two walks are independent, so pairing frame
    *i* with frame *i* samples the joint distribution -- provided both were
    produced on the same time base.

    \param[in] donor,acceptor positions, flat, three per frame, same length
    \param[in] R0,tau0,r_min2,kappa2_scale as above
    \param[out] out_view,n_out_view the rate at each frame, as above
*/
IMPBFFEXPORT void fret_rate_pair_trace_kernel(
        const std::vector<double>& donor,
        const std::vector<double>& acceptor,
        double R0,
        double tau0,
        double r_min2,
        double kappa2_scale,
        double** out_view, int* n_out_view
);

//! The largest physically meaningful orientation factor: collinear dipoles.
extern IMPBFFEXPORT const double MAX_KAPPA2;

//! Turn an orientation factor into a multiplier on \f$R_0^6\f$.
/*!
    A published Forster radius is quoted *at* the isotropic 2/3, so the transfer
    rate carries the ratio \f$\kappa^2/(2/3)\f$ -- the \f$1.5\kappa^2\f$ of
    the literature -- and passing 2/3 leaves the rate as it was.

    **The diffusing-sphere model does not compute \f$\kappa^2\f$, and cannot**:
    its dye is a structureless point in a volume, with no transition dipole to
    orient. This makes the *assumption* explicit and adjustable instead of
    silently baked into \f$R_0\f$. kappa2_from_dipoles() does compute it, but it
    takes the two dipole vectors, which is exactly what a point has not got.

    There is **no sentinel for "unspecified"**, and deliberately: the isotropic
    default is the value #IMP::bff::kappa2_isotropic() names, so a caller with
    nothing better to say passes that. NaN is an *error* here rather than a
    stand-in for it -- a NaN reaching this from a failed orientation calculation
    would otherwise be silently answered with the isotropic average.

    \throw ValueException for NaN, or outside [0, MAX_KAPPA2]
*/
IMPBFFEXPORT double kappa2_scale(double kappa2);

// The two trace functions below take \p kappa2 without a default: it precedes
// the output view, and C++ will not have a defaulted parameter before an
// undefaulted one. Pass kappa2_isotropic() for "no orientation information".

//! Per-frame FRET rate (1/ns) along a donor trajectory, acceptor averaged.
/*!
    \f$k(t) = (1/\tau_0)\,\langle (R_0/|r_D(t) - r_A|)^6\rangle_A\,
    (\kappa^2/(2/3))\f$

    This resolves the **donor** in time -- the point of simulating its diffusion
    at all -- while treating the acceptor as sampling its own volume quickly
    compared with the donor's excited-state lifetime. That is the *fast-acceptor
    limit*, and it is not valid for an immobilised acceptor; use
    fret_rate_pair_trace() then.

    \param[in] trajectory donor positions, flat, three per frame
    \param[in] acceptor_points acceptor cloud, flat, three per point
    \param[in] r_min closest approach of the two dye centres, A. Both clouds are
               centre positions and can overlap, where \f$(R_0/r)^6\f$ would
               diverge; two dye spheres cannot interpenetrate.
    \param[in] max_acceptor_points cap on the points averaged over, taken by an
               even stride. The average converges long before the full cloud,
               and the cost is `n_frames * n_acceptor`.
    \throw ValueException for an empty cloud, a non-positive \p tau0, or a
           \p kappa2 kappa2_scale() refuses
*/
IMPBFFEXPORT void fret_rate_trace(
        const std::vector<double>& trajectory,
        const std::vector<double>& acceptor_points,
        double R0,
        double tau0,
        double r_min,
        int max_acceptor_points,
        double kappa2,
        double** out_view, int* n_out_view);

//! Per-frame FRET rate (1/ns) from *two* trajectories, both resolved in time.
/*!
    The general case. fret_rate_trace() averages the acceptor over its volume
    instead, which is the fast-acceptor limit. The two do not agree, and the
    reason is **occupancy**, not the averaging order: a cloud weights every
    accessible voxel the same, while a trajectory is weighted by where the dye
    actually spends its time, and a sticky dye dwells near the surface. Which
    way that moves the efficiency depends on the geometry.

    **Unequal lengths are refused, not truncated.** A trajectory is a
    *concatenation of one walk per excitation*, so an arbitrary cut splits a
    walk and pairs one excitation's tail against another's head. Unequal lengths
    are the normal case -- the two dyes sit at different sites and are quenched
    differently -- and the lever is the number of *walks*, not `t_max` or
    `t_step`.

    \throw ValueException for an empty trajectory, a non-positive \p tau0, or
           a frame-count mismatch
*/
IMPBFFEXPORT void fret_rate_pair_trace(
        const std::vector<double>& donor_trajectory,
        const std::vector<double>& acceptor_trajectory,
        double R0,
        double tau0,
        double r_min,
        double kappa2,
        double** out_view, int* n_out_view);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_FRETRATETRACE_H
