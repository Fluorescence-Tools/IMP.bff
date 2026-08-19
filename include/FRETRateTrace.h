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

IMPBFF_END_NAMESPACE

#endif //IMPBFF_FRETRATETRACE_H
