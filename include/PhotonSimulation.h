/**
 *  \file IMP/bff/PhotonSimulation.h
 *  \brief Excited-state kinetics along a trajectory: photons, and the curve.
 *
 * Two ways to ask the same question of a quenching-rate trace. #photon_trace
 * samples one excitation at a time and answers *did this photon get out, and
 * when* -- shot noise included, which is what a burst measurement has.
 * #quenched_decay propagates the excited-state population instead and returns
 * the noise-free curve directly, far cheaper when only the curve is wanted.
 *
 * Neither convolves with an instrument response: that is experiment-side and
 * lives outside this package.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_PHOTONSIMULATION_H
#define IMPBFF_PHOTONSIMULATION_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Monte-Carlo photon trace against a per-frame quenching rate.
/*!
    Each excitation draws an intrinsic delay time \f$-\tau_0 \ln(1-u)\f$ and is
    then raced against quenching frame by frame from a random starting point on
    the trajectory. The draw is written as \f$1-u\f$ deliberately: a generator
    uniform on [0, 1) makes \f$1-u\f$ land in (0, 1], so the logarithm is finite
    and never positive. The ancestor of this code used
    \f$\ln(1/(u+\epsilon))\f$, which goes *negative* for \f$u > 1-\epsilon\f$ --
    a photon emitted before it was excited, 2.4e-4 of all photons, each then
    silently dropped by any histogram starting at zero.

    Parallel **and** reproducible: each thread gets its own generator seeded from
    \p seed and the photon index, so the result does not depend on how the
    scheduler distributed the work. The numba this replaces had to choose
    between the two.

    \param[in] n_ph excitation events to simulate
    \param[in] k_quench per-frame quenching rate, 1/ns
    \param[in] t_step trajectory time step, ns
    \param[in] tau0 intrinsic lifetime, ns
    \param[in] seed reproducible when non-negative
    \param[out] emitted per event: 1 if a photon got out, 0 if quenched first
    \return delay time of each event, ns; 0 for a quenched event
*/
IMPBFFEXPORT std::vector<double> photon_trace(
        int n_ph, const std::vector<double>& k_quench,
        double t_step, double tau0, int seed,
        std::vector<int>& emitted);

//! Trajectory-driven fluorescence decay, without shot noise.
/*!
    For each virtual trajectory a random starting frame is chosen and the
    excited-state population is propagated under
    \f$1/\tau_0 + k_q(\text{frame})\f$, the emitted intensity accumulating into
    TAC bins.

    Accumulation is blocked, one private histogram per block reduced in block
    order. That is a correctness requirement, not a tuning knob: the bin index
    is data-dependent, so accumulating into a shared histogram from parallel
    iterations loses updates. Measured on the numba ancestor with 8 threads and
    200 curves, sums scattered over 196.73-200.61 against a stable single-thread
    200.612864255 -- up to 2 % of the emitted intensity dropped, differently
    each run.

    \param[in] n_curves virtual trajectories to average
    \param[in] n_bins TAC bins
    \param[in] dt_tac TAC bin width, ns
    \param[in] k_quench per-frame quenching rate, 1/ns
    \param[in] t_step trajectory time step, ns
    \param[in] tau0 intrinsic lifetime, ns
    \param[in] seed reproducible when non-negative
    \return the decay, \p n_bins entries, to be added to the caller's histogram
*/
IMPBFFEXPORT std::vector<double> quenched_decay(
        int n_curves, int n_bins, double dt_tac,
        const std::vector<double>& k_quench,
        double t_step, double tau0, int seed);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_PHOTONSIMULATION_H
