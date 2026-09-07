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

    Parallel **and** reproducible: each thread gets its own generator seeded
    from \p seed and the photon index, so the result does not depend on how the
    scheduler distributed the work. Seeding per thread instead would force a
    choice between the two.

    \param[in] n_ph excitation events to simulate
    \param[in] k_quench per-frame quenching rate, 1/ns
    \param[in] t_step trajectory time step, ns
    \param[in] tau0 intrinsic lifetime, ns
    \param[in] seed reproducible when non-negative
    \return **two** values per event: the delay time in ns (0 for a quenched
            event), then 1.0 if a photon got out and 0.0 if quenching won.

    Interleaved rather than returned alongside an out-parameter. SWIG turns a
    returned `std::vector` into a Python tuple, which costs ~35-40 ns per element to build and walk back;
    an out-parameter stays a wrapper object that numpy walks one `__getitem__`
    at a time, about 340 ns per element against 34 ns for an *input* array of
    the same size. On a 40 000-photon trace that flag cost 13.7 ms of 152 ms.

    The input `k_quench` is the remaining marshalling cost here -- 36 ns per
    element, so ~180 ms for a five-million-frame trace. #quenched_donor_photons
    avoids it entirely by never letting the trace out of C++.
*/
IMPBFFEXPORT void photon_trace(
        int n_ph, const std::vector<double>& k_quench,
        double t_step, double tau0, int seed,
        double** out_view, int* n_out_view);

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

    The caller's \p decay histogram is accumulated **in place**, so several
    virtual-trajectory groups can fill one experiment.

    \param[in] n_curves virtual trajectories to average
    \param[in,out] decay,n_decay the histogram, \p n_decay entries, float64
    \param[in] dt_tac TAC bin width, ns
    \param[in] k_quench per-frame quenching rate, 1/ns
    \param[in] t_step trajectory time step, ns
    \param[in] tau0 intrinsic lifetime, ns
    \param[in] random_seed reproducible when non-negative; drawn freely otherwise
*/
IMPBFFEXPORT void simulate_quenched_decay(
        int n_curves, double* decay, int n_decay, double dt_tac,
        const std::vector<double>& k_quench, double t_step, double tau0,
        int random_seed = -1);

//! A Monte-Carlo photon trace, split into delay and emitted-flag views.
/*! The sampling is #photon_trace's; this is the shape its callers asked for: a
    delay time per event and a separate emitted-flag. \p out_delays holds
    \p n_ph doubles (0.0 for a quenched event) and \p out_emitted holds \p n_ph
    bytes, 1 or 0. `random_seed < 0` draws freely. */
IMPBFFEXPORT void simulate_photon_trace(
        int n_ph, const std::vector<double>& k_quench, double t_step = 0.01,
        double tau0 = 0.25, int random_seed = -1, double** out_delays = NULL,
        int* n_delays = NULL, unsigned char** out_emitted = NULL,
        int* n_emitted = NULL);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_PHOTONSIMULATION_H
