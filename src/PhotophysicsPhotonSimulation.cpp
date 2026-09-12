/**
 * \file PhotophysicsPhotonSimulation.cpp
 * \brief Excited-state kinetics along a trajectory: photons, and the curve.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/PhotophysicsPhotonSimulation.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/PhotonRace.h>

#include <cmath>
#include <random>

#ifdef _OPENMP
#include <omp.h>
#endif

IMPBFF_BEGIN_NAMESPACE

// The per-photon race lives in internal/PhotonRace.h, shared with the fused
// walk->rate->photons kernel so both run the same draws.

namespace {
std::vector<double> photon_trace_impl(
        int n_ph, const std::vector<double>& k_quench,
        double t_step, double tau0, int seed) {
    const int n = n_ph > 0 ? n_ph : 0;
    // Two per photon: delay, emitted. See the header on why not an out-param.
    std::vector<double> out(static_cast<std::size_t>(n) * 2, 0.0);
    if (k_quench.empty() || n == 0) return out;

#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        bool got = false;
        out[2 * i + 0] =
                internal::race_one_photon(k_quench, t_step, tau0, seed, i, got);
        out[2 * i + 1] = got ? 1.0 : 0.0;
    }
    return out;
}
}  // namespace

void photon_trace(int n_ph, const std::vector<double>& k_quench,
        double t_step, double tau0, int seed, double** out_view, int* n_out_view) {
    internal::copy_to_view(photon_trace_impl(n_ph, k_quench, t_step, tau0, seed),
                           out_view, n_out_view);
}

void simulate_photon_trace(int n_ph, const std::vector<double>& k_quench,
                           double t_step, double tau0, int random_seed,
                           double** out_delays, int* n_delays,
                           unsigned char** out_emitted, int* n_emitted) {
    const std::vector<double> flat =
            photon_trace_impl(n_ph, k_quench, t_step, tau0, random_seed);
    const int n = n_ph > 0 ? n_ph : 0;
    if (out_delays != NULL && n_delays != NULL) {
        double* delays = internal::new_double_view(n, out_delays, n_delays);
        if (delays != NULL) {
            for (int i = 0; i < n; ++i) delays[i] = flat[2 * i + 0];
        }
    }
    if (out_emitted != NULL && n_emitted != NULL) {
        unsigned char* emitted =
                internal::new_uchar_view(n, out_emitted, n_emitted);
        if (emitted != NULL) {
            for (int i = 0; i < n; ++i)
                emitted[i] = flat[2 * i + 1] != 0.0 ? 1 : 0;
        }
    }
}

void simulate_quenched_decay(int n_curves, double* decay, int n_decay,
                             double dt_tac, const std::vector<double>& k_quench,
                             double t_step, double tau0, int random_seed) {
    const std::vector<double> curve = quenched_decay(
            n_curves, n_decay, dt_tac, k_quench, t_step, tau0, random_seed);
    for (int i = 0; i < n_decay && i < static_cast<int>(curve.size()); ++i)
        decay[i] += curve[i];
}

std::vector<double> quenched_decay(
        int n_curves, int n_bins, double dt_tac,
        const std::vector<double>& k_quench,
        double t_step, double tau0, int seed) {
    std::vector<double> decay(n_bins > 0 ? n_bins : 0, 0.0);
    const long long n_frames = static_cast<long long>(k_quench.size());
    if (n_frames == 0 || n_bins <= 0 || n_curves <= 0) return decay;

    const double intrinsic_rate = tau0 > 0.0 ? 1.0 / tau0 : 0.0;
    const double max_time = dt_tac * n_bins;

    int n_blocks = 1;
#ifdef _OPENMP
    n_blocks = omp_get_max_threads();
#endif
    if (n_blocks > n_curves) n_blocks = n_curves;
    if (n_blocks < 1) n_blocks = 1;

    // One private histogram per block, reduced in block order below: the bin
    // index is data-dependent, so a shared accumulator loses updates.
    std::vector<double> partial(static_cast<std::size_t>(n_blocks) * n_bins, 0.0);

#pragma omp parallel for schedule(static)
    for (int block = 0; block < n_blocks; ++block) {
        const int lo = static_cast<int>(static_cast<long long>(block) * n_curves / n_blocks);
        const int hi = static_cast<int>(static_cast<long long>(block + 1) * n_curves / n_blocks);
        double* row = &partial[static_cast<std::size_t>(block) * n_bins];
        for (int i = lo; i < hi; ++i) {
            std::mt19937_64 rng = internal::photon_stream(seed, i);
            std::uniform_real_distribution<double> uni(0.0, 1.0);
            const long long span = n_frames > 1 ? n_frames / 2 : 1;
            const long long shift = span > 0
                    ? static_cast<long long>(uni(rng) * static_cast<double>(span))
                    : 0;
            double intensity = 1.0, t = 0.0;
            long long frame = 0;
            while (t < max_time && intensity > 1e-6 && frame < n_frames) {
                const std::size_t idx =
                        static_cast<std::size_t>((shift + frame) % n_frames);
                const double total_rate = intrinsic_rate + k_quench[idx];
                double dt = t_step;
                if (t + dt > max_time) dt = max_time - t;
                if (dt <= 0.0) break;
                const int bin = static_cast<int>(t / dt_tac);
                if (bin >= 0 && bin < n_bins) row[bin] += intensity * total_rate * dt;
                intensity *= std::exp(-total_rate * dt);
                t += dt;
                ++frame;
            }
        }
    }
    for (int block = 0; block < n_blocks; ++block) {
        const double* row = &partial[static_cast<std::size_t>(block) * n_bins];
        for (int bin = 0; bin < n_bins; ++bin) decay[bin] += row[bin];
    }
    return decay;
}

IMPBFF_END_NAMESPACE
