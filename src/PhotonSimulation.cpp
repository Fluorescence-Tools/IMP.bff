/**
 * \file PhotonSimulation.cpp
 * \brief Excited-state kinetics along a trajectory: photons, and the curve.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/PhotonSimulation.h>

#include <cmath>
#include <random>

#ifdef _OPENMP
#include <omp.h>
#endif

IMPBFF_BEGIN_NAMESPACE

namespace {
//! A generator seeded from (seed, index) so the stream does not depend on
//! which thread happened to pick up the work.
inline std::mt19937_64 stream_for(int seed, long long index) {
    if (seed < 0) {
        static std::random_device rd;
        return std::mt19937_64(((static_cast<std::uint64_t>(rd()) << 32) ^ rd()) +
                               static_cast<std::uint64_t>(index) * 0x9E3779B97F4A7C15ULL);
    }
    return std::mt19937_64(static_cast<std::uint64_t>(seed) +
                           static_cast<std::uint64_t>(index) * 0x9E3779B97F4A7C15ULL);
}
}  // namespace

std::vector<double> photon_trace(
        int n_ph, const std::vector<double>& k_quench,
        double t_step, double tau0, int seed, std::vector<int>& emitted) {
    const int n = n_ph > 0 ? n_ph : 0;
    std::vector<double> dts(n, 0.0);
    emitted.assign(n, 0);
    const long long n_frames = static_cast<long long>(k_quench.size());
    if (n_frames == 0 || n == 0) return dts;

#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        std::mt19937_64 rng = stream_for(seed, i);
        std::uniform_real_distribution<double> uni(0.0, 1.0);
        const long long shift = n_frames > 1
                ? static_cast<long long>(uni(rng) * static_cast<double>(n_frames))
                : 0;
        // 1 - u lands in (0, 1]: finite logarithm, never a negative delay.
        const double dt = -std::log(1.0 - uni(rng)) * tau0;
        const long long n_step = static_cast<long long>(dt / t_step);
        bool quenched = false;
        for (long long j = shift; j < shift + n_step; ++j) {
            if (uni(rng) < k_quench[static_cast<std::size_t>(j % n_frames)] * t_step) {
                quenched = true;
                break;
            }
        }
        if (!quenched) {
            dts[i] = dt;
            emitted[i] = 1;
        }
    }
    return dts;
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
            std::mt19937_64 rng = stream_for(seed, i);
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
