/**
 *  \file IMP/bff/internal/PhotonRace.h
 *  \brief One excitation, raced against a rate trace. Header-only.
 *
 * Templated on the trace's element type so a caller holding `float` and a
 * caller holding `double` run the *same* race. That is not a nicety: at a
 * 20 000 000-step trajectory the race makes tens of millions of random reads
 * into the trace, and whether that array is 80 MB or 160 MB is the difference
 * between fitting in cache and not.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#ifndef IMPBFF_INTERNAL_PHOTONRACE_H
#define IMPBFF_INTERNAL_PHOTONRACE_H

#include <IMP/bff/bff_config.h>

#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

IMPBFF_BEGIN_INTERNAL_NAMESPACE

//! A generator seeded from (seed, index), so the stream does not depend on
//! which thread happened to pick up the work.
inline std::mt19937_64 photon_stream(int seed, long long index) {
    if (seed < 0) {
        static std::random_device rd;
        return std::mt19937_64(((static_cast<std::uint64_t>(rd()) << 32) ^ rd()) +
                               static_cast<std::uint64_t>(index) * 0x9E3779B97F4A7C15ULL);
    }
    return std::mt19937_64(static_cast<std::uint64_t>(seed) +
                           static_cast<std::uint64_t>(index) * 0x9E3779B97F4A7C15ULL);
}

//! Race one excitation; returns its delay time, or 0 if quenching won.
template <class Rate>
inline double race_one_photon(const std::vector<Rate>& k_quench, double t_step,
                              double tau0, int seed, long long index,
                              bool& emitted) {
    const long long n_frames = static_cast<long long>(k_quench.size());
    emitted = false;
    if (n_frames == 0) return 0.0;

    std::mt19937_64 rng = photon_stream(seed, index);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    const long long shift = n_frames > 1
            ? static_cast<long long>(uni(rng) * static_cast<double>(n_frames))
            : 0;
    // 1 - u lands in (0, 1]: finite logarithm, never a negative delay.
    const double dt = -std::log(1.0 - uni(rng)) * tau0;
    const long long n_step = static_cast<long long>(dt / t_step);
    for (long long j = shift; j < shift + n_step; ++j) {
        const double k = static_cast<double>(
                k_quench[static_cast<std::size_t>(j % n_frames)]);
        if (uni(rng) < k * t_step) return 0.0;
    }
    emitted = true;
    return dt;
}

IMPBFF_END_INTERNAL_NAMESPACE

#endif //IMPBFF_INTERNAL_PHOTONRACE_H
