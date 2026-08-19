/**
 * \file QuenchedDecay.cpp
 * \brief The fused walk -> rate -> photons path.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/QuenchedDecay.h>
#include <IMP/bff/internal/OutputView.h>

#include <IMP/bff/internal/PhotonRace.h>
#include <IMP/bff/internal/RandomWalk.h>

#include <algorithm>
#include <cmath>

IMPBFF_BEGIN_NAMESPACE

namespace {
std::vector<double> quenched_donor_photons_impl(
        int* occupancy, int n_occupancy,
        double* mobility, int n_mobility,
        double* rate_map, int n_rate_map,
        int ng, double dg, double t_max, double t_step,
        double diffusion_coefficient,
        const std::vector<int>& walk_seeds,
        double tau0, int n_photons, int photon_seed,
        std::vector<double>& stats) {
    stats.assign(5, 0.0);
    const int n_steps = static_cast<int>(t_max / t_step);
    const std::size_t n = static_cast<std::size_t>(ng);
    const bool has_rates = static_cast<std::size_t>(n_rate_map) == n * n * n;

    // The whole point: this vector is the only large thing allocated, it is one
    // value per step rather than four, and it never crosses into Python.
    //
    // `float`, not `double`, for two reasons. It halves the memory the photon
    // race walks -- and that race is the bottleneck once the trajectory is long,
    // making tens of millions of random reads into this array. It also matches
    // `sample_grid`, which casts the rate map to float32, so the fused path and
    // the three-call path see bit-identical rates rather than rates that agree
    // to 2.6e-8 and could in principle disagree about a photon.
    std::vector<float> k_quench;
    k_quench.reserve(static_cast<std::size_t>(std::max(0, n_steps)) *
                     std::max<std::size_t>(1, walk_seeds.size()));

    long long n_accepted = 0, n_rejected = 0;
    bool any = false;
    for (std::size_t s = 0; s < walk_seeds.size(); ++s) {
        int acc = 0, rej = 0;
        const std::size_t before = k_quench.size();
        const bool ok = internal::run_walk(
                occupancy, static_cast<std::size_t>(n_occupancy),
                mobility, static_cast<std::size_t>(n_mobility),
                ng, t_step, diffusion_coefficient, dg,
                walk_seeds[s], n_steps, acc, rej,
                [&](int /*i*/, double /*px*/, double /*py*/, double /*pz*/,
                    bool /*accepted*/, std::size_t voxel) {
                    k_quench.push_back(has_rates ? static_cast<float>(rate_map[voxel]) : 0.0f);
                });
        if (!ok) {
            k_quench.resize(before);   // a walk that never started contributes nothing
            continue;
        }
        any = true;
        n_accepted += acc;
        n_rejected += rej;
    }
    if (!any) return std::vector<double>();

    // Accumulated in float, like numpy's float32 mean, so the reported average
    // matches the three-call path exactly rather than nearly.
    float sum = 0.0f;
    long long n_contact = 0;
    for (std::size_t i = 0; i < k_quench.size(); ++i) {
        sum += k_quench[i];
        if (k_quench[i] > 0.0f) ++n_contact;
    }
    const double frames = static_cast<double>(k_quench.size());
    stats[0] = frames;
    stats[1] = static_cast<double>(n_accepted);
    stats[2] = static_cast<double>(n_rejected);
    stats[3] = frames > 0.0 ? static_cast<double>(sum / static_cast<float>(frames)) : 0.0;
    stats[4] = frames > 0.0 ? static_cast<double>(n_contact) / frames : 0.0;

    // Interleaved so one returned vector carries delay and flag, for the same
    // reason the walk packs its accept flag: a SWIG out-parameter costs ~480 ns
    // per element to convert, while a returned vector becomes a tuple and costs
    // nothing.
    const int n_ph = n_photons > 0 ? n_photons : 0;
    std::vector<double> out(static_cast<std::size_t>(n_ph) * 2, 0.0);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n_ph; ++i) {
        bool got = false;
        out[2 * i + 0] = internal::race_one_photon(
                k_quench, t_step, tau0, photon_seed, i, got);
        out[2 * i + 1] = got ? 1.0 : 0.0;
    }
    return out;
}
}  // namespace

void quenched_donor_photons(int* occupancy, int n_occupancy,
        double* mobility, int n_mobility,
        double* rate_map, int n_rate_map,
        int ng, double dg, double t_max, double t_step,
        double diffusion_coefficient,
        const std::vector<int>& walk_seeds,
        double tau0, int n_photons, int photon_seed,
        std::vector<double>& stats, double** out_view, int* n_out_view) {
    internal::copy_to_view(quenched_donor_photons_impl(occupancy, n_occupancy, mobility, n_mobility, rate_map, n_rate_map, ng, dg, t_max, t_step, diffusion_coefficient, walk_seeds, tau0, n_photons, photon_seed, stats),
                           out_view, n_out_view);
}

IMPBFF_END_NAMESPACE
