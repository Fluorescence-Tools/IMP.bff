/**
 * \file BrownianWalk.cpp
 * \brief A rejection-sampled Brownian walk inside an occupancy grid.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/BrownianWalk.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/RandomWalk.h>

#include <cmath>
#include <random>

IMPBFF_BEGIN_NAMESPACE

void brownian_walk_in_volume(
        int* occupancy, int n_occupancy, double* mobility, int n_mobility,
        int ng, double dg, double t_max, double t_step,
        double diffusion_coefficient, int seed, std::vector<int>& counts,
        double** out_view, int* n_out_view) {
    counts.assign(2, 0);
    const int n_steps = static_cast<int>(t_max / t_step);
    if (n_steps <= 0 || ng <= 0) {
        internal::new_double_view(0, out_view, n_out_view);
        return;
    }
    // Four per step: x, y, z, accepted. See the header for why the flag rides
    // here rather than in a separate out-parameter.
    double* xyz = internal::new_double_view(
            static_cast<std::size_t>(n_steps) * 4, out_view, n_out_view);
    if (xyz == nullptr) return;

    // The offset is integral, matching grids.grid_center_index(): the float
    // corner (ng - 1) / 2 disagrees with it on every even ng.
    const double half = (ng - 1) / 2;
    int n_acc = 0, n_rej = 0;
    const bool ok = internal::run_walk(
            occupancy, static_cast<std::size_t>(n_occupancy),
            mobility, static_cast<std::size_t>(n_mobility),
            ng, t_step, diffusion_coefficient, dg, seed, n_steps, n_acc, n_rej,
            [&](int i, double px, double py, double pz, bool accepted,
                std::size_t /*voxel*/) {
                xyz[4 * i + 0] = (px - half) * dg;
                xyz[4 * i + 1] = (py - half) * dg;
                xyz[4 * i + 2] = (pz - half) * dg;
                xyz[4 * i + 3] = accepted ? 1.0 : 0.0;
            });
    if (!ok) {
        // No accessible starting voxel: an empty view, and the buffer just
        // published is replaced rather than leaked.
        std::free(xyz);
        internal::new_double_view(0, out_view, n_out_view);
        return;
    }
    counts[0] = n_acc;
    counts[1] = n_rej;
}

IMPBFF_END_NAMESPACE
