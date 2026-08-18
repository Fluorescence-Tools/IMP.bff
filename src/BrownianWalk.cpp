/**
 * \file BrownianWalk.cpp
 * \brief A rejection-sampled Brownian walk inside an occupancy grid.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/BrownianWalk.h>
#include <IMP/bff/internal/RandomWalk.h>

#include <cmath>
#include <random>

IMPBFF_BEGIN_NAMESPACE

std::vector<double> brownian_walk_in_volume(
        const std::vector<int>& occupancy, const std::vector<double>& mobility,
        int ng, double dg, double t_max, double t_step,
        double diffusion_coefficient, int seed, std::vector<int>& counts) {
    counts.assign(2, 0);
    const int n_steps = static_cast<int>(t_max / t_step);
    std::vector<double> xyz;
    if (n_steps <= 0 || ng <= 0) return xyz;
    // Four per step: x, y, z, accepted. See the header for why the flag rides
    // here rather than in an out-parameter.
    xyz.assign(static_cast<std::size_t>(n_steps) * 4, 0.0);

    // The offset is integral, matching grids.grid_center_index(): the float
    // corner (ng - 1) / 2 disagrees with it on every even ng.
    const double half = (ng - 1) / 2;
    int n_acc = 0, n_rej = 0;
    const bool ok = internal::run_walk(
            occupancy, mobility, ng, t_step, diffusion_coefficient, dg, seed,
            n_steps, n_acc, n_rej,
            [&](int i, double px, double py, double pz, bool accepted,
                std::size_t /*voxel*/) {
                xyz[4 * i + 0] = (px - half) * dg;
                xyz[4 * i + 1] = (py - half) * dg;
                xyz[4 * i + 2] = (pz - half) * dg;
                xyz[4 * i + 3] = accepted ? 1.0 : 0.0;
            });
    if (!ok) return std::vector<double>();
    counts[0] = n_acc;
    counts[1] = n_rej;
    return xyz;
}

IMPBFF_END_NAMESPACE
