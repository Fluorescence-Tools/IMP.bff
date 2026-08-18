/**
 * \file BrownianWalk.cpp
 * \brief A rejection-sampled Brownian walk inside an occupancy grid.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/BrownianWalk.h>

#include <cmath>
#include <random>

IMPBFF_BEGIN_NAMESPACE

std::vector<double> brownian_walk_in_volume(
        const std::vector<int>& occupancy, const std::vector<double>& mobility,
        int ng, double dg, double t_max, double t_step,
        double diffusion_coefficient, int seed,
        std::vector<int>& accepted, std::vector<int>& counts) {
    counts.assign(2, 0);
    const int n_steps = static_cast<int>(t_max / t_step);
    std::vector<double> xyz;
    accepted.assign(n_steps > 0 ? n_steps : 0, 0);
    if (n_steps <= 0 || ng <= 0) return xyz;
    xyz.assign(static_cast<std::size_t>(n_steps) * 3, 0.0);

    std::mt19937_64 rng(seed >= 0 ? static_cast<std::uint64_t>(seed)
                                  : std::random_device{}());
    // In voxel units. 2 D dt is the variance of ONE Cartesian component; the
    // total three-dimensional MSD is 6 D dt and belongs to no single axis.
    // This read sqrt(2 * D * 3 * dt) until 2026-08-18, i.e. it used the 3-D
    // total as one component's width, so the walk diffused at 3D while
    // GridDiffusionSolver -- the same dynamics as a density -- gave exactly 2Dt.
    const double sigma = std::sqrt(2.0 * diffusion_coefficient * t_step) / dg;
    std::normal_distribution<double> gauss(0.0, sigma);
    std::uniform_int_distribution<int> pick(0, ng - 1);

    const std::size_t n = static_cast<std::size_t>(ng);
    const bool has_mobility = mobility.size() == n * n * n;

    double px = 0.0, py = 0.0, pz = 0.0;
    bool found = false;
    for (int attempt = 0; attempt < 1000; ++attempt) {
        const int ix = pick(rng), iy = pick(rng), iz = pick(rng);
        if (occupancy[(static_cast<std::size_t>(ix) * n + iy) * n + iz] > 0) {
            px = ix; py = iy; pz = iz;
            found = true;
            break;
        }
    }
    if (!found) return std::vector<double>();

    xyz[0] = px; xyz[1] = py; xyz[2] = pz;
    accepted[0] = 1;
    counts[0] = 1;
    // The offset is integral, matching grids.grid_center_index(): the float
    // corner (ng - 1) / 2 disagrees with it on every even ng.
    const double half = (ng - 1) / 2;
    for (int i = 1; i < n_steps; ++i) {
        const std::size_t k =
                (static_cast<std::size_t>(static_cast<int>(px)) * n +
                 static_cast<int>(py)) * n + static_cast<int>(pz);
        double scale = 1.0;
        if (has_mobility) scale = std::sqrt(mobility[k]);

        const double nx = px + gauss(rng) * scale;
        const double ny = py + gauss(rng) * scale;
        const double nz = pz + gauss(rng) * scale;

        const int ix = static_cast<int>(nx);
        const int iy = static_cast<int>(ny);
        const int iz = static_cast<int>(nz);
        if (ix >= 0 && ix < ng && iy >= 0 && iy < ng && iz >= 0 && iz < ng &&
            occupancy[(static_cast<std::size_t>(ix) * n + iy) * n + iz] > 0) {
            px = nx; py = ny; pz = nz;
            accepted[i] = 1;
            ++counts[0];
        } else {
            // Rejected: stay put, but still emit the frame.
            ++counts[1];
        }
        xyz[3 * i + 0] = px;
        xyz[3 * i + 1] = py;
        xyz[3 * i + 2] = pz;
    }
    for (std::size_t i = 0; i < xyz.size(); ++i) xyz[i] = (xyz[i] - half) * dg;
    return xyz;
}

IMPBFF_END_NAMESPACE
