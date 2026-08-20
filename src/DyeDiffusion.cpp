/**
 * \file DyeDiffusion.cpp
 * \brief The dye's Brownian walk in its accessible volume, as an object.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/DyeDiffusion.h>
#include <IMP/bff/BrownianWalk.h>
#include <IMP/bff/internal/GridShape.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/exception.h>

#include <cmath>
#include <cstdlib>
#include <random>
#include <thread>

IMPBFF_BEGIN_NAMESPACE

const int MAX_PARALLEL_TRAJECTORIES = 8;

int default_trajectory_count() {
    const unsigned n = std::thread::hardware_concurrency();
    const unsigned cap = static_cast<unsigned>(MAX_PARALLEL_TRAJECTORIES);
    const unsigned capped = n == 0 ? 1u : (n < cap ? n : cap);
    return static_cast<int>(capped);
}

int resolve_trajectory_count(int requested) {
    if (requested < 0) return default_trajectory_count();
    if (requested < 1) return 1;
    return requested < MAX_PARALLEL_TRAJECTORIES ? requested
                                                 : MAX_PARALLEL_TRAJECTORIES;
}

std::vector<int> trajectory_seeds(int random_seed, int n_trajectories) {
    const int n = n_trajectories > 0 ? n_trajectories : 1;
    std::vector<int> seeds;
    seeds.reserve(n);
    if (random_seed < 0) {
        // One walk gets -1 -- "draw freely" -- rather than a seed of its own,
        // so an unseeded single walk stays as unconstrained as it was.
        if (n <= 1) { seeds.push_back(-1); return seeds; }
        std::random_device source;
        for (int i = 0; i < n; ++i) {
            seeds.push_back(static_cast<int>(source() & 0x7fffffffu));
        }
        return seeds;
    }
    // A large prime stride, so trajectories from one base seed do not share the
    // low-order pattern `base + i` would give them.
    const long long max_seed = (1LL << 31) - 1;
    for (int i = 0; i < n; ++i) {
        seeds.push_back(static_cast<int>(
                (static_cast<long long>(random_seed) + i * 104729LL) % max_seed));
    }
    return seeds;
}

DyeDiffusionSimulation::DyeDiffusionSimulation(
        const std::vector<int>& density, double dg,
        const std::vector<double>& x0,
        const std::vector<int>& slow_density,
        const std::vector<double>& slow_factor_map,
        const std::vector<double>& quenching_rate_map)
        : density_(density), slow_density_(slow_density),
          slow_factor_map_(slow_factor_map),
          quenching_rate_map_(quenching_rate_map),
          x0_(x0.empty() ? std::vector<double>(3, 0.0) : x0),
          dg_(dg), t_step_(0.0), ng_(internal::cube_side(density.size())),
          n_accepted_(0), n_rejected_(0) {
    if (!density_.empty() && ng_ == 0) {
        IMP_THROW("an occupancy grid must be cubic; " << density_.size()
                          << " values are not a whole cube",
                  IMP::ValueException);
    }
    if (x0_.size() != 3) {
        IMP_THROW("the grid anchor is three coordinates", IMP::ValueException);
    }
}

int DyeDiffusionSimulation::run(double D, double slow_fact, double t_step,
                                double t_max, int n_trajectories, int seed) {
    t_step_ = t_step;
    if (n_trajectories < 0) n_trajectories = default_trajectory_count();
    if (n_trajectories < 1) n_trajectories = 1;

    // Mobility is a *field*: one scaling per voxel. The scalar-plus-mask form
    // is one way to build it, and building it here is what collapsed two
    // near-identical kernels into one.
    std::vector<double> mobility;
    if (!slow_factor_map_.empty()) {
        mobility = slow_factor_map_;
    } else if (slow_fact != 1.0 && !slow_density_.empty()) {
        mobility.resize(density_.size(), 1.0);
        for (std::size_t i = 0; i < mobility.size() && i < slow_density_.size(); ++i) {
            if (slow_density_[i]) mobility[i] = slow_fact;
        }
    }

    // Distinct seeds per walk, derived from one. Reusing the base seed would
    // run n identical trajectories and concatenate them, which looks like n
    // times the sampling and is none of it.
    //
    // The stride is a large prime, so trajectories from one base do not share
    // the low-order pattern that `base + i` would give them. This is exactly
    // the derivation the Python used, kept digit for digit: it is what lets a
    // walk seeded here reproduce a walk seeded there, so the port could be
    // gated bit-for-bit instead of distributionally.
    std::vector<int> seeds(n_trajectories);
    if (seed < 0) {
        std::random_device rd;
        for (int i = 0; i < n_trajectories; ++i) {
            seeds[i] = static_cast<int>(rd() & 0x7fffffff);
        }
    } else {
        const long long max_seed = (1LL << 31) - 1;
        for (int i = 0; i < n_trajectories; ++i) {
            seeds[i] = static_cast<int>(
                    (static_cast<long long>(seed) + i * 104729LL) % max_seed);
        }
    }

    std::vector<std::vector<double>> parts(n_trajectories);
    std::vector<std::vector<int>> counts(n_trajectories);

    auto one = [&](int i) {
        double* buf = nullptr;
        int n = 0;
        counts[i].assign(2, 0);
        brownian_walk_in_volume(
                density_.empty() ? nullptr : const_cast<int*>(density_.data()),
                static_cast<int>(density_.size()),
                mobility.empty() ? nullptr : mobility.data(),
                static_cast<int>(mobility.size()),
                ng_, dg_, t_max, t_step, D, seeds[i], counts[i], &buf, &n);
        if (buf != nullptr) {
            parts[i].assign(buf, buf + n);
            std::free(buf);
        }
    };

    if (n_trajectories == 1) {
        one(0);
    } else {
        // std::thread, not OpenMP: this build leaves OpenMP_CXX_FLAGS empty, so
        // every `#pragma omp` in the module is inert. The Python this replaces
        // used a ThreadPoolExecutor and got real parallelism because the kernel
        // releases the GIL; losing that in the move would be a regression.
        std::vector<std::thread> pool;
        pool.reserve(n_trajectories);
        for (int i = 0; i < n_trajectories; ++i) pool.emplace_back(one, i);
        for (auto& t : pool) t.join();
    }

    trajectory_.clear();
    n_accepted_ = n_rejected_ = 0;
    for (int i = 0; i < n_trajectories; ++i) {
        if (counts[i].size() < 2 || counts[i][0] <= 0) continue;   // no start found
        n_accepted_ += counts[i][0];
        n_rejected_ += counts[i][1];
        // The walk reports (x, y, z, accepted) per step; only the coordinates
        // are kept, and shifted into the structure's frame here so that every
        // consumer reads absolute positions.
        const std::size_t n_steps = parts[i].size() / 4;
        const std::size_t base = trajectory_.size();
        trajectory_.resize(base + n_steps * 3);
        for (std::size_t s = 0; s < n_steps; ++s) {
            trajectory_[base + 3 * s + 0] = parts[i][4 * s + 0] + x0_[0];
            trajectory_[base + 3 * s + 1] = parts[i][4 * s + 1] + x0_[1];
            trajectory_[base + 3 * s + 2] = parts[i][4 * s + 2] + x0_[2];
        }
    }
    return get_n_frames();
}

void DyeDiffusionSimulation::get_trajectory(double** out_view,
                                            int* n_out_view) const {
    internal::copy_to_view(trajectory_, out_view, n_out_view);
}

void DyeDiffusionSimulation::get_x0(double** out_view, int* n_out_view) const {
    internal::copy_to_view(x0_, out_view, n_out_view);
}

void DyeDiffusionSimulation::get_mean_position(double** out_view,
                                               int* n_out_view) const {
    const int n = get_n_frames();
    if (n == 0) {
        IMP_THROW("Run the simulation first.", IMP::ValueException);
    }
    std::vector<double> m(3, 0.0);
    for (int i = 0; i < n; ++i) {
        m[0] += trajectory_[3 * i + 0];
        m[1] += trajectory_[3 * i + 1];
        m[2] += trajectory_[3 * i + 2];
    }
    for (int k = 0; k < 3; ++k) m[k] /= n;
    internal::copy_to_view(m, out_view, n_out_view);
}

void DyeDiffusionSimulation::sample_grid(const std::vector<double>& grid, int ng,
                                         double** out_view,
                                         int* n_out_view) const {
    const int n = get_n_frames();
    if (n == 0) {
        IMP_THROW("Run the simulation first.", IMP::ValueException);
    }
    double* out = internal::new_double_view(static_cast<std::size_t>(n),
                                            out_view, n_out_view);
    if (out == nullptr) return;
    const int centre = grid_center_index(ng);
    for (int f = 0; f < n; ++f) {
        bool inside = true;
        std::size_t flat = 0;
        for (int axis = 0; axis < 3; ++axis) {
            // floor, not trunc: trunc maps [-1, 0) to 0, so a position up to one
            // voxel below the grid would be read as voxel 0.
            const double u = (trajectory_[3 * f + axis] - x0_[axis]) / dg_ + centre;
            const long idx = static_cast<long>(std::floor(u));
            if (idx < 0 || idx >= ng) { inside = false; break; }
            flat = flat * static_cast<std::size_t>(ng) + static_cast<std::size_t>(idx);
        }
        out[f] = (inside && flat < grid.size()) ? grid[flat] : 0.0;
    }
}

void DyeDiffusionSimulation::get_density(double** out_view,
                                         int* n_out_view) const {
    std::vector<double> d(density_.begin(), density_.end());
    internal::copy_to_view(d, out_view, n_out_view);
}

void DyeDiffusionSimulation::set_density(const std::vector<int>& density) {
    const int ng = internal::cube_side(density.size());
    if (!density.empty() && ng == 0) {
        IMP_THROW("an occupancy grid must be cubic; " << density.size()
                          << " values are not a whole cube",
                  IMP::ValueException);
    }
    density_ = density;
    ng_ = ng;
    trajectory_.clear();
    n_accepted_ = n_rejected_ = 0;
}

void DyeDiffusionSimulation::get_slow_factor_map(double** out_view,
                                                 int* n_out_view) const {
    internal::copy_to_view(slow_factor_map_, out_view, n_out_view);
}

void DyeDiffusionSimulation::get_slow_density(double** out_view,
                                              int* n_out_view) const {
    std::vector<double> d(slow_density_.begin(), slow_density_.end());
    internal::copy_to_view(d, out_view, n_out_view);
}

void DyeDiffusionSimulation::get_quenching_rate_map(double** out_view,
                                                    int* n_out_view) const {
    internal::copy_to_view(quenching_rate_map_, out_view, n_out_view);
}

void DyeDiffusionSimulation::get_k_quench(double** out_view,
                                          int* n_out_view) const {
    const int n = get_n_frames();
    if (n == 0) {
        IMP_THROW("Run the simulation first.", IMP::ValueException);
    }
    if (quenching_rate_map_.empty()) {
        internal::new_double_view(static_cast<std::size_t>(n), out_view, n_out_view);
        return;                                   // calloc already zeroed it
    }
    sample_grid(quenching_rate_map_, internal::cube_side(quenching_rate_map_.size()),
                out_view, n_out_view);
    // Rounded through `float`, and deliberately. The fused kernel holds its own
    // trace as `std::vector<float>` -- it halves the memory the photon race
    // walks, and that race makes tens of millions of random reads into it -- so
    // the two paths see bit-identical rates rather than rates that agree to
    // 2.6e-8 and could in principle disagree about a photon. The precision costs
    // nothing real: a PET rate constant is a transferable starting value known
    // to perhaps two significant figures, and float carries seven.
    //
    // Here rather than in the Python property that used to do it, because both
    // the C++ and the Python consumers need the same rates for the equality to
    // hold, and only one of them went through Python.
    if (*out_view != nullptr) {
        for (int i = 0; i < *n_out_view; ++i) {
            (*out_view)[i] = static_cast<float>((*out_view)[i]);
        }
    }
}

double DyeDiffusionSimulation::get_collision_fraction() const {
    const int n = get_n_frames();
    if (n == 0) return 0.0;
    double* k = nullptr;
    int nk = 0;
    get_k_quench(&k, &nk);
    if (k == nullptr) return 0.0;
    int hit = 0;
    for (int i = 0; i < nk; ++i) if (k[i] > 0.0) ++hit;
    std::free(k);
    return static_cast<double>(hit) / n;
}

IMPBFF_END_NAMESPACE
