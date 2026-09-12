/**
 * \file ProbeDiffusionSimulation.cpp
 * \brief Probe diffusion with state retained in the simulation object.
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/ProbeDiffusionSimulation.h>
#include <IMP/bff/ProbeSampling.h>
#include <IMP/bff/BrownianWalk.h>
#include <IMP/bff/internal/GridShape.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/json.h>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <random>
#include <thread>

IMPBFF_BEGIN_NAMESPACE

ProbeDiffusionSimulation::ProbeDiffusionSimulation(
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
          n_accepted_(0), n_rejected_(0),
      diffusion_coefficient_(40.0), slow_fact_(0.01), t_max_(10000.0),
      n_trajectories_(1), random_seed_(-1) {
    if (!density_.empty() && ng_ == 0) {
        IMP_THROW("an occupancy grid must be cubic; " << density_.size()
                          << " values are not a whole cube",
                  IMP::ValueException);
    }
    if (x0_.size() != 3) {
        IMP_THROW("the grid anchor is three coordinates", IMP::ValueException);
    }
}

int ProbeDiffusionSimulation::simulate(double D, double slow_fact, double t_step,
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
    // a fixed derivation, digit for digit, so a walk seeded here reproduces a
    // walk seeded anywhere else that follows it.
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
        // every `#pragma omp` in the module is inert and a trajectory fan-out
        // written with one would silently run on a single core.
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

void ProbeDiffusionSimulation::get_trajectory(double** out_view,
                                            int* n_out_view) const {
    internal::copy_to_view(trajectory_, out_view, n_out_view);
}

void ProbeDiffusionSimulation::get_x0(double** out_view, int* n_out_view) const {
    internal::copy_to_view(x0_, out_view, n_out_view);
}

void ProbeDiffusionSimulation::get_mean_position(double** out_view,
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

void ProbeDiffusionSimulation::sample_grid(const std::vector<double>& grid, int ng,
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

void ProbeDiffusionSimulation::get_density(double** out_view,
                                         int* n_out_view) const {
    std::vector<double> d(density_.begin(), density_.end());
    internal::copy_to_view(d, out_view, n_out_view);
}

void ProbeDiffusionSimulation::set_density(const std::vector<int>& density) {
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

void ProbeDiffusionSimulation::get_slow_factor_map(double** out_view,
                                                 int* n_out_view) const {
    internal::copy_to_view(slow_factor_map_, out_view, n_out_view);
}

void ProbeDiffusionSimulation::get_slow_density(double** out_view,
                                              int* n_out_view) const {
    std::vector<double> d(slow_density_.begin(), slow_density_.end());
    internal::copy_to_view(d, out_view, n_out_view);
}

void ProbeDiffusionSimulation::get_quenching_rate_map(double** out_view,
                                                    int* n_out_view) const {
    internal::copy_to_view(quenching_rate_map_, out_view, n_out_view);
}

void ProbeDiffusionSimulation::get_k_quench(double** out_view,
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
    // Here rather than in a caller, so every consumer sees the same rates and
    // the equality holds for all of them.
    if (*out_view != nullptr) {
        for (int i = 0; i < *n_out_view; ++i) {
            (*out_view)[i] = static_cast<float>((*out_view)[i]);
        }
    }
}

double ProbeDiffusionSimulation::get_collision_fraction() const {
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

// ---- the shared simulation interface -----------------------------------

std::string ProbeDiffusionSimulation::get_parameters() const {
    nlohmann::json j;
    j["diffusion_coefficient"] = diffusion_coefficient_;
    j["slow_factor"] = slow_fact_;
    j["t_step"] = t_step_;
    j["t_max"] = t_max_;
    j["n_trajectories"] = n_trajectories_;
    j["seed"] = random_seed_;
    j["voxel_edge"] = dg_;
    return j.dump();
}

void ProbeDiffusionSimulation::set_parameters(const std::string& json) {
    if (json.empty()) return;
    nlohmann::json j;
    try {
        j = nlohmann::json::parse(json);
    } catch (const std::exception& e) {
        IMP_THROW("the parameters are not JSON: " << e.what(), ValueException);
    }
    if (!j.is_object()) {
        IMP_THROW("the parameters must be a JSON object, not " << j.type_name(),
                  ValueException);
    }
    for (nlohmann::json::const_iterator it = j.begin(); it != j.end(); ++it) {
        const std::string& k = it.key();
        if (k == "diffusion_coefficient") diffusion_coefficient_ = it.value().get<double>();
        else if (k == "slow_factor") slow_fact_ = it.value().get<double>();
        else if (k == "t_step") t_step_ = it.value().get<double>();
        else if (k == "t_max") t_max_ = it.value().get<double>();
        else if (k == "n_trajectories") n_trajectories_ = it.value().get<int>();
        else if (k == "seed") random_seed_ = it.value().get<int>();
        else if (k == "voxel_edge") dg_ = it.value().get<double>();
        else {
            IMP_THROW("a grid diffusion simulation has no parameter '" << k
                      << "'; it has diffusion_coefficient, slow_factor, t_step,"
                         " t_max, n_trajectories, seed and voxel_edge",
                      ValueException);
        }
    }
    if (t_step_ <= 0.0) {
        IMP_THROW("t_step is above zero, not " << t_step_, ValueException);
    }
}

void ProbeDiffusionSimulation::get_positions(double** out_view,
                                             int* n_out_view) const {
    // the walker as it stands: the last frame, or where it started
    if (trajectory_.size() >= 3) {
        const std::size_t last = trajectory_.size() - 3;
        std::vector<double> p(trajectory_.begin() + last, trajectory_.end());
        internal::copy_to_view(p, out_view, n_out_view);
        return;
    }
    internal::copy_to_view(x0_, out_view, n_out_view);
}

void ProbeDiffusionSimulation::set_positions(const std::vector<double>& xyz) {
    if (xyz.size() != 3) {
        IMP_THROW("a walk has one position, three coordinates, not "
                  << xyz.size(), ValueException);
    }
    // the walk starts from x0_, so putting the walker somewhere *is* moving
    // the start: there is no other state to move
    x0_ = xyz;
    trajectory_.clear();
}

void ProbeDiffusionSimulation::step(int n_steps) {
    if (n_steps <= 0) return;
    simulate(diffusion_coefficient_, slow_fact_, t_step_,
             t_step_ * static_cast<double>(n_steps), 1, random_seed_);
}

ProbeSimulationTrajectory ProbeDiffusionSimulation::run(int n_steps, int write_every) {
    if (write_every < 1) write_every = 1;
    simulate(diffusion_coefficient_, slow_fact_, t_step_,
             t_step_ * static_cast<double>(n_steps), n_trajectories_,
             random_seed_);
    ProbeSimulationTrajectory out;
    out.integrator = "grid-walk";
    out.temperature = std::numeric_limits<double>::quiet_NaN();
    // t_step is ns here and the record is femtoseconds, as every trajectory is
    out.timestep_fs = t_step_ * 1e6;
    out.n_atoms = 1;
    const int n = get_n_frames();
    for (int f = 0; f < n; f += write_every) {
        out.coordinates.push_back(trajectory_[f * 3]);
        out.coordinates.push_back(trajectory_[f * 3 + 1]);
        out.coordinates.push_back(trajectory_[f * 3 + 2]);
        out.times_fs.push_back(out.timestep_fs * f);
    }
    out.n_frames = static_cast<int>(out.times_fs.size());
    return out;
}

IMPBFF_END_NAMESPACE
