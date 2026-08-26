/**
 * \file GridDiffusionSolver.cpp
 * \brief The field picture of a tethered dye: an occupancy density, propagated.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/GridDiffusionSolver.h>
#include <IMP/bff/internal/GridShape.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/exception.h>

#include <algorithm>
#include <string>
#include <utility>
#include <cmath>
#include <cstdlib>

IMPBFF_BEGIN_NAMESPACE

GridDiffusionSolver::GridDiffusionSolver(
        const std::vector<double>& diffusion_map,
        const std::vector<double>& bounds, const std::vector<double>& density,
        const std::vector<double>& rate_map, double dg, double t_step,
        const std::string& flux_form, bool check_stability)
        : density_(density), bounds_(bounds),
          diffusion_map_(diffusion_map), rate_map_(rate_map), dg_(dg),
          t_step_(t_step),
          ng_(internal::cube_side(bounds.size())),
          check_stability_(check_stability), n_iterations_(0) {
    if (flux_form != "smoluchowski" && flux_form != "ito") {
        IMP_THROW("flux_form must be 'smoluchowski' or 'ito', not '"
                          << flux_form << "'", IMP::ValueException);
    }
    flux_form_ = flux_form == "ito" ? FLUX_ITO : FLUX_SMOLUCHOWSKI;
    if (!bounds_.empty() && ng_ == 0) {
        IMP_THROW("the bounds grid must be cubic; " << bounds_.size()
                          << " values are not a whole cube",
                  IMP::ValueException);
    }
    // The mask is used as a multiplier everywhere, so it is stored as one:
    // a caller may hand in a count, a probability or a bool array, and only
    // "inside or not" is meant.
    for (std::size_t i = 0; i < bounds_.size(); ++i) {
        bounds_[i] = bounds_[i] > 0.0 ? 1.0 : 0.0;
    }
    if (rate_map_.empty()) rate_map_.assign(density_.size(), 0.0);
}

void GridDiffusionSolver::validate() const {
    // The 7-point stencil cannot be evaluated on the outer shell, so a domain
    // touching it would lose population there with no warning.
    const std::size_t n = static_cast<std::size_t>(ng_);
    for (int ix = 0; ix < ng_; ++ix) {
        for (int iy = 0; iy < ng_; ++iy) {
            for (int iz = 0; iz < ng_; ++iz) {
                const bool edge = ix == 0 || iy == 0 || iz == 0 ||
                                  ix == ng_ - 1 || iy == ng_ - 1 || iz == ng_ - 1;
                if (!edge) continue;
                if (bounds_[(static_cast<std::size_t>(ix) * n + iy) * n + iz] != 0.0) {
                    IMP_THROW(
                            "The domain reaches the outer shell of the grid, "
                            "where the 7-point stencil cannot be evaluated and "
                            "population is discarded. Enlarge the grid so the "
                            "accessible volume is surrounded by at least one "
                            "empty voxel.",
                            IMP::ValueException);
                }
            }
        }
    }
    if (!check_stability_) return;
    double d_max = 0.0;
    for (std::size_t i = 0; i < diffusion_map_.size(); ++i) {
        d_max = std::max(d_max, diffusion_map_[i]);
    }
    // **Without k_max, deliberately.** The rate is carried as exp(-k dt),
    // applied as a factor -- the exact solution of dp/dt = -k p over the step
    // -- so it contributes no stability constraint at all here. It is the
    // `1 - k dt` form that diverges once k dt > 1, and this scheme does not
    // use it. Passing k_max would reject steps that are perfectly stable:
    // on a site with k_max = 96 1/ns it cuts the allowed step by a third.
    const double limit = diffusion_stability_limit(d_max, dg_);
    if (t_step_ > limit) {
        IMP_THROW("t_step " << t_step_ << " ns exceeds the explicit-scheme "
                            << "stability limit dg^2 / (6 D_max) = "
                            << limit << " ns for dg = " << dg_
                            << " A and D_max = " << d_max
                            << " A^2/ns. The scheme would diverge rather "
                            << "than lose accuracy.",
                  IMP::ValueException);
    }
}

std::vector<double> GridDiffusionSolver::start_density(
        const std::vector<double>& p) const {
    std::vector<double> cur(p.size());
    double total = 0.0;
    for (std::size_t i = 0; i < p.size(); ++i) {
        cur[i] = p[i] * (i < bounds_.size() ? bounds_[i] : 0.0);
        total += cur[i];
    }
    if (total > 0.0) {
        for (std::size_t i = 0; i < cur.size(); ++i) cur[i] /= total;
    }
    return cur;
}

GridDiffusionResult GridDiffusionSolver::run(int n_steps, int n_out) {
    validate();
    if (n_out < 1) n_out = 1;
    const std::vector<double> cur = start_density(density_);

    // Fold dt into the coefficients so the inner loop is pure arithmetic, and
    // carry the rate as exp(-k dt) -- the exact solution of dp/dt = -k p over
    // the step, so it contributes no stability constraint of its own.
    std::vector<double> d(diffusion_map_.size());
    const double scale = t_step_ / (dg_ * dg_);
    for (std::size_t i = 0; i < d.size(); ++i) d[i] = diffusion_map_[i] * scale;
    std::vector<double> decay(rate_map_.size());
    for (std::size_t i = 0; i < decay.size(); ++i) {
        decay[i] = std::exp(-rate_map_[i] * t_step_);
    }

    double* fluo_buf = nullptr;
    int n_fluo_buf = 0;
    double* final_density = nullptr;
    int n_final = 0;
    diffusion_propagate(cur, d, decay, bounds_, ng_, flux_form_, n_steps, n_out,
                        &fluo_buf, &n_fluo_buf, &final_density, &n_final);
    // Both outputs are malloc'ed views the caller owns; this one is read into
    // a vector and released here, the density below.
    std::vector<double> fluorescence;
    if (fluo_buf != nullptr) {
        fluorescence.assign(fluo_buf, fluo_buf + std::max(0, n_fluo_buf));
        std::free(fluo_buf);
    }
    n_iterations_ += n_steps;

    const int n_reports = n_steps / n_out + 1;
    const std::size_t n_fluo =
            std::min<std::size_t>(fluorescence.size(),
                                  static_cast<std::size_t>(std::max(0, n_reports)));

    // The time axis is the caller's report grid, in ns.
    std::vector<double> time;
    time.reserve(n_fluo);
    for (std::size_t i = 0; i < n_fluo; ++i) {
        time.push_back(static_cast<double>(i) * t_step_ * n_out);
    }
    std::vector<double> final(final_density,
                              final_density + std::max(n_final, 0));
    if (final_density != nullptr) {
        density_.assign(final_density, final_density + n_final);
        std::free(final_density);
    }
    fluorescence.resize(n_fluo);
    return GridDiffusionResult(std::move(time), std::move(fluorescence),
                               std::move(final));
}

GridDiffusionGradient GridDiffusionSolver::gradient(
        const std::vector<double>& dL_dF, int n_steps, int n_out,
        const std::vector<double>& density) {
    validate();
    if (n_out < 1) n_out = 1;
    const int n_reports = n_steps / n_out + 1;
    if (static_cast<int>(dL_dF.size()) != n_reports) {
        IMP_THROW("dL_dF has " << dL_dF.size() << " entries; run(n_steps="
                               << n_steps << ", n_out=" << n_out
                               << ") reports " << n_reports << " points",
                  IMP::ValueException);
    }
    const std::vector<double>& p = density.empty() ? density_ : density;
    double total = 0.0;
    for (std::size_t i = 0; i < p.size(); ++i) {
        total += p[i] * (i < bounds_.size() ? bounds_[i] : 0.0);
    }
    const std::vector<double> cur = start_density(p);

    std::vector<double> d(diffusion_map_.size());
    const double scale = t_step_ / (dg_ * dg_);
    for (std::size_t i = 0; i < d.size(); ++i) d[i] = diffusion_map_[i] * scale;
    std::vector<double> decay(rate_map_.size());
    for (std::size_t i = 0; i < decay.size(); ++i) {
        decay[i] = std::exp(-rate_map_[i] * t_step_);
    }

    double* g = nullptr;
    int ng3 = 0;
    diffusion_propagate_adjoint(cur, d, decay, bounds_, ng_, flux_form_,
                                n_steps, n_out, dL_dF, &g, &ng3);
    const std::size_t nv = density_.size();
    std::vector<double> gd(nv, 0.0), gk(nv, 0.0), gp(nv, 0.0);
    if (g != nullptr) {
        // The chain rule through the folding run() applies:
        //   d = D dt/dg^2      -> dL/dD  = dL/dd * dt/dg^2
        //   decay = exp(-k dt) -> dL/dk  = dL/ddecay * (-dt) * decay
        //   cur = p b / sum(p b) -> dL/dp = b (g - <g, cur>) / total
        double dot = 0.0;
        for (std::size_t i = 0; i < nv; ++i) dot += g[2 * nv + i] * cur[i];
        for (std::size_t i = 0; i < nv; ++i) {
            gd[i] = g[i] * scale;
            gk[i] = g[nv + i] * (-t_step_) * decay[i];
            gp[i] = total > 0.0 ? bounds_[i] * (g[2 * nv + i] - dot) / total : 0.0;
        }
        std::free(g);
    }
    return GridDiffusionGradient(std::move(gd), std::move(gk), std::move(gp));
}

void GridDiffusionSolver::equilibrium(int n_steps, double tolerance,
                                      int n_check, double** out_view,
                                      int* n_out_view) {
    validate();
    std::vector<double> cur = start_density(density_);
    std::vector<double> d(diffusion_map_.size());
    const double scale = t_step_ / (dg_ * dg_);
    for (std::size_t i = 0; i < d.size(); ++i) d[i] = diffusion_map_[i] * scale;
    const std::vector<double> one(cur.size(), 1.0);   // exp(-0 dt): no decay

    std::vector<double> previous = cur;
    int remaining = n_steps;
    const int chunk = std::max(1, n_check);
    while (remaining > 0) {
        const int take = std::min(chunk, remaining);
        double* fluo = nullptr;
        int n_fluo = 0;
        double* next = nullptr;
        int n_next = 0;
        diffusion_propagate(cur, d, one, bounds_, ng_, flux_form_, take, take,
                            &fluo, &n_fluo, &next, &n_next);
        std::free(fluo);
        if (next == nullptr) break;
        cur.assign(next, next + n_next);
        std::free(next);
        n_iterations_ += take;
        remaining -= take;
        double drift = 0.0;
        for (std::size_t i = 0; i < cur.size(); ++i) {
            drift += std::fabs(cur[i] - previous[i]);
        }
        if (drift < tolerance) break;
        previous = cur;
    }
    double total = 0.0;
    for (std::size_t i = 0; i < cur.size(); ++i) total += cur[i];
    if (total > 0.0) {
        for (std::size_t i = 0; i < cur.size(); ++i) cur[i] /= total;
    }
    density_ = cur;
    internal::copy_to_view(cur, out_view, n_out_view);
}

void GridDiffusionSolver::get_density(double** o, int* n) const {
    internal::copy_to_view(density_, o, n);
}
void GridDiffusionSolver::get_bounds(double** o, int* n) const {
    internal::copy_to_view(bounds_, o, n);
}
void GridDiffusionSolver::get_diffusion_map(double** o, int* n) const {
    internal::copy_to_view(diffusion_map_, o, n);
}
void GridDiffusionSolver::get_rate_map(double** o, int* n) const {
    internal::copy_to_view(rate_map_, o, n);
}

void GridDiffusionResult::get_time(double** out_view, int* n_out_view) const {
    internal::copy_to_view(time_, out_view, n_out_view);
}
void GridDiffusionResult::get_fluorescence(double** out_view, int* n_out_view) const {
    internal::copy_to_view(fluorescence_, out_view, n_out_view);
}
void GridDiffusionResult::get_density(double** out_view, int* n_out_view) const {
    internal::copy_to_view(density_, out_view, n_out_view);
}

void GridDiffusionGradient::get_d_diffusion(double** out_view, int* n_out_view) const {
    internal::copy_to_view(d_diffusion_, out_view, n_out_view);
}
void GridDiffusionGradient::get_d_rate(double** out_view, int* n_out_view) const {
    internal::copy_to_view(d_rate_, out_view, n_out_view);
}
void GridDiffusionGradient::get_d_density(double** out_view, int* n_out_view) const {
    internal::copy_to_view(d_density_, out_view, n_out_view);
}

IMPBFF_END_NAMESPACE
