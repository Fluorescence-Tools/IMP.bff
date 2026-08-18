/**
 * \file DiffusionSolver.cpp
 * \brief Explicit propagation of an excited-state density on an AV grid.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/DiffusionSolver.h>

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

IMPBFF_BEGIN_NAMESPACE

namespace {

//! One sweep, writing into an existing buffer. Shared by both entry points.
/*!
    The outer loop runs over x-slabs, so each thread writes a disjoint region
    and no synchronisation is needed. The outer shell is left at zero: the
    7-point stencil cannot be evaluated there, and leaving it *alone* instead
    would keep two-steps-ago state alive in a ping-pong buffer -- stale data
    that never decays and never diffuses, silently added into every population
    sum.
*/
void sweep(const std::vector<double>& cur, const std::vector<double>& d,
           const std::vector<double>& decay, const std::vector<double>& bounds,
           std::size_t n, bool smoluchowski, std::vector<double>& nxt) {
    std::fill(nxt.begin(), nxt.end(), 0.0);
    const long ngl = static_cast<long>(n);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long ixl = 1; ixl < ngl - 1; ++ixl) {
        const std::size_t ix = static_cast<std::size_t>(ixl);
        for (std::size_t iy = 1; iy + 1 < n; ++iy) {
            for (std::size_t iz = 1; iz + 1 < n; ++iz) {
                const std::size_t c = (ix * n + iy) * n + iz;
                if (bounds[c] == 0.0) continue;
                const double p0 = cur[c];
                const double d0 = d[c];
                const std::size_t nb[6] = {
                    ((ix - 1) * n + iy) * n + iz, ((ix + 1) * n + iy) * n + iz,
                    (ix * n + iy - 1) * n + iz,   (ix * n + iy + 1) * n + iz,
                    (ix * n + iy) * n + iz - 1,   (ix * n + iy) * n + iz + 1};
                double flux = 0.0;
                for (int q = 0; q < 6; ++q) {
                    const std::size_t m = nb[q];
                    // A zero bound is a reflecting wall: the dye cannot enter
                    // the protein, so flux to that neighbour is dropped.
                    if (smoluchowski) {
                        flux += 0.5 * (d0 + d[m]) * (p0 - cur[m]) * bounds[m];
                    } else {
                        flux += (d0 * p0 - d[m] * cur[m]) * bounds[m];
                    }
                }
                // The rate is applied as a factor: exp(-k dt) is the exact
                // solution of dp/dt = -k p over the step, so it contributes no
                // stability constraint. Subtracting k dt diverges once k dt > 1.
                nxt[c] = (p0 - flux) * decay[c];
            }
        }
    }
}

double population(const std::vector<double>& p) {
    double total = 0.0;
    for (double v : p) total += v;
    return total;
}

}  // namespace

std::vector<double> diffusion_step(
        const std::vector<double>& cur, const std::vector<double>& d,
        const std::vector<double>& decay, const std::vector<double>& bounds,
        int ng, int flux_form) {
    std::vector<double> nxt(cur.size(), 0.0);
    sweep(cur, d, decay, bounds, static_cast<std::size_t>(ng),
          flux_form == FLUX_SMOLUCHOWSKI, nxt);
    return nxt;
}

std::vector<double> diffusion_propagate(
        const std::vector<double>& cur_in, const std::vector<double>& d,
        const std::vector<double>& decay, const std::vector<double>& bounds,
        int ng, int flux_form, int n_steps, int n_out,
        std::vector<double>& fluorescence) {
    const bool smoluchowski = (flux_form == FLUX_SMOLUCHOWSKI);
    const std::size_t n = static_cast<std::size_t>(ng);
    std::vector<double> cur = cur_in;
    std::vector<double> nxt(cur.size(), 0.0);
    if (n_out < 1) n_out = 1;
    const int n_reports = n_steps / n_out + 1;
    fluorescence.assign(static_cast<std::size_t>(n_reports), 0.0);

    int i_out = 0;
    for (int step = 0; step < n_steps; ++step) {
        if (step % n_out == 0 && i_out < n_reports) {
            fluorescence[static_cast<std::size_t>(i_out++)] = population(cur);
        }
        sweep(cur, d, decay, bounds, n, smoluchowski, nxt);
        // Swap **every** step. Swapping only on odd steps while always
        // computing nxt <- f(cur) recomputes the previous step from a stale
        // buffer and throws away half the evolution.
        cur.swap(nxt);
    }
    while (i_out < n_reports) {
        fluorescence[static_cast<std::size_t>(i_out++)] = population(cur);
    }
    return cur;
}

IMPBFF_END_NAMESPACE
