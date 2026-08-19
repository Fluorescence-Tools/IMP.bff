/**
 * \file DiffusionSolver.cpp
 * \brief Explicit propagation of an excited-state density on an AV grid.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/DiffusionSolver.h>
#include <IMP/bff/internal/OutputView.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

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

//! One adjoint sweep: pbar_prev = A^T pbar, and the parameter gradients.
/*!
    Written as a *gather* over every voxel m, in one pass over its six
    neighbours: the contribution to \f$\bar p_m\f$ from m itself (the diagonal
    of A, if m is active -- interior and in the domain) and from each active
    neighbour c whose flux read \f$p_m\f$ (the off-diagonal), plus the parameter
    gradients, all from the same reads. Each output is owned by one iteration,
    so the loop parallelises like the forward one, and nothing is tabulated:
    the sweep is memory-bound, and a table is more memory. `cur` is the
    forward state the sweep was applied to (\f$p_n\f$), `pbar` the adjoint of
    its result (\f$\bar p_{n+1}\f$).

    Smoluchowski: \f$p'_c = \delta_c\,[\,p_c (1 - \sum_m w_{cm}) + \sum_m w_{cm} p_m\,]\f$
    with \f$w_{cm} = \tfrac12 (d_c + d_m) b_m\f$. Ito:
    \f$p'_c = \delta_c\,[\,p_c (1 - d_c \sum_m b_m) + \sum_m d_m b_m p_m\,]\f$.
    Both for interior c with \f$b_c \ne 0\f$; every other row of A is zero.
    Templated on the flux form so the inner loop carries no branch. A shell
    voxel has a zero row but, if its bound is non-zero, was read by an active
    neighbour -- the Python wrapper forbids that domain, the kernel stays exact.
*/
template <bool kSmoluchowski>
void adjoint_sweep(const std::vector<double>& cur, const std::vector<double>& d,
                   const std::vector<double>& decay, const std::vector<double>& bounds,
                   const std::vector<double>& pbar, std::size_t n,
                   std::vector<double>& pbar_prev, std::vector<double>& dbar,
                   std::vector<double>& decaybar) {
    const long ngl = static_cast<long>(n);
    const long stride[6] = {-ngl * ngl, ngl * ngl, -ngl, ngl, -1, 1};
    const double* P = cur.data();
    const double* D = d.data();
    const double* B = bounds.data();
    const double* DEC = decay.data();
    const double* PB = pbar.data();
    double* PBP = pbar_prev.data();
    double* DB = dbar.data();
    double* DECB = decaybar.data();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long ixl = 0; ixl < ngl; ++ixl) {
        for (long iy = 0; iy < ngl; ++iy) {
            for (long iz = 0; iz < ngl; ++iz) {
                const long m = (ixl * ngl + iy) * ngl + iz;
                const double bm = B[m];
                if (bm == 0.0) { PBP[m] = 0.0; continue; }  // out of the domain: zero row of A, and never read through a non-zero bound
                const bool xi = ixl >= 1 && ixl < ngl - 1, yi = iy >= 1 && iy < ngl - 1, zi = iz >= 1 && iz < ngl - 1;
                const bool m_interior = xi && yi && zi;
                const bool m_active = m_interior && bm != 0.0;
                // Is neighbour q interior (so possibly active)? Only the two
                // outermost layers can fail this.
                const bool k_int[6] = {ixl >= 2 && yi && zi, ixl < ngl - 2 && yi && zi,
                                       xi && iy >= 2 && zi, xi && iy < ngl - 2 && zi,
                                       xi && yi && iz >= 2, xi && yi && iz < ngl - 2};
                const double p0 = P[m], d0 = D[m], dec0 = DEC[m], pb0 = PB[m];
                double coef = 0.0, flux = 0.0, ddiag = 0.0;   // the diagonal (m active)
                double pb_off = 0.0, db_off = 0.0;             // gathered from active neighbours
                for (int q = 0; q < 6; ++q) {
                    if (!m_interior && !k_int[q]) continue;    // neither m's row nor k's exists
                    const long k = m + stride[q];
                    const double pk = P[k], dk = D[k], bk = B[k];
                    if (m_active) {
                        if (kSmoluchowski) {
                            const double w = 0.5 * (d0 + dk) * bk;
                            coef += w;
                            const double dp = p0 - pk;
                            flux += w * dp;
                            ddiag += 0.5 * dp * bk;
                        } else {
                            coef += d0 * bk;
                            flux += (d0 * p0 - dk * pk) * bk;
                            ddiag += p0 * bk;
                        }
                    }
                    if (k_int[q] && bk != 0.0) {   // k is active: it read p_m and d_m through b_m
                        const double pbk = PB[k], deck = DEC[k];
                        if (kSmoluchowski) {
                            pb_off += deck * 0.5 * (dk + d0) * bm * pbk;
                            db_off += -deck * pbk * 0.5 * (pk - p0) * bm;
                        } else {
                            pb_off += deck * d0 * bm * pbk;
                            db_off += deck * pbk * p0 * bm;
                        }
                    }
                }
                double pb = pb_off, db = db_off;
                if (m_active) {
                    pb += dec0 * (1.0 - coef) * pb0;
                    DECB[m] += pb0 * (p0 - flux);
                    db += -dec0 * pb0 * ddiag;
                }
                PBP[m] = pb;
                DB[m] += db;
            }
        }
    }
}

}  // namespace

void diffusion_step(
        const std::vector<double>& cur, const std::vector<double>& d,
        const std::vector<double>& decay, const std::vector<double>& bounds,
        int ng, int flux_form, double** out_view, int* n_out_view) {
    std::vector<double> nxt(cur.size(), 0.0);
    sweep(cur, d, decay, bounds, static_cast<std::size_t>(ng),
          flux_form == FLUX_SMOLUCHOWSKI, nxt);
    internal::copy_to_view(nxt, out_view, n_out_view);
}

void diffusion_propagate(
        const std::vector<double>& cur_in, const std::vector<double>& d,
        const std::vector<double>& decay, const std::vector<double>& bounds,
        int ng, int flux_form, int n_steps, int n_out,
        std::vector<double>& fluorescence,
        double** out_view, int* n_out_view) {
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
    // `cur` and `nxt` are swapped every step, so which buffer holds the answer
    // is only known here. One memcpy, against ng^3 x 66 ns of marshalling.
    internal::copy_to_view(cur, out_view, n_out_view);
}

void diffusion_propagate_adjoint(
        const std::vector<double>& cur_in, const std::vector<double>& d,
        const std::vector<double>& decay, const std::vector<double>& bounds,
        int ng, int flux_form, int n_steps, int n_out,
        const std::vector<double>& dL_dF,
        std::vector<double>& dL_dd, std::vector<double>& dL_ddecay,
        std::vector<double>& dL_dcur) {
    const bool smoluchowski = (flux_form == FLUX_SMOLUCHOWSKI);
    const std::size_t n = static_cast<std::size_t>(ng);
    const std::size_t nv = cur_in.size();
    if (n_out < 1) n_out = 1;
    if (n_steps < 0) n_steps = 0;
    const int n_reports = n_steps / n_out + 1;
    if (static_cast<int>(dL_dF.size()) != n_reports)
        throw std::invalid_argument("diffusion_propagate_adjoint: dL_dF must have n_steps / n_out + 1 entries");

    dL_dd.assign(nv, 0.0);
    dL_ddecay.assign(nv, 0.0);
    dL_dcur.assign(nv, 0.0);

    // Which report index each step's *pre-sweep* state feeds (or -1); the
    // reports diffusion_propagate() emits after the loop all read the final state.
    auto report_of = [&](int step) { return (step % n_out == 0 && step / n_out < n_reports) ? step / n_out : -1; };
    int first_final_report = n_reports;  // reports [first_final_report, n_reports) read p_{n_steps}
    for (int k = 0; k < n_reports; ++k)
        if (k * n_out >= n_steps) { first_final_report = k; break; }

    // Forward with checkpoints every K steps.
    const int K = std::max(1, static_cast<int>(std::ceil(std::sqrt(static_cast<double>(std::max(n_steps, 1))))));
    const int n_ck = n_steps / K + 1;
    std::vector<std::vector<double>> checkpoint(static_cast<std::size_t>(n_ck));
    {
        std::vector<double> cur = cur_in, nxt(nv, 0.0);
        for (int step = 0; step < n_steps; ++step) {
            if (step % K == 0) checkpoint[static_cast<std::size_t>(step / K)] = cur;
            sweep(cur, d, decay, bounds, n, smoluchowski, nxt);
            cur.swap(nxt);
        }
        if (n_steps % K == 0) checkpoint[static_cast<std::size_t>(n_steps / K)] = cur;
    }

    // Seed: every report of the final state.
    std::vector<double> pbar(nv, 0.0), pbar_prev(nv, 0.0);
    {
        double seed = 0.0;
        for (int k = first_final_report; k < n_reports; ++k) seed += dL_dF[static_cast<std::size_t>(k)];
        if (seed != 0.0) std::fill(pbar.begin(), pbar.end(), seed);
    }

    // Backward, segment by segment: re-run the forward within the segment to
    // recover p_s for every step in it, then sweep the adjoint back through.
    std::vector<std::vector<double>> seg;
    std::vector<double> nxt(nv, 0.0);
    for (int seg_start = ((n_steps - 1) / K) * K; seg_start >= 0 && n_steps > 0; seg_start -= K) {
        const int seg_end = std::min(seg_start + K, n_steps);  // steps [seg_start, seg_end)
        const int len = seg_end - seg_start;
        seg.resize(static_cast<std::size_t>(len));
        seg[0] = checkpoint[static_cast<std::size_t>(seg_start / K)];
        for (int i = 1; i < len; ++i) {
            seg[static_cast<std::size_t>(i)].assign(nv, 0.0);
            sweep(seg[static_cast<std::size_t>(i - 1)], d, decay, bounds, n, smoluchowski, seg[static_cast<std::size_t>(i)]);
        }
        for (int step = seg_end - 1; step >= seg_start; --step) {
            // pbar is the adjoint of p_{step+1}; produce the adjoint of p_step.
            if (smoluchowski)
                adjoint_sweep<true>(seg[static_cast<std::size_t>(step - seg_start)], d, decay, bounds, pbar, n,
                                    pbar_prev, dL_dd, dL_ddecay);
            else
                adjoint_sweep<false>(seg[static_cast<std::size_t>(step - seg_start)], d, decay, bounds, pbar, n,
                                     pbar_prev, dL_dd, dL_ddecay);
            pbar.swap(pbar_prev);
            const int r = report_of(step);
            if (r >= 0 && r < first_final_report) {
                const double f = dL_dF[static_cast<std::size_t>(r)];
                if (f != 0.0) for (double& v : pbar) v += f;
            }
        }
    }
    dL_dcur = pbar;
}

void diffusion_propagate_adjoint(
        const std::vector<double>& cur, const std::vector<double>& d,
        const std::vector<double>& decay, const std::vector<double>& bounds,
        int ng, int flux_form, int n_steps, int n_out,
        const std::vector<double>& dL_dF,
        double** out_view, int* n_out_view) {
    std::vector<double> gd, gdec, gcur;
    diffusion_propagate_adjoint(cur, d, decay, bounds, ng, flux_form, n_steps, n_out,
                                dL_dF, gd, gdec, gcur);
    const std::size_t nv = gd.size();
    double* buf = internal::new_double_view(3 * nv, out_view, n_out_view);
    if (buf == nullptr) return;
    std::copy(gd.begin(), gd.end(), buf);
    std::copy(gdec.begin(), gdec.end(), buf + nv);
    std::copy(gcur.begin(), gcur.end(), buf + 2 * nv);
}

IMPBFF_END_NAMESPACE
