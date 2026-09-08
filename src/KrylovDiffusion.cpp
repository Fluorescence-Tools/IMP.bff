// SPDX-License-Identifier: BSD-3-Clause
/**
 * The whole propagation as one matrix exponential.
 *
 * `diffusion_propagate()` takes `n_steps` explicit Euler sweeps because
 * forward Euler on a diffusion operator is stable only for
 * `dt <= dg^2/(6D)`. That bound is why the solver costs `dg^-5`. This file
 * does not take a better step -- it stops stepping.
 *
 * Two properties of *this* operator make that possible, and a general PDE has
 * neither:
 *
 *  * **It does not change over a run.** `d`, `decay` and `bounds` are fixed,
 *    so the propagation is `p(t) = exp(-L t) p0` for one `L`, not a sequence
 *    of different problems.
 *  * **It is symmetrisable.** The Smoluchowski flux form is symmetric as it
 *    stands. The Ito form is not, but `q = sqrt(d) p` makes it so: with
 *    `A[c,m] = -d_m b_m`, the similar matrix `D^(1/2) A D^(-1/2)` has entries
 *    `-sqrt(d_c d_m) b_m`, which is symmetric. So one Lanczos serves both.
 *
 * So: project `L` onto an `m`-dimensional Krylov space, take the exponential
 * of the resulting `m x m` tridiagonal, and read **every** reported time out
 * of the same projection. Measured at 64-128 operator applications where the
 * stepping route needs four thousand, at better accuracy than the stepping
 * route has (okf/validation/larger_time_steps.md).
 *
 * Three things this exploits that are worth naming, because each one removes
 * a cost that a textbook Lanczos would pay:
 *
 *  * **No reorthogonalisation.** It would cost `m^2 n` flops -- at m = 256,
 *    n = 92k that is more work than the four thousand steps being replaced --
 *    and it is measurably unnecessary here: identical to fourteen digits at
 *    every `m` tested. For `exp(-tL)v` the ghost copies belong to
 *    already-converged eigenvalues.
 *  * **No stored basis for the trace.** The observable is a plain sum, so all
 *    a Krylov vector has to leave behind is one scalar, its (weighted) column
 *    sum. Three vectors of memory for any `m`.
 *  * **A second pass instead of storage for the density.** Storing the basis
 *    is 544 MB at ng = 81, m = 128. Re-running the recursion with the
 *    coefficients already known costs `m` more matvecs and no memory, so the
 *    density is `2m` rather than `m` -- still 15-30x fewer than stepping.
 *
 * What this is *not*: a reproduction of `diffusion_propagate()`. It solves the
 * semi-discrete equation, where the stepping route adds a forward-Euler and a
 * Lie-splitting error of about 2e-6. This is the more accurate answer, and
 * therefore a different one.
 */
#include <IMP/bff/DiffusionSolver.h>


#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

namespace {

/// Which voxels the operator acts on: interior, and inside the domain.
/*! The outer shell is excluded for the same reason the sweep excludes it --
    the 7-point stencil cannot be evaluated there. */
inline bool active(std::size_t ix, std::size_t iy, std::size_t iz, std::size_t n,
                   const std::vector<double>& bounds, std::size_t c) {
    return ix >= 1 && iy >= 1 && iz >= 1 && ix + 1 < n && iy + 1 < n && iz + 1 < n
           && bounds[c] != 0.0;
}

/// y = L x, with L = A + K the generator of dp/dt = -L p.
/*!
    Smoluchowski: `(Ax)[c] = sum_q w_q (x_c - x_m)`, `w_q = 0.5(d_c+d_m) b_m`.
    Ito, already symmetrised by `q = sqrt(d) p`:
    `(Ax)[c] = d_c (sum_q b_m) x_c - sum_q sqrt(d_c d_m) b_m x_m`.
    Both plus `k_c x_c`. Rows outside the active set are zero.

    This is the same memory traffic as one `lattice_sweep`, so a matvec here
    and a step there cost the same -- which is what makes the count of them
    the whole story.
*/
template <bool kSmoluchowski>
void generator_apply(const double* x, const std::vector<double>& d,
                     const std::vector<double>& k, const std::vector<double>& bounds,
                     const std::vector<double>& rootd, std::size_t n, double* y) {
    // Zeroing the whole grid every matvec is a full pass of writes for the
    // sake of the shell: 1 GB over 256 matvecs at ng = 81. The interior loop
    // below writes every interior voxel, so only the shell needs it.
    const long ngl = static_cast<long>(n);
    for (std::size_t ix = 0; ix < n; ++ix)
        for (std::size_t iy = 0; iy < n; ++iy)
            for (std::size_t iz = 0; iz < n; ++iz)
                if (ix == 0 || iy == 0 || iz == 0 ||
                    ix + 1 == n || iy + 1 == n || iz + 1 == n)
                    y[(ix * n + iy) * n + iz] = 0.0;
    const long stride[6] = {-ngl * ngl, ngl * ngl, -ngl, ngl, -1, 1};
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long ixl = 1; ixl < ngl - 1; ++ixl) {
        const std::size_t ix = static_cast<std::size_t>(ixl);
        for (std::size_t iy = 1; iy + 1 < n; ++iy) {
            for (std::size_t iz = 1; iz + 1 < n; ++iz) {
                const std::size_t c = (ix * n + iy) * n + iz;
                if (bounds[c] == 0.0) { y[c] = 0.0; continue; }
                const double x0 = x[c], d0 = d[c];
                double diag = 0.0, off = 0.0;
                for (int q = 0; q < 6; ++q) {
                    const std::size_t m = static_cast<std::size_t>(
                        static_cast<long>(c) + stride[q]);
                    const double bm = bounds[m];
                    if (bm == 0.0) continue;              // reflecting wall
                    if (kSmoluchowski) {
                        const double w = 0.5 * (d0 + d[m]) * bm;
                        diag += w;
                        off += w * x[m];
                    } else {
                        diag += d0 * bm;
                        off += rootd[c] * rootd[m] * bm * x[m];
                    }
                }
                y[c] = (diag + k[c]) * x0 - off;
            }
        }
    }
}

/// Eigenvalues and eigenvectors of a symmetric tridiagonal, QL with implicit shifts.
/*!
    `m` is at most a few hundred, so this costs nothing beside the matvecs and
    saves a LAPACK dependency in a library that does not otherwise have one.
    `diag` and `off` are consumed; `z` comes in as the identity and leaves as
    the eigenvectors in columns.
*/
void tridiagonal_eigen(std::vector<double>& diag, std::vector<double>& off,
                       std::vector<double>& z, int m) {
    off.push_back(0.0);
    for (int l = 0; l < m; ++l) {
        for (int iter = 0; iter < 50; ++iter) {
            int mm = l;
            for (; mm < m - 1; ++mm) {
                const double dd = std::fabs(diag[mm]) + std::fabs(diag[mm + 1]);
                if (std::fabs(off[mm]) <= 1e-300 + 1e-16 * dd) break;
            }
            if (mm == l) break;
            double g = (diag[l + 1] - diag[l]) / (2.0 * off[l]);
            double r = std::hypot(g, 1.0);
            g = diag[mm] - diag[l] + off[l] / (g + (g >= 0 ? std::fabs(r) : -std::fabs(r)));
            double s = 1.0, c = 1.0, p = 0.0;
            int i = mm - 1;
            for (; i >= l; --i) {
                double f = s * off[i], b = c * off[i];
                r = std::hypot(f, g);
                off[i + 1] = r;
                if (r == 0.0) { diag[i + 1] -= p; off[mm] = 0.0; break; }
                s = f / r; c = g / r;
                g = diag[i + 1] - p;
                r = (diag[i] - g) * s + 2.0 * c * b;
                p = s * r;
                diag[i + 1] = g + p;
                g = c * r - b;
                for (int q = 0; q < m; ++q) {        // Eigenvektoren mitdrehen
                    f = z[q * m + i + 1];
                    z[q * m + i + 1] = s * z[q * m + i] + c * f;
                    z[q * m + i] = c * z[q * m + i] - s * f;
                }
            }
            if (r == 0.0 && i >= l) continue;
            diag[l] -= p; off[l] = g; off[mm] = 0.0;
        }
    }
}

}  // namespace

void diffusion_propagate_krylov(
        const std::vector<double>& cur, const std::vector<double>& d,
        const std::vector<double>& decay, const std::vector<double>& bounds,
        int ng, int flux_form, int n_steps, int n_out, int krylov_dim,
        double** out_fluorescence, int* n_out_fluorescence,
        double** out_view, int* n_out_view) {
    const std::size_t n = static_cast<std::size_t>(ng);
    const std::size_t nv = n * n * n;
    if (ng < 3) IMP_THROW("diffusion_propagate_krylov: ng must be >= 3", IMP::ValueException);
    if (cur.size() != nv || d.size() != nv || decay.size() != nv || bounds.size() != nv)
        IMP_THROW("diffusion_propagate_krylov: grids must be ng^3 long", IMP::ValueException);
    if (n_steps < 0 || n_out <= 0)
        IMP_THROW("diffusion_propagate_krylov: n_steps >= 0, n_out > 0", IMP::ValueException);
    const bool smol = (flux_form == FLUX_SMOLUCHOWSKI);

    // k from the decay factor: decay = exp(-k dt), and the generator wants k.
    std::vector<double> k(nv, 0.0), rootd(nv, 0.0), weight(nv, 0.0);
    for (std::size_t ix = 0; ix < n; ++ix)
        for (std::size_t iy = 0; iy < n; ++iy)
            for (std::size_t iz = 0; iz < n; ++iz) {
                const std::size_t c = (ix * n + iy) * n + iz;
                if (!active(ix, iy, iz, n, bounds, c)) continue;
                if (decay[c] <= 0.0)
                    IMP_THROW("diffusion_propagate_krylov: decay must be positive in the domain "
                        "(it is exp(-k dt))", IMP::ValueException);
                k[c] = -std::log(decay[c]);
                if (!smol) {
                    if (d[c] <= 0.0)
                        IMP_THROW("diffusion_propagate_krylov: the Ito flux form needs d > 0 "
                            "everywhere in the domain -- the change of variables that makes "
                            "it symmetric is q = sqrt(d) p. Use FLUX_SMOLUCHOWSKI, or "
                            "diffusion_propagate().", IMP::ValueException);
                    rootd[c] = std::sqrt(d[c]);
                }
                // What the population sum weighs this voxel by. In q-space the
                // density is q/sqrt(d), so the plain sum becomes a weighted one.
                weight[c] = smol ? 1.0 : 1.0 / rootd[c];
            }

    // The initial vector, in whichever space the operator is symmetric in,
    // masked to the active set so nothing outside can leak into the sums.
    std::vector<double> v0(nv, 0.0);
    for (std::size_t ix = 0; ix < n; ++ix)
        for (std::size_t iy = 0; iy < n; ++iy)
            for (std::size_t iz = 0; iz < n; ++iz) {
                const std::size_t c = (ix * n + iy) * n + iz;
                if (active(ix, iy, iz, n, bounds, c))
                    v0[c] = smol ? cur[c] : cur[c] * rootd[c];
            }

    const int n_rep = n_steps / n_out + 1;
    *n_out_fluorescence = n_rep;
    *out_fluorescence = static_cast<double*>(std::malloc(sizeof(double) * n_rep));
    const bool want_density = (out_view != nullptr && n_out_view != nullptr);
    if (want_density) {
        *n_out_view = static_cast<int>(nv);
        *out_view = static_cast<double*>(std::calloc(nv, sizeof(double)));
    }

    double beta0 = 0.0;
    for (std::size_t c = 0; c < nv; ++c) beta0 += v0[c] * v0[c];
    beta0 = std::sqrt(beta0);
    if (beta0 == 0.0) {                       // nothing to propagate
        for (int i = 0; i < n_rep; ++i) (*out_fluorescence)[i] = 0.0;
        return;
    }

    int m = krylov_dim > 0 ? krylov_dim : 128;
    m = std::min<int>(m, static_cast<int>(nv));

    // The Lanczos vectors are zero outside the domain, so every dot product,
    // axpy and norm walks a list of the active voxels rather than the cube --
    // 43% of it at ng = 81 -- and walks it in parallel. Serial full-grid
    // versions of these made eight threads worth nothing here: the generator
    // scaled and they did not, so they became the whole cost.
    std::vector<std::size_t> act;
    act.reserve(nv / 2);
    for (std::size_t ix = 0; ix < n; ++ix)
        for (std::size_t iy = 0; iy < n; ++iy)
            for (std::size_t iz = 0; iz < n; ++iz) {
                const std::size_t c = (ix * n + iy) * n + iz;
                if (active(ix, iy, iz, n, bounds, c)) act.push_back(c);
            }
    const long na = static_cast<long>(act.size());
    const std::size_t* AC = act.data();

    auto apply = [&](const double* x, double* y) {
        if (smol) generator_apply<true>(x, d, k, bounds, rootd, n, y);
        else      generator_apply<false>(x, d, k, bounds, rootd, n, y);
    };

    // --- pass one: alpha, beta, and one scalar per Krylov vector ------------
    std::vector<double> vj(nv), vp(nv, 0.0), w(nv);
    for (std::size_t c = 0; c < nv; ++c) vj[c] = v0[c] / beta0;
    std::vector<double> al, be, u;
    int mm = m;
    double* VJ = vj.data(); double* VP = vp.data(); double* W = w.data();
    const double* WT = weight.data();
    for (int j = 0; j < m; ++j) {
        apply(VJ, W);
        // One pass, two accumulators: the projection alpha_j and the scalar
        // this Krylov vector leaves behind for the trace.
        double a = 0.0, s = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : a, s) schedule(static)
#endif
        for (long t = 0; t < na; ++t) {
            const std::size_t c = AC[t];
            a += VJ[c] * W[c];
            s += WT[c] * VJ[c];
        }
        al.push_back(a); u.push_back(s);
        const double bprev = j ? be[j - 1] : 0.0;
        if (j + 1 >= m) break;
        double b2 = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : b2) schedule(static)
#endif
        for (long t = 0; t < na; ++t) {
            const std::size_t c = AC[t];
            const double wv = W[c] - a * VJ[c] - bprev * VP[c];
            W[c] = wv; b2 += wv * wv;
        }
        const double b = std::sqrt(b2);
        // An exhausted subspace is an exact answer, not a failure.
        if (b <= 1e-12 * (std::fabs(a) + 1.0)) { mm = j + 1; break; }
        be.push_back(b);
        std::swap(VP, VJ);
        const double inv = 1.0 / b;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (long t = 0; t < na; ++t) VJ[AC[t]] = W[AC[t]] * inv;
    }
    al.resize(mm); be.resize(std::max(mm - 1, 0)); u.resize(mm);

    // --- the m x m exponential ---------------------------------------------
    std::vector<double> lam(al), offd(be), z(static_cast<std::size_t>(mm) * mm, 0.0);
    for (int i = 0; i < mm; ++i) z[static_cast<std::size_t>(i) * mm + i] = 1.0;
    tridiagonal_eigen(lam, offd, z, mm);
    // L is positive semidefinite; a tiny negative eigenvalue is round-off, and
    // exp(-lam t) with t in the thousands would turn it into a large one.
    for (int i = 0; i < mm; ++i) if (lam[i] < 0.0) lam[i] = 0.0;

    // F(t) = beta0 * sum_i c_i exp(-lam_i t): the decay comes out as a sum of
    // exponentials with explicit rates and amplitudes, which is the form a
    // fluorescence decay is fitted with.
    std::vector<double> amp(mm, 0.0);
    for (int i = 0; i < mm; ++i) {
        double uq = 0.0;
        for (int j = 0; j < mm; ++j) uq += u[j] * z[static_cast<std::size_t>(j) * mm + i];
        amp[i] = beta0 * uq * z[i];             // z[0*mm + i] = first row
    }
    for (int r = 0; r < n_rep; ++r) {
        const double t = static_cast<double>(r) * n_out;
        double f = 0.0;
        for (int i = 0; i < mm; ++i) f += amp[i] * std::exp(-lam[i] * t);
        (*out_fluorescence)[r] = f;
    }

    if (!want_density) return;

    // --- pass two: the density at the final time, without storing the basis --
    const double tend = static_cast<double>(n_steps);
    std::vector<double> y(mm, 0.0);
    for (int j = 0; j < mm; ++j) {
        double s = 0.0;
        for (int i = 0; i < mm; ++i)
            s += z[static_cast<std::size_t>(j) * mm + i] * std::exp(-lam[i] * tend) * z[i];
        y[j] = beta0 * s;
    }
    std::vector<double> acc(nv, 0.0);
    double* ACC = acc.data();
    VJ = vj.data(); VP = vp.data(); W = w.data();
    std::fill(vp.begin(), vp.end(), 0.0);
    const double invb0 = 1.0 / beta0;
    for (std::size_t c = 0; c < nv; ++c) vj[c] = v0[c] * invb0;
    for (int j = 0; j < mm; ++j) {
        const double yj = y[j];
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (long t = 0; t < na; ++t) ACC[AC[t]] += yj * VJ[AC[t]];
        if (j + 1 >= mm) break;
        apply(VJ, W);
        const double aj = al[j], bp = j ? be[j - 1] : 0.0, inv = 1.0 / be[j];
        double* NV = VP;                      // the vector about to be overwritten
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (long t = 0; t < na; ++t) {
            const std::size_t c = AC[t];
            NV[c] = (W[c] - aj * VJ[c] - bp * VP[c]) * inv;
        }
        std::swap(VP, VJ);
    }
    // `weight` is 1 inside the domain for Smoluchowski and 1/sqrt(d) for Ito,
    // and zero outside either way, so it undoes the change of variables and
    // masks the shell in one multiply.
    for (std::size_t c = 0; c < nv; ++c) (*out_view)[c] = acc[c] * weight[c];
}

void diffusion_trace_krylov(
        const std::vector<double>& cur, const std::vector<double>& d,
        const std::vector<double>& decay, const std::vector<double>& bounds,
        int ng, int flux_form, int n_steps, int n_out, int krylov_dim,
        double** out_fluorescence, int* n_out_fluorescence) {
    diffusion_propagate_krylov(cur, d, decay, bounds, ng, flux_form, n_steps,
                               n_out, krylov_dim, out_fluorescence,
                               n_out_fluorescence, nullptr, nullptr);
}

IMPBFF_END_NAMESPACE
