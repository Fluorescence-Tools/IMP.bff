/**
 * \file FRETPair.cpp
 * \brief Distances, orientation factors and efficiencies over all pairs.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/FRETPair.h>
#include <IMP/bff/internal/OutputView.h>

#include <cmath>
#include <cstdlib>
#include <limits>

IMPBFF_BEGIN_NAMESPACE

void fret_pair_matrices(
        double* points1, int n_points1, double* points2, int n_points2,
        double* mu1, int n_mu1, double* mu2, int n_mu2, int n1, int n2,
        double** out_view, int* n_out_view) {
    (void)n_points1; (void)n_points2;
    const std::size_t total =
            static_cast<std::size_t>(std::max(0, n1)) * std::max(0, n2) * 2;
    // malloc, not new[]: numpy's ARGOUTVIEWM typemap takes ownership and
    // releases it with free(). Mixing the two is undefined behaviour, and the
    // pairing is invisible from this side of the boundary.
    double* out = static_cast<double*>(std::calloc(total ? total : 1,
                                                   sizeof(double)));
    if (out == nullptr) {
        // An empty view rather than a null pointer: numpy would wrap the null
        // and the caller would fault on first touch, a long way from here.
        *out_view = static_cast<double*>(std::calloc(1, sizeof(double)));
        *n_out_view = 0;
        return;
    }
    *out_view = out;
    *n_out_view = static_cast<int>(total);
    if (n1 <= 0 || n2 <= 0) return;
    const std::size_t block = static_cast<std::size_t>(n1) * n2;
    const bool oriented = n_mu1 == n1 * 3 && n_mu2 == n2 * 3;

    // Normalised acceptor dipoles, once. Recomputing them inside the inner loop
    // would repeat n1 square roots per acceptor for no reason.
    std::vector<double> a_hat;
    if (oriented) {
        a_hat.assign(static_cast<std::size_t>(n2) * 3, 0.0);
        for (int j = 0; j < n2; ++j) {
            const double x = mu2[3 * j], y = mu2[3 * j + 1], z = mu2[3 * j + 2];
            const double n = std::sqrt(x * x + y * y + z * z);
            if (n > 0.0) {
                a_hat[3 * j] = x / n;
                a_hat[3 * j + 1] = y / n;
                a_hat[3 * j + 2] = z / n;
            }
        }
    }

#pragma omp parallel for schedule(static)
    for (int i = 0; i < n1; ++i) {
        const double px = points1[3 * i], py = points1[3 * i + 1],
                     pz = points1[3 * i + 2];
        double dx_hat = 0.0, dy_hat = 0.0, dz_hat = 0.0;
        if (oriented) {
            const double x = mu1[3 * i], y = mu1[3 * i + 1], z = mu1[3 * i + 2];
            const double n = std::sqrt(x * x + y * y + z * z);
            if (n > 0.0) { dx_hat = x / n; dy_hat = y / n; dz_hat = z / n; }
        }
        for (int j = 0; j < n2; ++j) {
            const double rx = px - points2[3 * j];
            const double ry = py - points2[3 * j + 1];
            const double rz = pz - points2[3 * j + 2];
            const double d = std::sqrt(rx * rx + ry * ry + rz * rz);
            const std::size_t k = static_cast<std::size_t>(i) * n2 + j;
            out[k] = d;
            if (!oriented) {
                out[block + k] = 2.0 / 3.0;
                continue;
            }
            // A zero separation has no direction; the Python divided with
            // `where=r_norm > 0` and left the unit vector at zero, so the two
            // projection terms vanish and kappa^2 falls back to (mu_D . mu_A)^2.
            double ux = 0.0, uy = 0.0, uz = 0.0;
            if (d > 0.0) { ux = rx / d; uy = ry / d; uz = rz / d; }
            const double ax = a_hat[3 * j], ay = a_hat[3 * j + 1],
                         az = a_hat[3 * j + 2];
            const double cos_da = dx_hat * ax + dy_hat * ay + dz_hat * az;
            const double cos_dr = dx_hat * ux + dy_hat * uy + dz_hat * uz;
            const double cos_ar = ax * ux + ay * uy + az * uz;
            const double kappa = cos_da - 3.0 * cos_dr * cos_ar;
            out[block + k] = kappa * kappa;
        }
    }
}

void fret_pair_efficiency_matrices(
        double* r, int n_r, double* kappa2, int n_kappa2,
        double forster_radius, double** out_view, int* n_out_view) {
    const int n = n_r < n_kappa2 ? n_r : n_kappa2;
    double* out = internal::new_double_view(
            static_cast<std::size_t>(std::max(0, n)) * 2, out_view, n_out_view);
    if (out == nullptr || n <= 0) return;
    const std::size_t block = static_cast<std::size_t>(n);

#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        const double ratio = r[i] / forster_radius;
        const double r3 = ratio * ratio * ratio;
        const double ratio6 = r3 * r3;
        const double k2 = kappa2[i];

        // E = 1 / (1 + (2/3) ratio6 / kappa2). Both degenerate limits mean
        // complete transfer: kappa2 = 0 with ratio6 = 0 gives 0/0, and any
        // ratio6 with kappa2 = 0 gives +inf in the denominator's numerator.
        double e;
        if (k2 == 0.0) {
            e = ratio6 == 0.0 ? 1.0 : 0.0;   // nan -> 1, else 1/(1+inf) -> 0
        } else {
            const double denom = 1.0 + (2.0 / 3.0) * ratio6 / k2;
            e = std::isfinite(denom) ? 1.0 / denom : 0.0;
            if (std::isnan(e)) e = 1.0;
        }
        out[i] = e;
        // Left unsanitised on purpose: an infinite rate at zero separation is
        // the right answer, and the caller decides what to do with it.
        out[block + i] = ratio6 == 0.0
                                 ? std::numeric_limits<double>::infinity()
                                 : 1.5 * k2 / ratio6;
    }
}

void FRETPairGeometry::get_R(double** out_view, int* n_out_view) const {
    internal::copy_to_view(R, out_view, n_out_view);
}
void FRETPairGeometry::get_kappa2(double** out_view, int* n_out_view) const {
    internal::copy_to_view(kappa2, out_view, n_out_view);
}
void FRETPairGeometry::get_weight(double** out_view, int* n_out_view) const {
    internal::copy_to_view(weight, out_view, n_out_view);
}

void FRETPairEfficiencies::get_E(double** out_view, int* n_out_view) const {
    internal::copy_to_view(E, out_view, n_out_view);
}
void FRETPairEfficiencies::get_rate_ratio(double** out_view,
                                          int* n_out_view) const {
    internal::copy_to_view(rate_ratio, out_view, n_out_view);
}
void FRETPairEfficiencies::get_k_fret(double** out_view, int* n_out_view) const {
    internal::copy_to_view(k_fret, out_view, n_out_view);
}
void FRETPairEfficiencies::get_R(double** out_view, int* n_out_view) const {
    internal::copy_to_view(R, out_view, n_out_view);
}
void FRETPairEfficiencies::get_kappa2(double** out_view, int* n_out_view) const {
    internal::copy_to_view(kappa2, out_view, n_out_view);
}
void FRETPairEfficiencies::get_weight(double** out_view, int* n_out_view) const {
    internal::copy_to_view(weight, out_view, n_out_view);
}

FRETPairGeometry fret_pair_geometry(const std::vector<double>& points1,
                                    const std::vector<double>& weights1,
                                    const std::vector<double>& points2,
                                    const std::vector<double>& weights2,
                                    const std::vector<double>& mu1,
                                    const std::vector<double>& mu2) {
    FRETPairGeometry out;
    out.n1 = static_cast<int>(points1.size() / 3);
    out.n2 = static_cast<int>(points2.size() / 3);
    const std::size_t n = static_cast<std::size_t>(out.n1) * out.n2;
    if (n == 0) return out;

    double* packed = NULL;
    int n_packed = 0;
    // The kernel takes raw pointers because SWIG hands it numpy buffers; the
    // casts are the same buffers, const-stripped for that signature.
    fret_pair_matrices(const_cast<double*>(points1.data()),
                       static_cast<int>(points1.size()),
                       const_cast<double*>(points2.data()),
                       static_cast<int>(points2.size()),
                       mu1.empty() ? NULL : const_cast<double*>(mu1.data()),
                       static_cast<int>(mu1.size()),
                       mu2.empty() ? NULL : const_cast<double*>(mu2.data()),
                       static_cast<int>(mu2.size()),
                       out.n1, out.n2, &packed, &n_packed);
    out.R.assign(packed, packed + n);
    out.kappa2.assign(packed + n, packed + 2 * n);
    std::free(packed);

    out.weight.resize(n);
    double total = 0.0;
    for (int i = 0; i < out.n1; ++i) {
        for (int j = 0; j < out.n2; ++j) {
            const double w = weights1[i] * weights2[j];
            out.weight[static_cast<std::size_t>(i) * out.n2 + j] = w;
            total += w;
        }
    }
    if (total > 0.0) {
        for (std::size_t i = 0; i < n; ++i) out.weight[i] /= total;
    }
    out.kappa2_avg = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        out.kappa2_avg += out.kappa2[i] * out.weight[i];
    }
    return out;
}

FRETPairEfficiencies fret_pair_efficiencies(const FRETPairGeometry& geometry,
                                            double forster_radius,
                                            double tau0) {
    FRETPairEfficiencies out;
    out.n1 = geometry.n1;
    out.n2 = geometry.n2;
    out.kappa2_avg = geometry.kappa2_avg;
    out.forster_radius = forster_radius;
    out.R = geometry.R;
    out.kappa2 = geometry.kappa2;
    out.weight = geometry.weight;
    const std::size_t n = geometry.R.size();
    if (n == 0) return out;

    double* packed = NULL;
    int n_packed = 0;
    fret_pair_efficiency_matrices(
            const_cast<double*>(geometry.R.data()),
            static_cast<int>(geometry.R.size()),
            const_cast<double*>(geometry.kappa2.data()),
            static_cast<int>(geometry.kappa2.size()), forster_radius, &packed,
            &n_packed);
    out.E.assign(packed, packed + n);
    out.rate_ratio.assign(packed + n, packed + 2 * n);
    std::free(packed);

    double e_static = 0.0, e_dyn1 = 0.0, rate_avg = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        e_static += out.E[i] * geometry.weight[i];
        // `dynamic1` uses the *ensemble* kappa^2 rather than the per-pair one,
        // so it is a different expression from `E` and is formed here.
        const double x = geometry.R[i] / forster_radius;
        const double x3 = x * x * x;
        const double ratio6 = x3 * x3;
        double d1 = 1.0;
        if (out.kappa2_avg > 0.0 && std::isfinite(ratio6)) {
            d1 = 1.0 / (1.0 + (2.0 / 3.0) / out.kappa2_avg * ratio6);
        }
        e_dyn1 += d1 * geometry.weight[i];
        if (std::isfinite(out.rate_ratio[i])) {
            rate_avg += out.rate_ratio[i] * geometry.weight[i];
        }
    }
    out.static_efficiency = e_static;
    out.dynamic1 = e_dyn1;
    out.dynamic2 = rate_avg / (rate_avg + 1.0);
    if (tau0 > 0.0) {
        out.k_fret.resize(n);
        for (std::size_t i = 0; i < n; ++i) out.k_fret[i] = out.rate_ratio[i] / tau0;
    }
    return out;
}

IMPBFF_END_NAMESPACE
