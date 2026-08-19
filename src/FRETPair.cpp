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

IMPBFF_END_NAMESPACE
