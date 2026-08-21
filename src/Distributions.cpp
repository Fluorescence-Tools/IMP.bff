/**
 * \file Distributions.cpp
 * \brief Probability distributions used by the dye and linker models.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/Distributions.h>
#include <IMP/bff/internal/Normalize.h>
#include <IMP/bff/internal/OutputView.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

IMPBFF_BEGIN_NAMESPACE

//! The vector kernel behind the public (void, view-publishing) wrapper.
/*! Kept private so the skew and two-Gaussian kernels can reuse the normal
    density without allocating a second managed view. */
std::vector<double> normal_density(const std::vector<double>& x, double loc,
                                   double scale) {
    std::vector<double> y(x.size(), 0.0);
    if (scale == 0.0) return y;
    const double a = 1.0 / (std::sqrt(2.0 * M_PI) * scale);
    const double two_s2 = 2.0 * scale * scale;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double d = x[i] - loc;
        y[i] = a * std::exp(-(d * d) / two_s2);
    }
    return y;
}

void poisson_0toN(double lam, int n, double** out_view, int* n_out_view) {
    std::vector<double> p;
    if (n > 0) {
        p.assign(static_cast<std::size_t>(n), 0.0);
        p[0] = std::exp(-lam);
        for (std::size_t i = 1; i < p.size(); ++i) {
            p[i] = p[i - 1] * lam / static_cast<double>(i);
        }
    }
    internal::copy_to_view(p, out_view, n_out_view);
}

void normal_distribution(const std::vector<double>& x, double loc,
                         double scale, bool norm,
                         double** out_view, int* n_out_view) {
    std::vector<double> y(x.size(), 0.0);
    if (scale != 0.0) {
        const double a = 1.0 / (std::sqrt(2.0 * M_PI) * scale);
        const double two_s2 = 2.0 * scale * scale;
        for (std::size_t i = 0; i < x.size(); ++i) {
            const double d = x[i] - loc;
            y[i] = a * std::exp(-(d * d) / two_s2);
        }
        if (norm) internal::normalize_sum(y);
    }
    internal::copy_to_view(y, out_view, n_out_view);
}

void generalized_normal_distribution(const std::vector<double>& x, double loc,
                                     double scale, double shape, bool norm,
                                     double** out_view, int* n_out_view) {
    std::vector<double> z(x.size(), 0.0);
    if (shape == 0.0) {
        if (scale != 0.0) {
            for (std::size_t i = 0; i < x.size(); ++i) z[i] = (x[i] - loc) / scale;
        }
    } else {
        if (scale != 0.0) {
            // `np.spacing(1)` is the step to the next double above 1.0, which is
            // what the Python clamps a non-positive transform argument to. Using
            // zero instead would send the logarithm to -inf.
            const double tiny = std::nextafter(1.0, 2.0) - 1.0;
            for (std::size_t i = 0; i < x.size(); ++i) {
                double t = 1.0 - shape * (x[i] - loc) / scale;
                if (t < 0.0) t = tiny;
                z[i] = -std::log(t) / shape;
            }
        }
    }
    // The Python evaluates the *standard* normal at z -- loc and scale have
    // already been folded into the transform above.
    std::vector<double> n;
    if (scale != 0.0) {
        n = normal_density(z, 0.0, 1.0);
        if (norm) internal::normalize_sum(n);
    } else {
        n.assign(x.size(), 0.0);
    }
    internal::copy_to_view(n, out_view, n_out_view);
}

void distance_between_gaussian(const std::vector<double>& distances,
                               double separation_distance, double sigma,
                               bool normalize,
                               double** out_view, int* n_out_view) {
    std::vector<double> pr(distances.size(), 0.0);
    if (sigma != 0.0) {
        if (separation_distance != 0.0) {
            const std::vector<double> a =
                    normal_density(distances, separation_distance, sigma);
            const std::vector<double> b =
                    normal_density(distances, -separation_distance, sigma);
            for (std::size_t i = 0; i < distances.size(); ++i) {
                pr[i] = distances[i] / separation_distance * (a[i] - b[i]);
            }
        } else {
            const std::vector<double> n =
                    normal_density(distances, 0.0, sigma);
            const double inv_s2 = 1.0 / (sigma * sigma);
            for (std::size_t i = 0; i < distances.size(); ++i) {
                pr[i] = 2.0 * distances[i] * distances[i] * inv_s2 * n[i];
            }
        }
        if (normalize) internal::normalize_sum(pr);
    }
    internal::copy_to_view(pr, out_view, n_out_view);
}

IMPBFF_END_NAMESPACE
