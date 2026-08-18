/**
 * \file Distributions.cpp
 * \brief Probability distributions used by the dye and linker models.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/Distributions.h>

#include <algorithm>
#include <cmath>
#include <limits>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! Divide by the sum, leaving an all-zero vector alone.
void normalize_in_place(std::vector<double>& y) {
    double total = 0.0;
    for (double v : y) total += v;
    if (total > 0.0) {
        for (double& v : y) v /= total;
    }
}

}  // namespace

std::vector<double> poisson_0toN(double lam, int n) {
    if (n <= 0) return {};
    std::vector<double> p(static_cast<std::size_t>(n));
    p[0] = std::exp(-lam);
    for (std::size_t i = 1; i < p.size(); ++i) {
        p[i] = p[i - 1] * lam / static_cast<double>(i);
    }
    return p;
}

std::vector<double> normal_distribution(
        const std::vector<double>& x, double loc, double scale, bool norm) {
    std::vector<double> y(x.size(), 0.0);
    if (scale == 0.0) return y;
    const double a = 1.0 / (std::sqrt(2.0 * M_PI) * scale);
    const double two_s2 = 2.0 * scale * scale;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double d = x[i] - loc;
        y[i] = a * std::exp(-(d * d) / two_s2);
    }
    if (norm) normalize_in_place(y);
    return y;
}

std::vector<double> generalized_normal_distribution(
        const std::vector<double>& x, double loc, double scale, double shape,
        bool norm) {
    std::vector<double> z(x.size(), 0.0);
    if (shape == 0.0) {
        if (scale == 0.0) return std::vector<double>(x.size(), 0.0);
        for (std::size_t i = 0; i < x.size(); ++i) z[i] = (x[i] - loc) / scale;
    } else {
        if (scale == 0.0) return std::vector<double>(x.size(), 0.0);
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
    // The Python evaluates the *standard* normal at z -- loc and scale have
    // already been folded into the transform above.
    std::vector<double> y = normal_distribution(z, 0.0, 1.0, false);
    if (norm) normalize_in_place(y);
    return y;
}

std::vector<double> distance_between_gaussian(
        const std::vector<double>& distances, double separation_distance,
        double sigma, bool normalize) {
    std::vector<double> pr(distances.size(), 0.0);
    if (sigma == 0.0) return pr;
    if (separation_distance != 0.0) {
        const std::vector<double> a =
                normal_distribution(distances, separation_distance, sigma, false);
        const std::vector<double> b =
                normal_distribution(distances, -separation_distance, sigma, false);
        for (std::size_t i = 0; i < distances.size(); ++i) {
            pr[i] = distances[i] / separation_distance * (a[i] - b[i]);
        }
    } else {
        const std::vector<double> n =
                normal_distribution(distances, 0.0, sigma, false);
        const double inv_s2 = 1.0 / (sigma * sigma);
        for (std::size_t i = 0; i < distances.size(); ++i) {
            pr[i] = 2.0 * distances[i] * distances[i] * inv_s2 * n[i];
        }
    }
    if (normalize) normalize_in_place(pr);
    return pr;
}

IMPBFF_END_NAMESPACE
