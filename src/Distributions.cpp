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

std::vector<double> generalized_normal_density_impl(
        const std::vector<double>& x, double loc, double scale, double shape,
        bool norm) {
    std::vector<double> z(x.size(), 0.0);
    if (shape == 0.0) {
        if (scale != 0.0) {
            for (std::size_t i = 0; i < x.size(); ++i) z[i] = (x[i] - loc) / scale;
        }
    } else {
        if (scale != 0.0) {
            // A non-positive transform argument clamps to the step above 1.0; zero
            // instead would send the logarithm to -inf.
            const double tiny = std::nextafter(1.0, 2.0) - 1.0;
            for (std::size_t i = 0; i < x.size(); ++i) {
                double t = 1.0 - shape * (x[i] - loc) / scale;
                if (t < 0.0) t = tiny;
                z[i] = -std::log(t) / shape;
            }
        }
    }
    // The *standard* normal at z: loc and scale are already folded into the
    // transform above.
    std::vector<double> n;
    if (scale != 0.0) {
        n = normal_density(z, 0.0, 1.0);
        if (norm) internal::normalize_sum(n);
    } else {
        n.assign(x.size(), 0.0);
    }
    return n;
}

void generalized_normal_distribution(const std::vector<double>& x, double loc,
                                     double scale, double shape, bool norm,
                                     double** out_view, int* n_out_view) {
    internal::copy_to_view(
            generalized_normal_density_impl(x, loc, scale, shape, norm),
            out_view, n_out_view);
}

std::vector<double> distance_between_gaussian_impl(
        const std::vector<double>& distances, double separation_distance,
        double sigma) {
    std::vector<double> pr(distances.size(), 0.0);
    if (sigma != 0.0) {
        // `> 0` rather than `!= 0`: the reference selects the coincident branch
        // for every non-positive separation, and the separated branch divides
        // by it. A negative separation is not physical, but the two forms are
        // not the same function and this one is the reference's.
        if (separation_distance > 0.0) {
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
    }
    return pr;
}

void distance_between_gaussian(const std::vector<double>& distances,
                               double separation_distance, double sigma,
                               bool normalize,
                               double** out_view, int* n_out_view) {
    std::vector<double> pr =
            distance_between_gaussian_impl(distances, separation_distance, sigma);
    if (sigma != 0.0 && normalize) internal::normalize_sum(pr);
    internal::copy_to_view(pr, out_view, n_out_view);
}


std::vector<double> gaussian_distance_mixture_impl(
        const std::vector<double>& axis, const std::vector<double>& means,
        const std::vector<double>& sigmas, const std::vector<double>& shapes,
        const std::vector<double>& amplitudes, int kernel,
        bool normalize_components, bool normalize) {
    std::vector<double> density(axis.size(), 0.0);
    for (std::size_t c = 0; c < means.size(); ++c) {
        const double sigma = c < sigmas.size() ? sigmas[c] : 0.0;
        const double amplitude = c < amplitudes.size() ? amplitudes[c] : 0.0;
        std::vector<double> component;
        if (kernel == GAUSSIAN_MIXTURE_DISTANCE_BETWEEN_GAUSSIANS) {
            component = distance_between_gaussian_impl(axis, means[c], sigma);
            if (normalize_components) internal::normalize_sum(component);
        } else if (kernel == GAUSSIAN_MIXTURE_NORMAL) {
            component = sigma > 0.0 ? normal_density(axis, means[c], sigma)
                                    : std::vector<double>(axis.size(), 0.0);
            if (normalize_components) internal::normalize_sum(component);
        } else {
            const double shape = c < shapes.size() ? shapes[c] : 0.0;
            component = generalized_normal_density_impl(axis, means[c], sigma,
                                                        shape,
                                                        normalize_components);
        }
        for (std::size_t j = 0; j < axis.size(); ++j) {
            density[j] += amplitude * component[j];
        }
    }
    if (normalize) internal::normalize_sum(density);
    return density;
}

void gaussian_distance_mixture(
        const std::vector<double>& axis, const std::vector<double>& means,
        const std::vector<double>& sigmas, const std::vector<double>& shapes,
        const std::vector<double>& amplitudes, int kernel,
        bool normalize_components, bool normalize,
        double** out_view, int* n_out_view) {
    internal::copy_to_view(
            gaussian_distance_mixture_impl(axis, means, sigmas, shapes,
                                           amplitudes, kernel,
                                           normalize_components, normalize),
            out_view, n_out_view);
}

IMPBFF_END_NAMESPACE
