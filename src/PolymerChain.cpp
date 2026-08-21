/**
 * \file PolymerChain.cpp
 * \brief End-to-end distance distributions of ideal and worm-like chains.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/PolymerChain.h>
#include <IMP/bff/Distributions.h>
#include <IMP/bff/internal/Normalize.h>
#include <IMP/bff/internal/OutputView.h>

#include <algorithm>
#include <cmath>
#include <cstring>

IMPBFF_BEGIN_NAMESPACE

//! The ideal-chain kernel as a buffer, for the view-publishing wrapper.
std::vector<double> gaussian_chain_impl(const std::vector<double>& distances,
                                        double segment_length,
                                        int number_of_segments) {
    const double ree = gaussian_chain_ree(segment_length, number_of_segments);
    const double r2_mean = ree * ree;
    std::vector<double> out(distances.size(), 0.0);
    if (r2_mean == 0.0) return out;
    const double denom = std::pow(2.0 / 3.0 * M_PI * r2_mean, 1.5);
    for (std::size_t i = 0; i < distances.size(); ++i) {
        const double r = distances[i];
        out[i] = 4.0 * M_PI * r * r / denom * std::exp(-1.5 * r * r / r2_mean);
    }
    return out;
}

//! The worm-like-chain kernel as a buffer (the wrapper publishes a view).
std::vector<double> worm_like_chain_impl(const std::vector<double>& distances,
                                         double kappa, double chain_length,
                                         bool normalize, bool distance) {
    std::vector<double> pr(distances.size(), 0.0);
    if (distances.empty()) return pr;
    if (chain_length == 0.0) {
        chain_length = *std::max_element(distances.begin(), distances.end());
    }
    if (chain_length == 0.0) return pr;

    const double a = 14.054;
    const double b = 0.473;
    const double c =
            1.0 - std::pow(1.0 + std::pow(0.38 * std::pow(kappa, -0.95), -5.0), -0.2);
    // The branch is the paper's: below kappa = 0.125 the correction is linear,
    // above it the fitted form takes over. Written with `exp(0.783*log(x))`
    // rather than `pow(x, 0.783)` to match the Python term for term.
    const double d = (kappa < 0.125)
            ? kappa + 1.0
            : 1.0 - 1.0 / (0.177 / (kappa - 0.111)
                           + 6.4 * std::exp(0.783 * std::log(kappa - 0.111)));

    for (std::size_t i = 0; i < distances.size(); ++i) {
        const double r = distances[i];
        if (r >= chain_length) continue;   // the chain cannot exceed its contour
        const double r_n = r / chain_length;
        const double r_n2 = r_n * r_n;
        double pri = std::pow((1.0 - c * r_n2) / (1.0 - r_n2), 2.5);
        pri *= std::exp(-d * kappa * a * b * (1.0 + b)
                        / (1.0 - (b * r_n) * (b * r_n)) * r_n2);
        const double r_n4 = r_n2 * r_n2;
        const double r_n6 = r_n4 * r_n2;
        const double g = ((-0.75) / kappa - 0.5) * r_n2
                       + ((-0.359375) / kappa + 1.0625) * r_n4
                       + ((-0.109375) / kappa - 0.5625) * r_n6;
        pri *= std::exp(g / (1.0 - r_n2));
        pri *= std::exp(-d * kappa * a * (1.0 + b) * r_n
                        / (1.0 - (b * r_n) * (b * r_n)));
        pr[i] = pri;
    }
    if (distance) {
        for (std::size_t i = 0; i < pr.size(); ++i) pr[i] *= distances[i] * distances[i];
    }
    if (normalize) internal::normalize_sum(pr);
    return pr;
}

double gaussian_chain_ree(double segment_length, int number_of_segments) {
    return segment_length * std::sqrt(static_cast<double>(number_of_segments));
}

void gaussian_chain(const std::vector<double>& distances,
                    double segment_length, int number_of_segments,
                    double** out_view, int* n_out_view) {
    internal::copy_to_view(
            gaussian_chain_impl(distances, segment_length, number_of_segments),
            out_view, n_out_view);
}

void worm_like_chain(const std::vector<double>& distances, double kappa,
                     double chain_length, bool normalize, bool distance,
                     double** out_view, int* n_out_view) {
    internal::copy_to_view(
            worm_like_chain_impl(distances, kappa, chain_length, normalize, distance),
            out_view, n_out_view);
}

void worm_like_chain_linker(const std::vector<double>& distances, double kappa,
                            double chain_length, double sigma, bool normalize,
                            double** out_view, int* n_out_view) {
    std::vector<double> pn(distances.size(), 0.0);
    if (sigma != 0.0) {
        const std::vector<double> pr =
                worm_like_chain_impl(distances, kappa, chain_length, normalize, true);
        for (std::size_t i = 0; i < distances.size(); ++i) {
            if (pr[i] == 0.0) continue;
            const std::vector<double> broad =
                    normal_density(distances, distances[i], sigma);
            for (std::size_t j = 0; j < pn.size(); ++j) pn[j] += pr[i] * broad[j];
        }
        if (normalize) internal::normalize_sum(pn);
    }
    internal::copy_to_view(pn, out_view, n_out_view);
}

IMPBFF_END_NAMESPACE
