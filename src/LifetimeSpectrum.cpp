/**
 * \file LifetimeSpectrum.cpp
 * \brief The experiment-neutral output: (amplitude, rate constant) pairs.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/LifetimeSpectrum.h>

#include <algorithm>
#include <cmath>

IMPBFF_BEGIN_NAMESPACE

std::vector<double> lifetime_spectrum_decay(
        const std::vector<double>& amplitudes,
        const std::vector<double>& rate_constants,
        const std::vector<double>& time) {
    std::vector<double> out(time.size(), 0.0);
    const std::size_t n = std::min(amplitudes.size(), rate_constants.size());
#pragma omp parallel for schedule(static)
    for (long long j = 0; j < static_cast<long long>(time.size()); ++j) {
        const double t = time[static_cast<std::size_t>(j)];
        double f = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            f += amplitudes[i] * std::exp(-rate_constants[i] * t);
        }
        out[static_cast<std::size_t>(j)] = f;
    }
    return out;
}

std::vector<double> lifetime_spectrum_coarse_grain(
        const std::vector<double>& amplitudes,
        const std::vector<double>& rate_constants,
        int n_bins, std::vector<double>& out_rates) {
    out_rates.clear();
    std::vector<double> out_amplitudes;
    const std::size_t n = std::min(amplitudes.size(), rate_constants.size());
    if (n == 0 || n_bins < 1) return out_amplitudes;

    double lo = rate_constants[0], hi = rate_constants[0];
    for (std::size_t i = 1; i < n; ++i) {
        lo = std::min(lo, rate_constants[i]);
        hi = std::max(hi, rate_constants[i]);
    }
    // Every species shares one rate: there is nothing to reduce, and a zero
    // width would divide by zero below.
    if (!(hi > lo)) {
        double a_sum = 0.0;
        for (std::size_t i = 0; i < n; ++i) a_sum += amplitudes[i];
        out_rates.push_back(lo);
        out_amplitudes.push_back(a_sum);
        return out_amplitudes;
    }

    const double width = (hi - lo) / n_bins;
    std::vector<double> a_sum(n_bins, 0.0), ak_sum(n_bins, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t b = static_cast<std::size_t>((rate_constants[i] - lo) / width);
        if (b >= static_cast<std::size_t>(n_bins)) b = n_bins - 1;  // closed right edge
        a_sum[b] += amplitudes[i];
        ak_sum[b] += amplitudes[i] * rate_constants[i];
    }
    for (int b = 0; b < n_bins; ++b) {
        if (a_sum[b] == 0.0) continue;   // an empty bin is not a species
        out_amplitudes.push_back(a_sum[b]);
        // The amplitude-weighted mean, so sum(a) and sum(a*k) both survive.
        out_rates.push_back(ak_sum[b] / a_sum[b]);
    }
    return out_amplitudes;
}

std::vector<double> outer_product_histogram(
        const std::vector<double>& a, const std::vector<double>& weights_a,
        const std::vector<double>& b, const std::vector<double>& weights_b,
        int n_bins, double lo, double hi) {
    std::vector<double> hist(n_bins > 0 ? n_bins : 0, 0.0);
    if (n_bins <= 0 || !(hi > lo)) return hist;
    const std::size_t na = std::min(a.size(), weights_a.size());
    const std::size_t nb = std::min(b.size(), weights_b.size());
    const double inv_width = n_bins / (hi - lo);
    for (std::size_t i = 0; i < na; ++i) {
        const double ai = a[i], wi = weights_a[i];
        for (std::size_t j = 0; j < nb; ++j) {
            const double v = ai * b[j];
            if (v < lo || v > hi) continue;
            std::size_t k = static_cast<std::size_t>((v - lo) * inv_width);
            if (k >= hist.size()) k = hist.size() - 1;   // the closed top edge
            hist[k] += wi * weights_b[j];
        }
    }
    return hist;
}

IMPBFF_END_NAMESPACE
