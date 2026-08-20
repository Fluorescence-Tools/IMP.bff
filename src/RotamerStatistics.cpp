/**
 * \file RotamerStatistics.cpp
 * \brief Weighted averaging over a rotamer ensemble.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/RotamerStatistics.h>
#include <IMP/bff/internal/OutputView.h>

#include <cmath>
#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

void rotamer_frame_weights(double* z_values, int n_frames, int n_pair,
                           double** out_view, int* n_out_view) {
    if (n_pair != 2) {
        throw std::invalid_argument("Z must have shape (n_frames, 2)");
    }
    double* out = internal::new_double_view(n_frames < 0 ? 0 : n_frames,
                                           out_view, n_out_view);
    if (out == nullptr) return;

    double total = 0.0;
    for (int i = 0; i < n_frames; ++i) {
        out[i] = z_values[2 * i] * z_values[2 * i + 1];
        total += out[i];
    }
    if (total == 0.0) {
        // uniform, not undefined: a frame where neither dye has an accessible
        // conformer says nothing about the others
        const double u = n_frames > 0 ? 1.0 / (double) n_frames : 0.0;
        for (int i = 0; i < n_frames; ++i) out[i] = u;
        return;
    }
    for (int i = 0; i < n_frames; ++i) out[i] /= total;
}

std::vector<double> weighted_average_sd_se(double* values, int n_values,
                                           double* weights, int n_weights) {
    if (n_values != n_weights) {
        throw std::invalid_argument("values and weights must be the same length");
    }
    std::vector<double> v, w;
    double total = 0.0;
    for (int i = 0; i < n_values; ++i) {
        // a frame where the dye could not be placed contributes nothing rather
        // than poisoning the mean
        if (!std::isfinite(values[i])) continue;
        v.push_back(values[i]);
        w.push_back(weights[i]);
        total += weights[i];
    }
    std::vector<double> out(3, std::nan(""));
    if (v.empty()) return out;

    double mean = 0.0;
    for (size_t i = 0; i < v.size(); ++i) mean += v[i] * (w[i] / total);
    double variance = 0.0;
    for (size_t i = 0; i < v.size(); ++i) {
        const double d = v[i] - mean;
        variance += d * d * (w[i] / total);
    }
    out[0] = mean;
    out[1] = std::sqrt(variance);
    // by the surviving frame count, not the effective count
    out[2] = std::sqrt(variance / (double) v.size());
    return out;
}

double effective_frame_fraction(double* weights, int n_weights) {
    std::vector<double> w;
    for (int i = 0; i < n_weights; ++i)
        if (weights[i] != 0.0) w.push_back(weights[i]);
    if (w.empty()) return 0.0;

    const double uniform = 1.0 / (double) w.size();
    double entropy = 0.0;
    for (size_t i = 0; i < w.size(); ++i)
        entropy -= w[i] * std::log(w[i] / uniform);
    return std::exp(entropy);
}

IMPBFF_END_NAMESPACE
