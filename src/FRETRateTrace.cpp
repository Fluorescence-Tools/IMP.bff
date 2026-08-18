/**
 * \file FRETRateTrace.cpp
 * \brief FRET rate along a dye trajectory.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/FRETRateTrace.h>

#include <cmath>

IMPBFF_BEGIN_NAMESPACE

std::vector<double> fret_rate_trace_kernel(
        const std::vector<double>& trajectory,
        const std::vector<double>& acceptor_points, double R0, double tau0,
        double r_min2, double kappa2_scale) {
    const std::size_t n_frames = trajectory.size() / 3;
    const std::size_t n_acceptor = acceptor_points.size() / 3;
    std::vector<double> rates(n_frames, 0.0);
    if (n_acceptor == 0 || n_frames == 0) return rates;
    const double R0_6 = std::pow(R0, 6) * kappa2_scale;
    const double inv_tau0 = 1.0 / tau0;
    for (std::size_t f = 0; f < n_frames; ++f) {
        const double px = trajectory[3 * f + 0];
        const double py = trajectory[3 * f + 1];
        const double pz = trajectory[3 * f + 2];
        double total = 0.0;
        for (std::size_t j = 0; j < n_acceptor; ++j) {
            const double dx = px - acceptor_points[3 * j + 0];
            const double dy = py - acceptor_points[3 * j + 1];
            const double dz = pz - acceptor_points[3 * j + 2];
            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < r_min2) r2 = r_min2;   // 1/r^6 diverges; dyes cannot overlap
            total += R0_6 / (r2 * r2 * r2);
        }
        // Arithmetic mean of rates: the fast-exchange limit.
        rates[f] = inv_tau0 * total / static_cast<double>(n_acceptor);
    }
    return rates;
}

std::vector<double> fret_rate_pair_trace_kernel(
        const std::vector<double>& donor, const std::vector<double>& acceptor,
        double R0, double tau0, double r_min2, double kappa2_scale) {
    const std::size_t n_frames = donor.size() / 3;
    std::vector<double> rates(n_frames, 0.0);
    const double R0_6 = std::pow(R0, 6) * kappa2_scale;
    const double inv_tau0 = 1.0 / tau0;
    for (std::size_t f = 0; f < n_frames; ++f) {
        const double dx = donor[3 * f + 0] - acceptor[3 * f + 0];
        const double dy = donor[3 * f + 1] - acceptor[3 * f + 1];
        const double dz = donor[3 * f + 2] - acceptor[3 * f + 2];
        double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 < r_min2) r2 = r_min2;
        rates[f] = inv_tau0 * R0_6 / (r2 * r2 * r2);
    }
    return rates;
}

IMPBFF_END_NAMESPACE
