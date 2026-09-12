/**
 * \file FRETExchange.cpp
 * \brief FRET when the labels exchange between states during the lifetime.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/FRETExchange.h>

#include <IMP/bff/IMPCompatibility.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! The stationary distribution: the eigenvector of P^T at eigenvalue 1.
Eigen::VectorXd stationary_distribution(const Eigen::MatrixXd& p) {
    Eigen::EigenSolver<Eigen::MatrixXd> solver(p.transpose());
    const Eigen::VectorXcd values = solver.eigenvalues();
    int best = 0;
    double best_distance = std::abs(values[0] - std::complex<double>(1.0, 0.0));
    for (int i = 1; i < values.size(); ++i) {
        const double d = std::abs(values[i] - std::complex<double>(1.0, 0.0));
        if (d < best_distance) { best_distance = d; best = i; }
    }
    Eigen::VectorXd w = solver.eigenvectors().col(best).real();
    const double total = w.sum();
    // A degenerate eigenvector can come back summing to zero, or negated; the
    // sign is arbitrary and the scale is the normalisation, so both are fixed
    // here rather than propagating as a negative population.
    if (total != 0.0) w /= total;
    return w;
}

}  // namespace

double fret_efficiency_exact_kinetic(const std::vector<double>& p_matrix,
                                     const std::vector<double>& fret_rates,
                                     double tau0, double dt,
                                     const std::vector<double>& weights) {
    const int n = static_cast<int>(fret_rates.size());
    if (n == 0) return 0.0;
    if (p_matrix.size() != static_cast<std::size_t>(n) * n) {
        IMP_THROW("the transition matrix must be (n, n) for n = " << n
                          << " states, not " << p_matrix.size() << " values",
                  ValueException);
    }
    if (!(tau0 > 0.0)) {
        IMP_THROW("tau0 must be > 0, not " << tau0, ValueException);
    }
    if (!(dt > 0.0)) {
        IMP_THROW("dt must be > 0, not " << dt, ValueException);
    }

    Eigen::MatrixXd p(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) p(i, j) = p_matrix[i * n + j];
    }

    // P = exp(M dt) ~ I + M dt for a small step, so M = (P - I)/dt.
    const Eigen::MatrixXd m =
            (p - Eigen::MatrixXd::Identity(n, n)) / dt;

    Eigen::VectorXd w(n);
    if (weights.empty()) {
        w = stationary_distribution(p);
    } else {
        if (weights.size() != static_cast<std::size_t>(n)) {
            IMP_THROW("one weight per state: " << weights.size() << " against "
                                               << n,
                      ValueException);
        }
        for (int i = 0; i < n; ++i) w[i] = weights[i];
    }

    const double k_rad = 1.0 / tau0;
    Eigen::MatrixXd a = -m;
    for (int i = 0; i < n; ++i) a(i, i) += k_rad + fret_rates[i];

    // The **transposed** system. `p_matrix[i][j]` is the probability of i -> j,
    // so populations evolve as a row vector and G solves A^T G = w. Solving
    // A G = w is right only for a symmetric M, and for a real transition matrix
    // it gives the wrong fast-exchange limit.
    const Eigen::VectorXd g = a.transpose().fullPivLu().solve(w);

    double efficiency = 0.0;
    for (int i = 0; i < n; ++i) efficiency += fret_rates[i] * g[i];
    return efficiency;
}

double fret_efficiency_exact_kinetic_pair(
        const std::vector<double>& dist_matrix,
        const std::vector<double>& kappa2_matrix,
        const std::vector<double>& p_d, const std::vector<double>& p_a,
        const std::vector<double>& weights_d,
        const std::vector<double>& weights_a, double forster_radius,
        double tau0, double dt) {
    const int nd = static_cast<int>(weights_d.size());
    const int na = static_cast<int>(weights_a.size());
    if (nd == 0 || na == 0) return 0.0;
    if (dist_matrix.size() != static_cast<std::size_t>(nd) * na ||
        kappa2_matrix.size() != static_cast<std::size_t>(nd) * na) {
        IMP_THROW("the distance and kappa2 matrices must be (" << nd << ", "
                                                              << na << ")",
                  ValueException);
    }

    // The joint kinetics is the Kronecker product, which is (nd*na)^2 -- the
    // reason this is not the default for a large library.
    const int n = nd * na;
    std::vector<double> p_total(static_cast<std::size_t>(n) * n, 0.0);
    for (int i = 0; i < nd; ++i) {
        for (int j = 0; j < na; ++j) {
            for (int k = 0; k < nd; ++k) {
                for (int l = 0; l < na; ++l) {
                    p_total[static_cast<std::size_t>(i * na + j) * n + k * na + l] =
                            p_d[i * nd + k] * p_a[j * na + l];
                }
            }
        }
    }

    const double k_rad = 1.0 / tau0;
    std::vector<double> w(n), rates(n);
    for (int i = 0; i < nd; ++i) {
        for (int j = 0; j < na; ++j) {
            const std::size_t at = static_cast<std::size_t>(i) * na + j;
            w[at] = weights_d[i] * weights_a[j];
            const double x = forster_radius / dist_matrix[at];
            const double x3 = x * x * x;
            rates[at] = k_rad * x3 * x3 * 1.5 * kappa2_matrix[at];
        }
    }
    return fret_efficiency_exact_kinetic(p_total, rates, tau0, dt, w);
}

FRETRegimes fret_efficiency_regimes(const std::vector<double>& dist_matrix,
                                    const std::vector<double>& kappa2_matrix,
                                    const std::vector<double>& weights_d,
                                    const std::vector<double>& weights_a,
                                    double forster_radius) {
    FRETRegimes out;
    const int nd = static_cast<int>(weights_d.size());
    const int na = static_cast<int>(weights_a.size());
    if (nd == 0 || na == 0) return out;
    if (dist_matrix.size() != static_cast<std::size_t>(nd) * na ||
        kappa2_matrix.size() != static_cast<std::size_t>(nd) * na) {
        IMP_THROW("the distance and kappa2 matrices must be (" << nd << ", "
                                                              << na << ")",
                  ValueException);
    }

    double e_static = 0.0, k2_avg = 0.0, rate_avg = 0.0;
    for (int i = 0; i < nd; ++i) {
        for (int j = 0; j < na; ++j) {
            const std::size_t at = static_cast<std::size_t>(i) * na + j;
            const double w = weights_d[i] * weights_a[j];
            const double x = forster_radius / dist_matrix[at];
            const double x3 = x * x * x;
            // 1.5 kappa2 normalises kappa2 against the isotropic 2/3.
            const double ratio = x3 * x3 * 1.5 * kappa2_matrix[at];
            e_static += w * ratio / (1.0 + ratio);
            k2_avg += w * kappa2_matrix[at];
            rate_avg += w * ratio;
        }
    }

    double e_dynamic = 0.0;
    for (int i = 0; i < nd; ++i) {
        for (int j = 0; j < na; ++j) {
            const std::size_t at = static_cast<std::size_t>(i) * na + j;
            const double w = weights_d[i] * weights_a[j];
            const double x = forster_radius / dist_matrix[at];
            const double x3 = x * x * x;
            const double ratio = x3 * x3 * 1.5 * k2_avg;
            e_dynamic += w * ratio / (1.0 + ratio);
        }
    }

    out.static_efficiency = e_static;
    out.dynamic = e_dynamic;
    out.dynamic_plus = rate_avg / (1.0 + rate_avg);
    out.kappa2_avg = k2_avg;
    return out;
}

IMPBFF_END_NAMESPACE
