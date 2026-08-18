/**
 * \file GreedyOlga.cpp
 * \brief Greedy Olga -- which FRET pair to measure next.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/GreedyOlga.h>

#include <cmath>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! Regularized *lower* incomplete gamma by its series. Converges for x < a+1.
double gamma_p_series(double a, double x) {
    if (x <= 0.0) return 0.0;
    double ap = a;
    double del = 1.0 / a;
    double sum = del;
    for (int i = 0; i < 1000; ++i) {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if (std::fabs(del) < std::fabs(sum) * 1e-16) break;
    }
    return sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
}

//! Regularized upper incomplete gamma by its continued fraction (Lentz).
/*! Converges for x >= a+1; the series above covers the rest. Taking only this
    branch, as a first cut of this file did, returns 1.0 for every x below a+1 --
    which at ndof = 201 is every chi-squared the algorithm ever sees. */
double gamma_q_continued_fraction(double a, double x) {
    const double tiny = 1e-300;
    const double eps = 1e-15;
    double b = x + 1.0 - a;
    double c = 1.0 / tiny;
    double d = 1.0 / b;
    double h = d;
    for (int i = 1; i < 400; ++i) {
        const double an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (std::fabs(d) < tiny) d = tiny;
        c = b + an / c;
        if (std::fabs(c) < tiny) c = tiny;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (std::fabs(del - 1.0) < eps) break;
    }
    return std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
}

//! Q(ndof/2, x) for one value, with the branch chosen by the parity of ndof.
double chi2_tail_one(double x, int ndof) {
    const double a = 0.5 * ndof;
    if (a > 100.0) {
        // The series would need ~a terms here. Reached only past 200 selected
        // pairs, so the general function is affordable exactly where the cheap
        // one stops being cheap.
        if (x <= 0.0) return 1.0;
        // Series below a+1, continued fraction above: neither converges over
        // the whole range, and the crossover is where both do.
        return x < a + 1.0 ? 1.0 - gamma_p_series(a, x)
                           : gamma_q_continued_fraction(a, x);
    }
    if (ndof % 2 == 0) {
        // a is an integer: Q(a, x) = exp(-x) sum_{n=0}^{a-1} x^n / n!
        double term = std::exp(-x);
        double total = term;
        const int n_terms = static_cast<int>(a);
        for (int n = 1; n < n_terms; ++n) {
            term = term * x / n;
            total += term;
        }
        return total;
    }
    // a is a half-integer:
    //   Q(a, x) = erfc(sqrt(x)) + exp(-x)/sqrt(pi x) * series
    double total = std::erfc(std::sqrt(x));
    if (a > 1.0 && x > 0.0) {
        // The loop bound is `n < a` with a half-integer a. Writing it as an
        // integer cast drops one term on every odd ndof -- the defect this
        // module's docstring records, which was very nearly a factor of two.
        double term = std::exp(-x) / std::sqrt(M_PI * x) * x / 0.5;
        double series = term;
        for (int n = 2; n < a; ++n) {
            term = term / (n - 0.5) * x;
            series += term;
        }
        total += series;
    }
    return total;
}

}  // namespace

std::vector<double> chi2_right_tail(const std::vector<double>& chisq, int ndof) {
    std::vector<double> out(chisq.size(), 0.0);
#pragma omp parallel for schedule(static)
    for (long long i = 0; i < static_cast<long long>(chisq.size()); ++i) {
        const std::size_t k = static_cast<std::size_t>(i);
        // Q(a, 0) = 1 exactly, and the whole diagonal is zero on every call.
        out[k] = chisq[k] <= 0.0 ? 1.0 : chi2_tail_one(0.5 * chisq[k], ndof);
    }
    return out;
}

std::vector<double> expected_rmsd_after_adding(
        const std::vector<double>& rmsds, const std::vector<double>& chi2,
        const std::vector<double>& e_add, double inv_err_sq, int ndof,
        double diag_weight, int n_frames, int n_candidates) {
    std::vector<double> out(static_cast<std::size_t>(std::max(0, n_candidates)), 0.0);
    if (n_frames <= 0 || n_candidates <= 0) return out;
    const std::size_t n = static_cast<std::size_t>(n_frames);

#pragma omp parallel for schedule(static)
    for (int k = 0; k < n_candidates; ++k) {
        const double* e = &e_add[static_cast<std::size_t>(k) * n];
        // Column accumulators only: the (candidates, n, n) array the vectorised
        // form builds is never created.
        std::vector<double> sum_w(n, 0.0), prod(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            const double ei = e[i];
            const double* chi2_row = &chi2[i * n];
            const double* rmsd_row = &rmsds[i * n];
            for (std::size_t j = 0; j < n; ++j) {
                const double de = ei - e[j];
                const double c = chi2_row[j] + de * de * inv_err_sq;
                const double w = c <= 0.0 ? 1.0 : chi2_tail_one(0.5 * c, ndof);
                sum_w[j] += w;
                prod[j] += w * rmsd_row[j];
            }
        }
        double acc = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            acc += prod[j] / (sum_w[j] - 1.0 + diag_weight);
        }
        out[k] = acc / static_cast<double>(n);
    }
    return out;
}

IMPBFF_END_NAMESPACE
