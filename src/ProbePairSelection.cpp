/**
 * \file ProbePairSelection.cpp
 * \brief Greedy probe selection -- which probe pair to measure next.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/ProbePairSelection.h>
#include <IMP/bff/internal/OutputView.h>

#include <IMP/bff/IMPCompatibility.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

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

//! Expected mean RMSD per candidate, written into `out`.
void score_candidates(const double* rmsds, const double* chi2,
                      const double* e_add, double inv_err_sq, int ndof,
                      double diag_weight, std::size_t n, int n_candidates,
                      double* out) {
#pragma omp parallel for schedule(static)
    for (int k = 0; k < n_candidates; ++k) {
        const double* e = e_add + static_cast<std::size_t>(k) * n;
        // Column accumulators only: the (candidates, n, n) array the vectorised
        // form builds is never created.
        std::vector<double> sum_w(n, 0.0), prod(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            const double ei = e[i];
            const double* chi2_row = chi2 + i * n;
            const double* rmsd_row = rmsds + i * n;
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
}

//! Mean over columns of the weight-averaged RMSD, weights being Q(ndof/2, chi2/2).
double weighted_column_mean(const double* rmsds, const double* chi2, int ndof,
                            double diag_weight, std::size_t n) {
    std::vector<double> sum_w(n, 0.0), prod(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double* chi2_row = chi2 + i * n;
        const double* rmsd_row = rmsds + i * n;
        for (std::size_t j = 0; j < n; ++j) {
            const double c = chi2_row[j];
            const double w = c <= 0.0 ? 1.0 : chi2_tail_one(0.5 * c, ndof);
            sum_w[j] += w;
            prod[j] += w * rmsd_row[j];
        }
    }
    double acc = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
        acc += prod[j] / (sum_w[j] - 1.0 + diag_weight);
    }
    return acc / static_cast<double>(n);
}

//! chi2 += (e_i - e_j)^2 / sigma^2 for the pair just selected.
void add_pair_to_chi2(std::vector<double>& chi2, const double* e,
                      double inv_err_sq, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) {
        const double ei = e[i];
        double* row = &chi2[i * n];
        for (std::size_t j = 0; j < n; ++j) {
            const double de = ei - e[j];
            row[j] += de * de * inv_err_sq;
        }
    }
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
    score_candidates(rmsds.data(), chi2.data(), e_add.data(), inv_err_sq, ndof,
                     diag_weight, static_cast<std::size_t>(n_frames),
                     n_candidates, out.data());
    return out;
}

double expected_rmsd(const std::vector<double>& rmsds,
                     const std::vector<double>& chi2, int ndof,
                     double diag_weight, int n_frames) {
    if (n_frames <= 0) return 0.0;
    return weighted_column_mean(rmsds.data(), chi2.data(), ndof, diag_weight,
                                static_cast<std::size_t>(n_frames));
}

void select_probe_pairs(
        double* predicted_measurements, int n_frames, int n_pairs,
        double* rmsds, int n_rmsd_rows, int n_rmsd_cols,
        double measurement_error, int max_pairs, bool unique_only, double diag_weight,
        int** out_pairs, int* n_out_pairs,
        double** out_decay, int* n_out_decay) {
    if (n_rmsd_rows != n_rmsd_cols || n_rmsd_rows != n_frames) {
        IMP_THROW("rmsds must be square and match the frame count of predicted_measurements",
                  ValueException);
    }
    const int n_take = std::min(max_pairs, n_pairs);
    std::vector<int> pairs;
    std::vector<double> decay;
    if (n_frames > 0 && n_pairs > 0 && n_take > 0) {
        const std::size_t n = static_cast<std::size_t>(n_frames);
        const std::size_t m = static_cast<std::size_t>(n_pairs);
        const double inv_err_sq = 1.0 / (measurement_error * measurement_error);

        // Candidate-major, transposed once: the scorer reads one candidate's
        // frames contiguously, and every step would otherwise stride through
        // `predicted_measurements`.
        std::vector<double> e_t(m * n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t k = 0; k < m; ++k) e_t[k * n + i] = predicted_measurements[i * m + k];
        }

        std::vector<double> chi2(n * n, 0.0);
        std::vector<double> scores(m, 0.0);
        std::vector<char> taken(m, 0);
        pairs.reserve(static_cast<std::size_t>(n_take));

        for (int step = 0; step < n_take; ++step) {
            // ndof is clamped at the first two steps: zero degrees of freedom
            // is not a chi-squared distribution.
            const int ndof = std::max(step - 1, 1);
            score_candidates(rmsds, chi2.data(), e_t.data(), inv_err_sq, ndof,
                             diag_weight, n, n_pairs, scores.data());
            int best = -1;
            for (std::size_t k = 0; k < m; ++k) {
                if (unique_only && taken[k]) continue;
                if (best < 0 ||
                    scores[k] < scores[static_cast<std::size_t>(best)]) {
                    best = static_cast<int>(k);
                }
            }
            if (best < 0) break;   // unique_only, and every candidate is spent
            pairs.push_back(best);
            taken[static_cast<std::size_t>(best)] = 1;
            add_pair_to_chi2(chi2, &e_t[static_cast<std::size_t>(best) * n],
                             inv_err_sq, n);
        }

        // The decay is a second pass from a zeroed chi-squared: it reports
        // what each step left behind, and its `ndof` counts the pairs added
        // so far rather than the selector's clamped one.
        std::vector<double> decay_chi2(n * n, 0.0);
        decay.reserve(pairs.size());
        for (std::size_t s = 0; s < pairs.size(); ++s) {
            add_pair_to_chi2(decay_chi2,
                             &e_t[static_cast<std::size_t>(pairs[s]) * n],
                             inv_err_sq, n);
            decay.push_back(weighted_column_mean(
                    rmsds, decay_chi2.data(), static_cast<int>(s) + 1,
                    diag_weight, n));
        }
    }

    // Publish both views; the kernel's vectors are copied out exactly once.
    // A C++ caller may take the header's NULL defaults -- the allocators
    // decline to publish then (scratch), and the dereferences below must
    // not happen either.
    const std::size_t n_sel = pairs.size();
    if (out_pairs == nullptr || n_out_pairs == nullptr ||
        out_decay == nullptr || n_out_decay == nullptr) {
        int* s1 = internal::new_int_view(n_sel, nullptr, nullptr);
        double* s2 = internal::new_double_view(n_sel, nullptr, nullptr);
        std::free(s1);
        std::free(s2);
        return;
    }
    internal::new_int_view(n_sel, out_pairs, n_out_pairs);
    if (*out_pairs != nullptr && !pairs.empty()) {
        std::memcpy(*out_pairs, pairs.data(), n_sel * sizeof(int));
    }
    internal::new_double_view(n_sel, out_decay, n_out_decay);
    if (*out_decay != nullptr && !decay.empty()) {
        std::memcpy(*out_decay, decay.data(), n_sel * sizeof(double));
    }
}

void select_probe_positions(
        double* predicted_measurements, int n_frames, int n_pairs,
        double* rmsds, int n_rmsd_rows, int n_rmsd_cols,
        int* pair_sites, int n_site_rows, int n_site_cols,
        double measurement_error, int max_sites, double diag_weight,
        int** out_sites, int* n_out_sites,
        double** out_decay, int* n_out_decay) {
    if (n_rmsd_rows != n_rmsd_cols || n_rmsd_rows != n_frames) {
        IMP_THROW("rmsds must be square and match the frame count of predicted_measurements",
                  ValueException);
    }
    if (n_site_rows != n_pairs) {
        IMP_THROW("pair_sites must have one row per pair", ValueException);
    }
    if (n_site_cols != 2) {
        IMP_THROW("pair_sites must hold two site indices per pair",
                  ValueException);
    }
    for (int p = 0; p < n_pairs; ++p) {
        if (pair_sites[2 * p] < 0 || pair_sites[2 * p + 1] < 0) {
            IMP_THROW("pair_sites entries must be site indices >= 0",
                      ValueException);
        }
    }

    // Candidate sites are the ones at least one pair connects; a site no pair
    // mentions can never add a measurement and is not a candidate.
    int n_sites = 0;
    for (int p = 0; p < n_pairs; ++p) {
        n_sites = std::max(n_sites, pair_sites[2 * p] + 1);
        n_sites = std::max(n_sites, pair_sites[2 * p + 1] + 1);
    }
    const int n_take = std::min(max_sites, n_sites);
    std::vector<int> sites;
    std::vector<double> decay;
    if (n_frames > 0 && n_pairs > 0 && n_take > 0) {
        const std::size_t n = static_cast<std::size_t>(n_frames);
        const std::size_t m = static_cast<std::size_t>(n_pairs);
        const double inv_err_sq = 1.0 / (measurement_error * measurement_error);

        // Candidate-major predictions, as in the pair selector: the scorer
        // reads one pair's frames contiguously.
        std::vector<double> e_t(m * n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t k = 0; k < m; ++k) e_t[k * n + i] = predicted_measurements[i * m + k];
        }

        // Pairs per site, and the running bookkeeping: a pair is implied once
        // both of its sites are chosen.
        std::vector<std::vector<int>> pairs_of(static_cast<std::size_t>(n_sites));
        std::vector<char> appears(static_cast<std::size_t>(n_sites), 0);
        for (std::size_t p = 0; p < m; ++p) {
            pairs_of[static_cast<std::size_t>(pair_sites[2 * p])].push_back(
                    static_cast<int>(p));
            appears[static_cast<std::size_t>(pair_sites[2 * p])] = 1;
            if (pair_sites[2 * p + 1] != pair_sites[2 * p]) {
                pairs_of[static_cast<std::size_t>(pair_sites[2 * p + 1])]
                        .push_back(static_cast<int>(p));
            }
            appears[static_cast<std::size_t>(pair_sites[2 * p + 1])] = 1;
        }
        std::vector<char> chosen(static_cast<std::size_t>(n_sites), 0);
        std::vector<char> active(m, 0);
        std::vector<double> chi2(n * n, 0.0);
        std::vector<double> scores(static_cast<std::size_t>(n_sites), 0.0);
        sites.reserve(static_cast<std::size_t>(n_take));

        int n_measurements = 0;
        for (int step = 0; step < n_take; ++step) {
            // Score every site that is not yet chosen by what the pairs it
            // would activate leave behind. The scratch chi-squared is per
            // candidate, so the loop parallelises without reductions.
#pragma omp parallel for schedule(static)
            for (int i = 0; i < n_sites; ++i) {
                const std::size_t si = static_cast<std::size_t>(i);
                if (chosen[si] || !appears[si]) {
                    continue;   // a residue cannot be mutated twice
                }
                std::vector<double> scratch(chi2);
                int n_new = 0;
                for (int p : pairs_of[si]) {
                    const std::size_t sp = static_cast<std::size_t>(p);
                    const int other = pair_sites[2 * sp] == i
                                      ? pair_sites[2 * sp + 1] : pair_sites[2 * sp];
                    if (active[sp] || (other != i && !chosen[static_cast<std::size_t>(other)])) {
                        continue;
                    }
                    add_pair_to_chi2(scratch, &e_t[sp * n], inv_err_sq, n);
                    ++n_new;
                }
                // The pair selector scores with ndof clamped two below the
                // measurement count it is heading for; the site selector holds
                // the same rule with the count summed over the implied pairs.
                scores[si] = weighted_column_mean(
                        rmsds, scratch.data(),
                        std::max(n_measurements + n_new - 2, 1), diag_weight, n);
            }
            int best = -1;
            for (int i = 0; i < n_sites; ++i) {
                if (chosen[static_cast<std::size_t>(i)] ||
                    !appears[static_cast<std::size_t>(i)]) {
                    continue;
                }
                if (best < 0 ||
                    scores[static_cast<std::size_t>(i)] <
                            scores[static_cast<std::size_t>(best)]) {
                    best = i;
                }
            }
            if (best < 0) break;   // every candidate site is chosen
            const std::size_t sb = static_cast<std::size_t>(best);
            chosen[sb] = 1;
            sites.push_back(best);
            for (int p : pairs_of[sb]) {
                const std::size_t sp = static_cast<std::size_t>(p);
                const int other = pair_sites[2 * sp] == best
                                  ? pair_sites[2 * sp + 1] : pair_sites[2 * sp];
                if (active[sp] || (other != best &&
                                   !chosen[static_cast<std::size_t>(other)])) {
                    continue;
                }
                add_pair_to_chi2(chi2, &e_t[sp * n], inv_err_sq, n);
                active[sp] = 1;
                ++n_measurements;
            }
        }

        // The decay reports what each site left behind, from a second pass
        // over a zeroed chi-squared -- the same shape as the pair selector's
        // decay, so its ndof counts the measurements made so far rather than
        // the selector's clamped one. Walking the sites in selection order, a
        // pair joins when its second site arrives, exactly as during
        // selection; the first site of a set without homotypic pairs measures
        // nothing yet, and its ndof is clamped at one.
        std::fill(chi2.begin(), chi2.end(), 0.0);
        std::fill(active.begin(), active.end(), 0);
        std::fill(chosen.begin(), chosen.end(), 0);
        int n_so_far = 0;
        decay.reserve(sites.size());
        for (int s : sites) {
            chosen[static_cast<std::size_t>(s)] = 1;
            for (int p : pairs_of[static_cast<std::size_t>(s)]) {
                const std::size_t sp = static_cast<std::size_t>(p);
                const int other = pair_sites[2 * sp] == s
                                  ? pair_sites[2 * sp + 1] : pair_sites[2 * sp];
                if (active[sp] || (other != s &&
                                   !chosen[static_cast<std::size_t>(other)])) {
                    continue;
                }
                add_pair_to_chi2(chi2, &e_t[sp * n], inv_err_sq, n);
                active[sp] = 1;
                ++n_so_far;
            }
            decay.push_back(weighted_column_mean(
                    rmsds, chi2.data(), std::max(n_so_far, 1), diag_weight, n));
        }
    }

    const std::size_t n_sel = sites.size();
    if (out_sites == nullptr || n_out_sites == nullptr ||
        out_decay == nullptr || n_out_decay == nullptr) {
        int* s1 = internal::new_int_view(n_sel, nullptr, nullptr);
        double* s2 = internal::new_double_view(n_sel, nullptr, nullptr);
        std::free(s1);
        std::free(s2);
        return;
    }
    internal::new_int_view(n_sel, out_sites, n_out_sites);
    if (*out_sites != nullptr && !sites.empty()) {
        std::memcpy(*out_sites, sites.data(), n_sel * sizeof(int));
    }
    internal::new_double_view(n_sel, out_decay, n_out_decay);
    if (*out_decay != nullptr && !decay.empty()) {
        std::memcpy(*out_decay, decay.data(), n_sel * sizeof(double));
    }
}

IMPBFF_END_NAMESPACE
