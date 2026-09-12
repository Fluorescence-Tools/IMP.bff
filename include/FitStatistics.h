/**
 * \file IMP/bff/FitStatistics.h
 * \brief Whether a fit is good: the deviance a count model earns, and whether
 *        its residuals are structured.
 *
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_FITSTATISTICS_H
#define IMPBFF_FITSTATISTICS_H

#include <IMP/bff/IMPCompatibility.h>
#include <cmath>
#include <cstddef>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! \name Fit statistics
//! @{

/**
 * \brief Poisson deviance, `2 sum [ m - y + y log(y/m) ]`.
 *
 * **Why not chi-square.** For counting data the natural goodness-of-fit
 * statistic is the likelihood-ratio one, and dividing by an estimated variance
 * is a Gaussian approximation that fails where it matters -- in the tail,
 * where the counts are few and the long lifetimes live. Weighting by the
 * OBSERVED counts (`sigma = sqrt(y)`, Neyman) biases the fit low there;
 * weighting by the model (Pearson) does not, but is not the likelihood. The
 * deviance is the likelihood ratio against a model that fits every bin
 * exactly, so it needs no weights at all, and it is what "chi-square" should
 * mean for photon counting.
 *
 * A bin with `y = 0` contributes `2m`, which is the limit of `y log(y/m)` as
 * `y -> 0` and not a special case to be skipped.
 */
inline double poisson_deviance(const double* y, const double* m, std::size_t n) {
  double d = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double mi = (m[i] > 1e-300) ? m[i] : 1e-300;
    d += 2.0 * (mi - y[i]);
    if (y[i] > 0.0) d += 2.0 * y[i] * std::log(y[i] / mi);
  }
  return d;
}

//! Signed square roots of the per-bin deviance -- residuals whose sum of
//! squares IS the deviance, unlike `(y - m)/sqrt(y)`.
inline void deviance_residuals(const double* y, const double* m, std::size_t n, double* r) {
  for (std::size_t i = 0; i < n; ++i) {
    const double mi = (m[i] > 1e-300) ? m[i] : 1e-300;
    double d = 2.0 * (mi - y[i]);
    if (y[i] > 0.0) d += 2.0 * y[i] * std::log(y[i] / mi);
    if (d < 0.0) d = 0.0;
    r[i] = (y[i] >= mi ? 1.0 : -1.0) * std::sqrt(d);
  }
}

//! The outcome of a Wald-Wolfowitz runs test.
struct RunsTest {
  std::size_t n_runs = 0;      //!< runs of equal sign
  std::size_t n_above = 0;     //!< values above the cutoff
  std::size_t n_below = 0;     //!< values below it
  double z = 0.0;              //!< standardised deviation from the expected count
  double p_value = 1.0;        //!< two-sided
};

/**
 * \brief Are the residuals arranged in runs, or do their signs alternate like
 *        coin flips?
 *
 * **What it catches that a fit statistic does not.** A deviance per degree of
 * freedom near one says the residuals are the right SIZE. It says nothing
 * about their ORDER, and a model that is systematically low over one stretch
 * of the axis and high over another can have an excellent deviance. The runs
 * test asks whether the number of sign changes is what independent noise would
 * give: too few runs means the residuals are correlated along the axis, which
 * is what a missing component or a wrong instrument response looks like.
 *
 * Under the null the number of runs `R` among `n1` values above and `n2` below
 * has mean `2 n1 n2 / (n1 + n2) + 1` and variance
 * `2 n1 n2 (2 n1 n2 - n1 - n2) / ((n1+n2)^2 (n1+n2-1))`, and `R` is
 * asymptotically normal (Wald & Wolfowitz, *Ann. Math. Statist.* 11:147,
 * 1940).
 *
 * Values exactly at the cutoff are counted as above, which is the convention
 * `statsmodels.sandbox.stats.runs.runstest_1samp` uses; pass the residuals
 * themselves with `cutoff = 0` for the usual sign test.
 */
inline RunsTest runs_test(const double* x, std::size_t n, double cutoff = 0.0,
                          bool continuity_correction = false) {
  RunsTest r;
  if (n == 0) return r;
  std::vector<bool> above(n);
  for (std::size_t i = 0; i < n; ++i) {
    above[i] = (x[i] >= cutoff);
    if (above[i]) ++r.n_above; else ++r.n_below;
  }
  r.n_runs = 1;
  for (std::size_t i = 1; i < n; ++i) if (above[i] != above[i - 1]) ++r.n_runs;
  const double n1 = double(r.n_above), n2 = double(r.n_below), N = n1 + n2;
  if (n1 == 0.0 || n2 == 0.0 || N < 2.0) { r.z = 0.0; r.p_value = 1.0; return r; }
  const double mean = 2.0 * n1 * n2 / N + 1.0;
  const double var = 2.0 * n1 * n2 * (2.0 * n1 * n2 - N) / (N * N * (N - 1.0));
  if (!(var > 0.0)) { r.z = 0.0; r.p_value = 1.0; return r; }
  double diff = double(r.n_runs) - mean;
  if (continuity_correction) diff -= (diff > 0.0 ? 0.5 : -0.5);
  r.z = diff / std::sqrt(var);
  r.p_value = std::erfc(std::fabs(r.z) / std::sqrt(2.0));
  return r;
}

//! @}

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_FITSTATISTICS_H
