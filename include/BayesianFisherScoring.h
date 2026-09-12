/**
 * \file IMP/bff/BayesianFisherScoring.h
 * \brief The gradient and the information matrix of a Poisson count model.
 *
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_BAYESIANFISHERSCORING_H
#define IMPBFF_BAYESIANFISHERSCORING_H

#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! \name Fisher scoring for count models
//! @{

//! Which information matrix to build.
enum BayesianInformationKind {
  BAYESIAN_INFORMATION_EXPECTED = 0,  //!< Fisher scoring: `J' diag(w/m) J`
  BAYESIAN_INFORMATION_OBSERVED       //!< Newton-Raphson: `J' diag(w y / m^2) J`
};

/**
 * \brief `grad = J' w (y/m - 1)` and `A = J' W J` for a Poisson likelihood.
 *
 * For `log L = sum_b w_b [ y_b log m_b - m_b ]` with `m = m(theta)`, the
 * gradient is `J' w (y/m - 1)` and the two standard curvature matrices differ
 * only in `W`: the EXPECTED information puts `w/m` there and the OBSERVED one
 * `w y / m^2`. They agree in expectation because `E[y] = m`, and the choice
 * between them is Fisher scoring against Newton-Raphson (McCullagh & Nelder,
 * *Generalized Linear Models*, 2nd ed., 2.5).
 *
 * **Neither is an approximation to the objective.** The mode is where the
 * gradient vanishes either way; the choice only changes how fast the iteration
 * walks there, and which one is faster depends on how far from the mode you
 * are. Both drop the term `(y/m - 1) d2m/dtheta2`, which is what an exact
 * Newton step would restore -- and near the mode that term is small and
 * random-signed, which is why scoring converges at all.
 *
 * `w` is a per-bin weight, typically a 0/1 mask restricting the likelihood to
 * part of a histogram. Pass null for all ones.
 *
 * \param J `(n_bin, n_par)` row-major, `dm/dtheta`
 * \param grad `n_par` out
 * \param A `n_par x n_par` row-major out, the information matrix
 */
inline void bayesian_poisson_score(const double* y, const double* m, const double* J,
                          const double* w, std::size_t n_bin, std::size_t n_par,
                          BayesianInformationKind kind, double* grad, double* A,
                          double floor = 1e-12, std::size_t n_threads = 1) {
  std::vector<double> u(n_bin), Wd(n_bin);
  for (std::size_t b = 0; b < n_bin; ++b) {
    const double mb = m[b] > floor ? m[b] : floor;
    const double wb = w ? w[b] : 1.0;
    u[b] = wb * (y[b] / mb - 1.0);
    Wd[b] = (kind == BAYESIAN_INFORMATION_EXPECTED) ? wb / mb : wb * y[b] / (mb * mb);
  }
  for (std::size_t p = 0; p < n_par; ++p) {
    double s = 0.0;
    for (std::size_t b = 0; b < n_bin; ++b) s += J[b * n_par + p] * u[b];
    grad[p] = s;
  }
  //  `A = J' W J` accumulated as rank-one updates over the BINS rather than
  //  as a dot product per entry. Same arithmetic, same result -- but both
  //  reads then run along a contiguous row of `J`, where the obvious form
  //  strides by `n_par` down a column and misses cache on nearly every access.
  //  Measured at 114 parameters over 1808 bins: **19.1 ms the obvious way,
  //  1.94 ms this way**, and 2.06 ms for a blocked, threaded GEMM
  //  (`tttrlib::Mat`) doing the same product. The order is worth ten times
  //  what the library is, which is worth knowing before reaching for one.
  //  **And it partitions over the bins when asked.** Each thread owns a full
  //  `n_par x n_par` accumulator -- 104 kB at 114 parameters, which is nothing
  //  beside the Jacobian it is reading -- and the reduction is one pass at the
  //  end. Partitioning over PARAMETERS instead would have every thread stream
  //  the whole of `J`; partitioning over bins has each read a contiguous slice
  //  of it once. The serial path is kept exactly as it was, bit for bit, so a
  //  caller that does not ask for threads gets the same numbers as before.
  auto accumulate = [&](std::size_t b0, std::size_t b1, double* Aout) {
    for (std::size_t i = 0; i < n_par * n_par; ++i) Aout[i] = 0.0;
    for (std::size_t b = b0; b < b1; ++b) {
      const double* Jb = &J[b * n_par];
      const double wb = Wd[b];
      if (wb == 0.0) continue;
      for (std::size_t p = 0; p < n_par; ++p) {
        const double wj = wb * Jb[p];
        if (wj == 0.0) continue;
        double* Ap = &Aout[p * n_par];
        for (std::size_t q = p; q < n_par; ++q) Ap[q] += wj * Jb[q];
      }
    }
  };
  if (n_threads <= 1) {
    accumulate(0, n_bin, A);
  } else {
    std::vector<std::vector<double>> part(n_threads, std::vector<double>(n_par * n_par));
    std::vector<std::thread> th;
    for (std::size_t t = 0; t < n_threads; ++t)
      th.emplace_back(accumulate, n_bin * t / n_threads, n_bin * (t + 1) / n_threads,
                      part[t].data());
    for (auto& x : th) x.join();
    for (std::size_t i = 0; i < n_par * n_par; ++i) {
      double s = 0.0;
      for (std::size_t t = 0; t < n_threads; ++t) s += part[t][i];
      A[i] = s;
    }
  }
  for (std::size_t p = 0; p < n_par; ++p)
    for (std::size_t q = p + 1; q < n_par; ++q) A[q * n_par + p] = A[p * n_par + q];
}

/**
 * \brief `s log(1 + exp(x/s))`: a mean that stays positive without a kink.
 *
 * A Poisson mean must be positive -- the likelihood has `log m` in it -- but a
 * model built as a sum of terms can go slightly negative where there is no
 * signal. Clamping at zero puts a kink in the objective exactly where the
 * optimiser works, and a zero derivative beyond it, so the fit stops seeing the
 * data there. Bending instead costs nothing above a few `s` and never reaches
 * zero. Its derivative is `sigmoid(x/s)`, which the chain rule needs and which
 * is why the floor has to be smooth rather than clamped.
 */
inline double bayesian_soft_positive(double x, double s = 0.05) {
  const double z = x / s;
  return s * (z > 0.0 ? z + std::log1p(std::exp(-z)) : std::log1p(std::exp(z)));
}

//! `d bayesian_soft_positive / dx`.
inline double bayesian_soft_positive_derivative(double x, double s = 0.05) {
  const double z = x / s;
  return z > 0.0 ? 1.0 / (1.0 + std::exp(-z)) : std::exp(z) / (1.0 + std::exp(z));
}

//! @}

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_BAYESIANFISHERSCORING_H
