/**
 * \file IMP/bff/Transforms.h
 * \brief Bijections from a constrained parameter to the coordinate an
 *        optimiser or a sampler should actually move, with their Jacobians.
 *
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_TRANSFORMS_H
#define IMPBFF_TRANSFORMS_H

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! \name Parameter transforms
//! @{

/**
 * \brief Why a transform and not a bound.
 *
 * A positive quantity fitted with a lower bound at zero has a boundary the
 * optimiser can reach and the curvature is wrong there; a sampler proposes
 * across it and rejects. Fitting `z = log x` instead removes the boundary,
 * and the price is one term: a density stated on `x` becomes a density on `z`
 * only after multiplying by `|dx/dz|`, so `log p(z) = log p(x(z)) +
 * log|dx/dz|`. That term is `log_abs_det` below. Forgetting it does not make
 * a fit fail -- it silently changes the prior, which is why every transform
 * here carries it whether or not the caller happens to need it.
 *
 * Two of these carry no Jacobian for a reason rather than an oversight.
 * `ALR` is used with a density stated on the unconstrained coordinate already
 * (a logistic-normal), and `SumToZero` is an orthonormal map onto a subspace,
 * whose determinant on that subspace is one.
 */

//! `x = z`. The base case, so a caller can be generic.
template <typename T = double>
struct IdentityTransform {
  void to_constrained(const T* z, std::size_t n, T* x) const {
    for (std::size_t i = 0; i < n; ++i) x[i] = z[i];
  }
  void to_unconstrained(const T* x, std::size_t n, T* z) const {
    for (std::size_t i = 0; i < n; ++i) z[i] = x[i];
  }
  T log_abs_det(const T*, std::size_t) const { return T(0.0); }
};

//! Positive quantities: `x = exp(z)`, `log|dx/dz| = sum z`.
template <typename T = double>
struct LogTransform {
  void to_constrained(const T* z, std::size_t n, T* x) const {
    using std::exp;
    for (std::size_t i = 0; i < n; ++i) x[i] = exp(z[i]);
  }
  void to_unconstrained(const T* x, std::size_t n, T* z) const {
    using std::log;
    for (std::size_t i = 0; i < n; ++i) z[i] = log(x[i]);
  }
  T log_abs_det(const T* z, std::size_t n) const {
    T s = T(0.0);
    for (std::size_t i = 0; i < n; ++i) s += z[i];
    return s;
  }
};

/**
 * \brief Bounded quantities: `x = lo + (hi - lo) sigmoid(z)`.
 *
 * `log|dx/dz| = log(hi - lo) + log sigmoid(z) + log sigmoid(-z)`, written
 * through `logsigmoid` rather than as `log(s (1 - s))` because the latter
 * underflows to `-inf` for `|z|` beyond about 37 while the former stays
 * accurate: `log sigmoid(z) = -log1p(exp(-z))` for positive `z` and
 * `z - log1p(exp(z))` for negative.
 */
template <typename T = double>
struct LogitTransform {
  double lo = 0.0, hi = 1.0;
  LogitTransform() = default;
  LogitTransform(double lo_, double hi_) : lo(lo_), hi(hi_) {
    if (!(hi_ > lo_)) throw std::invalid_argument("LogitTransform: hi must exceed lo");
  }
  static T log_sigmoid(const T& z) {
    using std::exp; using std::log1p;
    return (z > T(0.0)) ? -log1p(exp(-z)) : z - log1p(exp(z));
  }
  void to_constrained(const T* z, std::size_t n, T* x) const {
    using std::exp;
    for (std::size_t i = 0; i < n; ++i) {
      const T s = (z[i] > T(0.0)) ? T(1.0) / (T(1.0) + exp(-z[i]))
                                  : exp(z[i]) / (T(1.0) + exp(z[i]));
      x[i] = T(lo) + T(hi - lo) * s;
    }
  }
  void to_unconstrained(const T* x, std::size_t n, T* z) const {
    using std::log;
    for (std::size_t i = 0; i < n; ++i) {
      T p = (x[i] - T(lo)) / T(hi - lo);
      if (p < T(1e-12)) p = T(1e-12);
      if (p > T(1.0 - 1e-12)) p = T(1.0 - 1e-12);
      z[i] = log(p) - log(T(1.0) - p);
    }
  }
  T log_abs_det(const T* z, std::size_t n) const {
    using std::log;
    T s = T(0.0);
    for (std::size_t i = 0; i < n; ++i)
      s += T(std::log(hi - lo)) + log_sigmoid(z[i]) + log_sigmoid(-z[i]);
    return s;
  }
};

/**
 * \brief A simplex of `n` weights from `n - 1` free numbers:
 *        `w = softmax([z, 0])` (the additive log-ratio transform).
 *
 * No Jacobian is added: this is used with a density stated on `z` itself -- a
 * logistic-normal -- so the density already lives in the coordinate. A caller
 * who instead states a density on the simplex must supply the Jacobian.
 */
template <typename T = double>
struct ALRTransform {
  std::size_t n = 0;                       //!< number of weights (z has n - 1)
  ALRTransform() = default;
  explicit ALRTransform(std::size_t n_) : n(n_) {
    if (n_ < 2) throw std::invalid_argument("ALRTransform: need at least two weights");
  }
  //! `x` has `n` entries, `z` has `n - 1`.
  void to_constrained(const T* z, std::size_t, T* x) const {
    using std::exp; using std::log;
    T m = T(0.0);                                   // max([z, 0]) for stability
    for (std::size_t i = 0; i + 1 < n; ++i) if (z[i] > m) m = z[i];
    T s = T(0.0);
    for (std::size_t i = 0; i + 1 < n; ++i) { x[i] = exp(z[i] - m); s += x[i]; }
    x[n - 1] = exp(-m); s += x[n - 1];
    for (std::size_t i = 0; i < n; ++i) x[i] /= s;
  }
  void to_unconstrained(const T* x, std::size_t, T* z) const {
    using std::log;
    for (std::size_t i = 0; i + 1 < n; ++i) z[i] = log(x[i]) - log(x[n - 1]);
  }
  T log_abs_det(const T*, std::size_t) const { return T(0.0); }
};

/**
 * \brief `n` coefficients summing to zero, from `n - 1` free numbers:
 *        `c = Q z` with `Q` an orthonormal basis of the sum-zero subspace.
 *
 * **The basis must be canonical, and this is not a detail.** The obvious
 * construction -- the left singular vectors of the centring projector
 * `I - 11'/n` -- is wrong for anything that persists a coordinate, because
 * that projector's eigenvalue-one subspace is `(n-1)`-fold degenerate and
 * LAPACK is free to return any rotation of it. It returns a different one on
 * macOS Accelerate than on Linux OpenBLAS, so a fit saved on one machine
 * restores on the other with the same numbers meaning a different function.
 * That happened in the project this came from (2026-09-07): a saved fit came
 * back with the right donor-only fraction and the wrong distribution.
 *
 * The Helmert contrasts are a fixed orthonormal basis of the same subspace,
 * computed identically everywhere with no eigen-decomposition: column `j` is
 * `(1, ..., 1, -j, 0, ..., 0) / sqrt(j (j + 1))`.
 *
 * `|det|` on the subspace is one, so `log_abs_det` is zero.
 */
template <typename T = double>
struct SumToZeroTransform {
  std::size_t n = 0;                       //!< length of the constrained vector
  std::vector<double> Q;                   //!< n x (n - 1), row-major
  SumToZeroTransform() = default;
  explicit SumToZeroTransform(std::size_t n_) : n(n_), Q(helmert_basis(n_)) {}

  //! The Helmert contrasts, `n x (n - 1)` row-major.
  static std::vector<double> helmert_basis(std::size_t n) {
    if (n < 2) throw std::invalid_argument("helmert_basis: need n >= 2");
    std::vector<double> q(n * (n - 1), 0.0);
    for (std::size_t j = 1; j < n; ++j) {
      const double s = std::sqrt(double(j) * (double(j) + 1.0));
      for (std::size_t i = 0; i < j; ++i) q[i * (n - 1) + (j - 1)] = 1.0 / s;
      q[j * (n - 1) + (j - 1)] = -double(j) / s;
    }
    return q;
  }
  //! `c` has `n` entries, `z` has `n - 1`.
  void to_constrained(const T* z, std::size_t, T* c) const {
    for (std::size_t i = 0; i < n; ++i) {
      T s = T(0.0);
      for (std::size_t j = 0; j + 1 < n; ++j) s += T(Q[i * (n - 1) + j]) * z[j];
      c[i] = s;
    }
  }
  void to_unconstrained(const T* c, std::size_t, T* z) const {
    for (std::size_t j = 0; j + 1 < n; ++j) {
      T s = T(0.0);
      for (std::size_t i = 0; i < n; ++i) s += T(Q[i * (n - 1) + j]) * c[i];
      z[j] = s;
    }
  }
  T log_abs_det(const T*, std::size_t) const { return T(0.0); }
};

//! @}

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_TRANSFORMS_H
