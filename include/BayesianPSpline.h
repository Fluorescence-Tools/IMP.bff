/**
 * \file IMP/bff/BayesianPSpline.h
 * \brief The Bayesian P-spline prior: a difference penalty on spline
 *        coefficients, read as a density.
 *
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_BAYESIANPSPLINE_H
#define IMPBFF_BAYESIANPSPLINE_H

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! \name Bayesian P-splines
//! @{

/**
 * \brief `p(c | lambda)`: a Gaussian (or Student-t) on the d-th differences of
 *        spline coefficients, with the penalty weight as a parameter of the
 *        density rather than a knob.
 *
 * **The idea.** Fit a smooth function with many more basis functions than the
 * data can support, and control the wiggliness with a penalty on the
 * differences of neighbouring coefficients instead of by choosing knots
 * (Eilers & Marx, *Statist. Sci.* 11:89, 1996). Written as a density rather
 * than as a penalty it is a prior, and its weight `lambda` becomes a parameter
 * that can be inferred (Lang & Brezger, *J. Comput. Graph. Statist.* 13:183,
 * 2004).
 *
 * **What the ORDER means, which is the part worth understanding.** The
 * difference matrix `D` of order `d` annihilates polynomials of degree
 * `d - 1`. So the order is a statement about the null hypothesis the prior
 * shrinks toward, not a smoothness dial:
 *
 *   - `d = 1` shrinks toward a constant,
 *   - `d = 2` toward a straight line (the usual choice),
 *   - `d = 3` toward a parabola -- which, for a log-density, means one
 *     Gaussian population is free at any penalty weight.
 *
 * **The normaliser matters when lambda is inferred.** `D` has a null space, so
 * `exp(-lambda/2 ||D c||^2)` is improper along it and its normalising constant
 * is `(rank/2) log lambda` rather than `(n/2) log lambda`. Get that wrong and
 * the marginal posterior of `lambda` is wrong -- it is the term that stops the
 * evidence from preferring an arbitrarily large penalty. `rank` is `n - d`
 * here: `D` has `n - d` rows and full row rank, and it already annihilates the
 * constant, so projecting onto the sum-to-zero subspace removes nothing
 * further.
 *
 * **The null space is given a weak proper prior** rather than left improper:
 * an N(0, `tilt_sd`^2) on the linear direction, and under order 3 an
 * N(0, `quad_sd`^2) on the quadratic one. Both directions are unit vectors
 * orthogonal to the constant, so this is a statement about the function's tilt
 * and curvature in units of the coefficient, not a ridge on everything.
 *
 * **The Student-t family** replaces the Gaussian on each difference with a
 * t of `nu` degrees of freedom and scale `1/sqrt(lambda)` -- the locally
 * adaptive P-spline with its local variances integrated out
 * (Baladandayuthapani, Mallick & Carroll, *Biometrics* 61:64, 2005). A few
 * large jumps then cost little while the rest stays smooth, which is what a
 * distribution with well separated populations needs: under the Gaussian
 * family the gap between two populations is expensive and the fit would rather
 * remove a population than pay for it.
 *
 * **Templated on the scalar** so the same code evaluates at `double` and at a
 * forward-mode dual number (`tttrlib::Dual<GradVec<N>>`), which makes the
 * closed-form gradient and Hessian below checkable against the very expression
 * they differentiate, in C++, with nothing else in the loop.
 */
template <typename T = double>
class BayesianPSplinePrior {
 public:
  /**
   * \param n number of spline coefficients
   * \param order difference order `d`; 1, 2 or 3
   * \param tilt_sd prior sd of the linear direction of the null space
   * \param quad_sd prior sd of the quadratic direction (order 3 only)
   * \param student Student-t rather than Gaussian on the differences
   * \param nu degrees of freedom of the Student-t
   */
  BayesianPSplinePrior(std::size_t n, int order = 2, double tilt_sd = 3.0,
               double quad_sd = 30.0, bool student = false, double nu = 3.0)
      : n_(n), order_(order), tilt_sd_(tilt_sd), quad_sd_(quad_sd),
        student_(student), nu_(nu) {
    if (order < 1 || order > 3) throw std::invalid_argument("BayesianPSplinePrior: order must be 1, 2 or 3");
    if (n <= static_cast<std::size_t>(order)) throw std::invalid_argument("BayesianPSplinePrior: n must exceed the order");
    build_difference();
    build_null_space();
    rank_ = static_cast<int>(n_) - order_;
  }

  std::size_t size() const { return n_; }
  int order() const { return order_; }
  //! Rank of the penalty, and so the power of `lambda` in the normaliser.
  int rank() const { return rank_; }
  //! The difference matrix, `(n - order)` rows of `n`, row-major.
  const std::vector<double>& difference_matrix() const { return D_; }
  //! Unit vector along the linear direction of the null space.
  const std::vector<double>& tilt_direction() const { return v_tilt_; }
  //! Unit vector along the quadratic direction, orthogonal to the linear one.
  const std::vector<double>& quadratic_direction() const { return v_quad_; }

  //! `log p(c | lambda)` up to a constant independent of `c` and `lambda`.
  T log_prob(const T* c, const T& lambda) const {
    std::vector<T> d(n_rows());
    differences(c, d.data());
    T rough = T(0.0);
    if (student_) {
      for (std::size_t k = 0; k < n_rows(); ++k)
        rough -= T(0.5) * T(nu_ + 1.0) * log1p_(lambda * d[k] * d[k] / T(nu_));
    } else {
      T q = T(0.0);
      for (std::size_t k = 0; k < n_rows(); ++k) q += d[k] * d[k];
      rough = T(-0.5) * lambda * q;
    }
    const T tilt = dot(v_tilt_.data(), c);
    T lp = T(0.5 * rank_) * log_(lambda) + rough - T(0.5) * (tilt / T(tilt_sd_)) * (tilt / T(tilt_sd_));
    if (order_ == 3) {
      const T qd = dot(v_quad_.data(), c);
      lp -= T(0.5) * (qd / T(quad_sd_)) * (qd / T(quad_sd_));
    }
    return lp;
  }

  //! `d log p / dc`, in closed form. `out` has `n` entries.
  void gradient(const T* c, const T& lambda, T* out) const {
    std::vector<T> d(n_rows());
    differences(c, d.data());
    std::vector<T> w(n_rows());
    for (std::size_t k = 0; k < n_rows(); ++k) {
      if (student_) {
        //  d/dx of -0.5 (nu+1) log1p(lambda x^2 / nu)
        const T den = T(1.0) + lambda * d[k] * d[k] / T(nu_);
        w[k] = -T(nu_ + 1.0) * lambda * d[k] / T(nu_) / den;
      } else {
        w[k] = -lambda * d[k];
      }
    }
    for (std::size_t i = 0; i < n_; ++i) out[i] = T(0.0);
    for (std::size_t k = 0; k < n_rows(); ++k) {
      const double* row = &D_[k * n_];
      for (std::size_t i = 0; i < n_; ++i)
        if (row[i] != 0.0) out[i] += w[k] * T(row[i]);
    }
    const T tilt = dot(v_tilt_.data(), c) / T(tilt_sd_ * tilt_sd_);
    for (std::size_t i = 0; i < n_; ++i) out[i] -= tilt * T(v_tilt_[i]);
    if (order_ == 3) {
      const T qd = dot(v_quad_.data(), c) / T(quad_sd_ * quad_sd_);
      for (std::size_t i = 0; i < n_; ++i) out[i] -= qd * T(v_quad_[i]);
    }
  }

  /**
   * \brief `d2 log p / dc2`, `n * n` row-major.
   *
   * Under the Gaussian family this does not depend on `c` at all -- the
   * density is a quadratic form -- so a caller that refreshes curvature every
   * few steps can hoist it out of the loop entirely. Under the Student-t it
   * does depend on `c`, which is the price of local adaptivity.
   */
  void hessian(const T* c, const T& lambda, T* out) const {
    for (std::size_t i = 0; i < n_ * n_; ++i) out[i] = T(0.0);
    std::vector<T> d(n_rows());
    if (student_) differences(c, d.data());
    for (std::size_t k = 0; k < n_rows(); ++k) {
      T s;
      if (student_) {
        const T u = lambda * d[k] * d[k] / T(nu_);
        const T den = T(1.0) + u;
        //  d2/dx2 of -0.5 (nu+1) log1p(lambda x^2/nu)
        s = -T(nu_ + 1.0) * lambda / T(nu_) * (T(1.0) - u) / (den * den);
      } else {
        s = -lambda;
      }
      const double* row = &D_[k * n_];
      for (std::size_t i = 0; i < n_; ++i) {
        if (row[i] == 0.0) continue;
        for (std::size_t j = 0; j < n_; ++j)
          if (row[j] != 0.0) out[i * n_ + j] += s * T(row[i] * row[j]);
      }
    }
    add_outer(out, v_tilt_, 1.0 / (tilt_sd_ * tilt_sd_));
    if (order_ == 3) add_outer(out, v_quad_, 1.0 / (quad_sd_ * quad_sd_));
  }

  //! The d-th differences of `c`; `out` has `n - order` entries.
  void differences(const T* c, T* out) const {
    for (std::size_t k = 0; k < n_rows(); ++k) {
      const double* row = &D_[k * n_];
      T s = T(0.0);
      for (std::size_t i = 0; i < n_; ++i)
        if (row[i] != 0.0) s += T(row[i]) * c[i];
      out[k] = s;
    }
  }

  std::size_t n_rows() const { return n_ - static_cast<std::size_t>(order_); }

 private:
  static T log_(const T& x) { using std::log; return log(x); }
  static T log1p_(const T& x) { using std::log1p; return log1p(x); }

  T dot(const double* v, const T* c) const {
    T s = T(0.0);
    for (std::size_t i = 0; i < n_; ++i) s += T(v[i]) * c[i];
    return s;
  }

  void add_outer(T* out, const std::vector<double>& v, double scale) const {
    for (std::size_t i = 0; i < n_; ++i)
      for (std::size_t j = 0; j < n_; ++j) out[i * n_ + j] -= T(scale * v[i] * v[j]);
  }

  //! `D = diff(I_n, order)`: rows of binomial coefficients with alternating signs.
  void build_difference() {
    D_.assign(n_rows() * n_, 0.0);
    std::vector<double> coef(static_cast<std::size_t>(order_) + 1, 0.0);
    double b = 1.0;
    for (int j = 0; j <= order_; ++j) {
      coef[static_cast<std::size_t>(j)] = ((j % 2) ? -b : b);
      b = b * double(order_ - j) / double(j + 1);
    }
    //  numpy's diff puts the newest sample first: row k is
    //  sum_j (-1)^j C(d, j) c[k + d - j]
    for (std::size_t k = 0; k < n_rows(); ++k)
      for (int j = 0; j <= order_; ++j)
        D_[k * n_ + k + static_cast<std::size_t>(order_ - j)] = coef[static_cast<std::size_t>(j)];
  }

  //! Unit vectors along the linear and quadratic directions, both orthogonal
  //! to the constant, and to each other.
  void build_null_space() {
    std::vector<double> v(n_);
    const double mid = 0.5 * (double(n_) - 1.0);
    double nv = 0.0;
    for (std::size_t i = 0; i < n_; ++i) { v[i] = double(i) - mid; nv += v[i] * v[i]; }
    nv = std::sqrt(nv);
    v_tilt_.resize(n_);
    for (std::size_t i = 0; i < n_; ++i) v_tilt_[i] = v[i] / nv;

    std::vector<double> q(n_);
    double mean = 0.0;
    for (std::size_t i = 0; i < n_; ++i) { q[i] = v[i] * v[i]; mean += q[i]; }
    mean /= double(n_);
    double proj = 0.0;
    for (std::size_t i = 0; i < n_; ++i) { q[i] -= mean; proj += q[i] * v_tilt_[i]; }
    double nq = 0.0;
    for (std::size_t i = 0; i < n_; ++i) { q[i] -= proj * v_tilt_[i]; nq += q[i] * q[i]; }
    nq = std::sqrt(nq);
    v_quad_.resize(n_);
    for (std::size_t i = 0; i < n_; ++i) v_quad_[i] = (nq > 0.0 ? q[i] / nq : 0.0);
  }

  std::size_t n_;
  int order_;
  double tilt_sd_, quad_sd_;
  bool student_;
  double nu_;
  int rank_;
  std::vector<double> D_, v_tilt_, v_quad_;
};

//! @}

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_BAYESIANPSPLINE_H
