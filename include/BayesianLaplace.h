/**
 * \file IMP/bff/BayesianLaplace.h
 * \brief The Gaussian approximation at a posterior mode, its evidence, and a
 *        grid of them mixed by that evidence.
 *
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_BAYESIANLAPLACE_H
#define IMPBFF_BAYESIANLAPLACE_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! \name The Laplace approximation
//! @{

/**
 * \brief `log det A` for a symmetric positive-definite `A`, by Cholesky.
 *
 * Returns false when `A` is not positive definite, which at a mode means the
 * point is not one. Computing the determinant as a sum of logarithms rather
 * than the logarithm of a product is not fastidiousness: the determinant of a
 * hundred-dimensional information matrix overflows a double long before its
 * logarithm is large.
 */
inline bool bayesian_log_det_spd(const double* A, std::size_t n, double* out) {
  std::vector<double> L(n * n, 0.0);
  double s_log = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      double s = A[i * n + j];
      for (std::size_t k = 0; k < j; ++k) s -= L[i * n + k] * L[j * n + k];
      if (i == j) {
        if (!(s > 0.0) || !std::isfinite(s)) return false;
        L[i * n + j] = std::sqrt(s);
        s_log += std::log(L[i * n + j]);
      } else {
        L[i * n + j] = s / L[j * n + j];
      }
    }
  }
  *out = 2.0 * s_log;
  return true;
}

/**
 * \brief `log p(y)` from a mode and the curvature there.
 *
 * Expand `log p(y, theta)` to second order about its maximum and integrate the
 * Gaussian that results:
 *
 *     log p(y) ~ log p(y, theta_hat) + (d/2) log(2 pi) - (1/2) log det A
 *
 * with `A = -d2 log p / dtheta2` at the mode. This is the marginal likelihood
 * -- the number that compares models, or, here, the nodes of a hyperparameter
 * grid -- and the `-(1/2) log det A` is the whole content of Occam's razor in
 * this approximation: a model whose posterior is sharply peaked in many
 * directions pays for each of them.
 *
 * \param log_post_at_mode `log p(y, theta_hat)`, unnormalised in `theta`
 * \param A minus the curvature at the mode, `n x n` row-major
 */
inline bool bayesian_laplace_log_evidence(double log_post_at_mode, const double* A,
                                 std::size_t n, double* out) {
  double ld = 0.0;
  if (!bayesian_log_det_spd(A, n, &ld)) return false;
  *out = log_post_at_mode + 0.5 * double(n) * std::log(2.0 * M_PI) - 0.5 * ld;
  return true;
}

/**
 * \brief A hyperparameter grid, its nodes weighted by their evidence.
 *
 * **Why a grid and not a joint mode.** A smoothing parameter cannot simply be
 * maximised over jointly with what it smooths: the penalty's own normaliser
 * rewards a large weight whenever the coefficients are smooth, so the joint
 * maximum sits at the smoothest possible answer -- no structure at all. The
 * remedy is to treat it as what it is, a parameter to be integrated out:
 * approximate the posterior at each node of a grid, weight the nodes by their
 * evidence, and report the mixture. That is the integrated-nested-Laplace
 * construction (Rue, Martino & Chopin, *J. R. Statist. Soc. B* 71:319, 2009;
 * Wood, *Biometrika* 98:53, 2011, for smoothing parameters specifically), and
 * it needs no sampler: the weights are evidences and the components are
 * Gaussians.
 *
 * **What this class does with them.** For any scalar summary the caller can
 * evaluate at each node -- a distribution's mean, one of its bins, an
 * amplitude -- it combines the per-node mean and variance by the law of total
 * variance:
 *
 *     E[g]   = sum_i w_i m_i
 *     Var[g] = sum_i w_i (v_i + m_i^2) - E[g]^2
 *
 * The second term is what a single node cannot give: the spread BETWEEN nodes,
 * which is the uncertainty in the hyperparameter itself. Reporting the best
 * node's interval instead is the common mistake, and it is too narrow by
 * exactly that term.
 *
 * Weights are formed by subtracting the largest log evidence before
 * exponentiating, because these are log evidences of real data sets and
 * differences of tens of thousands of nats are ordinary.
 */
class BayesianEvidenceMixture {
 public:
  //! Add one node: its log evidence and a summary's mean and variance there.
  void add(double log_evidence, double mean, double variance) {
    lz_.push_back(log_evidence);
    m_.push_back(mean);
    v_.push_back(variance);
  }
  //! Add one node whose summary is a vector; all nodes must agree in length.
  void add_vector(double log_evidence, const std::vector<double>& mean,
                  const std::vector<double>& variance) {
    lz_.push_back(log_evidence);
    mv_.push_back(mean);
    vv_.push_back(variance);
  }

  std::size_t size() const { return lz_.size(); }

  //! Normalised weights, `exp(lz - max)` renormalised.
  std::vector<double> weights() const {
    std::vector<double> w(lz_.size(), 0.0);
    if (lz_.empty()) return w;
    double mx = -std::numeric_limits<double>::infinity();
    for (double z : lz_) if (std::isfinite(z) && z > mx) mx = z;
    double s = 0.0;
    for (std::size_t i = 0; i < lz_.size(); ++i) {
      w[i] = std::isfinite(lz_[i]) ? std::exp(lz_[i] - mx) : 0.0;
      s += w[i];
    }
    if (s > 0.0) for (double& x : w) x /= s;
    return w;
  }

  //! The mixture's mean and standard deviation of a scalar summary.
  void moments(double* mean, double* sd) const {
    const auto w = weights();
    double m = 0.0, s2 = 0.0;
    for (std::size_t i = 0; i < w.size(); ++i) m += w[i] * m_[i];
    for (std::size_t i = 0; i < w.size(); ++i) s2 += w[i] * (v_[i] + m_[i] * m_[i]);
    s2 -= m * m;
    *mean = m;
    *sd = std::sqrt(s2 > 0.0 ? s2 : 0.0);
  }

  //! The same, elementwise, for a vector summary.
  void moments_vector(std::vector<double>* mean, std::vector<double>* sd) const {
    const auto w = weights();
    const std::size_t n = mv_.empty() ? 0 : mv_[0].size();
    mean->assign(n, 0.0);
    sd->assign(n, 0.0);
    for (std::size_t j = 0; j < n; ++j) {
      double m = 0.0, s2 = 0.0;
      for (std::size_t i = 0; i < w.size(); ++i) m += w[i] * mv_[i][j];
      for (std::size_t i = 0; i < w.size(); ++i) s2 += w[i] * (vv_[i][j] + mv_[i][j] * mv_[i][j]);
      s2 -= m * m;
      (*mean)[j] = m;
      (*sd)[j] = std::sqrt(s2 > 0.0 ? s2 : 0.0);
    }
  }

 private:
  std::vector<double> lz_, m_, v_;
  std::vector<std::vector<double>> mv_, vv_;
};

//! @}

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_BAYESIANLAPLACE_H
