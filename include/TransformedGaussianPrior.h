/**
 * \file IMP/bff/TransformedGaussianPrior.h
 * \brief Many small transformed priors, recognised as one Gaussian in the
 *        coordinate the optimiser moves.
 *
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_TRANSFORMEDGAUSSIANPRIOR_H
#define IMPBFF_TRANSFORMEDGAUSSIANPRIOR_H

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <IMP/bff/Transforms.h>

IMPBFF_BEGIN_NAMESPACE

//! \name Transformed priors
//! @{

//! The prior families this recognises.
enum PriorFamily {
  PRIOR_GAUSSIAN = 0,     //!< on the constrained value
  PRIOR_LOGNORMAL,        //!< median and the sd of the natural log
  PRIOR_UNIFORM,          //!< flat between two bounds
  PRIOR_GAUSSIAN_ON_Z     //!< already stated on the unconstrained coordinate
};

//! The transforms this recognises.
enum TransformKind {
  TRANSFORM_IDENTITY = 0,
  TRANSFORM_LOG,
  TRANSFORM_LOGIT,
  TRANSFORM_OTHER         //!< anything else: not reducible, the caller keeps it
};

/**
 * \brief A model's small priors, collapsed into one quadratic form.
 *
 * **The observation.** A model with many nuisance parameters declares its
 * priors one at a time and in the units each parameter is thought about --
 * a lognormal on a positive scale factor, a Gaussian on an offset, a uniform
 * on a bounded efficiency. Evaluated that way the prior costs one small
 * object per parameter, and its gradient and Hessian cost one differentiation
 * each. In the project this came from that was a Python loop over 45
 * variables: 0.88 ms for the value and 50 ms for the Hessian, refreshed every
 * few steps of every fit, and none of it arithmetic -- all of it dispatch.
 *
 * **The reduction.** Most of those declarations are the SAME distribution once
 * you look at them in the coordinate the optimiser actually moves. The one
 * worth writing out is the common case:
 *
 *     x = exp(z),  x ~ LogNormal(median m, sd s)
 *     log p(z) = log p(x) + log|dx/dz|
 *              = [-log x - log s - ((log x - log m)/s)^2 / 2 - log sqrt(2 pi)] + z
 *              = -((z - log m)/s)^2 / 2 - log s - log sqrt(2 pi)
 *
 * The density's `-log x` and the change of variable's `+z` are the same term
 * with opposite signs. A lognormal on a positive quantity IS a Gaussian on its
 * logarithm -- which is the reason the logarithm was the right coordinate.
 * A Gaussian under no transform is trivially one, and a prior already stated
 * on `z` is one by declaration. Three families become a single mean vector and
 * a single standard deviation vector, and the whole prior is one quadratic
 * form with a diagonal Hessian that does not depend on the parameters at all.
 *
 * A uniform under a logit keeps only the transform's Jacobian, since the
 * density itself is flat. Its contribution is
 * `log sigmoid(z) + log sigmoid(-z)` per element, whose second derivative is
 * `-2 sigmoid(z) (1 - sigmoid(z))` -- so it is the one part that does depend
 * on the parameters.
 *
 * **What is not recognised is not silently mis-grouped.** `add` returns false
 * for a pair it cannot reduce, and the caller keeps evaluating that variable
 * itself. A new prior family costs speed, never correctness.
 */
template <typename T = double>
class TransformedGaussianPrior {
 public:
  /**
   * \brief Declare one variable's prior.
   *
   * \param offset index of the variable's first coordinate in the parameter vector
   * \param size number of coordinates
   * \param family the prior family
   * \param transform the transform between the constrained value and `z`
   * \param a family parameter: mean (Gaussian), median (lognormal), lower bound (uniform)
   * \param b family parameter: sd (Gaussian, lognormal), upper bound (uniform)
   * \return true if it was reduced; false if the caller must keep it
   */
  bool add(std::size_t offset, std::size_t size, PriorFamily family,
           TransformKind transform, double a, double b) {
    const bool gaussian_in_z =
        (family == PRIOR_GAUSSIAN_ON_Z) ||
        (family == PRIOR_GAUSSIAN && transform == TRANSFORM_IDENTITY) ||
        (family == PRIOR_LOGNORMAL && transform == TRANSFORM_LOG);
    if (gaussian_in_z) {
      if (!(b > 0.0)) throw std::invalid_argument("TransformedGaussianPrior: sd must be positive");
      const double mu = (family == PRIOR_LOGNORMAL) ? std::log(a) : a;
      for (std::size_t i = 0; i < size; ++i) {
        idx_.push_back(offset + i);
        mu_.push_back(mu);
        sd_.push_back(b);
      }
      return true;
    }
    if (family == PRIOR_UNIFORM && transform == TRANSFORM_LOGIT) {
      //  The flat density contributes -log(hi - lo) per element and the
      //  logit's Jacobian contributes +log(hi - lo) per element, so when the
      //  transform's bounds ARE the prior's bounds -- which is the only case
      //  this overload describes, since it takes one pair -- the two cancel
      //  exactly and what survives is the shape term alone. Adding either one
      //  without the other is the mistake this comment exists to prevent.
      if (!(b > a)) throw std::invalid_argument("TransformedGaussianPrior: need hi > lo");
      for (std::size_t i = 0; i < size; ++i) logit_.push_back(offset + i);
      return true;
    }
    return false;
  }

  //! Number of coordinates recognised as Gaussian in `z`.
  std::size_t n_gaussian() const { return idx_.size(); }
  //! Number of coordinates contributing only a logit Jacobian.
  std::size_t n_logit() const { return logit_.size(); }

  //! `log p(z)`, including the constants, so evidences are comparable.
  T log_prob(const T* z) const {
    T s = T(0.0);
    for (std::size_t k = 0; k < idx_.size(); ++k) {
      const T d = (z[idx_[k]] - T(mu_[k])) / T(sd_[k]);
      s -= T(0.5) * d * d;
    }
    for (std::size_t k = 0; k < logit_.size(); ++k) {
      const T& zz = z[logit_[k]];
      s += LogitTransform<T>::log_sigmoid(zz) + LogitTransform<T>::log_sigmoid(-zz);
    }
    return s + T(normaliser());
  }

  //! `d log p / dz`, added into `out` (which the caller has sized and zeroed).
  void add_gradient(const T* z, T* out) const {
    for (std::size_t k = 0; k < idx_.size(); ++k)
      out[idx_[k]] -= (z[idx_[k]] - T(mu_[k])) / T(sd_[k] * sd_[k]);
    for (std::size_t k = 0; k < logit_.size(); ++k) {
      const T s = sigmoid(z[logit_[k]]);
      out[logit_[k]] += T(1.0) - T(2.0) * s;
    }
  }

  /**
   * \brief `d2 log p / dz2`, added into the diagonal of `out` (`n x n`
   *        row-major).
   *
   * It is diagonal, and for everything except the logit variables it does not
   * depend on `z` -- so a caller refreshing curvature in a loop can compute
   * the Gaussian part once and never again.
   */
  void add_hessian(const T* z, std::size_t n, T* out) const {
    for (std::size_t k = 0; k < idx_.size(); ++k)
      out[idx_[k] * n + idx_[k]] -= T(1.0 / (sd_[k] * sd_[k]));
    for (std::size_t k = 0; k < logit_.size(); ++k) {
      const T s = sigmoid(z[logit_[k]]);
      out[logit_[k] * n + logit_[k]] -= T(2.0) * s * (T(1.0) - s);
    }
  }

 private:
  static T sigmoid(const T& z) {
    using std::exp;
    return (z > T(0.0)) ? T(1.0) / (T(1.0) + exp(-z)) : exp(z) / (T(1.0) + exp(z));
  }
  double normaliser() const {
    double c = const_;
    for (std::size_t k = 0; k < sd_.size(); ++k)
      c += -std::log(sd_[k]) - 0.5 * std::log(2.0 * M_PI);
    return c;
  }

  std::vector<std::size_t> idx_, logit_;
  std::vector<double> mu_, sd_;
  double const_ = 0.0;
};

//! @}

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_TRANSFORMEDGAUSSIANPRIOR_H
