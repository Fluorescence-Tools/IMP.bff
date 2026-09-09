/**
 * \file IMP/bff/DeltaMethod.h
 * \brief The variance of a function of parameters, from the parameters' own.
 *
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_DELTAMETHOD_H
#define IMPBFF_DELTAMETHOD_H

#include <cmath>
#include <cstddef>

IMPBFF_BEGIN_NAMESPACE

//! \name The delta method
//! @{

/**
 * \brief `diag(J Sigma J')`: the variance of each component of a summary.
 *
 * A fit gives a covariance over its parameters, and what anyone wants is the
 * uncertainty of something else -- a mean distance, one bin of a distribution,
 * a ratio of amplitudes. Linearise the summary about the estimate and push the
 * covariance through it: `Var[g] ~ J Sigma J'` with `J = dg/dtheta`.
 *
 * **Two things it is easy to forget.** The approximation is only as good as the
 * linearisation, so for a summary that is strongly curved over the posterior's
 * width it understates the spread -- which is the case for a quantity near a
 * boundary, and the reason an interval formed this way should be clipped to
 * the quantity's support rather than allowed to cross it. And it is a variance,
 * not an interval: turning it into one assumes normality on whatever scale it
 * is applied, so a positive quantity's band is better formed in `p` and
 * clipped than formed in `log p` and exponentiated, unless the posterior really
 * is lognormal there.
 *
 * \param J `(m, n)` row-major, the summary's Jacobian
 * \param Sigma `(n, n)` row-major, the parameters' covariance
 * \param out `m` variances
 */
inline void delta_variance(const double* J, const double* Sigma,
                           std::size_t m, std::size_t n, double* out) {
  for (std::size_t i = 0; i < m; ++i) {
    const double* Ji = &J[i * n];
    double s = 0.0;
    for (std::size_t a = 0; a < n; ++a) {
      if (Ji[a] == 0.0) continue;
      double t = 0.0;
      for (std::size_t b = 0; b < n; ++b) t += Sigma[a * n + b] * Ji[b];
      s += Ji[a] * t;
    }
    out[i] = s > 0.0 ? s : 0.0;
  }
}

//! @}

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_DELTAMETHOD_H
