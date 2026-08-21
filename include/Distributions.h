/**
 *  \file IMP/bff/Distributions.h
 *  \brief Probability distributions used by the dye and linker models.
 *
 * Ported from Python by PRD-113: numba is a prototyping tool in this package,
 * not a runtime dependency, so every numerical kernel is C++.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_DISTRIBUTIONS_H
#define IMPBFF_DISTRIBUTIONS_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Poisson probabilities for k = 0 .. n-1.
/*!
    Uses the recursion \f$p_0 = e^{-\lambda},\; p_k = p_{k-1}\lambda/k\f$ rather
    than evaluating a factorial, which overflows well before the probabilities
    become negligible.

    \param[in] lam rate parameter
    \param[in] n number of terms
    \param[out] out_view,n_out_view the n probabilities as a managed view
*/
IMPBFFEXPORT void poisson_0toN(double lam, int n,
                              double** out_view, int* n_out_view);

//! Normal probability density on a given axis.
/*!
    \param[in] x the axis
    \param[in] loc mean
    \param[in] scale standard deviation; must be > 0
    \param[in] norm divide by the sum, so the discretised density sums to one
    \param[out] out_view,n_out_view the density as a managed view, one per x
*/
IMPBFFEXPORT void normal_distribution(
        const std::vector<double>& x,
        double loc = 0.0,
        double scale = 1.0,
        bool norm = false,
        double** out_view = 0, int* n_out_view = 0
);

//! Generalized normal density with a **skew** parameter.
/*!
    Not the exponential-power family: the shape parameter skews by transforming
    the axis, \f$z = -\log(1 - \kappa (x - \mu)/\sigma)/\kappa\f$, and the
    standard normal density is evaluated at \f$z\f$. ``shape = 0`` is the
    untransformed normal; positive skews left, negative right. Points where the
    transform argument would go non-positive are clamped to the smallest
    representable step rather than producing a NaN.

    \param[in] x the axis
    \param[in] loc location
    \param[in] scale scale
    \param[in] shape skewness; 0 gives a normal density
    \param[in] norm divide by the sum
    \param[out] out_view,n_out_view the density as a managed view, one per x
*/
IMPBFFEXPORT void generalized_normal_distribution(
        const std::vector<double>& x,
        double loc = 0.0,
        double scale = 1.0,
        double shape = 0.0,
        bool norm = true,
        double** out_view = 0, int* n_out_view = 0
);

//! Distance distribution between two isotropic 3-D Gaussians.
/*!
    For a non-zero separation \f$d\f$ this is
    \f$p(r) = (r/d)\,[\,N(r; d, \sigma) - N(r; -d, \sigma)\,]\f$; at \f$d = 0\f$
    it degenerates to the Maxwell form \f$2 r^2/\sigma^2 \cdot N(r; 0, \sigma)\f$,
    which the separate branch below exists to avoid dividing by.

    \param[in] distances the distance axis
    \param[in] separation_distance between the two means
    \param[in] sigma per-component width, shared by both Gaussians
    \param[in] normalize divide by the sum
    \param[out] out_view,n_out_view the distance distribution as a managed view
*/
IMPBFFEXPORT void distance_between_gaussian(
        const std::vector<double>& distances,
        double separation_distance,
        double sigma,
        bool normalize = false,
        double** out_view = 0, int* n_out_view = 0
);

//! The normal density at a mean and width, as a bare buffer (no view).
/*! Shared by the skew and two-Gaussian kernels, which build on the normal
    density without each allocating a second managed view. Not part of the
    public surface -- callers that handed the result back used to wrap it in
    `np.asarray`, and this helper keeps that to C++. */
IMPBFFEXPORT std::vector<double> normal_density(
        const std::vector<double>& x, double loc, double scale);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_DISTRIBUTIONS_H
