/**
 *  \file IMP/bff/Distributions.h
 *  \brief Probability distributions used by the dye and linker models.
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

//! distance_between_gaussian() as a bare buffer, for kernels that convolve.
/*! The worm-like-chain linker builds one of these per distance and reduces it
    immediately; publishing a managed view for each would allocate an ndarray
    per row of a convolution. */
IMPBFFEXPORT std::vector<double> distance_between_gaussian_impl(
        const std::vector<double>& distances, double separation_distance,
        double sigma);

//! The generalised normal as a bare buffer, for kernels that mix components.
IMPBFFEXPORT std::vector<double> generalized_normal_density_impl(
        const std::vector<double>& x, double loc, double scale,
        double shape, bool norm);

//! Which per-component density gaussian_distance_mixture() evaluates.
/*! Three, not two, and the third is easy to miss: with `norm = false` a
    generalised normal is **not** a normal. It evaluates the *standard* normal
    at `z = (x - loc)/scale`, so it omits the `1/scale` factor that
    normal_density() carries. Normalising each component hides the difference
    (the constant divides out); leaving them unnormalised does not, and then
    components of unequal width are weighted wrongly against each other. That
    is a 2e-2 error on a three-component mixture, found by measuring rather
    than by reading. */
enum GaussianMixtureKernel {
  //! `generalized_normal_distribution`: skewed, standard normal at z.
  GAUSSIAN_MIXTURE_GENERALIZED_NORMAL = 0,
  //! `distance_between_gaussian`: the distance between two Gaussian clouds.
  GAUSSIAN_MIXTURE_DISTANCE_BETWEEN_GAUSSIANS = 1,
  //! `normal_distribution`: the plain Gaussian density, `1/scale` included.
  GAUSSIAN_MIXTURE_NORMAL = 2
};

//! A weighted mixture of per-component distance distributions, in ONE call.
/*!
    The whole of a `Gaussians` distance distribution, so that a model with `k`
    components crosses the language boundary **once** rather than `k` times.
    The Python loop this replaces called back per component and summed the
    results in numpy, which is the half of the standing rule -- *the data stay
    where the computation is* -- that "it is already a C++ kernel" misses.

    The two booleans exist because the reference is **asymmetric between its
    own branches**, and hiding that in the kernel would be a silent behaviour
    change for one of the two callers:

    * the generalised-normal branch calls `generalized_normal_distribution`,
      whose `norm` defaults to **true**, so each component is normalised to
      unit sum before it is weighted;
    * the two-Gaussian branch calls `distance_between_gaussian`, whose
      `normalize` defaults to **false**, so each component is not.

    So the caller states which it wants rather than the kernel assuming.
    Normalising per component matters when a wide component falls partly off
    the axis: without it that component contributes less than its weight says.

    \param[in] axis the distance axis, shared by every component
    \param[in] means per-component centre; pass them already positive, as
        ChiSurf's getter does
    \param[in] sigmas per-component width
    \param[in] shapes per-component skew; empty or short is read as zero
    \param[in] amplitudes per-component weight
    \param[in] kernel which per-component density to evaluate, a
        GaussianMixtureKernel
    \param[in] normalize_components normalise each component before weighting
    \param[in] normalize normalise the combined density
    \param[out] out_view,n_out_view the mixture as a managed view
*/
IMPBFFEXPORT void gaussian_distance_mixture(
        const std::vector<double>& axis,
        const std::vector<double>& means,
        const std::vector<double>& sigmas,
        const std::vector<double>& shapes,
        const std::vector<double>& amplitudes,
        int kernel = GAUSSIAN_MIXTURE_GENERALIZED_NORMAL,
        bool normalize_components = true,
        bool normalize = true,
        double** out_view = 0, int* n_out_view = 0
);

//! gaussian_distance_mixture() as a bare buffer, for the node that owns one.
IMPBFFEXPORT std::vector<double> gaussian_distance_mixture_impl(
        const std::vector<double>& axis, const std::vector<double>& means,
        const std::vector<double>& sigmas, const std::vector<double>& shapes,
        const std::vector<double>& amplitudes, int kernel,
        bool normalize_components, bool normalize);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_DISTRIBUTIONS_H
