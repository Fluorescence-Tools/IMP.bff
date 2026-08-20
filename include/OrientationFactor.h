/**
 *  \file IMP/bff/OrientationFactor.h
 *  \brief The FRET orientation factor kappa^2, from geometry and from wobbling.
 *
 * Two questions that are usually confused. **Geometry**: two transition dipoles
 * and a separation vector give one number. **Distribution**: two dyes that
 * reorient during the donor's excited-state lifetime give a *distribution* of
 * numbers, whose width is what makes kappa^2 the dominant systematic in a FRET
 * distance -- and whose shape is set by the order parameters a time-resolved
 * anisotropy measures.
 *
 * The wobbling-in-a-cone form is eq. 9 of Sindbert et al., JACS 133, 2463
 * (2011).
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_ORIENTATIONFACTOR_H
#define IMPBFF_ORIENTATIONFACTOR_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! \f$\kappa^2\f$ for a wobbling-in-a-cone pair at given angles.
/*!
    \param[in] delta angle between the dyes' symmetry axes, radians
    \param[in] sD2,sA2 second-rank order parameters of donor and acceptor
    \param[in] beta1,beta2 angle of each symmetry axis to \f$R_{DA}\f$, radians
    \return the orientation factor; 2/3 when both order parameters vanish
*/
IMPBFFEXPORT double wobbling_kappa2(
        double delta, double sD2, double sA2, double beta1, double beta2);

//! Separation and \f$\kappa\f$ for two dipoles given by their end points.
/*!
    \param[in] d1,d2 the donor dipole's two ends, three coordinates each
    \param[in] a1,a2 the acceptor dipole's two ends
    \return two values: the centre-to-centre distance and \f$\kappa\f$ (not
            squared -- it is signed, and the sign is lost by squaring)
*/
IMPBFFEXPORT std::vector<double> dipole_kappa_distance(
        const std::vector<double>& d1, const std::vector<double>& d2,
        const std::vector<double>& a1, const std::vector<double>& a2);

//! \f$p(\kappa^2)\f$ for a known angle between the dyes' symmetry axes.
/*!
    Sweeps \f$\beta_1\f$ over \f$(0, \pi/2)\f$ and \f$\phi\f$ over
    \f$(0, 2\pi)\f$ on a regular grid, weighting each \f$\beta_1\f$ by
    \f$\sin\beta_1\f$ -- the solid-angle element, without which the poles are
    over-counted. Deterministic; no sampling.

    \param[in] delta angle between the symmetry axes, radians
    \param[in] sD2,sA2 second-rank order parameters
    \param[in] step angular step, degrees
    \param[in] n_bins histogram edges returned; the histogram has one fewer
    \param[in] k2_min,k2_max histogram range
    \param[out] k2_scale the \p n_bins bin edges
    \param[out] k2_hist the solid-angle-weighted histogram
    \return every computed \f$\kappa^2\f$, row-major over (beta1, phi)
*/
IMPBFFEXPORT std::vector<double> wobbling_kappa2_distribution_delta(
        double delta, double sD2, double sA2, double step,
        int n_bins, double k2_min, double k2_max,
        std::vector<double>& k2_scale, std::vector<double>& k2_hist);

//! \f$p(\kappa^2)\f$ over isotropically oriented dyes.
/*!
    Dipole directions are drawn **uniformly on the sphere**, by normalising
    three standard normals. The Python this replaces normalised three
    *uniform* variates instead, which fills only the positive octant of the
    cube and does so non-uniformly: it returned \f$\langle\kappa^2\rangle =
    0.333\f$ where the rigid isotropic limit is \f$2/3\f$, a factor of two.
    The function had no callers, so nothing downstream carried the error.

    \param[in] sD2,sA2 second-rank order parameters
    \param[in] n_bins histogram edges returned; the histogram has one fewer
    \param[in] k2_min,k2_max histogram range
    \param[in] n_samples orientation pairs to draw
    \param[in] seed for reproducibility
    \param[out] k2_scale the \p n_bins bin edges
    \param[out] k2_hist the histogram
    \return every sampled \f$\kappa^2\f$
*/
IMPBFFEXPORT std::vector<double> wobbling_kappa2_distribution(
        double sD2, double sA2, int n_bins, double k2_min, double k2_max,
        int n_samples, int seed,
        std::vector<double>& k2_scale, std::vector<double>& k2_hist);


//! \f$p(\kappa^2)\f$ for **diffusion with traps**, given a measured efficiency.
/*!
    Each dye is either freely diffusing -- rotating fast enough to average its
    orientation, contributing \f$\kappa^2 = 2/3\f$ -- or *trapped*, wobbling in
    a cone set by its order parameter. \p sD2 and \p sA2 are the trapped
    fractions, so the pair splits into four sub-populations: free/free,
    trapped/trapped, trapped/free and free/trapped. Their efficiencies are
    combined and the result inverted into the single \f$\kappa^2\f$ that would
    have produced it.

    That inversion is what separates this from wobbling_kappa2_distribution():
    it does not ask what \f$\kappa^2\f$ the geometry gives, but what
    \f$\kappa^2\f$ the geometry *and the measured efficiency* together imply.

    Per sample: two isotropic directions are drawn, \f$R_{DA}\f$ is taken along
    x, and the four sub-populations (donor free or trapped, acceptor free or
    trapped) are combined into one efficiency before being inverted back into a
    single \f$\kappa^2\f$.

    Drawing the directions with three standard normals is load-bearing: the
    normalised components of a *uniform* draw fill the cube's positive octant,
    not the sphere, and that halved \f$\langle\kappa^2\rangle\f$ to 0.333 in the
    isotropic limit where it must be exactly 2/3.

    \param[in] sD2,sA2 second-rank order parameters of donor and acceptor
    \param[in] fret_efficiency the measured efficiency
    \param[in] n_samples orientation pairs to draw
    \param[in] n_bins bins between \p k2_min and \p k2_max
    \param[in] k2_min,k2_max histogram range
    \param[in] seed reproducible sampling; negative draws freely
    \param[out] out_view,n_out_view the bin edges (\p n_bins of them), then the
                raw counts in the \p n_bins - 1 bins between them, then the
                \p n_samples sampled values -- the three things the Python
                returned as a tuple, concatenated, because a numpy view is one
                array and the shim splits it at known offsets
*/
IMPBFFEXPORT void sample_kappa2_diffusion_with_traps(
        double sD2, double sA2, double fret_efficiency,
        int n_samples, int n_bins, double k2_min, double k2_max, int seed,
        double** out_view, int* n_out_view);

//! \f$\kappa^2\f$ for every donor/acceptor dipole pair, given the separations.
/*!
    Flat in, flat out -- the shaped front door is
    ``IMP.bff.photophysics.kappa2_from_dipoles``, which reshapes to
    `(n_d, n_a)`. Two names because they are two things: this one cannot know
    the caller's intended shape from a 1-D view, and a Python function that
    only reshapes is exactly the thin shim this layer is supposed to be.

    \f$\kappa = \hat\mu_D\cdot\hat\mu_A - 3(\hat\mu_D\cdot\hat r)(\hat\mu_A\cdot\hat r)\f$,
    squared. A zero-length separation contributes zero rather than a division
    by zero: coincident states are reachable and are not an error.

    \param[in] mu_donor flat `(n_d, 3)` transition dipoles
    \param[in] mu_acceptor flat `(n_a, 3)`
    \param[in] r_vectors flat `(n_d, n_a, 3)` donor-to-acceptor separations
    \param[out] out_view,n_out_view the flat `(n_d, n_a)` matrix
*/
IMPBFFEXPORT void kappa2_dipole_matrix(
        const std::vector<double>& mu_donor,
        const std::vector<double>& mu_acceptor,
        const std::vector<double>& r_vectors,
        double** out_view, int* n_out_view);

//! \f$p(\kappa^2)\f$ for isotropically oriented, *static* dipoles.
/*!
    The closed form, not a sample: the classic two-branch expression, singular
    at \f$\kappa^2 = 1\f$ and zero above 4. Static because each molecule keeps
    its orientation for the whole excited-state lifetime -- the dynamic limit is
    the delta function at 2/3 instead.

    \param[in] k2 the abscissa
    \param[out] out_view,n_out_view the density at each \p k2
*/
IMPBFFEXPORT void isotropic_kappa2_density(
        const std::vector<double>& k2, double** out_view, int* n_out_view);


//! Turn a \f$\kappa^2\f$ distribution into the distance-ratio distribution.
/*!
    A FRET measurement does not see \f$\kappa^2\f$; it sees an *apparent*
    distance. The two are related by
    \f$R_{app}/R_{DA} = (\langle\kappa^2\rangle/\kappa^2)^{1/6}\f$ -- a low
    \f$\kappa^2\f$ reads as a longer distance -- and this performs that change
    of variable, Jacobian included.

    **The Jacobian is the point.** Transforming the abscissa and carrying the
    weights across unchanged would give a distribution that is wrong wherever
    the mapping is not linear, which is everywhere: \f$|dк^2/dr| = 6\langle
    \kappa^2\rangle / r^7\f$. Weights are re-normalised after it is applied.

    The transformed points are then resampled onto a uniform axis by linear
    interpolation, which is what makes the result convolvable with a distance
    distribution.

    \param[in] k2_amp amplitudes, non-negative, not required to be normalised
    \param[in] k2_val \f$\kappa^2\f$ values, strictly positive
    \param[in] n_bins points on the output axis
    \param[out] out_view,n_out_view the axis (\p n_bins), then the interpolated
                weights (\p n_bins), then \f$\langle\kappa^2\rangle\f$ as a
                single trailing value
    \throws IMP::ValueException on mismatched lengths, a non-positive
            \f$\kappa^2\f$, a negative amplitude, or zero total amplitude

    One flat buffer rather than three returns, because a numpy view is one
    array. ``IMP.bff.photophysics.kappa2_to_distance_ratio`` is the front door
    and splits it into the `(axis, weights, mean)` tuple callers expect.
*/
IMPBFFEXPORT void kappa2_distance_ratio_transform(
        const std::vector<double>& k2_amp, const std::vector<double>& k2_val,
        int n_bins, double** out_view, int* n_out_view);

//! \f$\langle\kappa^2\rangle\f$ for isotropically and rapidly reorienting
//! dipoles: 2/3.
/*! A named constant rather than a literal, because the literal is what makes an
    isotropic assumption invisible at a call site. */
IMPBFFEXPORT double kappa2_isotropic();

//! The order parameter of the inter-dye angle, from the residual anisotropies.
/*!
    Sindbert et al., JACS 133:2463 (2011), eq. 10.

    \param[in] s2_donor,s2_acceptor second-rank order parameters of the two dyes
    \param[in] r_inf_AD residual anisotropy of the FRET-sensitised emission
    \param[in] r_0 fundamental anisotropy
    \return \f$(S^2_\delta, \delta)\f$ — the order parameter of the angle
            between the two symmetry axes, and that angle in radians

    The name was `s2delta`, accurate domain shorthand. It was renamed to
    `kappa2_order_parameters`, which was wrong twice over: these are the *dyes'*
    order parameters, not \f$\kappa^2\f$'s, and it returns an angle as well.
*/
IMPBFFEXPORT std::vector<double> s2_delta_from_anisotropy(
        double s2_donor, double s2_acceptor, double r_inf_AD,
        double r_0 = 0.38);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_ORIENTATIONFACTOR_H
