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

IMPBFF_END_NAMESPACE

#endif //IMPBFF_ORIENTATIONFACTOR_H
