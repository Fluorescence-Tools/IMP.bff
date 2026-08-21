/**
 *  \file IMP/bff/DistanceCalibration.h
 *  \brief Empirical corrections from a computed distance to a measured one.
 *
 * A distance between the mean positions of two dye clouds is not the mean
 * dye-dye distance, and neither is what an experiment reports. These are the
 * fitted maps between them.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_DISTANCECALIBRATION_H
#define IMPBFF_DISTANCECALIBRATION_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Evaluate a polynomial by Horner's method.
/*!
    \param[in] x abscissa
    \param[in] coefficients **highest power first** -- what `np.polyfit`
               returns, and therefore what a fitted calibration carries. The
               opposite order was in use elsewhere and evaluated the same
               calibration to a different number (85.5 against 45.02 at
               `x = 45`, `c = [0, 1, 0.02]`); a caller holding ascending
               coefficients must reverse them.
    \return the polynomial's value; 0 for no coefficients
*/
IMPBFFEXPORT double polynomial_transfer(
        double x, const std::vector<double>& coefficients);

//! Evaluate a polynomial at many abscissae, as a managed view.
/*! The direct vector form; a caller with many abscissae uses this rather than
    dispatching a scalar across them. */
IMPBFFEXPORT void polynomial_transfer_vector(
        const std::vector<double>& x, const std::vector<double>& coefficients,
        double** out_view, int* n_out_view);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_DISTANCECALIBRATION_H
