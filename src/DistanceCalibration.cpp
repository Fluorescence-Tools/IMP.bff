/**
 * \file DistanceCalibration.cpp
 * \brief Empirical corrections from a computed distance to a measured one.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/DistanceCalibration.h>
#include <IMP/bff/internal/OutputView.h>

IMPBFF_BEGIN_NAMESPACE

double polynomial_transfer(double x, const std::vector<double>& coefficients) {
    double y = 0.0;
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        y = y * x + coefficients[i];
    }
    return y;
}

void polynomial_transfer_vector(
        const std::vector<double>& x, const std::vector<double>& coefficients,
        double** out_view, int* n_out_view) {
    std::vector<double> y(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        y[i] = polynomial_transfer(x[i], coefficients);
    }
    internal::copy_to_view(y, out_view, n_out_view);
}

IMPBFF_END_NAMESPACE
