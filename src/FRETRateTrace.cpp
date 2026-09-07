/**
 * \file FRETRateTrace.cpp
 * \brief FRET rate along a dye trajectory.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/FRETRateTrace.h>
#include <IMP/bff/OrientationFactor.h>
#include <IMP/bff/internal/OutputView.h>

#include <IMP/bff/Base.h>

#include <cmath>

IMPBFF_BEGIN_NAMESPACE

void fret_rate_trace_kernel(
        const std::vector<double>& trajectory,
        const std::vector<double>& acceptor_points, double R0, double tau0,
        double r_min2, double kappa2_scale,
        double** out_view, int* n_out_view) {
    const std::size_t n_frames = trajectory.size() / 3;
    const std::size_t n_acceptor = acceptor_points.size() / 3;
    double* rates = internal::new_double_view(n_frames, out_view, n_out_view);
    if (rates == nullptr || n_acceptor == 0 || n_frames == 0) return;
    const double R0_6 = std::pow(R0, 6) * kappa2_scale;
    const double inv_tau0 = 1.0 / tau0;
    for (std::size_t f = 0; f < n_frames; ++f) {
        const double px = trajectory[3 * f + 0];
        const double py = trajectory[3 * f + 1];
        const double pz = trajectory[3 * f + 2];
        double total = 0.0;
        for (std::size_t j = 0; j < n_acceptor; ++j) {
            const double dx = px - acceptor_points[3 * j + 0];
            const double dy = py - acceptor_points[3 * j + 1];
            const double dz = pz - acceptor_points[3 * j + 2];
            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < r_min2) r2 = r_min2;   // 1/r^6 diverges; dyes cannot overlap
            total += R0_6 / (r2 * r2 * r2);
        }
        // Arithmetic mean of rates: the fast-exchange limit.
        rates[f] = inv_tau0 * total / static_cast<double>(n_acceptor);
    }
}

void fret_rate_pair_trace_kernel(
        const std::vector<double>& donor, const std::vector<double>& acceptor,
        double R0, double tau0, double r_min2, double kappa2_scale,
        double** out_view, int* n_out_view) {
    const std::size_t n_frames = donor.size() / 3;
    double* rates = internal::new_double_view(n_frames, out_view, n_out_view);
    if (rates == nullptr) return;
    const double R0_6 = std::pow(R0, 6) * kappa2_scale;
    const double inv_tau0 = 1.0 / tau0;
    for (std::size_t f = 0; f < n_frames; ++f) {
        const double dx = donor[3 * f + 0] - acceptor[3 * f + 0];
        const double dy = donor[3 * f + 1] - acceptor[3 * f + 1];
        const double dz = donor[3 * f + 2] - acceptor[3 * f + 2];
        double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 < r_min2) r2 = r_min2;
        rates[f] = inv_tau0 * R0_6 / (r2 * r2 * r2);
    }
}

const double MAX_KAPPA2 = 4.0;

double kappa2_scale(double kappa2) {
    const double value = kappa2;
    // `>=` and `<=` negated rather than `<`/`>`, so NaN fails both: a NaN from a
    // failed orientation calculation must raise, not be answered with 2/3.
    if (!(value >= 0.0) || !(value <= MAX_KAPPA2)) {
        IMP_THROW("kappa2 must be in [0, " << MAX_KAPPA2
                          << "] (0 = perpendicular dipoles, 4 = collinear); got "
                          << kappa2,
                  ValueException);
    }
    return value / kappa2_isotropic();
}

void fret_rate_trace(const std::vector<double>& trajectory,
                     const std::vector<double>& acceptor_points, double R0,
                     double tau0, double r_min, int max_acceptor_points,
                     double kappa2, double** out_view, int* n_out_view) {
    const std::size_t n_points = acceptor_points.size() / 3;
    if (n_points == 0) {
        IMP_THROW("Acceptor accessible volume is empty; cannot compute FRET "
                  "rates.",
                  ValueException);
    }
    if (!(tau0 > 0.0)) {
        IMP_THROW("tau0 must be positive to derive a FRET rate, not " << tau0,
                  ValueException);
    }
    // The scale is validated before any striding, so a bad kappa2 raises
    // whether or not the cloud needed thinning.
    const double scale = kappa2_scale(kappa2);

    const std::vector<double>* points = &acceptor_points;
    std::vector<double> strided;
    if (max_acceptor_points > 0 &&
        n_points > static_cast<std::size_t>(max_acceptor_points)) {
        const std::size_t stride =
                (n_points + max_acceptor_points - 1) / max_acceptor_points;
        strided.reserve(3 * (n_points / stride + 1));
        for (std::size_t i = 0; i < n_points; i += stride) {
            strided.push_back(acceptor_points[3 * i + 0]);
            strided.push_back(acceptor_points[3 * i + 1]);
            strided.push_back(acceptor_points[3 * i + 2]);
        }
        points = &strided;
    }
    const double floor_r = r_min > 1e-6 ? r_min : 1e-6;
    fret_rate_trace_kernel(trajectory, *points, R0, tau0, floor_r * floor_r,
                           scale, out_view, n_out_view);
}

void fret_rate_pair_trace(const std::vector<double>& donor_trajectory,
                          const std::vector<double>& acceptor_trajectory,
                          double R0, double tau0, double r_min, double kappa2,
                          double** out_view, int* n_out_view) {
    const std::size_t nd = donor_trajectory.size() / 3;
    const std::size_t na = acceptor_trajectory.size() / 3;
    if (nd == 0 || na == 0) {
        IMP_THROW("Both trajectories must have frames to pair.", ValueException);
    }
    if (!(tau0 > 0.0)) {
        IMP_THROW("tau0 must be positive to derive a FRET rate, not " << tau0,
                  ValueException);
    }
    if (nd != na) {
        IMP_THROW(
                "Donor and acceptor trajectories must have the same number of "
                "frames to pair ("
                        << nd << " vs " << na
                        << "). This is expected whenever the two dyes are "
                           "quenched differently -- they sit at different sites "
                           "-- so they need different numbers of excitations for "
                           "the same photon count. Simulate the same number of "
                           "*walks* for both dyes; t_max and t_step are usually "
                           "already shared and are not the lever. Truncating "
                           "here is refused because a trajectory concatenates "
                           "one walk per excitation, so an arbitrary cut splits "
                           "a walk and pairs one excitation's tail against "
                           "another's head.",
                ValueException);
    }
    const double floor_r = r_min > 1e-6 ? r_min : 1e-6;
    fret_rate_pair_trace_kernel(donor_trajectory, acceptor_trajectory, R0, tau0,
                                floor_r * floor_r, kappa2_scale(kappa2),
                                out_view, n_out_view);
}

IMPBFF_END_NAMESPACE
