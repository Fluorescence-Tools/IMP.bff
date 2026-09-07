/**
 * \file FRET.cpp
 * \brief FRET between two labels: the pair values and the rate trace.
 *
 * Sections in the order of IMP/bff/FRET.h; each is marked with the file it
 * came from.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

// -------- from FRETPair.cpp --------
/**
 * (formerly FRETPair.cpp, now a section of this file)
 * \brief Distances, orientation factors and efficiencies over all pairs.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/FRET.h>
#include <IMP/bff/internal/OutputView.h>

#include <cmath>
#include <cstdlib>
#include <limits>

IMPBFF_BEGIN_NAMESPACE

void fret_pair_matrices(
        double* points1, int n_points1, double* points2, int n_points2,
        double* mu1, int n_mu1, double* mu2, int n_mu2, int n1, int n2,
        double** out_view, int* n_out_view) {
    (void)n_points1; (void)n_points2;
    const std::size_t total =
            static_cast<std::size_t>(std::max(0, n1)) * std::max(0, n2) * 2;
    // malloc, not new[]: numpy's ARGOUTVIEWM typemap takes ownership and
    // releases it with free(). Mixing the two is undefined behaviour, and the
    // pairing is invisible from this side of the boundary.
    double* out = static_cast<double*>(std::calloc(total ? total : 1,
                                                   sizeof(double)));
    if (out == nullptr) {
        // An empty view rather than a null pointer: numpy would wrap the null
        // and the caller would fault on first touch, a long way from here.
        *out_view = static_cast<double*>(std::calloc(1, sizeof(double)));
        *n_out_view = 0;
        return;
    }
    *out_view = out;
    *n_out_view = static_cast<int>(total);
    if (n1 <= 0 || n2 <= 0) return;
    const std::size_t block = static_cast<std::size_t>(n1) * n2;
    const bool oriented = n_mu1 == n1 * 3 && n_mu2 == n2 * 3;

    // Normalised acceptor dipoles, once. Recomputing them inside the inner loop
    // would repeat n1 square roots per acceptor for no reason.
    std::vector<double> a_hat;
    if (oriented) {
        a_hat.assign(static_cast<std::size_t>(n2) * 3, 0.0);
        for (int j = 0; j < n2; ++j) {
            const double x = mu2[3 * j], y = mu2[3 * j + 1], z = mu2[3 * j + 2];
            const double n = std::sqrt(x * x + y * y + z * z);
            if (n > 0.0) {
                a_hat[3 * j] = x / n;
                a_hat[3 * j + 1] = y / n;
                a_hat[3 * j + 2] = z / n;
            }
        }
    }

#pragma omp parallel for schedule(static)
    for (int i = 0; i < n1; ++i) {
        const double px = points1[3 * i], py = points1[3 * i + 1],
                     pz = points1[3 * i + 2];
        double dx_hat = 0.0, dy_hat = 0.0, dz_hat = 0.0;
        if (oriented) {
            const double x = mu1[3 * i], y = mu1[3 * i + 1], z = mu1[3 * i + 2];
            const double n = std::sqrt(x * x + y * y + z * z);
            if (n > 0.0) { dx_hat = x / n; dy_hat = y / n; dz_hat = z / n; }
        }
        for (int j = 0; j < n2; ++j) {
            const double rx = px - points2[3 * j];
            const double ry = py - points2[3 * j + 1];
            const double rz = pz - points2[3 * j + 2];
            const double d = std::sqrt(rx * rx + ry * ry + rz * rz);
            const std::size_t k = static_cast<std::size_t>(i) * n2 + j;
            out[k] = d;
            if (!oriented) {
                out[block + k] = 2.0 / 3.0;
                continue;
            }
            // A zero separation has no direction, so the unit vector stays zero: the
            // two projection terms vanish and kappa^2 falls back to (mu_D . mu_A)^2.
            double ux = 0.0, uy = 0.0, uz = 0.0;
            if (d > 0.0) { ux = rx / d; uy = ry / d; uz = rz / d; }
            const double ax = a_hat[3 * j], ay = a_hat[3 * j + 1],
                         az = a_hat[3 * j + 2];
            const double cos_da = dx_hat * ax + dy_hat * ay + dz_hat * az;
            const double cos_dr = dx_hat * ux + dy_hat * uy + dz_hat * uz;
            const double cos_ar = ax * ux + ay * uy + az * uz;
            const double kappa = cos_da - 3.0 * cos_dr * cos_ar;
            out[block + k] = kappa * kappa;
        }
    }
}

void fret_pair_efficiency_matrices(
        double* r, int n_r, double* kappa2, int n_kappa2,
        double forster_radius, double** out_view, int* n_out_view) {
    const int n = n_r < n_kappa2 ? n_r : n_kappa2;
    double* out = internal::new_double_view(
            static_cast<std::size_t>(std::max(0, n)) * 2, out_view, n_out_view);
    if (out == nullptr || n <= 0) return;
    const std::size_t block = static_cast<std::size_t>(n);

#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        const double ratio = r[i] / forster_radius;
        const double r3 = ratio * ratio * ratio;
        const double ratio6 = r3 * r3;
        const double k2 = kappa2[i];

        // E = 1 / (1 + (2/3) ratio6 / kappa2). Both degenerate limits mean
        // complete transfer: kappa2 = 0 with ratio6 = 0 gives 0/0, and any
        // ratio6 with kappa2 = 0 gives +inf in the denominator's numerator.
        double e;
        if (k2 == 0.0) {
            e = ratio6 == 0.0 ? 1.0 : 0.0;   // nan -> 1, else 1/(1+inf) -> 0
        } else {
            const double denom = 1.0 + (2.0 / 3.0) * ratio6 / k2;
            e = std::isfinite(denom) ? 1.0 / denom : 0.0;
            if (std::isnan(e)) e = 1.0;
        }
        out[i] = e;
        // Left unsanitised on purpose: an infinite rate at zero separation is
        // the right answer, and the caller decides what to do with it.
        out[block + i] = ratio6 == 0.0
                                 ? std::numeric_limits<double>::infinity()
                                 : 1.5 * k2 / ratio6;
    }
}

void FRETPairGeometry::get_R(double** out_view, int* n_out_view) const {
    internal::copy_to_view(R, out_view, n_out_view);
}
void FRETPairGeometry::get_kappa2(double** out_view, int* n_out_view) const {
    internal::copy_to_view(kappa2, out_view, n_out_view);
}
void FRETPairGeometry::get_weight(double** out_view, int* n_out_view) const {
    internal::copy_to_view(weight, out_view, n_out_view);
}

void FRETPairEfficiencies::get_E(double** out_view, int* n_out_view) const {
    internal::copy_to_view(E, out_view, n_out_view);
}
void FRETPairEfficiencies::get_rate_ratio(double** out_view,
                                          int* n_out_view) const {
    internal::copy_to_view(rate_ratio, out_view, n_out_view);
}
void FRETPairEfficiencies::get_k_fret(double** out_view, int* n_out_view) const {
    internal::copy_to_view(k_fret, out_view, n_out_view);
}
void FRETPairEfficiencies::get_R(double** out_view, int* n_out_view) const {
    internal::copy_to_view(R, out_view, n_out_view);
}
void FRETPairEfficiencies::get_kappa2(double** out_view, int* n_out_view) const {
    internal::copy_to_view(kappa2, out_view, n_out_view);
}
void FRETPairEfficiencies::get_weight(double** out_view, int* n_out_view) const {
    internal::copy_to_view(weight, out_view, n_out_view);
}

FRETPairGeometry fret_pair_geometry(const std::vector<double>& points1,
                                    const std::vector<double>& weights1,
                                    const std::vector<double>& points2,
                                    const std::vector<double>& weights2,
                                    const std::vector<double>& mu1,
                                    const std::vector<double>& mu2) {
    FRETPairGeometry out;
    out.n1 = static_cast<int>(points1.size() / 3);
    out.n2 = static_cast<int>(points2.size() / 3);
    const std::size_t n = static_cast<std::size_t>(out.n1) * out.n2;
    if (n == 0) return out;

    double* packed = NULL;
    int n_packed = 0;
    // The kernel takes raw pointers because SWIG hands it numpy buffers; the
    // casts are the same buffers, const-stripped for that signature.
    fret_pair_matrices(const_cast<double*>(points1.data()),
                       static_cast<int>(points1.size()),
                       const_cast<double*>(points2.data()),
                       static_cast<int>(points2.size()),
                       mu1.empty() ? NULL : const_cast<double*>(mu1.data()),
                       static_cast<int>(mu1.size()),
                       mu2.empty() ? NULL : const_cast<double*>(mu2.data()),
                       static_cast<int>(mu2.size()),
                       out.n1, out.n2, &packed, &n_packed);
    out.R.assign(packed, packed + n);
    out.kappa2.assign(packed + n, packed + 2 * n);
    std::free(packed);

    out.weight.resize(n);
    double total = 0.0;
    for (int i = 0; i < out.n1; ++i) {
        for (int j = 0; j < out.n2; ++j) {
            const double w = weights1[i] * weights2[j];
            out.weight[static_cast<std::size_t>(i) * out.n2 + j] = w;
            total += w;
        }
    }
    if (total > 0.0) {
        for (std::size_t i = 0; i < n; ++i) out.weight[i] /= total;
    }
    out.kappa2_avg = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        out.kappa2_avg += out.kappa2[i] * out.weight[i];
    }
    return out;
}

FRETPairEfficiencies fret_pair_efficiencies(const FRETPairGeometry& geometry,
                                            double forster_radius,
                                            double tau0) {
    FRETPairEfficiencies out;
    out.n1 = geometry.n1;
    out.n2 = geometry.n2;
    out.kappa2_avg = geometry.kappa2_avg;
    out.forster_radius = forster_radius;
    out.R = geometry.R;
    out.kappa2 = geometry.kappa2;
    out.weight = geometry.weight;
    const std::size_t n = geometry.R.size();
    if (n == 0) return out;

    double* packed = NULL;
    int n_packed = 0;
    fret_pair_efficiency_matrices(
            const_cast<double*>(geometry.R.data()),
            static_cast<int>(geometry.R.size()),
            const_cast<double*>(geometry.kappa2.data()),
            static_cast<int>(geometry.kappa2.size()), forster_radius, &packed,
            &n_packed);
    out.E.assign(packed, packed + n);
    out.rate_ratio.assign(packed + n, packed + 2 * n);
    std::free(packed);

    double e_static = 0.0, e_dyn1 = 0.0, rate_avg = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        e_static += out.E[i] * geometry.weight[i];
        // `dynamic1` uses the *ensemble* kappa^2 rather than the per-pair one,
        // so it is a different expression from `E` and is formed here.
        const double x = geometry.R[i] / forster_radius;
        const double x3 = x * x * x;
        const double ratio6 = x3 * x3;
        double d1 = 1.0;
        if (out.kappa2_avg > 0.0 && std::isfinite(ratio6)) {
            d1 = 1.0 / (1.0 + (2.0 / 3.0) / out.kappa2_avg * ratio6);
        }
        e_dyn1 += d1 * geometry.weight[i];
        if (std::isfinite(out.rate_ratio[i])) {
            rate_avg += out.rate_ratio[i] * geometry.weight[i];
        }
    }
    out.static_efficiency = e_static;
    out.dynamic1 = e_dyn1;
    out.dynamic2 = rate_avg / (rate_avg + 1.0);
    if (tau0 > 0.0) {
        out.k_fret.resize(n);
        for (std::size_t i = 0; i < n; ++i) out.k_fret[i] = out.rate_ratio[i] / tau0;
    }
    return out;
}

FRETPairEfficiencies fret_pair_distribution(
        const std::vector<double>& points1, const std::vector<double>& weights1,
        const std::vector<double>& points2, const std::vector<double>& weights2,
        double forster_radius, const std::vector<double>& mu1,
        const std::vector<double>& mu2, double tau0) {
    return fret_pair_efficiencies(fret_pair_geometry(points1, weights1,
                                                     points2, weights2, mu1,
                                                     mu2),
                                  forster_radius, tau0);
}

IMPBFF_END_NAMESPACE

// -------- from FRETRateTrace.cpp --------
/**
 * (formerly FRETRateTrace.cpp, now a section of this file)
 * \brief FRET rate along a dye trajectory.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/OrientationFactor.h>

#include <IMP/bff/Base.h>


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
