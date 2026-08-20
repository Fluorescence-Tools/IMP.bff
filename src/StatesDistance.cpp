/**
 * \file StatesDistance.cpp
 * \brief Distances between two labels, whatever represents them.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/StatesDistance.h>

#include <IMP/bff/AVBuilder.h>
#include <IMP/bff/AVDistance.h>
#include <IMP/bff/DistanceCalibration.h>
#include <IMP/bff/AVModel.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/exception.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <random>

IMPBFF_BEGIN_NAMESPACE

namespace sd {

//! A cloud as a flat (x, y, z, w) buffer this side owns.
std::vector<double> cloud(const States& s) {
    double* points = NULL;
    int n = 0;
    s.get_points(&points, &n);
    std::vector<double> out(points, points + n);
    std::free(points);
    return out;
}

std::vector<double> mean_position(const States& s) {
    double* p = NULL;
    int n = 0;
    s.get_mean_position(&p, &n);
    std::vector<double> out(p, p + n);
    std::free(p);
    out.resize(3, 0.0);
    return out;
}

//! `(distance, weight)` pairs, from the one sampler.
std::vector<double> sample(const States& s1, const States& s2, int n_samples) {
    if (s1.get_n_points() == 0 || s2.get_n_points() == 0) {
        IMP_THROW("cannot sample a distance: one or both labels have no states",
                  ValueException);
    }
    double* out = NULL;
    int n = 0;
    random_distances(cloud(s1), cloud(s2), n_samples, 0, &out, &n);
    std::vector<double> pairs(out, out + n);
    std::free(out);
    return pairs;
}

void split(const std::vector<double>& pairs, std::vector<double>* d,
           std::vector<double>* w) {
    d->reserve(pairs.size() / 2);
    w->reserve(pairs.size() / 2);
    for (std::size_t i = 0; i + 1 < pairs.size(); i += 2) {
        d->push_back(pairs[i]);
        w->push_back(pairs[i + 1]);
    }
}

}  // namespace sd

double states_average_distance(const States& s1, const States& s2,
                               int n_samples) {
    std::vector<double> d, w;
    sd::split(sd::sample(s1, s2, n_samples), &d, &w);
    double num = 0.0, den = 0.0;
    for (std::size_t i = 0; i < d.size(); ++i) {
        num += d[i] * w[i];
        den += w[i];
    }
    return den > 0.0 ? num / den : 0.0;
}

double states_mean_fret_distance(const States& s1, const States& s2,
                                 double forster_radius, int n_samples) {
    std::vector<double> d, w;
    sd::split(sd::sample(s1, s2, n_samples), &d, &w);
    return distance_sample_statistics(d, w, forster_radius)[1];
}

double distance_between_mean_positions(const States& s1, const States& s2) {
    const std::vector<double> a = sd::mean_position(s1);
    const std::vector<double> b = sd::mean_position(s2);
    const double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double standard_deviation_of_distances(const States& s1, const States& s2,
                                       int n_samples) {
    std::vector<double> d, w;
    sd::split(sd::sample(s1, s2, n_samples), &d, &w);
    return distance_sample_statistics(d, w)[3];
}

std::vector<double> av_pair_statistics(const States& s1, const States& s2,
                                       double forster_radius, int n_samples) {
    const double rmp = distance_between_mean_positions(s1, s2);
    std::vector<double> out(4, 0.0);
    out[0] = out[1] = out[2] = rmp;
    if (s1.get_n_points() == 0 || s2.get_n_points() == 0) return out;

    std::vector<double> d, w;
    sd::split(sd::sample(s1, s2, n_samples), &d, &w);
    double total = 0.0;
    for (std::size_t i = 0; i < w.size(); ++i) total += w[i];
    // A zero total weight is a cloud that exists but carries no occupancy; the
    // mean position is still the one thing it has.
    if (total <= 0.0) return out;

    const std::vector<double> stats =
            distance_sample_statistics(d, w, forster_radius);
    out[1] = stats[0];
    out[2] = stats[1];
    out[3] = stats[3];
    return out;
}

double model_distance(const States& s1, const States& s2,
                      const std::string& distance_type, double forster_radius,
                      int n_samples) {
    if (distance_type == "Rmp") return distance_between_mean_positions(s1, s2);
    if (distance_type == "RDAMean") {
        return states_average_distance(s1, s2, n_samples);
    }
    if (distance_type == "RDAMeanE") {
        return states_mean_fret_distance(s1, s2, forster_radius, n_samples);
    }
    IMP_THROW("Unknown distance type: " << distance_type, ValueException);
}

void histogram_rda(const States& s1, const States& s2,
                   const std::vector<double>& axis, int n_samples,
                   bool normalize, double** out_view, int* n_out_view) {
    const std::size_t n_bins = axis.size() > 1 ? axis.size() - 1 : 0;
    double* out = internal::new_double_view(n_bins, out_view, n_out_view);
    if (out == NULL || n_bins == 0) return;

    std::vector<double> d, w;
    sd::split(sd::sample(s1, s2, n_samples), &d, &w);

    double total = 0.0;
    for (std::size_t i = 0; i < d.size(); ++i) {
        const double r = d[i];
        if (r < axis.front() || r > axis.back()) continue;
        // Upper edge closed, as numpy's histogram has it.
        std::size_t b = 0;
        while (b + 1 < n_bins && r >= axis[b + 1]) ++b;
        out[b] += w[i];
        total += w[i];
    }
    if (normalize && total > 0.0) {
        for (std::size_t b = 0; b < n_bins; ++b) out[b] /= total;
    }
}

std::vector<double> fit_transfer_polynomial(const States& s1, const States& s2,
                                            const std::string& distance_type,
                                            double forster_radius, int degree,
                                            int n_samples) {
    // y = x: there is nothing to correct.
    std::vector<double> identity(degree + 1, 0.0);
    if (degree >= 1) identity[degree - 1] = 1.0;

    const double rmp = distance_between_mean_positions(s1, s2);
    if (s1.get_n_points() == 0 || s2.get_n_points() == 0 || rmp <= 1e-6) {
        return identity;
    }

    // Separation vectors rather than distances: the fit translates one cloud
    // along the join, and `|v + t u|` needs the vector.
    const std::vector<double> c1 = sd::cloud(s1), c2 = sd::cloud(s2);
    const std::size_t n1 = c1.size() / 4, n2 = c2.size() / 4;
    std::vector<double> v_dot_u(n_samples), v_sq(n_samples), w(n_samples);

    const std::vector<double> m1 = sd::mean_position(s1);
    const std::vector<double> m2 = sd::mean_position(s2);
    const double ux = (m2[0] - m1[0]) / rmp;
    const double uy = (m2[1] - m1[1]) / rmp;
    const double uz = (m2[2] - m1[2]) / rmp;

    // The same linear congruential draw the sampler uses, so a fit and a
    // distance drawn at the same seed see the same pairs.
    unsigned int state = 1u;
    double w_sum = 0.0;
    for (int i = 0; i < n_samples; ++i) {
        state = state * 1103515245u + 12345u;
        const std::size_t i1 = (state >> 16) % n1;
        state = state * 1103515245u + 12345u;
        const std::size_t i2 = (state >> 16) % n2;
        const double vx = c2[i2 * 4 + 0] - c1[i1 * 4 + 0];
        const double vy = c2[i2 * 4 + 1] - c1[i1 * 4 + 1];
        const double vz = c2[i2 * 4 + 2] - c1[i1 * 4 + 2];
        v_dot_u[i] = vx * ux + vy * uy + vz * uz;
        v_sq[i] = vx * vx + vy * vy + vz * vz;
        w[i] = c1[i1 * 4 + 3] * c2[i2 * 4 + 3];
        w_sum += w[i];
    }
    if (w_sum <= 0.0) return identity;

    const int n_points = 7;
    const double t_min = -std::min(rmp - 5.0, 15.0);
    const double t_max = 20.0;
    std::vector<double> xs(n_points), ys(n_points);
    for (int k = 0; k < n_points; ++k) {
        const double t = t_min + (t_max - t_min) * k / (n_points - 1);
        xs[k] = rmp + t;
        double num = 0.0, e_sum = 0.0;
        for (int i = 0; i < n_samples; ++i) {
            const double d = std::sqrt(std::max(
                    v_sq[i] + 2.0 * t * v_dot_u[i] + t * t, 1e-10));
            num += d * w[i];
            if (distance_type == "RDAMeanE") {
                const double x = d / forster_radius;
                const double x3 = x * x * x;
                e_sum += w[i] / (1.0 + x3 * x3);
            }
        }
        if (distance_type == "RDAMeanE") {
            const double mean_e = e_sum / w_sum;
            if (mean_e <= 0.0) ys[k] = num / w_sum;
            else if (mean_e >= 1.0) ys[k] = 0.0;
            else ys[k] = forster_radius * std::pow(1.0 / mean_e - 1.0, 1.0 / 6.0);
        } else {
            ys[k] = num / w_sum;
        }
    }

    // Least squares on the Vandermonde normal equations, highest power first --
    // the order `numpy.polyfit` returns and `polynomial_transfer` expects.
    const int m = degree + 1;
    std::vector<double> ata(m * m, 0.0), atb(m, 0.0);
    for (int k = 0; k < n_points; ++k) {
        std::vector<double> row(m);
        double p = 1.0;
        for (int j = m - 1; j >= 0; --j) { row[j] = p; p *= xs[k]; }
        for (int a = 0; a < m; ++a) {
            atb[a] += row[a] * ys[k];
            for (int b = 0; b < m; ++b) ata[a * m + b] += row[a] * row[b];
        }
    }
    // Gaussian elimination with partial pivoting; m is 4 in every caller.
    for (int col = 0; col < m; ++col) {
        int pivot = col;
        for (int r = col + 1; r < m; ++r) {
            if (std::fabs(ata[r * m + col]) > std::fabs(ata[pivot * m + col])) {
                pivot = r;
            }
        }
        if (std::fabs(ata[pivot * m + col]) < 1e-12) return identity;
        if (pivot != col) {
            for (int c = 0; c < m; ++c) {
                std::swap(ata[col * m + c], ata[pivot * m + c]);
            }
            std::swap(atb[col], atb[pivot]);
        }
        for (int r = col + 1; r < m; ++r) {
            const double f = ata[r * m + col] / ata[col * m + col];
            for (int c = col; c < m; ++c) ata[r * m + c] -= f * ata[col * m + c];
            atb[r] -= f * atb[col];
        }
    }
    std::vector<double> coeffs(m, 0.0);
    for (int r = m - 1; r >= 0; --r) {
        double acc = atb[r];
        for (int c = r + 1; c < m; ++c) acc -= ata[r * m + c] * coeffs[c];
        coeffs[r] = acc / ata[r * m + r];
    }
    return coeffs;
}

double gaussian_rmp_to_rda_mean(double rmp, double sigma) {
    if (!(rmp > 0.0)) return 0.0;
    return rmp + sigma * sigma / rmp;
}

double polynomial_transfer_ascending(double rmp,
                                     const std::vector<double>& coeffs) {
    // One evaluator underneath: reverse and delegate.
    std::vector<double> descending(coeffs.rbegin(), coeffs.rend());
    return polynomial_transfer(rmp, descending);
}

double mean_position_distance(const std::vector<double>& points_a,
                              const std::vector<double>& points_b,
                              const std::vector<double>& weights_a,
                              const std::vector<double>& weights_b) {
    const std::size_t na = points_a.size() / 3, nb = points_b.size() / 3;
    if (na == 0 || nb == 0) {
        IMP_THROW("both point clouds must be non-empty to define R_mp",
                  ValueException);
    }
    double ma[3] = {0, 0, 0}, mb[3] = {0, 0, 0}, wa = 0.0, wb = 0.0;
    for (std::size_t i = 0; i < na; ++i) {
        const double w = weights_a.empty() ? 1.0 : weights_a[i];
        for (int c = 0; c < 3; ++c) ma[c] += w * points_a[i * 3 + c];
        wa += w;
    }
    for (std::size_t i = 0; i < nb; ++i) {
        const double w = weights_b.empty() ? 1.0 : weights_b[i];
        for (int c = 0; c < 3; ++c) mb[c] += w * points_b[i * 3 + c];
        wb += w;
    }
    if (wa == 0.0 || wb == 0.0) {
        IMP_THROW("a cloud with zero total weight has no mean position",
                  ValueException);
    }
    double d2 = 0.0;
    for (int c = 0; c < 3; ++c) {
        const double d = ma[c] / wa - mb[c] / wb;
        d2 += d * d;
    }
    return std::sqrt(d2);
}

// --------------------------------------------------------------------------
// label distributions
// --------------------------------------------------------------------------

LabelDistribution::LabelDistribution(const std::string& simulation_type,
                                     const std::vector<double>& origin,
                                     double simulation_grid_resolution,
                                     const std::string& position_name)
    : simulation_type_(simulation_type), origin_(origin),
      simulation_grid_resolution_(simulation_grid_resolution),
      position_name_(position_name), computed_(false) {}

const AccessibleVolume& LabelDistribution::get_accessible_volume() const {
    if (!computed_) {
        do_compute();
        computed_ = true;
    }
    return av_;
}

void LabelDistribution::get_origin(double** out_view, int* n_out_view) const {
    internal::copy_to_view(origin_, out_view, n_out_view);
}

LabelDistributionAV::LabelDistributionAV(double* atoms_xyzr, int n_atoms,
                                         int n_cols,
                                         const std::vector<double>& source_xyz,
                                         double linker_length,
                                         double linker_width, double r1,
                                         double r2, double r3,
                                         double simulation_grid_resolution,
                                         const std::string& position_name)
    : LabelDistribution(r2 == 0.0 ? "AV1" : "AV3", source_xyz,
                        simulation_grid_resolution, position_name),
      source_xyz_(source_xyz), linker_length_(linker_length),
      linker_width_(linker_width), r1_(r1), r2_(r2), r3_(r3) {
    if (n_cols != 4) {
        IMP_THROW("obstacles must be (N, 4) of x, y, z, radius, not (" << n_atoms
                          << ", " << n_cols << ")",
                  ValueException);
    }
    atoms_xyzr_.assign(atoms_xyzr, atoms_xyzr + static_cast<std::size_t>(n_atoms) * 4);
}

void LabelDistributionAV::do_compute() const {
    av_ = compute_av(const_cast<double*>(atoms_xyzr_.data()),
                     static_cast<int>(atoms_xyzr_.size() / 4), 4, source_xyz_,
                     linker_length_, linker_width_, r1_, r2_, r3_,
                     simulation_grid_resolution_, DEFAULT_ALLOWED_SPHERE_RADIUS,
                     0);
    av_.set_position_name(position_name_);
}

DyeDistributionNormal::DyeDistributionNormal(const std::vector<double>& origin,
                                             double width, int n_points,
                                             int seed,
                                             const std::string& position_name)
    : LabelDistribution("Normal", origin, 1.0, position_name), width_(width),
      n_points_(n_points), seed_(seed) {
    if (origin.size() != 3) {
        IMP_THROW("the origin is three coordinates, not " << origin.size(),
                  ValueException);
    }
}

void DyeDistributionNormal::do_compute() const {
    // Box-Muller off a seeded engine. The Python drew from the *global*
    // `numpy.random` state, so two labels built in one process were not
    // reproducible and neither was a single one across runs.
    std::mt19937 engine(static_cast<unsigned int>(seed_));
    std::normal_distribution<double> normal(0.0, 1.0);

    std::vector<double> points(static_cast<std::size_t>(n_points_) * 4);
    double total = 0.0;
    for (int i = 0; i < n_points_; ++i) {
        double r2 = 0.0;
        for (int c = 0; c < 3; ++c) {
            const double d = normal(engine) * width_;
            points[i * 4 + c] = origin_[c] + d;
            r2 += d * d;
        }
        // The weight is the Gaussian height at the drawn point, which makes the
        // cloud an importance-weighted sample of the same distribution it was
        // drawn from rather than a uniform one.
        const double w = std::exp(-0.5 * r2 / (width_ * width_));
        points[i * 4 + 3] = w;
        total += w;
    }
    if (total > 0.0) {
        for (int i = 0; i < n_points_; ++i) points[i * 4 + 3] /= total;
    }
    av_ = AccessibleVolume(points, std::vector<double>(), std::vector<double>(),
                           1.0, position_name_, origin_);
}

// --------------------------------------------------------------------------
// the Gaussian pair converter
// --------------------------------------------------------------------------

FRETDistanceConverter::FRETDistanceConverter(double forster_radius,
                                             double sigma, double distance_min,
                                             double distance_max,
                                             int n_distances)
    : forster_radius_(forster_radius), sigma_(sigma) {
    distances_.resize(n_distances);
    for (int i = 0; i < n_distances; ++i) {
        distances_[i] = distance_min + (distance_max - distance_min) * i /
                                               std::max(1, n_distances - 1);
    }
    update_efficiencies();
    update_lookup();
}

void FRETDistanceConverter::update_efficiencies() {
    efficiencies_.resize(distances_.size());
    for (std::size_t i = 0; i < distances_.size(); ++i) {
        const double x = distances_[i] / forster_radius_;
        const double x3 = x * x * x;
        efficiencies_[i] = 1.0 / (1.0 + x3 * x3);
    }
}

void FRETDistanceConverter::update_lookup() {
    const std::size_t n = distances_.size();
    d_mean_.assign(n, 0.0);
    d_mean_fret_.assign(n, 0.0);
    e_mean_.assign(n, 0.0);

    const double norm = 1.0 / (sigma_ * std::sqrt(2.0 * M_PI));
    for (std::size_t k = 0; k < n; ++k) {
        const double dcc = distances_[k];
        // Two mirrored normals: the separation of two isotropic clouds is
        // symmetric about zero, and only its magnitude is observable.
        std::vector<double> p(n);
        double total = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double a = (distances_[i] + dcc) / sigma_;
            const double b = (distances_[i] - dcc) / sigma_;
            p[i] = norm * (std::exp(-0.5 * a * a) + std::exp(-0.5 * b * b));
            total += p[i];
        }
        if (total > 0.0) {
            for (std::size_t i = 0; i < n; ++i) p[i] /= total;
        }
        double dm = 0.0, em = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            dm += p[i] * distances_[i];
            em += p[i] * efficiencies_[i];
        }
        d_mean_[k] = dm;
        e_mean_[k] = em;
        d_mean_fret_[k] =
                em > 0.0 ? forster_radius_ * std::pow(1.0 / em - 1.0, 1.0 / 6.0)
                         : dm;
    }
}

void FRETDistanceConverter::set_forster_radius(double v) {
    forster_radius_ = v;
    update_efficiencies();
    update_lookup();
}

void FRETDistanceConverter::set_sigma(double v) {
    sigma_ = v;
    update_lookup();
}

namespace sd {

//! Linear interpolation on an ascending axis, clamped at both ends.
double interp(double x, const std::vector<double>& xs,
              const std::vector<double>& ys) {
    if (xs.empty()) return 0.0;
    if (x <= xs.front()) return ys.front();
    if (x >= xs.back()) return ys.back();
    std::size_t i = 0;
    while (i + 1 < xs.size() && xs[i + 1] < x) ++i;
    const double t = (x - xs[i]) / (xs[i + 1] - xs[i]);
    return ys[i] + t * (ys[i + 1] - ys[i]);
}

}  // namespace sd

double FRETDistanceConverter::get_distance_mean(double dcc) const {
    return sd::interp(dcc, distances_, d_mean_);
}

double FRETDistanceConverter::get_distance_mean_fret(double dcc) const {
    return sd::interp(dcc, distances_, d_mean_fret_);
}

double FRETDistanceConverter::get_fret_efficiency_mean(double dcc) const {
    return sd::interp(dcc, distances_, e_mean_);
}

double FRETDistanceConverter::get_effective_distance(double value,
                                                     int distance_type) const {
    if (distance_type == DYE_PAIR_DISTANCE_MEAN) return get_distance_mean(value);
    if (distance_type == DYE_PAIR_DISTANCE_E) return get_distance_mean_fret(value);
    if (distance_type == DYE_PAIR_DISTANCE_MP) return value;
    IMP_THROW("unknown dye-pair distance type " << distance_type,
              ValueException);
}

double effective_distance(double rmp,
                          const std::string& transfer_function_type,
                          double sigma_rda,
                          const std::vector<double>& coeffs) {
    if (transfer_function_type == "None") return rmp;
    if (transfer_function_type == "Gaussian") {
        return sigma_rda > 0.0 ? gaussian_rmp_to_rda_mean(rmp, sigma_rda) : rmp;
    }
    if (transfer_function_type == "Polynomial") {
        if (!coeffs.empty()) return polynomial_transfer_ascending(rmp, coeffs);
        // A calibration that names a polynomial and carries no coefficients has
        // not been fitted; sigma is the parametric stand-in.
        return sigma_rda > 0.0 ? gaussian_rmp_to_rda_mean(rmp, sigma_rda) : rmp;
    }
    return rmp;
}

IMPBFF_END_NAMESPACE
