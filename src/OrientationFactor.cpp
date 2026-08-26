/**
 * \file OrientationFactor.cpp
 * \brief The FRET orientation factor kappa^2, from geometry and from wobbling.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/OrientationFactor.h>
#include <IMP/exception.h>
#include <IMP/bff/internal/OutputView.h>

#include <algorithm>
#include <cmath>
#include <random>

IMPBFF_BEGIN_NAMESPACE

namespace {
//! Second-rank Legendre polynomial of a cosine, P2(cos x).
inline double p2_of_angle(double x) {
    const double c = std::cos(x);
    return (3.0 * c * c - 1.0) / 2.0;
}

//! Bin into edges that are uniform by construction; the last edge is inclusive.
void accumulate(std::vector<double>& hist, const std::vector<double>& scale,
                double value, double weight) {
    if (scale.size() < 2) return;
    const double lo = scale.front(), hi = scale.back();
    if (value < lo || value > hi) return;
    const double step = (hi - lo) / static_cast<double>(scale.size() - 1);
    std::size_t k = static_cast<std::size_t>((value - lo) / step);
    if (k >= hist.size()) k = hist.size() - 1;   // the closed right edge
    hist[k] += weight;
}

std::vector<double> bin_edges(int n_bins, double k2_min, double k2_max) {
    std::vector<double> scale;
    if (n_bins < 2) return scale;
    const double step = (k2_max - k2_min) / (n_bins - 1);
    scale.reserve(n_bins);
    for (int i = 0; i < n_bins; ++i) scale.push_back(k2_min + i * step);
    return scale;
}
}  // namespace

double wobbling_kappa2(double delta, double sD2, double sA2,
                       double beta1, double beta2) {
    const double s2delta = p2_of_angle(delta);
    const double s2beta1 = p2_of_angle(beta1);
    const double s2beta2 = p2_of_angle(beta2);
    return 2.0 / 3.0 * (
            1.0
            + sD2 * s2beta1
            + sA2 * s2beta2
            + sD2 * sA2 * (
                    s2delta
                    + 6.0 * s2beta1 * s2beta2
                    + 1.0
                    + 2.0 * s2beta1
                    + 2.0 * s2beta2
                    - 9.0 * std::cos(beta1) * std::cos(beta2) * std::cos(delta)
            )
    );
}

std::vector<double> dipole_kappa_distance(
        const std::vector<double>& d1, const std::vector<double>& d2,
        const std::vector<double>& a1, const std::vector<double>& a2) {
    double muD[3], muA[3], dM[3], aM[3], rda[3];
    double dD = 0.0, dA = 0.0;
    for (int i = 0; i < 3; ++i) {
        dD += (d1[i] - d2[i]) * (d1[i] - d2[i]);
        dA += (a1[i] - a2[i]) * (a1[i] - a2[i]);
    }
    dD = std::sqrt(dD);
    dA = std::sqrt(dA);
    for (int i = 0; i < 3; ++i) {
        muD[i] = (d2[i] - d1[i]) / dD;
        muA[i] = (a2[i] - a1[i]) / dA;
        dM[i] = d1[i] + dD * muD[i] / 2.0;
        aM[i] = a1[i] + dA * muA[i] / 2.0;
        rda[i] = dM[i] - aM[i];
    }
    const double dRDA = std::sqrt(rda[0] * rda[0] + rda[1] * rda[1] + rda[2] * rda[2]);
    double mu_dot = 0.0, d_dot_n = 0.0, a_dot_n = 0.0;
    for (int i = 0; i < 3; ++i) {
        const double n = rda[i] / dRDA;
        mu_dot += muD[i] * muA[i];
        d_dot_n += muD[i] * n;
        a_dot_n += muA[i] * n;
    }
    return {dRDA, mu_dot - 3.0 * d_dot_n * a_dot_n};
}

void Kappa2Distribution::get_values(double** out_view, int* n_out_view) const {
    internal::copy_to_view(values, out_view, n_out_view);
}

void Kappa2Distribution::get_scale(double** out_view, int* n_out_view) const {
    internal::copy_to_view(scale, out_view, n_out_view);
}

void Kappa2Distribution::get_hist(double** out_view, int* n_out_view) const {
    internal::copy_to_view(hist, out_view, n_out_view);
}

Kappa2Distribution wobbling_kappa2_distribution_delta(
        double delta, double sD2, double sA2, double step,
        int n_bins, double k2_min, double k2_max) {
    Kappa2Distribution out;
    std::vector<double>& k2_scale = out.scale;
    std::vector<double>& k2_hist = out.hist;
    k2_scale = bin_edges(n_bins, k2_min, k2_max);
    k2_hist.assign(k2_scale.empty() ? 0 : k2_scale.size() - 1, 0.0);

    const double d_rad = step * M_PI / 180.0;
    std::vector<double> beta1, phi;
    for (double b = 0.001; b < M_PI / 2.0; b += d_rad) beta1.push_back(b);
    for (double p = 0.001; p < 2.0 * M_PI; p += d_rad) phi.push_back(p);

    std::vector<double> k2(beta1.size() * phi.size(), 0.0);
    const double sd = std::sin(delta), cd = std::cos(delta);
    for (std::size_t i = 0; i < beta1.size(); ++i) {
        const double cb = std::cos(beta1[i]), sb = std::sin(beta1[i]);
        // d1 is the donor axis; n1, n2 span the cone that phi sweeps.
        const double d1[3] = {cb, 0.0, sb};
        const double n1[3] = {-sb, 0.0, cb};
        const double n2[3] = {0.0, 1.0, 0.0};
        const double weight = std::sin(beta1[i]);   // the solid-angle element
        for (std::size_t j = 0; j < phi.size(); ++j) {
            const double cp = std::cos(phi[j]), sp = std::sin(phi[j]);
            // R_DA is along x, so beta2 needs only the x component of d2.
            const double d2x = (n1[0] * cp + n2[0] * sp) * sd + d1[0] * cd;
            const double beta2 = std::acos(std::min(1.0, std::abs(d2x)));
            const double v = wobbling_kappa2(delta, sD2, sA2, beta1[i], beta2);
            k2[i * phi.size() + j] = v;
            accumulate(k2_hist, k2_scale, v, weight);
        }
    }
    out.values.swap(k2);
    return out;
}

Kappa2Distribution wobbling_kappa2_distribution(
        double sD2, double sA2, int n_bins, double k2_min, double k2_max,
        int n_samples, int seed) {
    Kappa2Distribution out;
    std::vector<double>& k2_scale = out.scale;
    std::vector<double>& k2_hist = out.hist;
    k2_scale = bin_edges(n_bins, k2_min, k2_max);
    k2_hist.assign(k2_scale.empty() ? 0 : k2_scale.size() - 1, 0.0);

    std::vector<double> k2(std::max(0, n_samples), 0.0);
    std::mt19937_64 rng(static_cast<std::uint64_t>(seed));
    // Three standard normals, normalised: the only cheap way to hit the sphere
    // uniformly. Normalising three *uniform* variates gives the cube's positive
    // octant instead, which is what the Python here used to do.
    std::normal_distribution<double> gauss(0.0, 1.0);
    for (int i = 0; i < n_samples; ++i) {
        double d1[3], d2[3], n1 = 0.0, n2 = 0.0;
        for (int c = 0; c < 3; ++c) {
            d1[c] = gauss(rng);
            d2[c] = gauss(rng);
            n1 += d1[c] * d1[c];
            n2 += d2[c] * d2[c];
        }
        n1 = std::sqrt(n1);
        n2 = std::sqrt(n2);
        double dot = 0.0;
        for (int c = 0; c < 3; ++c) dot += (d1[c] / n1) * (d2[c] / n2);
        // R_DA is along x by assumption, so beta is the angle to the x axis.
        const double delta = std::acos(std::max(-1.0, std::min(1.0, dot)));
        const double beta1 = std::acos(std::max(-1.0, std::min(1.0, d1[0] / n1)));
        const double beta2 = std::acos(std::max(-1.0, std::min(1.0, d2[0] / n2)));
        const double v = wobbling_kappa2(delta, sD2, sA2, beta1, beta2);
        k2[i] = v;
        accumulate(k2_hist, k2_scale, v, 1.0);
    }
    out.values.swap(k2);
    return out;
}


void sample_kappa2_diffusion_with_traps(
        double sD2, double sA2, double fret_efficiency,
        int n_samples, int n_bins, double k2_min, double k2_max, int seed,
        double** out_view, int* n_out_view) {
    if (n_samples < 0) n_samples = 0;
    if (n_bins < 2) n_bins = 2;
    const std::size_t n_edges = static_cast<std::size_t>(n_bins);
    const std::size_t n_counts = n_edges - 1;
    double* out = internal::new_double_view(
            n_edges + n_counts + static_cast<std::size_t>(n_samples),
            out_view, n_out_view);
    if (out == nullptr) return;

    const double step = (k2_max - k2_min) / (n_bins - 1);
    for (std::size_t i = 0; i < n_edges; ++i) out[i] = k2_min + step * i;
    double* counts = out + n_edges;
    double* samples = counts + n_counts;

    // Three standard normals per direction. Normalising a *uniform* draw fills
    // the cube's positive octant rather than the sphere, and that halved
    // <kappa^2> to 0.333 in the isotropic limit where it must be exactly 2/3.
    std::mt19937_64 rng(seed < 0 ? std::random_device{}()
                                 : static_cast<std::uint64_t>(seed));
    std::normal_distribution<double> gauss(0.0, 1.0);
    const double x = 1.0 / fret_efficiency - 1.0;

    for (int i = 0; i < n_samples; ++i) {
        double d[3], a[3];
        for (int k = 0; k < 3; ++k) { d[k] = gauss(rng); a[k] = gauss(rng); }
        const double nd = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
        const double na = std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
        const double dot = d[0]*a[0] + d[1]*a[1] + d[2]*a[2];
        // R_DA is taken along x, so each beta is the angle to that axis.
        const double delta = std::acos(dot / (nd * na));
        const double beta1 = std::acos(d[0] / nd);
        const double beta2 = std::acos(a[0] / na);

        const double k2_tf = wobbling_kappa2(delta, 1.0, 0.0, beta1, beta2);
        const double k2_ft = wobbling_kappa2(delta, 0.0, 1.0, beta1, beta2);
        const double k2_tt = wobbling_kappa2(delta, 1.0, 1.0, beta1, beta2);

        // The four sub-populations -- donor free or trapped against acceptor
        // free or trapped -- combined into one efficiency, then inverted back
        // into the single kappa^2 that would have produced it.
        const double two_thirds = 2.0 / 3.0;
        const double e = (1.0 - sD2) * (1.0 - sA2) / (1.0 + x)
                       + sD2 * sA2 / (1.0 + two_thirds / k2_tt * x)
                       + sD2 * (1.0 - sA2) / (1.0 + two_thirds / k2_tf * x)
                       + (1.0 - sD2) * sA2 / (1.0 + two_thirds / k2_ft * x);
        samples[i] = two_thirds * x / (1.0 / e - 1.0);
    }

    // numpy's histogram convention: bins are half-open except the last, which
    // is closed, so a sample exactly at k2_max lands in it rather than nowhere.
    for (int i = 0; i < n_samples; ++i) {
        const double v = samples[i];
        if (v < out[0] || v > out[n_edges - 1]) continue;
        std::size_t b = 0;
        while (b + 1 < n_counts && v >= out[b + 1]) ++b;
        counts[b] += 1.0;
    }
}

void kappa2_dipole_matrix(const std::vector<double>& mu_donor,
                         const std::vector<double>& mu_acceptor,
                         const std::vector<double>& r_vectors,
                         double** out_view, int* n_out_view) {
    const std::size_t nd = mu_donor.size() / 3;
    const std::size_t na = mu_acceptor.size() / 3;
    double* out = internal::new_double_view(nd * na, out_view, n_out_view);
    if (out == nullptr) return;
    for (std::size_t i = 0; i < nd; ++i) {
        const double* md = &mu_donor[3 * i];
        for (std::size_t j = 0; j < na; ++j) {
            const double* ma = &mu_acceptor[3 * j];
            const std::size_t base = 3 * (i * na + j);
            if (base + 2 >= r_vectors.size()) continue;
            const double rx = r_vectors[base + 0];
            const double ry = r_vectors[base + 1];
            const double rz = r_vectors[base + 2];
            const double rn = std::sqrt(rx * rx + ry * ry + rz * rz);
            // A zero separation contributes zero rather than a division by
            // zero: coincident states are reachable and are not an error.
            const double ux = rn > 0.0 ? rx / rn : 0.0;
            const double uy = rn > 0.0 ? ry / rn : 0.0;
            const double uz = rn > 0.0 ? rz / rn : 0.0;
            const double cos_da = md[0]*ma[0] + md[1]*ma[1] + md[2]*ma[2];
            const double cos_dr = md[0]*ux + md[1]*uy + md[2]*uz;
            const double cos_ar = ma[0]*ux + ma[1]*uy + ma[2]*uz;
            const double k = cos_da - 3.0 * cos_dr * cos_ar;
            out[i * na + j] = k * k;
        }
    }
}

void isotropic_kappa2_density(const std::vector<double>& k2,
                              double** out_view, int* n_out_view) {
    double* out = internal::new_double_view(k2.size(), out_view, n_out_view);
    if (out == nullptr) return;
    const double s3 = std::sqrt(3.0);
    for (std::size_t i = 0; i < k2.size(); ++i) {
        const double k = std::sqrt(k2[i]);
        if (k >= 0.0 && k <= 1.0) {
            out[i] = 0.5 / (s3 * k) * std::log(2.0 + s3);
        } else if (k > 1.0 && k <= 2.0) {
            out[i] = 0.5 / (s3 * k) * std::log((2.0 + s3) / (k + std::sqrt(k * k - 1.0)));
        } else {
            out[i] = 0.0;
        }
    }
}


void kappa2_distance_ratio_transform(const std::vector<double>& k2_amp,
                              const std::vector<double>& k2_val, int n_bins,
                              double** out_view, int* n_out_view) {
    if (k2_amp.size() != k2_val.size()) {
        IMP_THROW("k2_amp and k2_val must have the same shape.",
                  IMP::ValueException);
    }
    for (std::size_t i = 0; i < k2_val.size(); ++i) {
        if (k2_val[i] <= 0.0) {
            IMP_THROW("k2_val must be strictly positive (kappa^2 > 0).",
                      IMP::ValueException);
        }
        if (k2_amp[i] < 0.0) {
            IMP_THROW("k2_amp must be non-negative.", IMP::ValueException);
        }
    }
    double total = 0.0;
    for (std::size_t i = 0; i < k2_amp.size(); ++i) total += k2_amp[i];
    if (total <= 0.0) {
        IMP_THROW("Total amplitude must be positive.", IMP::ValueException);
    }
    if (n_bins < 1) n_bins = 1;

    const std::size_t n = k2_val.size();
    std::vector<double> w(n), r(n);
    double k2_mean = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        w[i] = k2_amp[i] / total;
        k2_mean += w[i] * k2_val[i];
    }
    // The change of variable, with its Jacobian. Carrying the weights across
    // unchanged would be wrong wherever the mapping is non-linear, which is
    // everywhere: |dk2/dr| = 6 <k2> / r^7.
    double wj_total = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        r[i] = std::pow(k2_mean / k2_val[i], 1.0 / 6.0);
        w[i] *= 6.0 * k2_mean / std::pow(r[i], 7.0);
        wj_total += w[i];
    }
    if (wj_total > 0.0) {
        for (std::size_t i = 0; i < n; ++i) w[i] /= wj_total;
    }

    std::vector<std::size_t> order(n);
    for (std::size_t i = 0; i < n; ++i) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&r](std::size_t a, std::size_t b) { return r[a] < r[b]; });
    std::vector<double> rs(n), ws(n);
    for (std::size_t i = 0; i < n; ++i) { rs[i] = r[order[i]]; ws[i] = w[order[i]]; }

    // A 5 % margin either side, clamped at zero: the transformed points are
    // not uniformly spaced, and an axis that stopped exactly at the extremes
    // would put half the end bins' mass outside it.
    double r_min = rs.front(), r_max = rs.back();
    const double span = r_max - r_min;
    r_min = std::max(0.0, r_min - 0.05 * span);
    r_max = r_max + 0.05 * span;

    double* out = internal::new_double_view(2 * static_cast<std::size_t>(n_bins) + 1,
                                            out_view, n_out_view);
    if (out == nullptr) return;
    double* axis = out;
    double* interp = out + n_bins;
    const double step = n_bins > 1 ? (r_max - r_min) / (n_bins - 1) : 0.0;
    for (int b = 0; b < n_bins; ++b) {
        const double x = n_bins > 1 ? r_min + step * b : r_min;
        axis[b] = x;
        // **Zero** outside the sampled range, not clamped to the end value.
        // The 5 % margin above puts the first and last bins outside it by
        // construction, so clamping instead would plant the extreme weight on
        // a ratio no kappa^2 in the input maps to.
        if (x < rs.front() || x > rs.back()) {
            interp[b] = 0.0;
        } else {
            const std::size_t hi = static_cast<std::size_t>(
                    std::upper_bound(rs.begin(), rs.end(), x) - rs.begin());
            const std::size_t lo = hi - 1;
            const double dx = rs[hi] - rs[lo];
            const double t = dx > 0.0 ? (x - rs[lo]) / dx : 0.0;
            interp[b] = ws[lo] + t * (ws[hi] - ws[lo]);
        }
    }
    // Renormalised *after* interpolation, not before: resampling onto a
    // uniform axis does not preserve the sum, and the caller is handed a
    // distribution.
    double interp_total = 0.0;
    for (int b = 0; b < n_bins; ++b) interp_total += interp[b];
    if (interp_total > 0.0) {
        for (int b = 0; b < n_bins; ++b) interp[b] /= interp_total;
    }
    out[2 * n_bins] = k2_mean;
}

double kappa2_isotropic() { return 2.0 / 3.0; }

std::vector<double> s2_delta_from_anisotropy(double s2_donor,
                                             double s2_acceptor,
                                             double r_inf_AD, double r_0) {
    const double s2_delta = r_inf_AD / (r_0 * s2_donor * s2_acceptor);
    std::vector<double> out(2);
    out[0] = s2_delta;
    out[1] = std::acos(std::sqrt((2.0 * s2_delta + 1.0) / 3.0));
    return out;
}

IMPBFF_END_NAMESPACE
