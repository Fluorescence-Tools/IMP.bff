/**
 * \file OrientationFactor.cpp
 * \brief The FRET orientation factor kappa^2, from geometry and from wobbling.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/OrientationFactor.h>

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

std::vector<double> wobbling_kappa2_distribution_delta(
        double delta, double sD2, double sA2, double step,
        int n_bins, double k2_min, double k2_max,
        std::vector<double>& k2_scale, std::vector<double>& k2_hist) {
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
    return k2;
}

std::vector<double> wobbling_kappa2_distribution(
        double sD2, double sA2, int n_bins, double k2_min, double k2_max,
        int n_samples, int seed,
        std::vector<double>& k2_scale, std::vector<double>& k2_hist) {
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
    return k2;
}

IMPBFF_END_NAMESPACE
