/**
 * \file PolymerChain.cpp
 * \brief End-to-end distance distributions of ideal and worm-like chains.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/PolymerChain.h>
#include <IMP/bff/SpecialFunctions.h>
#include <IMP/bff/Distributions.h>
#include <IMP/bff/internal/Normalize.h>
#include <IMP/bff/internal/OutputView.h>

#include <algorithm>
#include <cmath>
#include <cstring>

IMPBFF_BEGIN_NAMESPACE

//! The ideal-chain kernel as a buffer, for the view-publishing wrapper.
std::vector<double> gaussian_chain_impl(const std::vector<double>& distances,
                                        double segment_length,
                                        int number_of_segments) {
    const double ree = gaussian_chain_ree(segment_length, number_of_segments);
    const double r2_mean = ree * ree;
    std::vector<double> out(distances.size(), 0.0);
    if (r2_mean == 0.0) return out;
    const double denom = std::pow(2.0 / 3.0 * M_PI * r2_mean, 1.5);
    for (std::size_t i = 0; i < distances.size(); ++i) {
        const double r = distances[i];
        out[i] = 4.0 * M_PI * r * r / denom * std::exp(-1.5 * r * r / r2_mean);
    }
    return out;
}

//! The worm-like-chain kernel as a buffer (the wrapper publishes a view).
std::vector<double> worm_like_chain_impl(const std::vector<double>& distances,
                                         double kappa, double chain_length,
                                         bool normalize, bool distance) {
    std::vector<double> pr(distances.size(), 0.0);
    if (distances.empty()) return pr;
    if (chain_length == 0.0) {
        chain_length = *std::max_element(distances.begin(), distances.end());
    }
    if (chain_length == 0.0) return pr;

    const double a = 14.054;
    const double b = 0.473;
    const double c =
            1.0 - std::pow(1.0 + std::pow(0.38 * std::pow(kappa, -0.95), -5.0), -0.2);
    // The branch is the paper's: below kappa = 0.125 the correction is linear,
    // above it the fitted form takes over. Written with `exp(0.783*log(x))`
    // rather than `pow(x, 0.783)`, term for term as the paper writes it.
    const double d = (kappa < 0.125)
            ? kappa + 1.0
            : 1.0 - 1.0 / (0.177 / (kappa - 0.111)
                           + 6.4 * std::exp(0.783 * std::log(kappa - 0.111)));

    // The reference implementation `break`s at the first distance that reaches
    // the contour length, so everything past that index stays zero *whether or
    // not* those distances are themselves shorter. That is a PREFIX, not a
    // mask, and on an unsorted axis the two differ -- nothing requires
    // `distances` to be sorted. Reproduce the prefix.
    std::size_t limit = distances.size();
    for (std::size_t i = 0; i < distances.size(); ++i) {
        if (distances[i] >= chain_length) { limit = i; break; }
    }

    for (std::size_t i = 0; i < limit; ++i) {
        const double r = distances[i];
        const double r_n = r / chain_length;
        const double r_n2 = r_n * r_n;
        double pri = std::pow((1.0 - c * r_n2) / (1.0 - r_n2), 2.5);
        pri *= std::exp(-d * kappa * a * b * (1.0 + b)
                        / (1.0 - (b * r_n) * (b * r_n)) * r_n2);
        const double r_n4 = r_n2 * r_n2;
        const double r_n6 = r_n4 * r_n2;
        const double g = ((-0.75) / kappa - 0.5) * r_n2
                       + ((-0.359375) / kappa + 1.0625) * r_n4
                       + ((-0.109375) / kappa - 0.5625) * r_n6;
        pri *= std::exp(g / (1.0 - r_n2));
        // I0, NOT exp. This line read `std::exp(...)` until 2026-09-02 and the
        // two are nothing alike here: the argument is negative and I0 is even,
        // so I0(-x) GROWS where exp(-x) decays -- a factor of 4e3 at x = -5 and
        // 2e16 at x = -20. The symptom was qualitative rather than subtle: with
        // exp, <R^2> *fell* as the chain stiffened, where Kratky-Porod says it
        // must rise. Nothing outside this file's own tests called it, which is
        // why it survived.
        pri *= i0(-d * kappa * a * (1.0 + b) * r_n
                  / (1.0 - (b * r_n) * (b * r_n)));
        pr[i] = pri;
    }
    if (distance) {
        for (std::size_t i = 0; i < pr.size(); ++i) pr[i] *= distances[i] * distances[i];
    }
    if (normalize) internal::normalize_sum(pr);
    return pr;
}

double gaussian_chain_ree(double segment_length, int number_of_segments) {
    return segment_length * std::sqrt(static_cast<double>(number_of_segments));
}

void gaussian_chain(const std::vector<double>& distances,
                    double segment_length, int number_of_segments,
                    double** out_view, int* n_out_view) {
    internal::copy_to_view(
            gaussian_chain_impl(distances, segment_length, number_of_segments),
            out_view, n_out_view);
}

void worm_like_chain(const std::vector<double>& distances, double kappa,
                     double chain_length, bool normalize, bool distance,
                     double** out_view, int* n_out_view) {
    internal::copy_to_view(
            worm_like_chain_impl(distances, kappa, chain_length, normalize, distance),
            out_view, n_out_view);
}

void worm_like_chain_linker(const std::vector<double>& distances, double kappa,
                            double chain_length, double sigma, bool normalize,
                            double** out_view, int* n_out_view) {
    std::vector<double> pn(distances.size(), 0.0);
    if (sigma != 0.0) {
        // `distance` is false here: the reference broadens the bare chain
        // distribution, and applying r^2 first would weight the convolution.
        const std::vector<double> pr =
                worm_like_chain_impl(distances, kappa, chain_length, normalize, false);
        // The broadening kernel is distance_between_gaussian -- the distribution
        // of the distance between two points that are each Gaussian-distributed
        // -- NOT a plain normal density. This read `normal_density` until
        // 2026-09-02; a normal is the density of an offset, and what a linker
        // pair contributes is the distance between two clouds, which carries the
        // extra r/separation factor and the antisymmetric second term.
        for (std::size_t i = 0; i < distances.size(); ++i) {
            if (pr[i] == 0.0) continue;
            const std::vector<double> broad =
                    distance_between_gaussian_impl(distances, distances[i], sigma);
            for (std::size_t j = 0; j < pn.size(); ++j) pn[j] += pr[i] * broad[j];
        }
        if (normalize) internal::normalize_sum(pn);
    }
    internal::copy_to_view(pn, out_view, n_out_view);
}


void ising_chain(const std::vector<double>& distances, int number_of_residues,
                 double b_structured, double b_unstructured, double coupling,
                 double field, int n_k, double** out_view, int* n_out_view) {
    const std::size_t n_r = distances.size();
    std::vector<double> pr(n_r, 0.0);
    if (number_of_residues < 1 || n_r == 0 || n_k < 2) {
        internal::copy_to_view(pr, out_view, n_out_view);
        return;
    }
    const int n = number_of_residues;

    // Per-residue bond variance in the characteristic function, b^2/6.
    const double vS = b_structured * b_structured / 6.0;
    const double vU = b_unstructured * b_unstructured / 6.0;

    // Ising nearest-neighbour weights W(sigma, sigma'), states 0 = S, 1 = U,
    // for the energy -J*delta(sigma,sigma') - (h/2)*(is_S(sigma)+is_S(sigma')).
    const double w00 = std::exp(coupling + field);
    const double w01 = std::exp(-coupling + 0.5 * field);
    const double w10 = w01;
    const double w11 = std::exp(coupling);

    // The k grid runs to where phi has decayed, which the *smallest* bond
    // variance sets -- the stiffest state is the one still oscillating when
    // the other has died away.
    const double r_max = *std::max_element(distances.begin(), distances.end());
    const double scale = std::max(std::max(std::sqrt(std::min(vS, vU) * n),
                                           r_max / n), 1e-6);
    const double k_min = 1e-6;
    const double k_max = 30.0 / scale;
    const double step = (k_max - k_min) / (n_k - 1);

    std::vector<double> k(n_k), phi(n_k);
    for (int j = 0; j < n_k; ++j) k[j] = k_min + step * j;
    k[n_k - 1] = k_max;                      // linspace pins its own endpoint

    for (int j = 0; j < n_k; ++j) {
        const double kk = k[j];
        const double g0 = std::exp(-kk * kk * vS);
        const double g1 = std::exp(-kk * kk * vU);
        // M[s][s'] = W[s][s'] * g[s'], and v is a row vector: v <- v M.
        const double m00 = w00 * g0, m01 = w01 * g1;
        const double m10 = w10 * g0, m11 = w11 * g1;
        double v0 = 1.0, v1 = 1.0;
        for (int i = 0; i < n; ++i) {
            const double t0 = v0 * m00 + v1 * m10;
            const double t1 = v0 * m01 + v1 * m11;
            v0 = t0;
            v1 = t1;
        }
        phi[j] = v0 + v1;
    }
    const double phi0 = phi[0];
    if (phi0 != 0.0) {
        for (int j = 0; j < n_k; ++j) phi[j] /= phi0;
    }

    // P(R) = (2R/pi) * trapz_k[ k sin(kR) phi(k) ].
    for (std::size_t i = 0; i < n_r; ++i) {
        const double r = distances[i];
        double integral = 0.0;
        double f_prev = k[0] * std::sin(k[0] * r) * phi[0];
        for (int j = 1; j < n_k; ++j) {
            const double f = k[j] * std::sin(k[j] * r) * phi[j];
            integral += 0.5 * (f + f_prev) * (k[j] - k[j - 1]);
            f_prev = f;
        }
        const double value = 2.0 * r / M_PI * integral;
        // A non-finite phi (see the header on overflow) must not poison the
        // normalisation; a negative lobe of the transform is a truncation
        // artefact, not a probability.
        pr[i] = (std::isfinite(value) && value > 0.0) ? value : 0.0;
    }

    double area = 0.0;
    for (std::size_t i = 1; i < n_r; ++i) {
        area += 0.5 * (pr[i] + pr[i - 1]) * (distances[i] - distances[i - 1]);
    }
    if (area > 0.0) {
        for (std::size_t i = 0; i < n_r; ++i) pr[i] /= area;
    }
    internal::copy_to_view(pr, out_view, n_out_view);
}


void saw_nu(const std::vector<double>& distances, double r_rms, double nu,
            double gamma_exp, double** out_view, int* n_out_view) {
    std::vector<double> pr(distances.size(), 0.0);
    if (!(nu > 0.0 && nu < 1.0) || r_rms <= 0.0) {
        internal::copy_to_view(pr, out_view, n_out_view);
        return;
    }
    const double theta = (gamma_exp - 1.0) / nu;
    const double delta = 1.0 / (1.0 - nu);
    // <r^2> = r0^2 Gamma((5+theta)/delta) / Gamma((3+theta)/delta), which fixes
    // r0 from the requested RMS distance.
    const double ratio = std::tgamma((5.0 + theta) / delta)
                       / std::tgamma((3.0 + theta) / delta);
    const double r0 = r_rms / std::sqrt(ratio);
    const double norm = delta / (std::pow(r0, 3.0 + theta)
                                 * std::tgamma((3.0 + theta) / delta));
    for (std::size_t i = 0; i < distances.size(); ++i) {
        const double r = distances[i];
        const double value = norm * std::pow(r, 2.0 + theta)
                           * std::exp(-std::pow(r / r0, delta));
        // The reference squashes non-finite values to zero rather than letting
        // an overflow at large r/r0 propagate; delta > 1 makes that reachable.
        pr[i] = std::isfinite(value) ? value : 0.0;
    }
    internal::copy_to_view(pr, out_view, n_out_view);
}

IMPBFF_END_NAMESPACE
