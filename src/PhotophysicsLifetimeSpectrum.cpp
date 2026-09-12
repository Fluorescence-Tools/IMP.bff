/**
 * \file PhotophysicsLifetimeSpectrum.cpp
 * \brief The experiment-neutral output: (amplitude, rate constant) pairs.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/PhotophysicsLifetimeSpectrum.h>
#include <IMP/bff/IMPCompatibility.h>
#include <cstdlib>
#include <limits>
#include <IMP/bff/internal/OutputView.h>

#include <algorithm>
#include <cmath>

IMPBFF_BEGIN_NAMESPACE

namespace {
std::vector<double> lifetime_spectrum_decay_impl(
        const std::vector<double>& amplitudes,
        const std::vector<double>& rate_constants,
        const std::vector<double>& time) {
    std::vector<double> out(time.size(), 0.0);
    const std::size_t n = std::min(amplitudes.size(), rate_constants.size());
#pragma omp parallel for schedule(static)
    for (long long j = 0; j < static_cast<long long>(time.size()); ++j) {
        const double t = time[static_cast<std::size_t>(j)];
        double f = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            f += amplitudes[i] * std::exp(-rate_constants[i] * t);
        }
        out[static_cast<std::size_t>(j)] = f;
    }
    return out;
}
}  // namespace

void lifetime_spectrum_decay(const std::vector<double>& amplitudes,
        const std::vector<double>& rate_constants,
        const std::vector<double>& time, double** out_view, int* n_out_view) {
    internal::copy_to_view(lifetime_spectrum_decay_impl(amplitudes, rate_constants, time),
                           out_view, n_out_view);
}

std::vector<double> lifetime_spectrum_coarse_grain(
        const std::vector<double>& amplitudes,
        const std::vector<double>& rate_constants,
        int n_bins, std::vector<double>& out_rates) {
    out_rates.clear();
    std::vector<double> out_amplitudes;
    const std::size_t n = std::min(amplitudes.size(), rate_constants.size());
    if (n == 0 || n_bins < 1) return out_amplitudes;

    double lo = rate_constants[0], hi = rate_constants[0];
    for (std::size_t i = 1; i < n; ++i) {
        lo = std::min(lo, rate_constants[i]);
        hi = std::max(hi, rate_constants[i]);
    }
    // Every species shares one rate: there is nothing to reduce, and a zero
    // width would divide by zero below.
    if (!(hi > lo)) {
        double a_sum = 0.0;
        for (std::size_t i = 0; i < n; ++i) a_sum += amplitudes[i];
        out_rates.push_back(lo);
        out_amplitudes.push_back(a_sum);
        return out_amplitudes;
    }

    const double width = (hi - lo) / n_bins;
    std::vector<double> a_sum(n_bins, 0.0), ak_sum(n_bins, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t b = static_cast<std::size_t>((rate_constants[i] - lo) / width);
        if (b >= static_cast<std::size_t>(n_bins)) b = n_bins - 1;  // closed right edge
        a_sum[b] += amplitudes[i];
        ak_sum[b] += amplitudes[i] * rate_constants[i];
    }
    for (int b = 0; b < n_bins; ++b) {
        if (a_sum[b] == 0.0) continue;   // an empty bin is not a species
        out_amplitudes.push_back(a_sum[b]);
        // The amplitude-weighted mean, so sum(a) and sum(a*k) both survive.
        out_rates.push_back(ak_sum[b] / a_sum[b]);
    }
    return out_amplitudes;
}

std::vector<double> outer_product_histogram(
        const std::vector<double>& a, const std::vector<double>& weights_a,
        const std::vector<double>& b, const std::vector<double>& weights_b,
        int n_bins, double lo, double hi) {
    std::vector<double> hist(n_bins > 0 ? n_bins : 0, 0.0);
    if (n_bins <= 0 || !(hi > lo)) return hist;
    const std::size_t na = std::min(a.size(), weights_a.size());
    const std::size_t nb = std::min(b.size(), weights_b.size());
    const double inv_width = n_bins / (hi - lo);
    for (std::size_t i = 0; i < na; ++i) {
        const double ai = a[i], wi = weights_a[i];
        for (std::size_t j = 0; j < nb; ++j) {
            const double v = ai * b[j];
            if (v < lo || v > hi) continue;
            std::size_t k = static_cast<std::size_t>((v - lo) * inv_width);
            if (k >= hist.size()) k = hist.size() - 1;   // the closed top edge
            hist[k] += wi * weights_b[j];
        }
    }
    return hist;
}

void convolve_distance_with_k2_ratio(
        const std::vector<double>& r_da, const std::vector<double>& amp,
        const std::vector<double>& r_ratio,
        const std::vector<double>& weights_ratio, int n_bins,
        double** out_centres, int* n_centres, double** out_hist, int* n_hist) {
    if (r_da.empty() || r_ratio.empty() || n_bins <= 0) {
        internal::new_double_view(0, out_centres, n_centres);
        internal::new_double_view(0, out_hist, n_hist);
        return;
    }
    double lo = *std::min_element(r_da.begin(), r_da.end()) *
                *std::min_element(r_ratio.begin(), r_ratio.end());
    double hi = *std::max_element(r_da.begin(), r_da.end()) *
                *std::max_element(r_ratio.begin(), r_ratio.end());
    if (!(hi > lo)) {   // every product identical: a delta, like numpy's range
        lo -= 0.5;
        hi += 0.5;
    }
    const std::vector<double> hist = outer_product_histogram(
            r_da, amp, r_ratio, weights_ratio, n_bins, lo, hi);
    double peak = 0.0;
    for (std::size_t i = 0; i < hist.size(); ++i) peak = std::max(peak, hist[i]);

    double* centres = internal::new_double_view(n_bins, out_centres, n_centres);
    double* kept = internal::new_double_view(n_bins, out_hist, n_hist);
    const double step = (hi - lo) / n_bins;
    int n_kept = 0;
    for (int i = 0; i < n_bins; ++i) {
        if (hist[i] <= 1e-10 * peak) continue;
        centres[n_kept] = lo + (i + 0.5) * step;
        kept[n_kept] = hist[i];
        ++n_kept;
    }
    // The views were sized for n_bins; shrink the published length to n_kept.
    if (out_centres != NULL) *n_centres = n_kept;
    if (out_hist != NULL) *n_hist = n_kept;
}


// --------------------------------------------------------------------------
// PhotophysicsLifetimeSpectrum
// --------------------------------------------------------------------------

PhotophysicsLifetimeSpectrum::PhotophysicsLifetimeSpectrum(const std::vector<double>& amplitudes,
                                   const std::vector<double>& rate_constants,
                                   bool exact)
        : amplitudes_(amplitudes), rate_constants_(rate_constants), exact_(exact) {
    if (amplitudes_.size() != rate_constants_.size()) {
        IMP_THROW("a spectrum needs one rate per amplitude: "
                          << amplitudes_.size() << " against "
                          << rate_constants_.size(),
                  IMP::ValueException);
    }
    for (std::size_t i = 0; i < rate_constants_.size(); ++i) {
        if (rate_constants_[i] < 0.0) {
            IMP_THROW("a negative rate constant is a growing population, "
                      "not a decaying one", IMP::ValueException);
        }
    }
}

void PhotophysicsLifetimeSpectrum::get_amplitudes(double** out_view, int* n_out_view) const {
    internal::copy_to_view(amplitudes_, out_view, n_out_view);
}

void PhotophysicsLifetimeSpectrum::get_rate_constants(double** out_view, int* n_out_view) const {
    internal::copy_to_view(rate_constants_, out_view, n_out_view);
}

void PhotophysicsLifetimeSpectrum::get_lifetimes(double** out_view, int* n_out_view) const {
    const std::size_t n = rate_constants_.size();
    double* out = internal::new_double_view(n, out_view, n_out_view);
    if (out == nullptr) return;
    for (std::size_t i = 0; i < n; ++i) {
        // A zero rate is a species that never decays, and its lifetime is
        // infinite rather than undefined. Dividing would give the same value
        // with a signal that says "error"; this says "does not decay".
        out[i] = rate_constants_[i] > 0.0
                         ? 1.0 / rate_constants_[i]
                         : std::numeric_limits<double>::infinity();
    }
}

double PhotophysicsLifetimeSpectrum::get_total_amplitude() const {
    double s = 0.0;
    for (std::size_t i = 0; i < amplitudes_.size(); ++i) s += amplitudes_[i];
    return s;
}

double PhotophysicsLifetimeSpectrum::get_species_averaged_lifetime() const {
    const double total = get_total_amplitude();
    if (total == 0.0) return 0.0;
    double acc = 0.0;
    for (std::size_t i = 0; i < amplitudes_.size(); ++i) {
        // One non-decaying species makes the whole average infinite -- it holds
        // population forever. Skipping it instead would report the average of
        // the species that *do* decay, under a name that says otherwise.
        if (rate_constants_[i] <= 0.0) {
            return std::numeric_limits<double>::infinity();
        }
        acc += amplitudes_[i] / rate_constants_[i];
    }
    return acc / total;
}

double PhotophysicsLifetimeSpectrum::get_intensity_averaged_lifetime() const {
    double num = 0.0, den = 0.0;
    for (std::size_t i = 0; i < amplitudes_.size(); ++i) {
        if (rate_constants_[i] <= 0.0) {
            return std::numeric_limits<double>::infinity();
        }
        const double tau = 1.0 / rate_constants_[i];
        num += amplitudes_[i] * tau * tau;
        den += amplitudes_[i] * tau;
    }
    if (den == 0.0) return 0.0;
    return num / den;
}

PhotophysicsLifetimeSpectrum PhotophysicsLifetimeSpectrum::normalized() const {
    const double total = get_total_amplitude();
    if (total == 0.0) return *this;
    std::vector<double> a(amplitudes_.size());
    for (std::size_t i = 0; i < a.size(); ++i) a[i] = amplitudes_[i] / total;
    return PhotophysicsLifetimeSpectrum(a, rate_constants_, exact_);
}

void PhotophysicsLifetimeSpectrum::decay(const std::vector<double>& time,
                             double** out_view, int* n_out_view) const {
    lifetime_spectrum_decay(amplitudes_, rate_constants_, time,
                            out_view, n_out_view);
}

PhotophysicsLifetimeSpectrum PhotophysicsLifetimeSpectrum::coarse_grain(int n_bins) const {
    std::vector<double> rates;
    std::vector<double> amps = lifetime_spectrum_coarse_grain(
            amplitudes_, rate_constants_, n_bins, rates);
    // Exact only if nothing was actually merged. The reduction preserves the
    // population and the initial slope for any n_bins, but not the curvature,
    // so a spectrum that lost species is an approximation and says so.
    const bool still_exact =
            exact_ && amps.size() == amplitudes_.size();
    return PhotophysicsLifetimeSpectrum(amps, rates, still_exact);
}

double fret_efficiency_from_lifetimes(const PhotophysicsLifetimeSpectrum& donor_only,
                                      const PhotophysicsLifetimeSpectrum& donor_acceptor) {
    const double tau_d = donor_only.get_species_averaged_lifetime();
    if (tau_d <= 0.0 || !std::isfinite(tau_d)) return 0.0;
    return 1.0 - donor_acceptor.get_species_averaged_lifetime() / tau_d;
}

PhotophysicsLifetimeSpectrum lifetime_spectrum_from_rates(const std::vector<double>& rates,
                                              const std::vector<double>& weights,
                                              bool exact) {
    // An empty weight vector is "uniform", which is a different thing from a
    // length mismatch: the caller that omits weights means one per species.
    if (weights.empty()) {
        return PhotophysicsLifetimeSpectrum(std::vector<double>(rates.size(), 1.0), rates,
                                exact);
    }
    if (weights.size() != rates.size()) {
        IMP_THROW("one weight per rate: " << weights.size() << " against "
                                          << rates.size(), ValueException);
    }
    return PhotophysicsLifetimeSpectrum(weights, rates, exact);
}

PhotophysicsLifetimeSpectrum lifetime_spectrum_from_states(
        const InteractionTerms& terms, const States& first,
        const States& second, const std::vector<double>& weights, bool exact) {
    const std::vector<double> k = total_rate(terms, first, second);
    if (!weights.empty()) return lifetime_spectrum_from_rates(k, weights, exact);

    // The first participant's own weights are the populations: an accessible
    // volume's occupancy already is one, and a rotamer ensemble's are its
    // Boltzmann weights.
    double* points = nullptr;
    int n_points = 0;
    first.get_points(&points, &n_points);
    std::vector<double> w;
    w.reserve(static_cast<std::size_t>(n_points) / 4);
    for (int i = 3; i < n_points; i += 4) w.push_back(points[i]);
    std::free(points);
    return lifetime_spectrum_from_rates(k, w, exact);
}

void interleaved_scale_amplitudes(const std::vector<double>& spectrum,
                                  double n, std::vector<double>* out) {
  const std::size_t pairs = spectrum.size() / 2;
  for (std::size_t i = 0; i < pairs; ++i) {
    out->push_back(spectrum[2 * i] * n);
    out->push_back(spectrum[2 * i + 1]);
  }
}

std::vector<double> interleaved_product(const std::vector<double>& first,
                                        const std::vector<double>& second) {
  const std::size_t n1 = first.size() / 2, n2 = second.size() / 2;
  std::vector<double> out;
  out.reserve(2 * n1 * n2);
  for (std::size_t i = 0; i < n1; ++i) {
    for (std::size_t j = 0; j < n2; ++j) {
      out.push_back(first[2 * i] * second[2 * j]);
      // 1 / (1/t1 + 1/t2). Written this way, and not as the algebraically
      // equal t1 t2 / (t1 + t2), because that spelling is what ChiSurf
      // evaluates and the two do not round identically.
      out.push_back(1.0 / (1.0 / first[2 * i + 1] + 1.0 / second[2 * j + 1]));
    }
  }
  return out;
}

IMPBFF_END_NAMESPACE
