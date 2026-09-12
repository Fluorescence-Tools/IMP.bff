/**
 * \file TcspcDecay.cpp
 * \brief A time-correlated single-photon-counting decay, as a node.
 *
 * The kernels are tttrlib's, from the vendored copy of `DecayConvolution.h`:
 * `fconv_per_cs_ad<double>` for the periodic reconvolution and
 * `shift_lamp_ad<double>` for the timeshift. What is written here is the
 * *composition* -- the order the instrument model applies its terms in --
 * and the graph around it.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/TcspcDecay.h>

// The kernels. A byte-identical copy of tttrlib's
// `modules/spectroscopy/decay/include/DecayConvolution.h`, kept in step by
// `test/decay/test_decay_convolution_copy_is_identical.py` -- the same
// arrangement, and for the same reasons, as the expression engine (see
// `src/standalone/GraphExpression.cpp`). Only the header-only pieces are used;
// nothing here links tttrlib.
#include <IMP/bff/internal/DecayConvolution.h>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! ChiSurf's `rescale_w_bg`, which is **not** the photon library's.
/*! Same least-squares solution, three deliberate differences that change the
    answer: this one guards on `e > 0` *and* a finite weight (an empty channel
    can carry an infinite weight, which would otherwise poison the whole sum),
    adds no epsilon to the squared weight, and returns the factor rather than
    rescaling the model as a side effect. Reproduced rather than unified,
    because unifying them moves fitted amplitudes and that is a decision about
    the fit, not about where the code lives. */
double rescale_factor(const std::vector<double>& model,
                      const std::vector<double>& data,
                      const std::vector<double>& errors, double background,
                      int begin, int end) {
  double sum_nom = 0.0, sum_denom = 0.0;
  for (int i = begin; i < end; ++i) {
    const double e = data[static_cast<std::size_t>(i)];
    if (!(e > 0.0)) continue;
    const double error = errors[static_cast<std::size_t>(i)];
    // The weight is 1/error; a zero error is an infinite weight, which is
    // what the `isfinite` guard on the Python side removes.
    if (!(error > 0.0) || !std::isfinite(error)) continue;
    const double w2 = 1.0 / (error * error);
    if (!std::isfinite(w2)) continue;
    const double m = model[static_cast<std::size_t>(i)];
    sum_nom += m * (e - background) * w2;
    sum_denom += m * m * w2;
  }
  if (sum_denom == 0.0) return 0.0;
  return sum_nom / sum_denom;
}

}  // namespace

// Deliberately empty. A `GraphNode` that owns ports has to be owned by a
// `shared_ptr` first -- `add_port` reaches for `shared_from_this()` -- so a
// constructor cannot create any, and every Python-wrapped node is
// constructed before it is owned. The ports are therefore made by
// set_number_of_lifetimes(), which is the one call a caller cannot skip.
TcspcDecay::TcspcDecay(const std::string& name) : GraphNode(name) {}

void TcspcDecay::add_scalar_port(const std::string& key, double value,
                                 GraphPort** slot) {
  std::shared_ptr<GraphPort> port(new GraphPort(value));
  add_input_port(key, port);
  *slot = port.get();
}

void TcspcDecay::set_number_of_lifetimes(int n) {
  if (n < 0) {
    throw std::domain_error(
        "TcspcDecay::set_number_of_lifetimes: a negative component count");
  }
  // Ports cannot be removed from a node, so calling this with a *smaller*
  // count leaves the surplus `a`/`t` ports in the map. They feed nothing --
  // only the ports in `lifetime_ports_` are read -- but they are still
  // visible, so build the node once rather than reshaping it.
  if (scatter_port_ == nullptr) {
    // The interleaved spectrum, for a caller that computes it upstream. A
    // vector port, and one that does *not* sanitise: it is fit transport,
    // and a NaN lifetime has to reach `ChiSquared` and make the misfit
    // infinite rather than be floored to `tiny` and read as a good fit.
    std::shared_ptr<GraphPort> spectrum(new GraphPort(std::vector<double>(1, 0.0)));
    spectrum->set_sanitize(false);
    add_input_port(spectrum_port_key(), spectrum);
    spectrum_port_ = spectrum.get();

    add_scalar_port("scatter", 0.0, &scatter_port_);
    add_scalar_port("background", 0.0, &background_port_);
    add_scalar_port("n0", 1.0, &n0_port_);
    add_scalar_port("timeshift", 0.0, &timeshift_port_);
  }
  lifetime_ports_.clear();
  lifetime_ports_.reserve(static_cast<std::size_t>(2 * n));
  for (int i = 0; i < n; ++i) {
    std::ostringstream a, t;
    a << "a" << i;
    t << "t" << i;
    GraphPort* amplitude = nullptr;
    GraphPort* lifetime = nullptr;
    add_scalar_port(a.str(), 1.0, &amplitude);
    add_scalar_port(t.str(), 1.0, &lifetime);
    lifetime_ports_.push_back(amplitude);
    lifetime_ports_.push_back(lifetime);
  }
  n_lifetimes_ = n;
  spectrum_.assign(static_cast<std::size_t>(2 * n), 0.0);
  set_valid(false);
}

bool TcspcDecay::get_basis_is_jacobian() const {
  // Both of these put the amplitudes through something before the
  // contraction, so the curve depends on each of them by a route the basis
  // does not carry.
  if (normalize_amplitudes_ || autoscale_) return false;
  // |a| is differentiable away from zero and its derivative is the sign, so
  // the basis is the Jacobian up to a per-species sign -- and only a caller
  // who knows that can use it. The spectrum here is the one the last
  // evaluation built, which is where the signs are.
  //
  // The sign cannot be read back out of `spectrum_`: `fabs` is applied in
  // place, so by the time anything can look, every amplitude there is
  // non-negative whatever the caller wrote. `build_spectrum` records it.
  if (absolute_amplitudes_ && saw_negative_amplitude_) return false;
  return true;
}

void TcspcDecay::set_spectrum_from_port(bool v) {
  if (v && spectrum_port_ == nullptr) {
    throw std::domain_error(
        "TcspcDecay::set_spectrum_from_port: the ports do not exist yet; "
        "call set_number_of_lifetimes() first, which is what builds them");
  }
  spectrum_from_port_ = v;
  set_valid(false);
}

void TcspcDecay::set_response(const std::vector<double>& response) {
  if (response.empty()) {
    throw std::domain_error("TcspcDecay::set_response: the response is empty");
  }
  response_ = response;
  shifted_.assign(response_.size(), 0.0);
  ++response_epoch_;
  irf_valid_ = false;
  set_valid(false);
}

void TcspcDecay::set_data(const std::vector<double>& y,
                          const std::vector<double>& ey) {
  if (y.size() != ey.size()) {
    std::ostringstream m;
    m << "TcspcDecay::set_data: " << y.size() << " values and " << ey.size()
      << " errors";
    throw std::domain_error(m.str());
  }
  data_y_ = y;
  data_ey_ = ey;
  set_valid(false);
}

void TcspcDecay::set_timing(double dt, double period) {
  if (!(dt > 0.0)) {
    throw std::domain_error("TcspcDecay::set_timing: the channel width is not "
                            "positive");
  }
  if (!(period > 0.0)) {
    throw std::domain_error("TcspcDecay::set_timing: the excitation period is "
                            "not positive");
  }
  dt_ = dt;
  period_ = period;
  set_valid(false);
}

void TcspcDecay::set_convolution_range(int convolution_stop, int stop) {
  convolution_stop_ = convolution_stop;
  stop_ = stop;
  set_valid(false);
}

void TcspcDecay::set_scale_range(int start, int stop) {
  scale_start_ = start;
  scale_stop_ = stop;
  set_valid(false);
}

void TcspcDecay::set_linearization(const std::vector<double>& table) {
  lin_table_ = table;
  set_valid(false);
}

void TcspcDecay::set_linearization_array(double* in_table, int n_table) {
  set_linearization(std::vector<double>(in_table, in_table + n_table));
}

void TcspcDecay::set_response_array(double* in_response, int n_response) {
  set_response(std::vector<double>(in_response, in_response + n_response));
}

void TcspcDecay::set_data_arrays(double* in_data_y, int n_data_y,
                                 double* in_data_ey, int n_data_ey) {
  set_data(std::vector<double>(in_data_y, in_data_y + n_data_y),
           std::vector<double>(in_data_ey, in_data_ey + n_data_ey));
}

void TcspcDecay::build_spectrum() {
  if (spectrum_from_port_) {
    // Whatever the upstream node wrote, in the interleaved layout. An odd
    // length is a caller error rather than something to round down: the last
    // amplitude would silently lose its lifetime.
    const std::vector<double>& incoming = spectrum_port_->get_values_ref();
    if (incoming.size() < 2 || incoming.size() % 2 != 0) {
      std::ostringstream m;
      m << "TcspcDecay '" << get_name() << "': the lifetime spectrum has "
        << incoming.size()
        << " entries; it is interleaved (a0, t0, a1, t1, ...) so the count "
           "is even and at least two";
      throw std::domain_error(m.str());
    }
    spectrum_ = incoming;
    n_lifetimes_ = static_cast<int>(incoming.size() / 2);
  } else {
    spectrum_.resize(static_cast<std::size_t>(2 * n_lifetimes_));
    for (int i = 0; i < n_lifetimes_; ++i) {
      spectrum_[static_cast<std::size_t>(2 * i)] =
          lifetime_ports_[static_cast<std::size_t>(2 * i)]->get_value();
      spectrum_[static_cast<std::size_t>(2 * i + 1)] =
          lifetime_ports_[static_cast<std::size_t>(2 * i + 1)]->get_value();
    }
  }

  // `absolute_amplitudes` and `normalize_amplitudes` describe what the model
  // means by an amplitude, so they apply whichever way the pairs arrived.
  double sum = 0.0;
  saw_negative_amplitude_ = false;
  for (int i = 0; i < n_lifetimes_; ++i) {
    double amplitude = spectrum_[static_cast<std::size_t>(2 * i)];
    if (amplitude < 0.0) saw_negative_amplitude_ = true;
    // `fabs`, where ChiSurf writes `sqrt(a**2)`: the same number for every
    // finite amplitude, and it does not overflow on the way.
    if (absolute_amplitudes_) amplitude = std::fabs(amplitude);
    spectrum_[static_cast<std::size_t>(2 * i)] = amplitude;
    sum += amplitude;
    // Lifetimes are positive: a fit that walks one through zero would
    // otherwise put a growing exponential in the model. ChiSurf takes the
    // absolute value in the getter and writes it back; the write-back is the
    // application's business, the abs() is the model's.
    spectrum_[static_cast<std::size_t>(2 * i + 1)] =
        std::fabs(spectrum_[static_cast<std::size_t>(2 * i + 1)]);
  }
  if (normalize_amplitudes_) {
    const double scale = std::fabs(sum);
    for (int i = 0; i < n_lifetimes_; ++i) {
      spectrum_[static_cast<std::size_t>(2 * i)] /= scale;
    }
  }

  // Drop the species that cannot pay for themselves. The reconvolution below
  // is one serial recursion over every channel *per species*, so this is the
  // only knob on a decay whose spectrum came from a distribution: those have
  // as many species as the distribution has bins, whatever their weight.
  //
  // Compacted in place and `n_lifetimes_` moved down with it, because that
  // count is what the kernel is handed. At the default threshold of zero
  // this drops only exact zeros and is bit-exact -- and it still fires,
  // because a donor-only fraction of zero contributes a whole block of them.
  std::size_t kept = 0;
  double largest = 0.0;
  for (int i = 0; i < n_lifetimes_; ++i) {
    largest = std::max(largest,
                       std::fabs(spectrum_[static_cast<std::size_t>(2 * i)]));
  }
  // A relative threshold against nothing is meaningless: if every amplitude
  // is zero the spectrum is empty, and an empty spectrum is not something
  // the kernel can be handed. Keep it as it is and let the curve be zero.
  // Compaction is skipped when the basis is wanted, so that its columns line
  // up with the input spectrum one for one (see set_emit_basis).
  if (largest > 0.0 && !emit_basis_) {
    const double cutoff = amplitude_threshold_ * largest;
    for (int i = 0; i < n_lifetimes_; ++i) {
      const double amplitude = spectrum_[static_cast<std::size_t>(2 * i)];
      // `>` and not `>=`, so a threshold of exactly zero keeps everything a
      // non-zero amplitude contributes and drops only true zeros.
      if (std::fabs(amplitude) > cutoff) {
        spectrum_[2 * kept] = amplitude;
        spectrum_[2 * kept + 1] =
            spectrum_[static_cast<std::size_t>(2 * i + 1)];
        ++kept;
      }
    }
    // Every amplitude at or below the cutoff, with a positive largest: only
    // reachable when the threshold is 1 or more. Keep the largest species so
    // there is still a decay rather than an empty spectrum.
    if (kept == 0) {
      std::size_t best = 0;
      for (int i = 0; i < n_lifetimes_; ++i) {
        if (std::fabs(spectrum_[static_cast<std::size_t>(2 * i)]) == largest) {
          best = static_cast<std::size_t>(i);
          break;
        }
      }
      spectrum_[0] = spectrum_[2 * best];
      spectrum_[1] = spectrum_[2 * best + 1];
      kept = 1;
    }
    spectrum_.resize(2 * kept);
  }
  n_active_ = static_cast<int>(spectrum_.size() / 2);
}

void TcspcDecay::set_emit_basis(bool on) {
  if (on == emit_basis_) return;
  emit_basis_ = on;
  if (on && !get_output_port(basis_port_key())) {
    add_output_port(basis_port_key(), std::make_shared<GraphPort>(
                        std::vector<double>{0.0}));
  }
  set_valid(false);
}

void TcspcDecay::set_amplitude_threshold(double relative) {
  if (!(relative >= 0.0)) {
    throw std::domain_error(
        "TcspcDecay::set_amplitude_threshold: the threshold is relative to "
        "the largest amplitude and cannot be negative");
  }
  amplitude_threshold_ = relative;
  set_valid(false);
}

void TcspcDecay::evaluate() {
  if (response_.empty()) {
    throw std::domain_error("TcspcDecay '" + get_name() +
                            "' has no response function");
  }
  // Only when the pairs come from the scalar ports: reading them from the
  // spectrum port is how a node with *no* `a`/`t` ports drives this one, and
  // `build_spectrum` sets the count from what actually arrived.
  if (!spectrum_from_port_ && n_lifetimes_ <= 0) {
    throw std::domain_error("TcspcDecay '" + get_name() +
                            "' has no lifetime components");
  }
  build_spectrum();

  const int n_points = static_cast<int>(response_.size());

  // The shift, then the unit-sum normalisation -- that order, because the
  // shift zeroes what moves past either end and normalising first would
  // leave the model's area depending on the timeshift.
  const double timeshift = timeshift_port_->get_value();
  // The shift and the renormalisation depend on the response and the
  // timeshift and on nothing else, so an evaluation that moved neither --
  // every Jacobian column on an amplitude or a lifetime, which is most of
  // them -- can keep the copy it already has. A NaN timeshift compares
  // unequal to itself and rebuilds, which is the safe direction.
  const bool irf_is_current =
      irf_valid_ && irf_epoch_ == response_epoch_ &&
      irf_timeshift_ == timeshift &&
      irf_.size() == static_cast<std::size_t>(n_points);
  if (!irf_is_current) {
  const double* response = response_.data();
  if (timeshift != 0.0) {
    // Sign: ChiSurf's `shift_array(v, s)` is tttrlib's `shift_lamp(v, -s)`.
    // The two index in opposite directions *and* interpolate toward opposite
    // neighbours, and the two flips cancel exactly -- verified across
    // integer and fractional shifts of both signs. At `s == 0` they do
    // differ (`shift_lamp` drops the last sample), which is why this branch
    // exists rather than an unconditional call.
    shift_lamp_ad<double>(shifted_.data(), response_.data(), -timeshift,
                          n_points, 0.0);
    response = shifted_.data();
  }

  // Into a member buffer rather than a local: this runs once per objective
  // evaluation, and a local is a heap allocation plus a full copy each time.
  // The sum is taken from the source and folded into the copy, so the
  // response is walked twice rather than three times -- same arithmetic, in
  // the same order, which is what keeps the curve bit-identical.
  irf_.resize(static_cast<std::size_t>(n_points));
  double total = 0.0;
  for (int i = 0; i < n_points; ++i) total += response[i];
  if (total > 0.0) {
    for (int i = 0; i < n_points; ++i) {
      irf_[static_cast<std::size_t>(i)] = response[i] / total;
    }
  } else {
    std::copy(response, response + n_points, irf_.begin());
  }
  irf_epoch_ = response_epoch_;
  irf_timeshift_ = timeshift;
  irf_valid_ = true;
  }
  const std::vector<double>& irf = irf_;

  curve_.assign(static_cast<std::size_t>(n_points), 0.0);

  // tttrlib's periodic reconvolution. Its stop arguments are *inclusive*
  // indices, so the last valid one is `n_points - 1`; passing a length reads
  // and writes one element past both buffers.
  const int last = n_points - 1;
  int convolution_stop =
      convolution_stop_ < 0 ? last : std::min(convolution_stop_, last);
  int stop = stop_ < 0 ? last : std::min(stop_, last);
  convolution_stop = std::max(0, convolution_stop);
  stop = std::max(0, stop);
  fconv_per_cs_ad<double>(curve_.data(), spectrum_.data(), irf.data(),
                          n_active_, stop, n_points, period_,
                          convolution_stop, dt_);

  // The basis: each species reconvolved on its own, with unit amplitude, in
  // the order the input spectrum gave them.
  //
  // This is NOT the same work as the summed call above, though it recurses
  // over the same species and this comment claimed it was until it was
  // measured: 7.5x the curve at K = 33 over 1 563 channels. Each call here
  // passes `numexp = 1`, which is below FCONV_AD_BLOCK_MIN, so it takes the
  // serial recursion while the summed call takes the 8-way blocked body --
  // the basis loses the blocking entirely. Closing that needs a kernel that
  // writes K columns instead of accumulating them, which belongs in
  // tttrlib's DecayConvolution.h (this file's copy is vendored and pinned
  // byte-identical by test/test_vendored_headers.py), not here.
  if (emit_basis_) {
    const std::size_t n_species = static_cast<std::size_t>(n_active_);
    const std::size_t n_bins = static_cast<std::size_t>(n_points);
    // Species-major first: the kernel writes each column into its own
    // contiguous run, so there is no separate column buffer and no copy out
    // of one. Writing the port's bins x species layout directly would put
    // consecutive writes `n_species` doubles apart -- a different cache line
    // every time, and at 33 species over 1563 channels that is 51k of them
    // against a buffer far larger than L1.
    basis_columns_.assign(n_bins * n_species, 0.0);
    double single[2];
    for (std::size_t s = 0; s < n_species; ++s) {
      single[0] = 1.0;
      single[1] = spectrum_[2 * s + 1];
      fconv_per_cs_ad<double>(basis_columns_.data() + s * n_bins, single,
                              irf.data(), 1, stop, n_points, period_,
                              convolution_stop, dt_);
    }
    // Then one blocked transpose into the contract the port promises. Tiled
    // because the naive loop is strided on whichever side it does not walk,
    // which is the cost this is here to avoid.
    basis_.assign(n_bins * n_species, 0.0);
    const std::size_t tile = 32;
    for (std::size_t b0 = 0; b0 < n_bins; b0 += tile) {
      const std::size_t b1 = std::min(b0 + tile, n_bins);
      for (std::size_t s0 = 0; s0 < n_species; s0 += tile) {
        const std::size_t s1 = std::min(s0 + tile, n_species);
        for (std::size_t b = b0; b < b1; ++b) {
          for (std::size_t s = s0; s < s1; ++s) {
            basis_[b * n_species + s] = basis_columns_[s * n_bins + b];
          }
        }
      }
    }
    const std::shared_ptr<GraphPort> bp = get_output_port(basis_port_key());
    if (!bp) {
      throw std::domain_error(
          "TcspcDecay '" + get_name() +
          "' emits the basis but has no '" + basis_port_key() + "' port");
    }
    bp->set_sanitize(false);
    bp->set_value_vector(basis_);
  }

  const double scatter = scatter_port_->get_value();
  if (scatter != 0.0) {
    for (int i = 0; i < n_points; ++i) {
      curve_[static_cast<std::size_t>(i)] += scatter * irf[static_cast<std::size_t>(i)];
    }
  }

  // Coates pile-up, chisurf's order: after the scatter term, before the
  // (auto)scaling -- the correction rescales the model, so scaling has to
  // come after it. The factors are computed from the data.
  if (pile_up_) {
    if (data_y_.size() != response_.size()) {
      throw std::domain_error(
          "TcspcDecay '" + get_name() +
          "' corrects pile-up, which needs data as long as the response");
    }
    if (period_ <= 0.0) {
      throw std::domain_error("TcspcDecay '" + get_name() +
                              "' corrects pile-up, which needs a period");
    }
    // Lifetimes and the period are in ns, so the repetition rate in MHz is
    // 1000 / period -- exactly the inverse of how chisurf builds the
    // period from `convolve.rep_rate`.
    const double rep_rate_mhz = 1000.0 / period_;
    add_pile_up_to_model_ad<double>(
        curve_.data(), n_points, data_y_.data(),
        static_cast<int>(data_y_.size()), rep_rate_mhz,
        pile_up_dead_time_ns_, pile_up_measurement_time_s_, 0, -1);
  }

  const double background = background_port_->get_value();
  if (autoscale_) {
    if (data_y_.size() != response_.size()) {
      throw std::domain_error(
          "TcspcDecay '" + get_name() +
          "' autoscales, which needs data as long as the response");
    }
    int begin = std::max(0, scale_start_);
    int end = scale_stop_ < 0 ? n_points : std::min(scale_stop_, n_points);
    end = std::min(end, static_cast<int>(data_y_.size()));
    if (end < begin) end = begin;
    n0_ = rescale_factor(curve_, data_y_, data_ey_, background, begin, end);
    // Published, because a fit that autoscales still has to report the
    // amplitude it settled on -- the port is the only place a caller can
    // read it from without evaluating the model a second time.
    n0_port_->set_value(n0_);
  } else {
    n0_ = n0_port_->get_value();
  }
  for (double& v : curve_) v = v * n0_ + background;

  // DNL: the measured channel-width table multiplies the finished curve --
  // chisurf's order, after the scaling and the constant background, before
  // the non-negativity clamp.
  if (!lin_table_.empty()) {
    if (lin_table_.size() != curve_.size()) {
      throw std::domain_error(
          "TcspcDecay '" + get_name() +
          "' linearizes, which needs a table as long as the response");
    }
    for (std::size_t i = 0; i < curve_.size(); ++i) {
      curve_[i] *= lin_table_[i];
    }
  }

  // A negative expected count is not a decay. ChiSurf clamps, and the clamp
  // is part of the objective rather than a cosmetic step: without it a fit
  // can trade a negative channel against a positive one.
  //
  // `v < 0.0`, not `!(v > 0.0)`: the second spelling also catches NaN and
  // would turn it into a *zero*, which in a fit reads as a good fit near
  // zero -- the same trap the port sanitiser sets, one line earlier. NaN
  // compares false against everything, so this leaves it alone, which is
  // what `np.maximum` does too.
  for (double& v : curve_) {
    if (v < 0.0) v = 0.0;
  }

  const std::shared_ptr<GraphPort> out = get_output_port(get_name());
  if (!out) {
    throw std::domain_error(
        "TcspcDecay '" + get_name() +
        "' writes its curve to the output port keyed by its own name, which "
        "this node does not have");
  }
  // Fit transport: a NaN must survive rather than be floored to `tiny`,
  // which in a fit reads as a *good* fit near zero. See ChiSquared::update.
  out->set_sanitize(false);
  out->set_value_vector(curve_);
  set_valid(true);
}

std::string TcspcDecay::describe() const {
  std::ostringstream out;
  out << "components     : " << n_lifetimes_ << "\n"
      << "channels       : " << response_.size() << "\n"
      << "dt / period    : " << dt_ << " / " << period_ << "\n"
      << "autoscale      : " << (autoscale_ ? "yes" : "no") << "\n"
      << "n0             : " << n0_ << "\n";
  return out.str();
}

IMPBFF_END_NAMESPACE
