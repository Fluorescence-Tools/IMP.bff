/**
 * \file IMP/bff/TcspcDecay.h
 * \brief A time-correlated single-photon-counting decay, as a node.
 *
 * `LifetimeSpectrum` is the experiment-neutral answer: (amplitude, rate)
 * pairs, with nothing an instrument adds. This is the other half -- what a
 * TCSPC setup makes of such a spectrum: convolved with a measured response,
 * shifted against it, scaled to the data, plus scatter and a constant
 * background. As a `Node`, so a whole fit is one C++ graph:
 * `TcspcDecay -> ChiSquared -> Minimizer`, exactly the arrangement
 * `Expression` gives a parse model.
 *
 * \par The arithmetic is not here either
 * The periodic convolution is tttrlib's `fconv_per_cs_ad<double>` and the
 * timeshift is its `shift_lamp_ad<double>`, both taken from the vendored
 * copy of `DecayConvolution.h`. That is deliberate and it is the same
 * decision `Expression` records about the expression engine: a decay is a
 * *curve*, `AGENTS.md` puts curves in tttrlib, and a second implementation
 * of a convolution is how two libraries end up disagreeing about what a
 * lifetime is. What this class contributes is the graph -- ports, ordering,
 * invalidation -- and the composition of the instrument model around those
 * kernels.
 *
 * \par What is a port and what is data
 * A port is something an optimiser writes: the amplitudes and lifetimes,
 * the scatter fraction, the constant background, the amplitude `n0` and the
 * timeshift. Everything else is measured or configured and is set once --
 * the response function, the data (needed only for autoscaling), the
 * channel width, the excitation period and the windows. The split is the
 * one that matters for speed: nothing that is set once is re-read, and
 * nothing that is a port crosses the boundary during a fit.
 *
 * \par The response is handed over already processed
 * Background subtraction, truncation and normalisation of the measured
 * response are the application's business and depend on parameters that do
 * not move during a fit; the caller does them once and hands the result in.
 * The *timeshift* does move, so it is a port and is applied here -- to the
 * response as given, followed by the unit-sum normalisation, which is the
 * order ChiSurf uses and is not interchangeable with the other one.
 *
 * \see LifetimeSpectrum, ChiSquared, Expression, Minimizer, Node
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_TCSPCDECAY_H
#define IMPBFF_TCSPCDECAY_H

#include <IMP/bff/bff_config.h>

#include <string>
#include <vector>

#include <IMP/bff/Node.h>
#include <IMP/bff/Port.h>

IMPBFF_BEGIN_NAMESPACE

//! The model curve of a multi-exponential decay through a TCSPC instrument.
/**
 * Ports, all scalars, created by set_number_of_lifetimes():
 *
 * | port | what it is |
 * |---|---|
 * | `a0`, `t0`, `a1`, `t1`, ... | amplitude and lifetime of each component |
 * | `scatter` | fraction of the (normalised) response added to the curve |
 * | `background` | constant offset, added after scaling |
 * | `n0` | the curve's amplitude; **written** by the node when autoscaling |
 * | `timeshift` | shift of the response against the data, in samples |
 *
 * The curve is written to the output port keyed by the node's own name,
 * which is the protocol `ChiSquared` and `Expression` already use.
 */
class IMPBFFEXPORT TcspcDecay : public Node {
 public:
  explicit TcspcDecay(const std::string& name = "decay");

  //! Build the node's ports: `2 * n` lifetime ports named `a0`, `t0`,
  //! `a1`, `t1`, ..., plus the four scalars, the first time it is called.
  /** **This is the call that makes the node usable**, and it cannot be done
      in the constructor: a `Node` that owns ports must already be owned by a
      `shared_ptr` (`add_port` reaches for `shared_from_this()`), and a
      Python-wrapped node is constructed before it is owned. */
  void set_number_of_lifetimes(int n);
  int get_number_of_lifetimes() const { return n_lifetimes_; }

  //! The measured response function, already background-subtracted,
  //! truncated and scaled by the caller -- everything except the shift.
  void set_response(const std::vector<double>& response);
  const std::vector<double>& get_response() const { return response_; }

  //! The measured curve and its errors, used **only** for autoscaling.
  /** A node that does not autoscale needs neither; `ChiSquared` downstream
      holds its own copy of the data, and that is the one the misfit uses. */
  void set_data(const std::vector<double>& y, const std::vector<double>& ey);

  //! Channel width and excitation period, in the units of the lifetimes.
  void set_timing(double dt, double period);
  double get_dt() const { return dt_; }
  double get_period() const { return period_; }

  //! Last channel the convolution recursion runs to, and the last channel
  //! written. Both are *inclusive* indices, as tttrlib's kernel takes them.
  void set_convolution_range(int convolution_stop, int stop);

  //! Whether `n0` is computed from the data rather than read from its port.
  void set_autoscale(bool v) { autoscale_ = v; }
  bool get_autoscale() const { return autoscale_; }

  //! Coates pile-up on the model curve, chisurf's order: after the scatter
  //! term, before the (auto)scaling.
  /** The scaling factors are computed from the *data* (set_data), so
      enabling pile-up requires the data even when the node does not
      autoscale. The repetition rate is derived from the excitation period
      (`1000 / period` MHz for lifetimes in ns), so there is no second copy
      of it to disagree with set_timing. */
  void set_pile_up(bool v) { pile_up_ = v; }
  bool get_pile_up() const { return pile_up_; }
  //! Detection dead time (ns) and measurement time (s) of the correction.
  void set_pile_up_parameters(double dead_time_ns, double measurement_time_s) {
    pile_up_dead_time_ns_ = dead_time_ns;
    pile_up_measurement_time_s_ = measurement_time_s;
  }
  double get_pile_up_dead_time() const { return pile_up_dead_time_ns_; }
  double get_pile_up_measurement_time() const {
    return pile_up_measurement_time_s_;
  }

  //! Multiply the finished curve by a fixed linearization table (DNL).
  /** chisurf's order: after the scaling and the constant background, before
      the non-negativity clamp. The table is *measured* -- a smoothed ratio
      of a reference decay, describing each channel's effective width -- and
      does not move during a fit, so it is configuration rather than a port.
      An empty table (the default) disables the stage; a non-empty one must
      be as long as the response. */
  void set_linearization(const std::vector<double>& table);
  const std::vector<double>& get_linearization() const { return lin_table_; }
  //! Numpy spelling of set_linearization (the IN_ARRAY1 typemap).
  void set_linearization_array(double* in_table, int n_table);

  //! The half-open channel window the autoscale is computed over.
  void set_scale_range(int start, int stop);

  //! `|a|` on the amplitudes before use, as ChiSurf's
  //! `absolute_amplitudes` does (it spells it `sqrt(a**2)`).
  void set_absolute_amplitudes(bool v) { absolute_amplitudes_ = v; }
  bool get_absolute_amplitudes() const { return absolute_amplitudes_; }

  //! Divide the amplitudes by `|sum|`, as ChiSurf's `normalize_amplitudes`
  //! does. One amplitude then pins the scale and is not a free parameter.
  void set_normalize_amplitudes(bool v) { normalize_amplitudes_ = v; }
  bool get_normalize_amplitudes() const { return normalize_amplitudes_; }

  //! Drop species whose amplitude is negligible against the largest one.
  /*!
      **Why a decay model needs this at all.** A multi-exponential decay costs
      one serial recursion over every channel *per species*, so the
      reconvolution is linear in the species count -- and a model whose
      spectrum comes from a *distribution* has as many species as the
      distribution has bins, whatever their weight. A FRET decay over the
      96-point distance axis is 97 species, of which (measured, on a Gaussian
      at 45 A with sigma 6) only 44 carry a weight above `1e-12` of the
      largest and only 53 above `1e-14`. The rest cost a full recursion each
      to contribute nothing.

      The threshold is **relative to the largest amplitude** and is applied
      after `absolute_amplitudes` and `normalize_amplitudes`, because those
      describe what the model means by an amplitude and this describes which
      of them are worth computing.

      Measured on that FRET spectrum, against the unpruned curve:

      | threshold | species | max abs. curve change / peak | reconvolution |
      |---|---|---|---|
      | `0` (default) | 97 | 0 (exact) | 96.2 us |
      | `1e-14` | 53 | 1.1e-15 | 55.0 us |
      | `1e-12` | 44 | 1.1e-13 | 44.9 us |
      | `1e-10` | 37 | 2.0e-11 | 39.8 us |
      | `1e-6` | 26 | 2.1e-7 | 28.4 us |

      **The default is `0`, which drops only amplitudes that are exactly
      zero and is therefore bit-exact.** A library that silently drops terms
      is the failure this module spent 2026-09-01 fixing at a different
      layer; the caller that knows what its data are worth chooses the
      threshold. At `1e-14` the change to the curve is smaller than the
      change from summing the same species in a different order, which is why
      that is the value ChiSurf asks for.
   */
  void set_amplitude_threshold(double relative);
  double get_amplitude_threshold() const { return amplitude_threshold_; }

  //! Species kept by the last evaluation, after the threshold.
  /** Differs from `get_number_of_lifetimes()` whenever pruning fired, and it
      is the count the reconvolution actually paid for. */
  int get_number_of_active_lifetimes() const {
    return static_cast<int>(spectrum_.size() / 2);
  }

  //! Take the spectrum from the `lifetime_spectrum` port, not the scalars.
  /*!
      **This is what makes the node composable, and it is the whole reason a
      decay model other than a plain multi-exponential can be a graph.**

      Every TCSPC model in an application this serves has the *same*
      instrument model -- reconvolve, add scatter, scale to the data, add a
      background -- and differs only in how the (amplitude, lifetime) pairs
      are arrived at. A FRET model derives them from a distance, a
      distribution model from a distance distribution, a mixture from two
      spectra. Wiring each of those as its own set of `a`/`t` ports would put
      the *deriving* back in the caller, once per iteration, which is the one
      thing that measured as a regression.

      So the spectrum can instead arrive on an input port, which some
      upstream node writes: `<whatever computes a spectrum> -> TcspcDecay ->
      ChiSquared`. The node keeps its opinion about the instrument and holds
      none about the photophysics.

      The port carries the interleaved `(a0, t0, a1, t1, ...)` layout the
      kernels take, and `absolute_amplitudes` / `normalize_amplitudes` are
      applied to it exactly as they are to the scalar ports -- they describe
      what the *model* means by an amplitude, not where the number came from.
   */
  void set_spectrum_from_port(bool v);
  bool get_spectrum_from_port() const { return spectrum_from_port_; }

  //! The key of the input port carrying an interleaved lifetime spectrum.
  static const char* spectrum_port_key() { return "lifetime_spectrum"; }

  //! The interleaved `(a0, t0, a1, t1, ...)` spectrum the last evaluation
  //! built from the ports, after abs() and normalisation.
  const std::vector<double>& get_lifetime_spectrum() const {
    return spectrum_;
  }

  //! The amplitude the last evaluation used -- the autoscaled one when
  //! autoscaling, otherwise whatever the `n0` port held.
  double get_n0() const { return n0_; }

  //! The model curve from the last evaluation.
  const std::vector<double>& get_curve() const { return curve_; }

  //! Numpy in, for the arrays a caller holds as arrays.
  /** The `std::vector` overloads above make a Python caller build a list of
      every channel before a single one is stored, which on a 4096-channel
      decay costs more than the evaluation it serves. */
  void set_response_array(double* in_response, int n_response);
  void set_data_arrays(double* in_data_y, int n_data_y, double* in_data_ey,
                       int n_data_ey);

  //! Build the curve from the ports and write it out.
  void evaluate() override;

  std::string describe() const;

 private:
  std::vector<double> response_;
  std::vector<double> shifted_;
  //! The shifted, unit-sum response the last evaluation convolved with.
  /** A member and not a local, because a local is a heap allocation and
      a 512-double copy *per objective evaluation*, and an objective
      evaluation is the unit a fit spends its time in. */
  std::vector<double> irf_;
  std::vector<double> data_y_;
  std::vector<double> data_ey_;
  std::vector<double> lin_table_;
  std::vector<double> spectrum_;
  std::vector<double> curve_;
  //! The lifetime ports, in the order the spectrum interleaves them.
  std::vector<Port*> lifetime_ports_;
  Port* spectrum_port_ = nullptr;
  Port* scatter_port_ = nullptr;
  Port* background_port_ = nullptr;
  Port* n0_port_ = nullptr;
  Port* timeshift_port_ = nullptr;

  int n_lifetimes_ = 0;
  double dt_ = 1.0;
  double period_ = 0.0;
  int convolution_stop_ = -1;
  int stop_ = -1;
  int scale_start_ = 0;
  int scale_stop_ = -1;
  double n0_ = 1.0;
  bool autoscale_ = false;
  bool pile_up_ = false;
  double pile_up_dead_time_ns_ = 85.0;
  double pile_up_measurement_time_s_ = 0.0;
  bool spectrum_from_port_ = false;
  bool absolute_amplitudes_ = false;
  bool normalize_amplitudes_ = false;
  double amplitude_threshold_ = 0.0;
  //! Species surviving the threshold -- what the kernel is handed.
  /** Kept apart from `n_lifetimes_`, which is the number of `a`/`t`
      *ports* and must survive an evaluation: on the scalar path nothing
      restores it, so pruning it in place would make a component with a
      zero amplitude disappear permanently rather than for one curve. */
  int n_active_ = 0;

  void add_scalar_port(const std::string& key, double value, Port** slot);
  void build_spectrum();
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_TCSPCDECAY_H
