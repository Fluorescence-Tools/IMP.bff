/**
 * \file IMP/bff/SpectrumNode.h
 * \brief The nodes that *produce* an interleaved lifetime spectrum.
 *
 * `TcspcDecay` holds the instrument -- reconvolve, add scatter, scale to the
 * data, add a background -- and since it grew a `lifetime_spectrum` input
 * port it holds no opinion at all about where the (amplitude, lifetime)
 * pairs came from. This file is the other half: the small nodes that arrive
 * at those pairs, so a decay model that is *not* a plain multi-exponential
 * still becomes one C++ graph rather than a Python callback per iteration.
 *
 * \par Why the producers are separate nodes
 * Every TCSPC model in the applications this serves has the same instrument
 * model and differs only in the photophysics upstream of it. Wiring each
 * one's physics into `TcspcDecay` would make that class the union of every
 * model anyone ever fits; wiring it into the *caller* would put the deriving
 * back in Python once per iteration, which is the arrangement that measured
 * as a regression (`okf/log.md` 2026-09-01 (9)). A node per producer is the
 * third option and the only one that composes:
 *
 *     LifetimeSpectrumNode -> AnisotropySpectrum -> TcspcDecay -> ChiSquared
 *
 * \see TcspcDecay, LifetimeSpectrum, ChiSquared, Minimizer, Node
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_SPECTRUMNODE_H
#define IMPBFF_SPECTRUMNODE_H

#include <IMP/bff/bff_config.h>

#include <string>
#include <vector>

#include <IMP/bff/Node.h>
#include <IMP/bff/Port.h>

IMPBFF_BEGIN_NAMESPACE

//! An interleaved `(a0, t0, a1, t1, ...)` spectrum built from scalar ports.
/**
 * The smallest producer there is, and the one every other model starts from:
 * a fitted amplitude and lifetime per species, published as one vector the
 * downstream node reads. `TcspcDecay` can do this itself for the plain
 * multi-exponential case; this node exists for when something has to sit
 * *between* the parameters and the instrument -- a polarisation, a FRET
 * quenching -- and that something needs the donor spectrum as an input.
 *
 * Ports, created by set_number_of_lifetimes():
 *
 * | port | what it is |
 * |---|---|
 * | `a0`, `t0`, `a1`, `t1`, ... | amplitude and lifetime of each species |
 *
 * The spectrum is written to the output port keyed by the node's own name,
 * which is the protocol `TcspcDecay`, `ChiSquared` and `Expression` use.
 */
class IMPBFFEXPORT LifetimeSpectrumNode : public Node {
 public:
  explicit LifetimeSpectrumNode(const std::string& name = "lifetimes");

  //! Build `2 * n` ports named `a0`, `t0`, `a1`, `t1`, ...
  /** As with `TcspcDecay`, this cannot happen in the constructor: a `Node`
      that owns ports must already be held by a `shared_ptr`, and a
      Python-wrapped node is constructed before it is owned. */
  void set_number_of_lifetimes(int n);
  int get_number_of_lifetimes() const { return n_lifetimes_; }

  //! `|a|` on the amplitudes, as ChiSurf's `absolute_amplitudes` does.
  /** Invalidates, so a flag flipped after an evaluation is seen. It reads
      like configuration a caller sets once, but it changes the *value* the
      node publishes, and a setter that changes a value without invalidating
      hands the next reader a stale spectrum with nothing to say so. */
  void set_absolute_amplitudes(bool v) {
    absolute_amplitudes_ = v;
    set_valid(false);
  }
  bool get_absolute_amplitudes() const { return absolute_amplitudes_; }

  //! Divide the amplitudes by `|sum|`, as `normalize_amplitudes` does.
  void set_normalize_amplitudes(bool v) {
    normalize_amplitudes_ = v;
    set_valid(false);
  }
  bool get_normalize_amplitudes() const { return normalize_amplitudes_; }

  const std::vector<double>& get_spectrum() const { return spectrum_; }

  void evaluate() override;

 private:
  std::vector<double> spectrum_;
  std::vector<Port*> lifetime_ports_;
  int n_lifetimes_ = 0;
  bool absolute_amplitudes_ = false;
  bool normalize_amplitudes_ = false;
};

//! The polarised lifetime spectrum a VV or VH detector sees.
/**
 * A magic-angle decay carries no anisotropy, so a VM model needs none of
 * this. A polarisation-resolved one does: the measured decay is the
 * *product* of the fluorescence decay and the anisotropy decay
 * \f$r(t) = \sum_i b_i e^{-t/\rho_i}\f$, and because both are sums of
 * exponentials that product is again a sum of exponentials -- over the
 * Cartesian product of the two spectra, with the harmonic mean of the two
 * time constants. That is the whole reason this can be a spectrum transform
 * rather than a curve one, and it is why the instrument node downstream does
 * not have to know a polarisation happened.
 *
 * \f[
 *   f_{VV} = f\,(1 + (2 - 3 l_1)\,r), \qquad
 *   f_{VH} = f\,(1 - (1 - 3 l_2)\,r) / G
 * \f]
 *
 * \f$G = S_\parallel / S_\perp\f$ is the sensitivity ratio, so the
 * perpendicular channel records \f$1/G\f$ of what an equally sensitive one
 * would -- it **divides**. `l1` and `l2` enter the amplitudes in the
 * Schaffer/Eggeling parameterisation, *not* as a 2x2 mixing of an ideal
 * pair; the two spellings use the same symbols for different quantities and
 * only this one inverts back to the anisotropy it was built from.
 *
 * \par The rotational amplitudes are normalised to r0, and that is a model
 * decision this node reproduces rather than invents. The fitted \f$b_i\f$
 * are taken absolute, divided by their sum and multiplied by `r0`, so `r0`
 * alone sets \f$r(0)\f$ and the \f$b_i\f$ only divide it up. The
 * correlation times are taken absolute for the same reason lifetimes are: a
 * fit that walks one through zero would otherwise put a growing exponential
 * in the model.
 *
 * Ports:
 *
 * | port | what it is |
 * |---|---|
 * | `lifetime_spectrum` | the unpolarised interleaved spectrum, from upstream |
 * | `r0` | the fundamental anisotropy, \f$r(0)\f$ |
 * | `g` | \f$S_\parallel / S_\perp\f$ |
 * | `l1`, `l2` | the depolarisation mixing factors |
 * | `b0`, `rho0`, `b1`, `rho1`, ... | amplitude and correlation time of each rotation |
 */
class IMPBFFEXPORT AnisotropySpectrum : public Node {
 public:
  //! What the detector selects; anything else leaves the spectrum alone.
  enum Polarization { VM = 0, VV = 1, VH = 2, VV_VH = 3 };

  explicit AnisotropySpectrum(const std::string& name = "anisotropy");

  //! Build the ports: `2 * n` rotation ports plus `r0`, `g`, `l1`, `l2`.
  void set_number_of_rotations(int n);
  int get_number_of_rotations() const { return n_rotations_; }

  //! The key of the input port carrying the unpolarised spectrum.
  static const char* spectrum_port_key() { return "lifetime_spectrum"; }

  void set_polarization(Polarization p);
  Polarization get_polarization() const { return polarization_; }

  //! The same choice by the name the application spells it with.
  /** `"vm"`, `"vv"`, `"vh"` and `"vv/vh"`, case-insensitively. An unknown
      name is refused rather than silently treated as VM, which would return
      the spectrum unchanged and read as a model with no anisotropy at all. */
  void set_polarization_name(const std::string& name);
  std::string get_polarization_name() const;

  //! The rotation spectrum the last evaluation used, after normalisation.
  const std::vector<double>& get_rotation_spectrum() const {
    return rotation_;
  }

  //! The polarised spectrum from the last evaluation.
  const std::vector<double>& get_spectrum() const { return spectrum_; }

  void evaluate() override;

  std::string describe() const;

 private:
  std::vector<double> rotation_;
  std::vector<double> spectrum_;
  std::vector<Port*> rotation_ports_;
  Port* spectrum_port_ = nullptr;
  Port* r0_port_ = nullptr;
  Port* g_port_ = nullptr;
  Port* l1_port_ = nullptr;
  Port* l2_port_ = nullptr;
  int n_rotations_ = 0;
  Polarization polarization_ = VM;

  void build_rotation_spectrum();
};

//! A distance distribution built from a sum of generalised normals.
/**
 * The producer for the FRET model everyone actually fits: a distance
 * distribution as a weighted sum of (possibly skewed) normal distributions,
 * evaluated on a fixed distance axis. Downstream it becomes a rate
 * distribution and then a lifetime spectrum, so the whole chain from a
 * fitted mean distance to a model curve is
 *
 *     GaussianDistances -> FretSpectrum -> TcspcDecay -> ChiSquared
 *
 * and nothing in it returns to the caller.
 *
 * \par The axis is data, the components are ports
 * The distance axis is a setting -- its range and resolution do not move
 * during a fit -- so the caller hands it over once. The mean, width, skew
 * and weight of each component are what an optimiser writes, so they are
 * ports.
 *
 * \par Two normalisations, in this order, and they are not interchangeable
 * Each component is normalised to unit sum on the axis *before* it is
 * weighted, and the weighted sum is normalised again. Normalising only at
 * the end would let a wide component contribute less than its weight says
 * merely because more of it falls off the end of the axis.
 *
 * Ports, created by set_number_of_components():
 *
 * | port | what it is |
 * |---|---|
 * | `mean0`, `sigma0`, `shape0`, `amplitude0`, ... | one distance component |
 *
 * The output is interleaved `(p0, r0, p1, r1, ...)` -- weight and distance,
 * the layout `FretSpectrum` reads.
 */
//! A polymer-chain distance distribution as a producer node.
/**
 * The closed-form polymer models of chisurf's `FRETModel` subclasses, each a
 * distance distribution and nothing else (board `T-20260901-08`): the
 * worm-like chain (with or without dye-linker broadening), the SAW-nu
 * des Cloizeaux form, and the Ising two-state Gaussian chain. One node with
 * a mode rather than one node per model, because after the kernel ports of
 * `T-20260902-01/-14` every mode is a thin dispatch into `PolymerChain.h`
 * -- the cost lives in the kernel, and the kernels are shared with the
 * Python forwarders in `chisurf/core/math/functions/rdf.py`, so the graph
 * path and the numpy path cannot drift apart.
 *
 * Ports, created by set_mode():
 *
 * | mode | ports |
 * |---|---|
 * | `worm_like_chain` | `chain_length`, `persistence_length` |
 * | `worm_like_chain_linker` | the same plus `sigma_linker` |
 * | `saw_nu` | `r_rms`, `nu` |
 * | `ising_chain` | `n_residues`, `b_structured`, `b_unstructured`, `coupling`, `field` |
 *
 * The worm-like chain takes the *contour* and *persistence* lengths as its
 * ports and derives the kernel's dimensionless `kappa = lp / l` itself, the
 * same derivation chisurf's model property makes; `n_residues` is rounded to
 * the nearest integer as chisurf rounds it. The output is interleaved
 * `(p0, r0, p1, r1, ...)`, the layout `FretSpectrum` reads, with the
 * kernel's own normalisation -- the same array the Python property returns.
 */
class IMPBFFEXPORT PolymerDistances : public Node {
 public:
  explicit PolymerDistances(const std::string& name = "distances");

  //! Select the distribution and build its input ports (once).
  void set_mode(const std::string& mode);
  const std::string& get_mode() const { return mode_; }

  //! The distance axis the distribution is evaluated on, in Angstrom.
  void set_axis(const std::vector<double>& axis);
  const std::vector<double>& get_axis() const { return axis_; }
  //! Numpy in, for the axis a caller holds as an array.
  void set_axis_array(double* in_axis, int n_axis);

  //! k points of the Ising inverse transform (set-once configuration).
  void set_n_k(int n_k) { n_k_ = n_k; set_valid(false); }
  int get_n_k() const { return n_k_; }

  //! The interleaved `(p, r)` distribution from the last evaluation.
  const std::vector<double>& get_distribution() const { return spectrum_; }

  void evaluate() override;

 private:
  std::string mode_;
  std::vector<double> axis_;
  std::vector<double> spectrum_;
  std::vector<Port*> parameter_ports_;
  int n_k_ = 2000;
};

class IMPBFFEXPORT GaussianDistances : public Node {
 public:
  explicit GaussianDistances(const std::string& name = "distances");

  //! Build `4 * n` ports named `mean0`, `sigma0`, `shape0`, `amplitude0`,
  //! `mean1`, ... -- the index is a suffix with no separator.
  void set_number_of_components(int n);
  int get_number_of_components() const { return n_components_; }

  //! The distance axis the distribution is evaluated on, in Angstrom.
  void set_axis(const std::vector<double>& axis);
  const std::vector<double>& get_axis() const { return axis_; }
  //! Numpy in, for the axis a caller holds as an array.
  void set_axis_array(double* in_axis, int n_axis);

  //! The interleaved `(p, r)` distribution from the last evaluation.
  const std::vector<double>& get_distribution() const { return spectrum_; }

  //! Use the distance-between-two-Gaussians kernel instead of the skewed one.
  /*! These are different probability densities, not two parameterisations of
      one: the two-cloud form carries an extra `r/mean` factor and an
      antisymmetric second term. The builder used to refuse a model with this
      set, which left `GaussianModel` without a graph whenever it was on. */
  void set_distance_between_gaussians(bool value) {
    distance_between_gaussians_ = value;
    set_valid(false);
  }
  bool get_distance_between_gaussians() const {
    return distance_between_gaussians_;
  }

  void evaluate() override;

 private:
  std::vector<double> axis_;
  std::vector<double> density_;
  std::vector<double> spectrum_;
  std::vector<Port*> component_ports_;
  int n_components_ = 0;
  bool distance_between_gaussians_ = false;
};

//! The lifetime spectrum of a donor quenched by FRET over a distance
//! distribution.
/**
 * The piece that makes a FRET decay a graph. It takes two spectra and
 * returns one:
 *
 * * the **donor** lifetime spectrum -- what the donor would do with no
 *   acceptor present -- from upstream;
 * * a **distance distribution**, weight and distance interleaved, from
 *   whichever producer the model uses;
 *
 * and returns the lifetime spectrum of the mixture, ready for the
 * instrument node.
 *
 * \par The arithmetic, and why it is a spectrum transform
 * A distance becomes a transfer rate,
 * $k_{FRET} = rac{3}{2}\kappa^2 	au_0^{-1} (R_0/r)^6$, and a
 * quenched species decays at the *sum* of its own rate and the transfer
 * rate. Rates add, so the quenched spectrum is the Cartesian product of the
 * donor's rates and the transfer rates with the amplitudes multiplied --
 * which is exactly why this can be done on spectra rather than on curves.
 *
 * \par Donor-only is a mixture, not a term
 * A sample never labels perfectly. The fraction `x_donly` that carries no
 * acceptor decays at the *donor's* rates, so it enters as its own set of
 * components appended to the quenched ones, scaled by `x_donly` while the
 * quenched ones are scaled by `1 - x_donly`. They are appended
 * unconditionally, including at `x_donly = 0` where they carry zero
 * amplitude: the length of a spectrum is observable, and a node that
 * sometimes returns a shorter one is a second code path for no gain.
 *
 * Ports:
 *
 * | port | what it is |
 * |---|---|
 * | `donor_lifetime_spectrum` | interleaved `(a, tau)`, from upstream |
 * | `distance_distribution` | interleaved `(p, r)`, from upstream |
 * | `x_donly` | fraction of the sample carrying no acceptor |
 * | `forster_radius` | $R_0$, in the units of the distances |
 * | `tau0` | the donor lifetime $R_0$ was determined at |
 * | `kappa2` | the orientation factor |
 */
class IMPBFFEXPORT FretSpectrum : public Node {
 public:
  explicit FretSpectrum(const std::string& name = "fret");

  //! Build the node's ports. Cannot be done in the constructor.
  void build_ports();

  static const char* donor_port_key() { return "donor_lifetime_spectrum"; }
  static const char* distance_port_key() { return "distance_distribution"; }

  //! The interleaved `(a, tau)` spectrum from the last evaluation.
  const std::vector<double>& get_spectrum() const { return spectrum_; }

  void evaluate() override;

  std::string describe() const;

 private:
  std::vector<double> spectrum_;
  Port* donor_port_ = nullptr;
  Port* distance_port_ = nullptr;
  Port* x_donly_port_ = nullptr;
  Port* forster_radius_port_ = nullptr;
  Port* tau0_port_ = nullptr;
  Port* kappa2_port_ = nullptr;
};

//! Amplitudes times `n`, leaving the time constants alone (ChiSurf's `e1tn`).
/** Free rather than private because the FRET producers need the same
    algebra, and a second copy of it is how two spectra end up scaled
    differently. Appends to \p out. */
IMPBFFEXPORT void interleaved_scale_amplitudes(
    const std::vector<double>& spectrum, double n, std::vector<double>* out);

//! The Cartesian product of two interleaved spectra (ChiSurf's `elte2`).
/** Amplitudes multiply; time constants combine as the harmonic mean
    \f$1/(1/t_1 + 1/t_2)\f$, which is what the product of two exponentials
    is. First spectrum major, i.e. entry \f$k = i n_2 + j\f$, because the
    order is observable in the spectrum a caller reads back. */
IMPBFFEXPORT std::vector<double> interleaved_product(
    const std::vector<double>& first, const std::vector<double>& second);

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_SPECTRUMNODE_H
