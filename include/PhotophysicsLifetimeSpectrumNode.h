/**\file IMP/bff/PhotophysicsLifetimeSpectrumNode.h
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#ifndef IMPBFF_PHOTOPHYSICSLIFETIMESPECTRUMNODE_H
#define IMPBFF_PHOTOPHYSICSLIFETIMESPECTRUMNODE_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/GraphNode.h>
#include <IMP/bff/GraphPort.h>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! An interleaved `(a0, t0, a1, t1, ...)` spectrum built from scalar ports.
/**
 * The smallest producer there is, and the one every other model starts from:
 * a fitted amplitude and lifetime per species, published as one vector the
 * downstream node reads. `TCSPCDecay` can do this itself for the plain
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
 * which is the protocol `TCSPCDecay`, `FitChiSquared` and `GraphExpression` use.
 */
class IMPBFFEXPORT PhotophysicsLifetimeSpectrumNode : public GraphNode {
 public:
  explicit PhotophysicsLifetimeSpectrumNode(const std::string& name = "lifetimes");

  //! Build `2 * n` ports named `a0`, `t0`, `a1`, `t1`, ...
  /** As with `TCSPCDecay`, this cannot happen in the constructor: a `GraphNode`
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
  std::vector<GraphPort*> lifetime_ports_;
  int n_lifetimes_ = 0;
  bool absolute_amplitudes_ = false;
  bool normalize_amplitudes_ = false;
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_PHOTOPHYSICSLIFETIMESPECTRUMNODE_H
