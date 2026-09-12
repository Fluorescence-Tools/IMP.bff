/**\file IMP/bff/PolymerDistances.h
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#ifndef IMPBFF_POLYMERDISTANCES_H
#define IMPBFF_POLYMERDISTANCES_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/GraphNode.h>
#include <IMP/bff/GraphPort.h>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

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
 * `(p0, r0, p1, r1, ...)`, the layout `FRETSpectrumNode` reads, with the
 * kernel's own normalisation -- the same array the Python property returns.
 */
class IMPBFFEXPORT PolymerDistances : public GraphNode {
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
  std::vector<GraphPort*> parameter_ports_;
  int n_k_ = 2000;
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_POLYMERDISTANCES_H
