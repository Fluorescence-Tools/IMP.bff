/**\file IMP/bff/FRETSpectrumNode.h
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#ifndef IMPBFF_FRETSPECTRUMNODE_H
#define IMPBFF_FRETSPECTRUMNODE_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/GraphNode.h>
#include <IMP/bff/GraphPort.h>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

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
class IMPBFFEXPORT FRETSpectrumNode : public GraphNode {
 public:
  explicit FRETSpectrumNode(const std::string& name = "fret");

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
  GraphPort* donor_port_ = nullptr;
  GraphPort* distance_port_ = nullptr;
  GraphPort* x_donly_port_ = nullptr;
  GraphPort* forster_radius_port_ = nullptr;
  GraphPort* tau0_port_ = nullptr;
  GraphPort* kappa2_port_ = nullptr;
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_FRETSPECTRUMNODE_H
