/**\file IMP/bff/PhotophysicsAnisotropySpectrumNode.h
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#ifndef IMPBFF_PHOTOPHYSICSANISOTROPYSPECTRUMNODE_H
#define IMPBFF_PHOTOPHYSICSANISOTROPYSPECTRUMNODE_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/GraphNode.h>
#include <IMP/bff/GraphPort.h>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

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
class IMPBFFEXPORT PhotophysicsAnisotropySpectrumNode : public GraphNode {
 public:
  //! What the detector selects; anything else leaves the spectrum alone.
  enum Polarization { VM = 0, VV = 1, VH = 2, VV_VH = 3 };

  explicit PhotophysicsAnisotropySpectrumNode(const std::string& name = "anisotropy");

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
  std::vector<GraphPort*> rotation_ports_;
  GraphPort* spectrum_port_ = nullptr;
  GraphPort* r0_port_ = nullptr;
  GraphPort* g_port_ = nullptr;
  GraphPort* l1_port_ = nullptr;
  GraphPort* l2_port_ = nullptr;
  int n_rotations_ = 0;
  Polarization polarization_ = VM;

  void build_rotation_spectrum();
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_PHOTOPHYSICSANISOTROPYSPECTRUMNODE_H
