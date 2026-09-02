/**
 *  \file IMP/bff/FcsSaturationCurve.h
 *  \brief Saturated FCS forward model as a graph node.
 *
 *  The saturated FCS diffusion shape (see FcsSaturation.h) wrapped as a
 *  `Node`, so an FCS kinetics model in "full" mode joins the graph the way
 *  `FcsMdfCurve` does for the MDF mode.  The ports carry the quantities a
 *  fit varies (power, extinction, beam waists, diffusion, N, baseline b,
 *  background bg); the photokinetic scheme (dark matrix, excitation matrix,
 *  brightness) and the quadrature grids are configuration -- set once at
 *  build time, not re-set per evaluation -- so the expensive Hankel-matrix
 *  cache built inside `fcs_saturated_curve_shape` survives across the
 *  evaluations the optimiser drives rather than being rebuilt every time.
 *
 *  Port names follow the model's FittingParameter names (power, extinction,
 *  w0, z0, D, N, b, bg) so the minimizer's `_claim` matches them up with the
 *  free parameters.  The node reads them in the units the model exposes
 *  (power in mW, lengths in nm, D in um^2/s) and scales to SI internally,
 *  exactly as the Python wrapper `saturated_curve_shape` does.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#ifndef IMPBFF_FCS_SATURATION_CURVE_H
#define IMPBFF_FCS_SATURATION_CURVE_H

#include <IMP/bff/bff_config.h>

#include <string>
#include <vector>

#include <IMP/bff/Node.h>
#include <IMP/bff/Port.h>

IMPBFF_BEGIN_NAMESPACE

//! Saturated FCS diffusion shape as a graph node over `fcs_saturated_curve_shape`.
/*!
    Ports (all scalar, created by build_ports()):

    | port | what it is | unit |
    |---|---|---|
    | `power` | total measured excitation power | mW |
    | `extinction` | molar extinction coefficient | M^-1 cm^-1 |
    | `w0` | lateral 1/e^2 beam waist | nm |
    | `z0` | axial 1/e^2 beam waist | nm |
    | `D` | translational diffusion coefficient | um^2/s |
    | `N` | average number of molecules | count |
    | `b` | baseline offset | dimensionless |
    | `bg` | background count rate | kHz |

    Configuration (set once, not ports):

    | config | what it is |
    |---|---|
    | dark_matrix | flat row-major N x N dark rates (Hz) |
    | exc_matrix | flat row-major N x N excitation cross sections |
    | brightness | per-state brightness (N) |
    | n_r, n_z | radial/axial grid points |
    | wavelength_m | excitation wavelength (m) |
    | include_bunching | multiply by X(tau) |

    `evaluate()` reads the ports, scales to SI, and calls
    `fcs_saturated_curve_shape` to produce the shape G(tau) (amplitude
    V0/Veff, before the 1/N normalisation and baseline the caller applies).
    The shape is published as the output port keyed by the node's name.

    \see FcsSaturation.h, Node, Expression, ChiSquared
*/
class IMPBFFEXPORT FcsSaturationCurve : public Node {
 public:
  explicit FcsSaturationCurve(const std::string& name = "fcs_saturation");

  //! Build the scalar input ports.  Not in the constructor (see FcsMdfCurve).
  void build_ports();

  //! The lag axis, in **seconds**, set once from the data.
  void set_axis(const std::vector<double>& tau);
  //! Numpy in, for the axis a caller holds as an array.
  void set_axis_array(double* in_axis, int n_axis);
  const std::vector<double>& get_axis() const { return tau_; }

  //! The fixed photokinetic scheme: dark rates (Hz), excitation cross
  //! sections, and per-state brightness.  Set once at build time.
  void set_scheme(const std::vector<double>& dark_matrix,
                  const std::vector<double>& exc_matrix,
                  int n_states,
                  const std::vector<double>& brightness);

  //! Quadrature resolution.
  void set_quadrature(int n_r, int n_z);
  int get_n_r() const { return n_r_; }
  int get_n_z() const { return n_z_; }

  //! Excitation wavelength (m).
  void set_wavelength(double wavelength_m);
  double get_wavelength() const { return wavelength_m_; }

  //! Whether to include the bunching factor X(tau).
  void set_include_bunching(bool v);
  bool get_include_bunching() const { return include_bunching_; }

  //! The shape from the last evaluation.
  const std::vector<double>& get_curve() const { return curve_; }

  void evaluate() override;

  std::string describe() const;

 private:
  std::vector<double> tau_;
  std::vector<double> curve_;

  Port* power_port_ = nullptr;
  Port* extinction_port_ = nullptr;
  Port* w0_port_ = nullptr;
  Port* z0_port_ = nullptr;
  Port* d_port_ = nullptr;
  Port* n_port_ = nullptr;
  Port* b_port_ = nullptr;
  Port* bg_port_ = nullptr;

  // Fixed configuration (set once, not ports).
  std::vector<double> dark_matrix_;
  std::vector<double> exc_matrix_;
  std::vector<double> brightness_;
  int n_states_ = 0;
  int n_r_ = 120;
  int n_z_ = 40;
  double wavelength_m_ = 488e-9;
  bool include_bunching_ = true;

  //! The port values the published curve belongs to.
  bool curve_built_ = false;
  double curve_power_ = 0.0, curve_extinction_ = 0.0;
  double curve_w0_ = 0.0, curve_z0_ = 0.0, curve_d_ = 0.0;
};

IMPBFF_END_NAMESPACE

#endif /* IMPBFF_FCS_SATURATION_CURVE_H */
