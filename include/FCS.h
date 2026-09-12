#ifndef IMPBFF_FCS_H
#define IMPBFF_FCS_H

/**
 *  \file IMP/bff/FCS.h
 *  \brief FCS forward models beyond the closed forms: the MDF and saturation shapes.
 *
 * Three former headers, in the order a reader needs them:
 *
 * 1. **MDF** (formerly `FcsMdf.h`) -- Enderlein's molecule-detection-function
 *    diffusion shape, and #IMP::bff::FCSMdfCurve, the graph node that puts the
 *    `"mdf"` FCS mode next to the closed-form ones.
 * 2. **Saturation** (formerly `FcsSaturation.h`) -- steady-state photophysics
 *    on a grid, the kernel of the saturated diffusion shape.
 * 3. **Saturation as a node** (formerly `FCSSaturationCurve.h`) --
 *    #IMP::bff::FCSSaturationCurve, that kernel on the graph.
 *
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

// -------- from FcsMdf.h --------
/**
 *  (formerly IMP/bff/FcsMdf.h, now a section of this file)
 *  \brief Enderlein molecule-detection-function FCS forward model.
 *
 *  The Gauss--Lorentz MDF of Enderlein et al. (2005): a confocal detection
 *  volume built from a Gaussian excitation beam and a Gaussian-imaged
 *  pinhole, giving the effective volume (hence absolute concentrations) and
 *  the diffusion autocorrelation without the 3-D Gaussian approximation.
 *  Port of chisurf's `core/fluorescence/fcs/enderlein.py` (itself a port of
 *  Fretica's FEnderleinMDF / FVeffEnderlein / FGdiffEnderlein), moved here
 *  because it is a *forward model* -- parameters in, model curve out, no
 *  measured data touched -- which the placement rule (owner, 2026-09-02)
 *  sends to bff rather than to the photon library. The photon library's
 *  `SimGrid` PSF builders serve its simulator and share the same w(z) law;
 *  they are the cited twin, not a copy of this.
 *
 *  The lateral integral of the MDF--propagator convolution is analytic; the
 *  axial double integral is one trapezoid over a fixed z grid (which only
 *  has to resolve the smooth kappa/w profile) times a Gauss--Hermite
 *  quadrature in the *difference* direction (z' = z + sqrt(4 D tau) xi),
 *  which resolves the propagator at any lag -- a fixed z' grid fails once
 *  sqrt(4 D tau) drops below its spacing.
 *
 *  All lengths are micrometres, `tau` seconds, `diffusion` um^2/s. The
 *  pinhole enters as the back-projected radius a = pinhole/2/magnification.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/bff_config.h>

#include <string>
#include <vector>

#include <IMP/bff/GraphNode.h>
#include <IMP/bff/GraphPort.h>

IMPBFF_BEGIN_NAMESPACE

//! Gauss--Hermite nodes and weights for \f$\int e^{-x^2} f(x)\,dx\f$.
/*!
    Newton-refined roots of the physicists' Hermite polynomial over the
    orthonormal recurrence (no overflow at any practical n), agreeing with
    the usual numpy reference to machine precision. Exposed because the
    quadrature is useful beyond this file and the library deliberately
    carries no numerical-recipes dependency.

    \param[in] n number of nodes
    \param[out] nodes,weights filled with n entries each, ascending nodes
*/
IMPBFFEXPORT void hermgauss(int n, std::vector<double>& nodes,
                            std::vector<double>& weights);

//! Effective detection volume of the Enderlein MDF (um^3).
/*!
    \f$V_\mathrm{eff} = \pi (\int \kappa\,dz)^2 / \int (\kappa^2/w^2)\,dz\f$
    by axial trapezoidal integration over `n_grid` points spanning `span`
    Rayleigh ranges.

    \param[in] w0 lateral 1/e^2 excitation waist at focus (um)
    \param[in] r0 emission-beam waist parameter at focus (um)
    \param[in] excitation_wavelength lambda_ex (um)
    \param[in] emission_wavelength lambda_em (um)
    \param[in] refractive_index n of the immersion/sample medium
    \param[in] pinhole_radius back-projected pinhole radius a (um)
    \param[in] n_grid axial grid points
    \param[in] span axial half-range in Rayleigh ranges
*/
IMPBFFEXPORT double fcs_mdf_effective_volume(
        double w0, double r0,
        double excitation_wavelength = 0.485,
        double emission_wavelength = 0.520,
        double refractive_index = 1.33,
        double pinhole_radius = 25.0 / 60.0,
        int n_grid = 4001,
        double span = 60.0);

//! Diffusion (cross-)correlation shape of the Enderlein MDF.
/*!
    One value per lag. `normalize` true returns the shape with the
    small-lag plateau at one (an auto-correlation's g(0) = 1; a two-focus
    cross-correlation with `separation` > 0 starts below one, which is the
    point -- the known focus distance calibrates an absolute D); false
    divides by the effective volume so g(0) = 1/V_eff.

    `tau <= 0` takes the analytic small-lag plateau, exactly as the
    reference implementation does.

    \param[in] tau lag times (s)
    \param[in] w0 lateral 1/e^2 excitation waist at focus (um)
    \param[in] r0 emission-beam waist parameter at focus (um)
    \param[in] diffusion translational diffusion coefficient (um^2/s)
    \param[in] excitation_wavelength lambda_ex (um)
    \param[in] emission_wavelength lambda_em (um)
    \param[in] refractive_index n of the medium
    \param[in] pinhole_radius back-projected pinhole radius a (um)
    \param[in] n_grid axial grid points of the trapezoid
    \param[in] span axial half-range in Rayleigh ranges
    \param[in] normalize plateau at one (true) or at 1/V_eff (false)
    \param[in] n_herm Gauss--Hermite nodes of the difference integral
    \param[in] separation two-focus distance d (um); 0 is autocorrelation
    \param[out] out_view,n_out_view the correlation, one value per lag
*/
IMPBFFEXPORT void fcs_mdf_g_diff(
        const std::vector<double>& tau,
        double w0, double r0, double diffusion,
        double excitation_wavelength = 0.485,
        double emission_wavelength = 0.520,
        double refractive_index = 1.33,
        double pinhole_radius = 25.0 / 60.0,
        int n_grid = 201,
        double span = 40.0,
        bool normalize = true,
        int n_herm = 40,
        double separation = 0.0,
        double** out_view = 0, int* n_out_view = 0);

//! The axial grid and the kappa / w^2 profiles evaluated on it.
/*!
    The half of `fcs_mdf_g_diff` that does not depend on the lag or on the
    diffusion coefficient, split out so a caller evaluating the same optics
    at many lags -- or at many diffusion coefficients, which is what an
    optimiser does -- builds it once. Exported rather than kept private
    because `FCSMdfCurve` is that caller and a second transcription of the
    grid is exactly how two code paths start disagreeing.

    \param[in] w0 lateral 1/e^2 excitation waist at focus (um)
    \param[in] r0 emission-beam waist parameter at focus (um)
    \param[in] excitation_wavelength lambda_ex (um)
    \param[in] emission_wavelength lambda_em (um)
    \param[in] refractive_index n of the medium
    \param[in] pinhole_radius back-projected pinhole radius a (um)
    \param[in] n_grid axial grid points
    \param[in] span axial half-range in Rayleigh ranges
    \param[out] z the grid, \p h its (uniform) spacing
    \param[out] kappa_z,w2_z kappa(z) and w(z)^2 on that grid
*/
IMPBFFEXPORT void fcs_mdf_axial_profiles(
        double w0, double r0, double excitation_wavelength,
        double emission_wavelength, double refractive_index,
        double pinhole_radius, int n_grid, double span,
        std::vector<double>& z, double& h, std::vector<double>& kappa_z,
        std::vector<double>& w2_z);

//! The raw (unnormalised) correlation at one lag, over precomputed grids.
/*!
    The inner loop of `fcs_mdf_g_diff`: the z trapezoid of kappa(z) times the
    Gauss--Hermite sum over the difference coordinate, with the analytic
    lateral factor. `fcs_mdf_g_diff` and `FCSMdfCurve` both call this, so
    there is one implementation of the integral and the two cannot drift.

    \param[in] tau lag time (s); clamped up to 1e-18
    \param[in] separation_squared d^2 (um^2); 0 is an autocorrelation
    \param[in] diffusion translational diffusion coefficient (um^2/s)
    \param[in] z,h the axial grid and its spacing
    \param[in] kappa_z,w2_z kappa(z) and w(z)^2 on that grid
    \param[in] herm_nodes,herm_weights the Gauss--Hermite quadrature
    \param[in] w0,r0 the waists (um), for the off-grid z' evaluations
    \param[in] excitation_wavelength,emission_wavelength wavelengths (um)
    \param[in] refractive_index,pinhole_radius the remaining optics (um)
*/
IMPBFFEXPORT double fcs_mdf_g_raw(
        double tau, double separation_squared, double diffusion,
        const std::vector<double>& z, double h,
        const std::vector<double>& kappa_z, const std::vector<double>& w2_z,
        const std::vector<double>& herm_nodes,
        const std::vector<double>& herm_weights,
        double w0, double r0, double excitation_wavelength,
        double emission_wavelength, double refractive_index,
        double pinhole_radius);

//! The Enderlein MDF diffusion shape as a node in a model graph.
/**
 * The `"mdf"` diffusion mode of an FCS model is a *numerical kernel*, not a
 * formula, so unlike the 3-D-Gaussian modes it cannot be handed to
 * `GraphExpression` as a string. This is the third option, and the same one the
 * decay path takes for its photophysics (\see SpectrumNode.h): a producer
 * node that publishes the shape, with an `GraphExpression` downstream
 * multiplying the terms that *are* formulas onto it --
 *
 *     FCSMdfCurve -> GraphExpression -> FitChiSquared -> FitMinimizer
 *
 * -- so a fit crosses into Python once per `run()` rather than once per
 * iteration.
 *
 * \par What it publishes, and what it deliberately does not
 * The **shape** g(tau), plateau at one. Not `b + g/N`: the bunching and
 * anticorrelation factors multiply *g* and the baseline `b` is added after
 * them, so a node that folded `b` in could not have those terms downstream
 * of it. `N`, `b` and the background factor stay variables of the
 * expression, which is where they compose correctly.
 *
 * \par The quadrature and the axial profiles are cached, and that is the
 * whole point
 * Rebuilding the Gauss--Hermite quadrature and the kappa/w^2 profiles on
 * every evaluation would make the graph route *slower* than the Python
 * director it replaces -- the kernel is already 96% native, so what is left
 * to win is per-iteration overhead, and this is where it would be spent
 * instead. The quadrature depends only on `n_herm` (set-once configuration)
 * and is built once per node; the profiles depend on the waists and the
 * optics and are rebuilt only when a waist actually moves, which during a
 * Levenberg-Marquardt Jacobian is 2 columns out of however many the model
 * has. `get_quadrature_builds()` / `get_profile_builds()` /
 * `get_kernel_evaluations()` report the counts, so a caller can *measure*
 * that the cache survives rather than assume it.
 *
 * Ports, created by build_ports():
 *
 * | port | what it is |
 * |---|---|
 * | `w0` | lateral 1/e^2 excitation waist at focus |
 * | `wem` | emission-beam waist parameter at focus (Enderlein's R0) |
 * | `D` | translational diffusion coefficient (um^2/s) |
 * | `diam` | two-focus separation; 0 is an autocorrelation |
 *
 * `w0`, `wem` and `diam` are read in the caller's length unit and scaled by
 * set_length_scale() -- an application that fits waists in nanometres keeps
 * fitting them in nanometres, exactly as `PolymerDistances` takes the
 * contour and persistence lengths a model fits rather than the kernel's
 * dimensionless ratio.
 *
 * \see GraphExpression, FitChiSquared, FitMinimizer, GraphNode
 */
class IMPBFFEXPORT FCSMdfCurve : public GraphNode {
 public:
  explicit FCSMdfCurve(const std::string& name = "fcs_mdf");

  //! Build the `w0`, `wem`, `D` and `diam` ports. Not in the constructor.
  /** A `GraphNode` that owns ports must already be held by a `shared_ptr`, and a
      Python-wrapped node is constructed before it is owned -- the same
      reason `FRETSpectrumNode::build_ports()` exists. */
  void build_ports();

  //! The lag axis, in **seconds**, set once from the data.
  void set_axis(const std::vector<double>& tau);
  //! Numpy in, for the axis a caller holds as an array.
  void set_axis_array(double* in_axis, int n_axis);
  const std::vector<double>& get_axis() const { return tau_; }

  //! The fixed optics: wavelengths and pinhole radius in micrometres.
  void set_optics(double excitation_wavelength, double emission_wavelength,
                  double refractive_index, double pinhole_radius);

  //! Quadrature resolution: axial points, axial half-range, Hermite nodes.
  void set_quadrature(int n_grid, double span, int n_herm);
  int get_n_grid() const { return n_grid_; }
  double get_span() const { return span_; }
  int get_n_herm() const { return n_herm_; }

  //! Micrometres per unit of the `w0` / `wem` / `diam` ports (1e-3 for nm).
  void set_length_scale(double micrometres_per_unit);
  double get_length_scale() const { return length_scale_; }

  //! Plateau at one (true, the default) or at 1/V_eff (false).
  void set_normalize(bool v);
  bool get_normalize() const { return normalize_; }

  //! The shape from the last evaluation.
  const std::vector<double>& get_curve() const { return curve_; }

  //! How many times the Gauss--Hermite quadrature has been built.
  int get_quadrature_builds() const { return quadrature_builds_; }
  //! How many times the axial kappa / w^2 profiles have been rebuilt.
  int get_profile_builds() const { return profile_builds_; }
  //! How many times the double integral has actually run.
  /** Lower than the number of evaluations whenever the ports came back to
      values the node has already integrated -- a Levenberg-Marquardt
      Jacobian differences one parameter at a time, so the columns belonging
      to `N` or to a bunching term leave every port of this node untouched
      and re-publish the cached shape. */
  int get_kernel_evaluations() const { return kernel_evaluations_; }

  void evaluate() override;

  std::string describe() const;

 private:
  //! Bring the cached quadrature up to date with `n_herm_`.
  void refresh_quadrature();
  //! Bring the cached axial profiles up to date with the waists and optics.
  void refresh_profiles(double w0, double r0);

  std::vector<double> tau_;
  std::vector<double> curve_;

  GraphPort* w0_port_ = nullptr;
  GraphPort* wem_port_ = nullptr;
  GraphPort* d_port_ = nullptr;
  GraphPort* diam_port_ = nullptr;

  double excitation_wavelength_ = 0.485;
  double emission_wavelength_ = 0.520;
  double refractive_index_ = 1.33;
  double pinhole_radius_ = 25.0 / 60.0;
  int n_grid_ = 201;
  double span_ = 40.0;
  int n_herm_ = 40;
  double length_scale_ = 1.0;
  bool normalize_ = true;

  //! The cached quadrature, and the order it was built for.
  std::vector<double> herm_nodes_, herm_weights_;
  int quadrature_n_ = 0;
  //! The cached axial profiles, and the arguments they were built from.
  std::vector<double> z_, kappa_z_, w2_z_;
  double h_ = 0.0;
  bool profiles_built_ = false;
  double profile_w0_ = 0.0, profile_r0_ = 0.0;
  double profile_lambda_ex_ = 0.0, profile_lambda_em_ = 0.0;
  double profile_n_ = 0.0, profile_a_ = 0.0;
  int profile_n_grid_ = 0;
  double profile_span_ = 0.0;
  //! The port values the published curve belongs to.
  bool curve_built_ = false;
  double curve_w0_ = 0.0, curve_r0_ = 0.0, curve_d_ = 0.0, curve_sep_ = 0.0;

  int quadrature_builds_ = 0;
  int profile_builds_ = 0;
  int kernel_evaluations_ = 0;
};

IMPBFF_END_NAMESPACE

 /* IMPBFF_FCS_MDF_H */

// -------- from FcsSaturation.h --------
/**
 *  (formerly IMP/bff/FcsSaturation.h, now a section of this file)
 *  \brief Saturated FCS forward model: steady-state photophysics on a grid.
 *
 *  Port of chisurf's `core/fluorescence/fcs/saturation.py` compute core —
 *  the numerically integrated, power-dependent FCS saturation shape. Moved
 *  here on the owner's placement ruling (2026-09-02): a forward model
 *  belongs in bff regardless of how fast its Python is; the measurement
 *  (1.0 ms/evaluation, BLAS-bound) set the port's priority, never its
 *  placement.
 *
 *  The pipeline: a Gaussian excitation profile scaled to the measured laser
 *  power becomes a spatial excitation rate; an N-state master-equation
 *  scheme (dark rates + excitation cross sections, `K[target, source]`,
 *  diagonals derived so columns sum to zero) is solved for its steady state
 *  at every grid point; the brightness-weighted populations form the
 *  emission profile; its spatial autocorrelation under free diffusion is
 *  evaluated in reciprocal space (0th-order Hankel transform along r — the
 *  quadrature matrix is cached on its grids, exactly as the Python's
 *  lru_cache did — and a real DFT along z, with the separable propagator
 *  factored so the exponentials are (n_kr + n_kz) * n_tau, not their
 *  product); and the photokinetic bunching factor X(tau) comes from the
 *  eigen-decomposition of the generator (Eigen::EigenSolver — the modes may
 *  be complex for a cyclic scheme, the sum is real).
 *
 *  All quantities SI: lengths m, times s, power W, rates Hz. The scheme
 *  matrices arrive flat row-major with their state count.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */


IMPBFF_BEGIN_NAMESPACE

//! Bessel J0, to machine precision, via the periodic trapezoid.
/*!
    \f$J_0(x) = \frac{1}{\pi}\int_0^\pi \cos(x\sin\theta)\,d\theta\f$; the
    integrand is smooth and periodic, so the trapezoid converges
    geometrically — 64 nodes give ~1e-15 for the |x| <= 40 this model
    reaches. Exposed because the library deliberately carries no special-
    functions dependency, and a truncated polynomial here would put a 1e-7
    floor under every curve parity test downstream.
*/
IMPBFFEXPORT double bessel_j0(double x);

//! Peak focal excitation rate \f$k_{exc}(0,0)\f$ in 1/s.
/*!
    \f$\sigma_{abs}\, 2\Phi/(\pi w_0^2)\f$ with
    \f$\Phi = P/(h c/\lambda)\f$ and
    \f$\sigma_{abs} = 10^{-4}\cdot 1000\ln(10)\,\epsilon/N_A\f$.
    \param[in] power_w total measured excitation power (W)
    \param[in] extinction molar extinction coefficient (1/M/cm)
    \param[in] w0 radial beam waist (1/e^2, m)
    \param[in] wavelength_m excitation wavelength (m)
*/
IMPBFFEXPORT double fcs_excitation_rate_peak(
        double power_w, double extinction, double w0,
        double wavelength_m = 488e-9);

//! Analytical 3D-Gaussian diffusion autocorrelation shape, G(0) = 1.
IMPBFFEXPORT void fcs_gaussian_g_diff(
        const std::vector<double>& tau,
        double w0, double z0, double diffusion,
        double** out_view = 0, int* n_out_view = 0);

//! Photokinetic state-relaxation (bunching) factor X(tau).
/*!
    \f$X(\tau) = q_a^T e^{K\tau} (q_b \circ p_{eq}) /
    ((q_a\cdot p_{eq})(q_b\cdot p_{eq}))\f$ for the generator
    \f$K = K_{dark} + k_{exc} K_{exc}\f$, via its eigen-decomposition —
    the eigenvalues ARE the relaxation rates a bunching fit reports.
    Degenerate schemes (no stationary state, dark scheme with no bright
    reachable state) return ones, as the reference does.

    \param[in] tau lag times (s)
    \param[in] k_exc_0 excitation rate to evaluate at (1/s)
    \param[in] dark_matrix flat row-major N x N dark rates (Hz),
               K[target, source], diagonal derived internally
    \param[in] exc_matrix flat row-major N x N excitation cross sections
    \param[in] n_states N
    \param[in] brightness per-state brightness of channel a (N)
    \param[in] brightness_b channel b; empty = autocorrelation (b = a)
    \param[out] out_view,n_out_view X(tau), one value per lag
*/
IMPBFFEXPORT void fcs_bunching_factor(
        const std::vector<double>& tau,
        double k_exc_0,
        const std::vector<double>& dark_matrix,
        const std::vector<double>& exc_matrix,
        int n_states,
        const std::vector<double>& brightness,
        const std::vector<double>& brightness_b = std::vector<double>(),
        double** out_view = 0, int* n_out_view = 0);

//! The numerically integrated saturated FCS diffusion shape G(tau).
/*!
    The steady-state emission profile on an (r, z) grid spanning five beam
    waists, spatially (cross-)correlated in reciprocal space, amplitude
    \f$V_0/V_{eff}\f$ (the saturation volume expansion), optionally times
    the bunching factor. At zero power the exact limit — the unsaturated
    Gaussian at amplitude one — is returned.

    \param[in] tau lag times (s)
    \param[in] power_w total measured excitation power (W)
    \param[in] extinction molar extinction coefficient (1/M/cm)
    \param[in] dark_matrix,exc_matrix,n_states the scheme, as in
               fcs_bunching_factor
    \param[in] brightness per-state brightness (N)
    \param[in] w0,z0 radial/axial beam waists (1/e^2, m)
    \param[in] diffusion diffusion coefficient (m^2/s)
    \param[in] include_bunching multiply by X(tau)
    \param[in] n_r,n_z radial/axial grid points
    \param[in] wavelength_m excitation wavelength (m)
    \param[in] brightness_b second channel; empty = autocorrelation
    \param[out] out_view,n_out_view G(tau), one value per lag
*/
IMPBFFEXPORT void fcs_saturated_curve_shape(
        const std::vector<double>& tau,
        double power_w,
        double extinction,
        const std::vector<double>& dark_matrix,
        const std::vector<double>& exc_matrix,
        int n_states,
        const std::vector<double>& brightness,
        double w0, double z0, double diffusion,
        bool include_bunching = true,
        int n_r = 120, int n_z = 40,
        double wavelength_m = 488e-9,
        const std::vector<double>& brightness_b = std::vector<double>(),
        double** out_view = 0, int* n_out_view = 0);

IMPBFF_END_NAMESPACE

 /* IMPBFF_FCS_SATURATION_H */

// -------- from FCSSaturationCurve.h --------
/**
 *  (formerly IMP/bff/FCSSaturationCurve.h, now a section of this file)
 *  \brief Saturated FCS forward model as a graph node.
 *
 *  The saturated FCS diffusion shape (the section above) wrapped as a
 *  `GraphNode`, so an FCS kinetics model in "full" mode joins the graph the way
 *  `FCSMdfCurve` does for the MDF mode.  The ports carry the quantities a
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

    \see the FcsSaturation section above, GraphNode, GraphExpression, FitChiSquared
*/
class IMPBFFEXPORT FCSSaturationCurve : public GraphNode {
 public:
  explicit FCSSaturationCurve(const std::string& name = "fcs_saturation");

  //! Build the scalar input ports.  Not in the constructor (see FCSMdfCurve).
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

  GraphPort* power_port_ = nullptr;
  GraphPort* extinction_port_ = nullptr;
  GraphPort* w0_port_ = nullptr;
  GraphPort* z0_port_ = nullptr;
  GraphPort* d_port_ = nullptr;
  GraphPort* n_port_ = nullptr;
  GraphPort* b_port_ = nullptr;
  GraphPort* bg_port_ = nullptr;

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

 /* IMPBFF_FCS_SATURATION_CURVE_H */

#endif  // IMPBFF_FCS_H
