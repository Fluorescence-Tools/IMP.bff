/**
 *  \file IMP/bff/FcsSaturation.h
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
#ifndef IMPBFF_FCS_SATURATION_H
#define IMPBFF_FCS_SATURATION_H

#include <IMP/bff/bff_config.h>

#include <vector>

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

#endif /* IMPBFF_FCS_SATURATION_H */
