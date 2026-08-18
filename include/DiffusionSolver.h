/**
 *  \file IMP/bff/DiffusionSolver.h
 *  \brief Explicit propagation of an excited-state density on an AV grid.
 *
 * Ported from Python by PRD-113: numba is a prototyping tool in this package,
 * not a runtime dependency, so every numerical kernel is C++.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_DIFFUSIONSOLVER_H
#define IMPBFF_DIFFUSIONSOLVER_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! How the flux between two voxels is discretised.
/*!
    The two differ exactly when the mobility varies in space, and the difference
    decides **where the dye sits at equilibrium**.
*/
enum FluxForm {
    //! \f$D_{ij}(p_i - p_j)\f$: equilibrium uniform, independent of D.
    /*! Equilibrium is thermodynamics and mobility is kinetics, so a dye slowed
        by friction with no attraction is still found everywhere it can reach.
        This is the structure the Haas-Steinberg equation is written in, so that
        its stationary state is the given p(r) for any D. */
    FLUX_SMOLUCHOWSKI = 0,
    //! \f$D_i p_i - D_j p_j\f$: equilibrium \f$p \propto 1/D\f$.
    /*! Inherited from ChiSurf. It makes a friction field act as an attractive
        potential -- halving the local mobility doubles the local population --
        which turned `slow_factor` into a disguised attraction to the quenchers.
        Kept so the old behaviour can be reproduced and compared. */
    FLUX_ITO = 1
};

//! One explicit Euler step of \f$\partial p/\partial t = \nabla\cdot(D\nabla p) - kp\f$.
/*!
    The rate term is applied as a **factor**, \p decay carrying
    \f$e^{-k\,dt}\f$, which is the exact solution of \f$dp/dt = -kp\f$ over the
    step. Subtracting \f$k\,dt\f$ instead goes negative and diverges once
    \f$k\,dt > 1\f$, which on a strongly quenched site it does -- and the
    divergence need not look like one, having once returned a smooth, finite,
    entirely plausible decay that was 2.6 % wrong.

    The outer shell is written to zero rather than left alone: the 7-point
    stencil cannot be evaluated there, and with ping-pong buffers those voxels
    would otherwise hold the state from two steps ago -- stale data that never
    decays and never diffuses, silently added into every population sum.

    \param[in] cur current density, flat ng^3
    \param[in] d mobility, already carrying \f$D\,dt/dg^2\f$
    \param[in] decay per-voxel \f$e^{-k\,dt}\f$
    \param[in] bounds non-zero inside the accessible volume
    \param[in] ng voxels per axis
    \param[in] flux_form FLUX_SMOLUCHOWSKI or FLUX_ITO
    \return the density after one step
*/
IMPBFFEXPORT std::vector<double> diffusion_step(
        const std::vector<double>& cur,
        const std::vector<double>& d,
        const std::vector<double>& decay,
        const std::vector<double>& bounds,
        int ng,
        int flux_form
);

//! Propagate many steps, reporting the surviving population periodically.
/*!
    **The loop belongs here, not in Python.** Calling diffusion_step() from a
    Python loop crosses the binding once per step and allocates a fresh grid
    each time; a solve is 1000-10000 steps, and doing that made the quenching
    suite forty times slower than the numba it replaced. Kept inside, the
    boundary is crossed once per solve and the buffers are reused.

    \param[in] cur initial density, flat ng^3
    \param[in] d mobility carrying $D\,dt/dg^2$
    \param[in] decay per-voxel $e^{-k\,dt}$
    \param[in] bounds non-zero inside the accessible volume
    \param[in] ng voxels per axis
    \param[in] flux_form FLUX_SMOLUCHOWSKI or FLUX_ITO
    \param[in] n_steps steps to take
    \param[in] n_out report the population every this many steps
    \param[out] fluorescence surviving population at each reported step
    eturn the final density
*/
IMPBFFEXPORT std::vector<double> diffusion_propagate(
        const std::vector<double>& cur,
        const std::vector<double>& d,
        const std::vector<double>& decay,
        const std::vector<double>& bounds,
        int ng,
        int flux_form,
        int n_steps,
        int n_out,
        std::vector<double>& fluorescence
);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_DIFFUSIONSOLVER_H
