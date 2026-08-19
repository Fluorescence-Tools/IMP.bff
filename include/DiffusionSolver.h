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
    \param[out] out_view,n_out_view the density after one step, as a numpy
                view over the kernel's buffer
*/
IMPBFFEXPORT void diffusion_step(
        const std::vector<double>& cur,
        const std::vector<double>& d,
        const std::vector<double>& decay,
        const std::vector<double>& bounds,
        int ng,
        int flux_form,
        double** out_view, int* n_out_view
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
    \param[out] out_view,n_out_view the final density
*/
IMPBFFEXPORT void diffusion_propagate(
        const std::vector<double>& cur,
        const std::vector<double>& d,
        const std::vector<double>& decay,
        const std::vector<double>& bounds,
        int ng,
        int flux_form,
        int n_steps,
        int n_out,
        std::vector<double>& fluorescence,
        double** out_view, int* n_out_view
);

//! Adjoint of diffusion_propagate(): the gradient of a loss on the fluorescence.
/*!
    For \f$L = \sum_k \bar F_k\, F(t_k)\f$, with \f$F(t_k)\f$ the populations
    diffusion_propagate() reports for the same \p cur, \p d, \p decay, \p bounds,
    \p ng, \p flux_form, \p n_steps, \p n_out, returns \f$\partial L/\partial d\f$,
    \f$\partial L/\partial \mathrm{decay}\f$ and \f$\partial L/\partial \mathrm{cur}\f$
    -- every voxel at once, for the cost of about two forward propagations.

    **Why this is possible without a tape.** One sweep is linear in the
    density, \f$p_{n+1} = A p_n\f$, with \f$A\f$ built from \p d, \p decay and
    \p bounds only. So the adjoint density \f$\bar p_n = \partial L/\partial p_n\f$
    obeys \f$\bar p_n = A^T \bar p_{n+1}\f$ (plus \f$\bar F_k\f$ on every voxel at a
    reported step, since \f$F\f$ is a plain sum) -- the transposed 7-point
    stencil, run backwards -- and the parameter gradients are local products
    of the forward state and the adjoint at each step. The Smoluchowski flux
    is symmetric in the interface coefficient, so its transpose is the same
    stencil up to the \p bounds mask and where \p decay multiplies; the Ito flux
    is not, and both are written out.

    **Memory.** The backward sweep needs the forward state at every step and
    the forward keeps none. Storing all of them is \f$n_\mathrm{steps}\, ng^3\f$
    doubles (2.6 GB at 4800 steps, ng = 41), so this checkpoints every
    \f$\lceil\sqrt{n_\mathrm{steps}}\rceil\f$ steps and re-runs each segment forward
    on the way back: about twice the forward cost, two segments' worth of
    memory.

    The gradients are of the **discrete** scheme actually run -- exact for it,
    which is what makes them agree with a central difference of
    diffusion_propagate() to roundoff (the dot-product identity is the test).
    They are with respect to the folded coefficients \p d (\f$D\,dt/dg^2\f$) and
    \p decay (\f$e^{-k\,dt}\f$); the chain rule to \f$D\f$ and \f$k\f$ is the
    caller's, one multiplication each.

    \param[in] cur,d,decay,bounds,ng,flux_form,n_steps,n_out as diffusion_propagate()
    \param[in] dL_dF sensitivity of the loss to each reported population,
               length n_steps / n_out + 1 (as diffusion_propagate() fills)
    \param[out] dL_dd \f$\partial L/\partial d\f$, flat ng^3
    \param[out] dL_ddecay \f$\partial L/\partial \mathrm{decay}\f$, flat ng^3
    \param[out] dL_dcur \f$\partial L/\partial \mathrm{cur}\f$ (the initial density), flat ng^3
*/
#ifndef SWIG
IMPBFFEXPORT void diffusion_propagate_adjoint(
        const std::vector<double>& cur,
        const std::vector<double>& d,
        const std::vector<double>& decay,
        const std::vector<double>& bounds,
        int ng,
        int flux_form,
        int n_steps,
        int n_out,
        const std::vector<double>& dL_dF,
        std::vector<double>& dL_dd,
        std::vector<double>& dL_ddecay,
        std::vector<double>& dL_dcur
);
#endif

//! The same, for the bindings: the three gradients concatenated in one array.
/*!
    Three ``ng^3`` out-parameters would be walked back element by element by
    SWIG (~340 ns each -- 70 ms for a 41^3 grid); one ``out_view`` is adopted
    by numpy without a copy. Layout: ``[dL/dd | dL/ddecay | dL/dcur]``, each
    ``ng^3`` long, so ``out.reshape(3, ng, ng, ng)`` separates them.
*/
IMPBFFEXPORT void diffusion_propagate_adjoint(
        const std::vector<double>& cur,
        const std::vector<double>& d,
        const std::vector<double>& decay,
        const std::vector<double>& bounds,
        int ng,
        int flux_form,
        int n_steps,
        int n_out,
        const std::vector<double>& dL_dF,
        double** out_view,
        int* n_out_view
);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_DIFFUSIONSOLVER_H
