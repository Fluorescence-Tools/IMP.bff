/**
 *  \file IMP/bff/BrownianWalk.h
 *  \brief A rejection-sampled Brownian walk inside an occupancy grid.
 *
 * The particle-based counterpart of #DiffusionSolver: the same dynamics, one
 * trajectory at a time instead of as a density. A walk resolves the dye's
 * *history* -- which is what a correlation function needs and what a
 * Fokker-Planck solve cannot give -- at the cost of sampling noise.
 *
 * Mobility enters as a **field**, one scaling per voxel, rather than as a
 * mask-and-scalar pair. Sticky patches, a smooth near-surface slow-down and a
 * uniform medium are then the same code with a different field, and the field
 * is what #slow_factor_grid already builds.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_BROWNIANWALK_H
#define IMPBFF_BROWNIANWALK_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! One Brownian trajectory of a point particle confined to an occupancy grid.
/*!
    Steps are isotropic Gaussians of standard deviation
    \f$\sqrt{2 D \Delta t}/dg\f$ **in voxel units** -- the variance of one
    Cartesian component, so that \f$\langle \Delta x^2\rangle = 2 D t\f$ and
    #DiffusionSolver, which is the same dynamics as a density, agrees at the
    same \f$D\f$. Scaled per step by
    \f$\sqrt{m}\f$ for the mobility \f$m\f$ at the particle's current voxel. A
    step that would leave the accessible region is **rejected**: the particle
    stays where it was and the frame is still emitted, so the trajectory has one
    entry per time step whether or not it moved. That is what makes a rejection
    walk sample the uniform equilibrium on the region -- dropping rejected
    frames would bias it away from the boundary.

    The starting voxel is drawn uniformly from the grid and retried up to 1000
    times until it lands inside the accessible region.

    \param[in] occupancy flat ng^3; nonzero marks accessible
    \param[in] mobility flat ng^3 scaling of the step variance, or empty for a
               uniform medium. Applied unconditionally -- a value above 1 speeds
               the particle up rather than being ignored, which the two Python
               kernels this replaces disagreed about.
    \param[in] ng voxels per axis
    \param[in] dg voxel edge, Angstrom
    \param[in] t_max total simulated time, ns
    \param[in] t_step time step, ns
    \param[in] diffusion_coefficient \f$D\f$, A^2/ns
    \param[in] seed reproducible when non-negative; drawn from the system
               otherwise
    \param[out] accepted per step: 1 if the step was taken, 0 if rejected
    \param[out] counts two entries, accepted and rejected
    \return the trajectory, three coordinates per step, in Angstrom **relative
            to the grid anchor** -- add the attachment point for the structure's
            frame. Empty if no accessible starting voxel was found.
*/
IMPBFFEXPORT std::vector<double> brownian_walk_in_volume(
        const std::vector<int>& occupancy,
        const std::vector<double>& mobility,
        int ng, double dg, double t_max, double t_step,
        double diffusion_coefficient, int seed,
        std::vector<int>& accepted,
        std::vector<int>& counts);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_BROWNIANWALK_H
