/**
 *  \file IMP/bff/QuenchedDecay.h
 *  \brief Walk, read the quenching rate along the walk, and race photons —
 *         without the trajectory ever leaving C++.
 *
 * The three steps of a PET-quenching prediction used to be three calls with a
 * Python array between each: the walk returns a trajectory, Python reads the
 * rate map along it, the photon simulator races against that trace. The
 * trajectory is the largest object in the chain and **nobody wants it** — a
 * 5 000 000-step walk is 20 million doubles crossing the boundary twice to
 * produce a few thousand photons.
 *
 * This is the same computation with the intermediates kept where they are made.
 * It is not an approximation and not a different model: given the same seeds it
 * reproduces the three-call path exactly, which is what
 * `test/quenching/test_fused_decay.py` asserts.
 *
 * The separate calls remain, and should: a trajectory is what a correlation
 * function or a visualisation needs, and the rate trace is worth looking at on
 * its own. This is the path for when only the decay is wanted, which is the
 * common case and the one inside every fitting loop.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_QUENCHEDDECAY_H
#define IMPBFF_QUENCHEDDECAY_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Photon trace of a dye diffusing in its accessible volume, fused.
/*!
    Runs \p walk_seeds.size() independent walks, concatenates the quenching rate
    each one saw, and races \p n_photons excitations against the result.
    Concatenated rather than averaged, because each photon draws its own start
    frame — a longer record is exactly what that wants.

    The rate is read at the voxel the walk is standing on, with no coordinate
    round trip. That is the same voxel the Python path reaches by
    `floor((xyz - x0)/dg + (ng-1)//2)`: the attachment point cancels, and the
    walk's own occupancy test already resolved the index.

    \param[in] occupancy flat ng^3, nonzero where the dye may be
    \param[in] mobility flat ng^3 step-variance scaling, or empty
    \param[in] rate_map flat ng^3 quenching rate, 1/ns
    \param[in] ng,dg grid size and voxel edge
    \param[in] t_max,t_step walk duration and step, ns
    \param[in] diffusion_coefficient A^2/ns, per-component variance 2 D dt
    \param[in] walk_seeds one per independent trajectory
    \param[in] tau0 intrinsic lifetime, ns
    \param[in] n_photons excitation events
    \param[in] photon_seed reproducible when non-negative
    \param[out] stats five values: frames, accepted, rejected, mean rate,
                collision fraction — everything the caller would have derived
                from the trajectory it no longer receives
    \return two values per photon: delay time in ns (0 if quenched), and 1 or 0
            for emitted. Empty if no walk found a starting voxel.
*/
IMPBFFEXPORT std::vector<double> quenched_donor_photons(
        const std::vector<int>& occupancy,
        const std::vector<double>& mobility,
        const std::vector<double>& rate_map,
        int ng, double dg, double t_max, double t_step,
        double diffusion_coefficient,
        const std::vector<int>& walk_seeds,
        double tau0, int n_photons, int photon_seed,
        std::vector<double>& stats);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_QUENCHEDDECAY_H
