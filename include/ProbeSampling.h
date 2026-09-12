/**
 * \file IMP/bff/ProbeSampling.h
 * \brief Probe sampling kernels, trajectories and reference-library helpers.
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#ifndef IMPBFF_PROBESAMPLING_H
#define IMPBFF_PROBESAMPLING_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/IMPCompatibility.h>
#include <IMP/bff/ProbeRotamerLibrary.h>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Where voxel \p i sits: `x0 + (i - grid_center_index(ng)) * dg`.
/*!
    The integer `(ng - 1) / 2`, not the float `(ng - 1) / 2.0`. They differ on
    every even edge length, and even is the normal case -- a grid stamped with
    one convention and read with the other returns rates from beside where the
    quenchers were placed. One function, so the two cannot drift.
*/
inline int grid_center_index(int ng) { return (ng - 1) / 2; }

//! How many walks to run when the caller does not say.
/*!
    Capped rather than one-per-core: the walks are independent and each holds
    its own trajectory, so the memory grows with the count while the useful
    sampling saturates long before a large machine's core count does.

    Not IMP::bff::openmp_thread_count(), which reports *OpenMP's* view and is 1 in
    this build because `OpenMP_CXX_FLAGS` is empty. These walks are
    `std::thread`s and run whatever the hardware has.
*/
IMPBFFEXPORT int default_trajectory_count();

//! The cap default_trajectory_count() and resolve_trajectory_count() obey.
extern IMPBFFEXPORT const int MAX_PARALLEL_TRAJECTORIES;

//! How many walks to run for a request of \p requested.
/*! Negative asks for the default; anything else is clamped to
    [1, MAX_PARALLEL_TRAJECTORIES]. */
IMPBFFEXPORT int resolve_trajectory_count(int requested);

//! One seed per trajectory, derived from a base seed or drawn freshly.
/*!
    Derived seeds step by a large prime rather than by one, so trajectories from
    a single base do not share the low-order pattern `base + i` would give them.

    \param[in] random_seed negative draws freshly, and a single trajectory then
               gets -1 -- "draw freely" -- rather than a seed of its own
    \param[in] n_trajectories how many seeds are wanted
*/
IMPBFFEXPORT std::vector<int> trajectory_seeds(int random_seed,
                                               int n_trajectories);

//! A Brownian trajectory of the dye centre, in Angstrom.
/*! The coordinates are relative to the grid anchor: add the attachment point to
    place them in the structure's frame. */
class IMPBFFEXPORT ProbeDiffusionTrajectory {
    std::vector<double> xyz_;      //!< flat, three per frame
    std::vector<int> accepted_;    //!< one per frame, 0 or 1
    int n_accepted_, n_rejected_;

public:
    ProbeDiffusionTrajectory() : n_accepted_(0), n_rejected_(0) {}
    ProbeDiffusionTrajectory(const std::vector<double>& xyz,
                           const std::vector<int>& accepted, int n_accepted,
                           int n_rejected)
        : xyz_(xyz), accepted_(accepted), n_accepted_(n_accepted),
          n_rejected_(n_rejected) {}

    void get_xyz(double** out_view, int* n_out_view) const;
    void get_accepted(int** out_view_i, int* n_out_view_i) const;

    int get_n_frames() const { return static_cast<int>(xyz_.size() / 3); }
    int get_n_accepted() const { return n_accepted_; }
    int get_n_rejected() const { return n_rejected_; }
    //! Accepted steps as a fraction of attempted; zero for an empty walk.
    double get_acceptance_ratio() const;

    IMP_SHOWABLE_INLINE(ProbeDiffusionTrajectory,
                        out << "ProbeDiffusionTrajectory(" << get_n_frames()
                            << " frames, " << get_acceptance_ratio()
                            << " accepted)");
};
IMP_VALUES(ProbeDiffusionTrajectory, ProbeDiffusionTrajectories);

//! One Brownian trajectory of the dye centre in its accessible volume.
/*!
    Mobility is a **field**: one per-voxel scaling of the step variance. This
    function builds that field from the two cheaper spellings a caller has and
    runs the walk on it:

    - \p slow_fact with more than one entry is the field directly;
    - otherwise \p slow_fact (empty means 1.0) is a scalar, applied wherever
      \p slow_density is nonzero; an empty \p slow_density or a scalar of 1.0
      means a uniform medium.

    \warning The field is built here rather than passed in finished, so a
    caller cannot hand this a field built by a second, divergent spelling of
    the same construction.

    \param[in] density flat `ng^3` binary occupancy of the volume; `ng` is its
               cube root
    \param[in] slow_density flat `ng^3` binary contact grid, or empty for none;
               used only when \p slow_fact is not itself a field
    \param[in] dg voxel edge, A
    \param[in] t_max,t_step total time and step, ns
    \param[in] D diffusion coefficient, A^2/ns, in the standard convention: the
               per-Cartesian-component step variance is `2 D dt`, so
               `\langle dx^2\rangle = 2Dt`. The width was
               `sqrt(2D * 3 * dt)` per component until 2026-08-18 -- the
               three-dimensional MSD used as one component's width, so the walk
               diffused at 3D. No calibration of D exists anywhere in this
               stack; the documents call these uncalibrated transferable
               starting values.
    \param[in] slow_fact scalar or per-voxel field, see the mobility rule above
    \param[in] random_seed negative draws freely
    \throw ValueException when the density has no integer cube root
*/
IMPBFFEXPORT ProbeDiffusionTrajectory simulate_probe_diffusion(
        const std::vector<int>& density,
        const std::vector<int>& slow_density = std::vector<int>(), double dg = 0.5,
        double t_max = 10000.0, double t_step = 0.002, double D = 40.0,
        const std::vector<double>& slow_fact = std::vector<double>(),
        int random_seed = -1);

//! The stationary occupancy of the volume, in closed form.
/*!
    Which closed form depends on how the flux is discretised, and the two
    disagree completely:

    - `"smoluchowski"` (default, and the physics): \f$p_{eq} \propto 1\f$ on the
      accessible domain — **independent of D**.
    - `"ito"` (the inherited ChiSurf behaviour): \f$p_{eq} \propto 1/D\f$, so the
      dye piles up wherever it moves slowly.

    **The default is Smoluchowski because equilibrium is thermodynamics and
    mobility is kinetics.** A dye slowed by friction near a surface, with no
    attractive interaction, must still be found uniformly across its accessible
    volume at equilibrium — it merely takes longer to get around. Letting a
    friction field set the distribution asserts a potential that was never
    specified. On T4L site 132 the difference is a factor of ~75 between peak
    and mean occupancy.

    \throw ValueException for any other \p flux_form
*/
IMPBFFEXPORT void equilibrium_occupancy(
        const std::vector<double>& diffusion_map,
        const std::vector<int>& bounds,
        const std::string& flux_form = "smoluchowski", double** out_view = NULL,
        int* n_out_view = NULL);

//! Load a rotamer library from a PDB plus a trajectory.
/*!
    Atom names come from the PDB and coordinates from the trajectory, so this
    needs nothing beyond IMP.

    A weights file is one number per line. When it is shorter or longer than the
    trajectory the two are truncated to the shorter -- a library whose weights
    and frames disagree is a library whose extra frames have no weight, and
    dropping them is the only reading that does not invent one. Absent, or
    summing to zero, the frames are weighted uniformly.

    \throw IOException when the trajectory has no frames
*/
IMPBFFEXPORT ProbeRotamerLibrary load_rotamer_library_trajectory(
        const std::string& pdb_path, const std::string& trajectory_path,
        const std::string& weights_path = "", int max_frames = -1);

//! `(pdb, trajectory, weights)` for a reference dye+linker name.
/*! The weights entry is empty when there is no weights file. \throw
    IOException when the PDB or the trajectory is missing. */
IMPBFFEXPORT std::vector<std::string> get_reference_rotamer_files(
        const std::string& lib_dir, const std::string& probe_name,
        int cutoff = 30);

//! Draw an index in proportion to `weights`.
/*! Nothing about the draw is particular to a rotamer: it is the weighted
    choice any ensemble makes when it picks a member. */
/*! \param[in] weights need not be normalised
    \param[in] seed negative draws freely */
IMPBFFEXPORT int sample_weighted_index(const std::vector<double>& weights,
                                      int seed = -1);

//! Put a coordinate set onto a hierarchy's atoms, in place.
/*! The coordinates are one conformer's, flat and in the hierarchy's own atom
    order; the hierarchy may be a dye, a side chain or anything else whose
    atoms that order describes. */

// The two frame writers that were here are in `bin/imp_bff`, with the one
// command that attaches them. They are `IMP::OptimizerState`s, so they could be
// C++ -- but writing RMF needs `IMP.rmf`, and widening this module's dependency
// graph for a program's convenience is the wrong trade: `bin/` is where a
// program lives, and its Python already imports RMF.

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_PROBESAMPLING_H
