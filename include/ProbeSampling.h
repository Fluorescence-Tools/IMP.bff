#ifndef IMPBFF_PROBESAMPLING_H
#define IMPBFF_PROBESAMPLING_H

/**
 *  \file IMP/bff/ProbeSampling.h
 *  \brief Sampling a tethered dye: the walk, the photons, the frame writers --
 *         and the walk as an object.
 *
 *  Two sections: the **object** (formerly `ProbeDiffusion.h`:
 *  #IMP::bff::ProbeDiffusionSimulation, the dye's Brownian walk in its
 *  accessible volume with its fields kept on the C++ side) and the
 *  **kernels** (the former body of this file), described below.
 *
 * Three entry points over kernels that already exist, and the shapes their
 * callers want:
 *
 * - #IMP::bff::simulate_probe_diffusion — one rejection-sampled Brownian
 *   trajectory of the dye centre in its accessible volume;
 * - #IMP::bff::simulate_photon_trace and
 *   #IMP::bff::simulate_quenched_decay — what that trajectory emits, sampled
 *   photon by photon or integrated as a curve;
 * - #IMP::bff::equilibrium_occupancy — the stationary distribution in closed
 *   form, which is what the iterative solver spends tens of thousands of
 *   iterations converging to.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */

// -------- from ProbeDiffusion.h --------
/**
 *  (formerly IMP/bff/ProbeDiffusion.h, now a section of this file)
 *  \brief The dye's Brownian walk in its accessible volume, as an object.
 *
 * The particle picture of a tethered dye: a random walk confined to the
 * accessible volume, slowed where the dye touches the protein, reading a
 * per-voxel quenching rate as it goes. Its counterpart is the field picture in
 * #DiffusionSolver, which propagates an occupancy density instead of a
 * trajectory -- the two agree on equilibrium and on the mean squared
 * displacement, and only the particle picture can produce a correlation
 * function, because only it resolves the dye's history.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/bff_config.h>

#include <IMP/bff/Base.h>

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

    Not IMP::bff::parallel_threads(), which reports *OpenMP's* view and is 1 in
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

//! A dye's Brownian walk in its accessible volume, and the field it samples.
/*!
    \param density binary occupancy of the accessible volume, flat `ng^3`
    \param dg voxel edge, A
    \param x0 the grid anchor -- the attachment point
    \param slow_density binary contact grid, used with a scalar slow factor
    \param slow_factor_map per-voxel stickiness, used *instead* of the pair
           above when given
    \param quenching_rate_map per-voxel quenching rate, 1/ns, sampled along the
           trajectory to give k_quench()
*/
class IMPBFFEXPORT ProbeDiffusionSimulation {
    std::vector<int> density_;
    std::vector<int> slow_density_;
    std::vector<double> slow_factor_map_;
    std::vector<double> quenching_rate_map_;
    std::vector<double> trajectory_;     //!< flat, three per frame, structure frame
    std::vector<double> x0_;
    double dg_;
    double t_step_;
    int ng_;
    int n_accepted_;
    int n_rejected_;

public:
    ProbeDiffusionSimulation(
            const std::vector<int>& density = std::vector<int>(),
            double dg = 0.5,
            const std::vector<double>& x0 = std::vector<double>(),
            const std::vector<int>& slow_density = std::vector<int>(),
            const std::vector<double>& slow_factor_map = std::vector<double>(),
            const std::vector<double>& quenching_rate_map = std::vector<double>());

    //! Simulate the walk; the trajectory is in the structure's frame.
    /*!
        \p n_trajectories independent walks are run and **concatenated**, not
        averaged: the photon Monte-Carlo downstream draws a random start frame
        per photon, so a longer record is exactly what it wants. They run on
        `std::thread` rather than OpenMP, because this build's
        `OpenMP_CXX_FLAGS` is empty and every `#pragma omp` in the module is
        inert (see okf/validation/openmp_is_not_enabled.md).

        \return the number of frames, or 0 when no walk found a starting point
                -- an empty accessible volume
    */
    int run(double D = 40.0, double slow_fact = 0.01, double t_step = 0.002,
            double t_max = 10000.0, int n_trajectories = 1,
            int random_seed = -1);

    //! `(n_frames, 3)` positions in Angstrom, as a numpy view.
    void get_trajectory(double** out_view, int* n_out_view) const;
    void get_x0(double** out_view, int* n_out_view) const;
    //! Mean position over the trajectory, three values.
    void get_mean_position(double** out_view, int* n_out_view) const;

    int get_n_frames() const { return static_cast<int>(trajectory_.size() / 3); }
    int get_n_accepted() const { return n_accepted_; }
    int get_n_rejected() const { return n_rejected_; }
    double get_dg() const { return dg_; }
    double get_t_step() const { return t_step_; }
    //! The fused decay path builds a walk it does not run, and states the step
    //! itself -- the object has to accept that rather than only learn it from
    //! run().
    void set_t_step(double t) { t_step_ = t; }

    //! Record the outcome of a walk this object did not run.
    /*! The fused decay kernel walks and races photons in one pass, so no
        trajectory is ever materialised -- but it does report how many steps
        were accepted and rejected, and the walk record is where callers look
        for that. One named operation rather than three setters, because
        "adopt this outcome" is what is happening and "assign to n_accepted"
        is not. */
    void set_step_counts(int n_accepted, int n_rejected) {
        n_accepted_ = n_accepted;
        n_rejected_ = n_rejected;
    }

    //! The occupancy grid this walk is confined to.
    void get_density(double** out_view, int* n_out_view) const;

    //! Confine the walk to a different volume.
    /*! Any trajectory already recorded is discarded: it was a walk in the old
        volume, and keeping it would let a caller read positions the new
        occupancy forbids. */
    void set_density(const std::vector<int>& density);
    //! The stored mobility field; a zero-length view when none was given.
    void get_slow_factor_map(double** out_view, int* n_out_view) const;
    //! The stored contact grid, as doubles; zero-length when none was given.
    void get_slow_density(double** out_view, int* n_out_view) const;
    //! The stored quenching-rate grid; zero-length when none was given.
    void get_quenching_rate_map(double** out_view, int* n_out_view) const;
    int get_ng() const { return ng_; }

    //! Read a per-voxel field along the trajectory.
    /*!
        The index map here must be the one the grids were stamped with, or the
        walk reads rates from beside where the quenchers were placed. Two
        details carry that, and both were bugs once: the centre offset is the
        *integer* grid_center_index(), and the conversion is `floor`, not
        `trunc` -- `trunc` maps `[-1, 0)` to 0, so a position up to one voxel
        below the grid would be treated as inside it and read voxel 0's value.

        Frames outside the grid read zero.
    */
    void sample_grid(const std::vector<double>& grid, int ng,
                     double** out_view, int* n_out_view) const;

    //! The quenching rate the dye experiences, frame by frame, 1/ns.
    //! The quenching rate the dye sees, frame by frame, 1/ns.
    /*! Rounded through `float`: the fused kernel holds its own trace at that
        width, so the two paths see bit-identical rates. See the note in
        get_k_quench()'s body. */
    void get_k_quench(double** out_view, int* n_out_view) const;

    //! The fraction of frames spent in contact with a quencher.
    double get_collision_fraction() const;

    IMP_SHOWABLE_INLINE(ProbeDiffusionSimulation,
                        out << "ProbeDiffusionSimulation(" << get_n_frames()
                            << " frames, dg = " << dg_ << " A)");
};
IMP_VALUES(ProbeDiffusionSimulation, ProbeDiffusionSimulations);

IMPBFF_END_NAMESPACE

// -------- from ProbeSampling.h --------
#include <IMP/bff/RotamerLibrary.h>


#include <IMP/atom/Hierarchy.h>

#include <string>

IMPBFF_BEGIN_NAMESPACE

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
IMPBFFEXPORT RotamerLibrary load_rotamer_library_trajectory(
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
/*!
    The hierarchy's leaves must be the dye's atoms and \p coords three values
    per leaf. \throw ValueException on an atom-count mismatch.
*/
IMPBFFEXPORT void apply_coordinates(const IMP::atom::Hierarchy hierarchy,
                                            const std::vector<double>& coords);

// The two frame writers that were here are in `bin/imp_bff`, with the one
// command that attaches them. They are `IMP::OptimizerState`s, so they could be
// C++ -- but writing RMF needs `IMP.rmf`, and widening this module's dependency
// graph for a program's convenience is the wrong trade: `bin/` is where a
// program lives, and its Python already imports RMF.

IMPBFF_END_NAMESPACE


#endif  // IMPBFF_PROBESAMPLING_H
