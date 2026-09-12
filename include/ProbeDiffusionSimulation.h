/**
 *  \file IMP/bff/ProbeDiffusionSimulation.h
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
#ifndef IMPBFF_PROBE_DIFFUSION_SIMULATION_H
#define IMPBFF_PROBE_DIFFUSION_SIMULATION_H

#include <IMP/bff/ProbeSimulation.h>

IMPBFF_BEGIN_NAMESPACE

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
class IMPBFFEXPORT ProbeDiffusionSimulation : public ProbeSimulation {
    std::vector<int> density_;
    std::vector<int> slow_density_;
    std::vector<double> slow_factor_map_;
    std::vector<double> quenching_rate_map_;
    std::vector<double> trajectory_;     //!< flat, three per frame, structure frame
    std::vector<double> x0_;
    double dg_;
    double t_step_;
    //! What #run() drives #simulate() with; #set_parameters changes them.
    double diffusion_coefficient_, slow_fact_, t_max_;
    int n_trajectories_, random_seed_;
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
    int simulate(double D = 40.0, double slow_fact = 0.01, double t_step = 0.002,
                 double t_max = 10000.0, int n_trajectories = 1,
                 int random_seed = -1);

    // ---- IMP::bff::ProbeSimulation: the shape every simulation shares ----
    //
    // The walk is one particle on a grid, so `get_n_atoms()` is 1 and the
    // positions are its own. It has no energy and no minimiser, and says so
    // rather than returning a polite zero. `run()` drives #simulate() with
    // the parameters in force, taking `n_steps * t_step` as the time to
    // cover; `write_every` thins the record the walk already produces.

    std::string get_type() const override { return "grid-diffusion"; }
    std::string get_parameters() const override;
    void set_parameters(const std::string& json) override;
    int get_n_atoms() const override { return 1; }
    void get_positions(double** out_view, int* n_out_view) const override;
    void set_positions(const std::vector<double>& xyz) override;
    void step(int n_steps) override;
    ProbeSimulationTrajectory run(int n_steps, int write_every = 1) override;

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

#endif  // IMPBFF_PROBE_DIFFUSION_SIMULATION_H
