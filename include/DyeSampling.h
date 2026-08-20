/**
 *  \file IMP/bff/DyeSampling.h
 *  \brief Sampling a tethered dye: the walk, the photons, and the frame writers.
 *
 * Three entry points over kernels that already exist, and the shapes their
 * callers want:
 *
 * - #IMP::bff::simulate_dye_diffusion — one rejection-sampled Brownian
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
#ifndef IMPBFF_DYESAMPLING_H
#define IMPBFF_DYESAMPLING_H

#include <IMP/bff/bff_config.h>

#include <IMP/showable_macros.h>
#include <IMP/value_macros.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! A Brownian trajectory of the dye centre, in Angstrom.
/*! The coordinates are relative to the grid anchor: add the attachment point to
    place them in the structure's frame. */
class IMPBFFEXPORT DyeDiffusionTrajectory {
    std::vector<double> xyz_;      //!< flat, three per frame
    std::vector<int> accepted_;    //!< one per frame, 0 or 1
    int n_accepted_, n_rejected_;

public:
    DyeDiffusionTrajectory() : n_accepted_(0), n_rejected_(0) {}
    DyeDiffusionTrajectory(const std::vector<double>& xyz,
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

    IMP_SHOWABLE_INLINE(DyeDiffusionTrajectory,
                        out << "DyeDiffusionTrajectory(" << get_n_frames()
                            << " frames, " << get_acceptance_ratio()
                            << " accepted)");
};
IMP_VALUES(DyeDiffusionTrajectory, DyeDiffusionTrajectories);

//! One Brownian trajectory of the dye centre in its accessible volume.
/*!
    Mobility is a **field**: one scaling per voxel. The scalar-plus-mask form is
    one way to build it, and building it here is what collapsed two
    near-identical kernels into one.

    \param[in] density,n_density binary occupancy of the volume, flat `ng^3`
    \param[in] mobility,n_mobility per-voxel mobility, or length 0 for a uniform
               medium
    \param[in] ng grid side
    \param[in] dg voxel edge, A
    \param[in] t_max,t_step total time and step, ns
    \param[in] diffusion_coefficient \f$A^2/ns\f$, in the standard convention:
               the per-Cartesian-component step variance is \f$2Ddt\f$, so
               \f$\langle dx^2\rangle = 2Dt\f$ — the same D
               #IMP::bff::GridDiffusionSolver takes.

               The width was \f$\sqrt{2D \cdot 3 \cdot dt}\f$ per component
               until 2026-08-18 — the three-dimensional MSD used as one
               component's width, so the walk diffused at 3D. The tempting
               inference, that D had been calibrated around it, is false: no
               calibration of D exists anywhere in this stack, and every
               document calls these uncalibrated transferable starting values.
    \param[in] seed negative draws freely
*/
IMPBFFEXPORT DyeDiffusionTrajectory simulate_dye_diffusion(
        int* density, int n_density, double* mobility, int n_mobility, int ng,
        double dg = 0.5, double t_max = 10000.0, double t_step = 0.002,
        double diffusion_coefficient = 40.0, int seed = -1);

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
        double* diffusion_map, int n_diffusion_map, int* bounds, int n_bounds,
        const std::string& flux_form, double** out_view, int* n_out_view);

//! A reference rotamer library: conformer coordinates, weights and atom names.
struct IMPBFFEXPORT RotamerLibrary {
    //! Flat, `n_frames * n_atoms * 3`.
    std::vector<double> coords;
    //! One per frame, normalised to sum 1.
    std::vector<double> weights;
    std::vector<std::string> atom_names;
    int n_frames, n_atoms;

    RotamerLibrary() : n_frames(0), n_atoms(0) {}

    void get_coords(double** out_view, int* n_out_view) const;
    void get_weights(double** out_view, int* n_out_view) const;

    IMP_SHOWABLE_INLINE(RotamerLibrary,
                        out << "RotamerLibrary(" << n_frames << " x "
                            << n_atoms << " atoms)");
};
IMP_VALUES(RotamerLibrary, RotamerLibraries);

//! Load a reference rotamer library from a PDB plus a trajectory.
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
IMPBFFEXPORT RotamerLibrary load_rotamer_library(
        const std::string& pdb_path, const std::string& trajectory_path,
        const std::string& weights_path = "", int max_frames = -1);

//! `(pdb, trajectory, weights)` for a reference dye+linker name.
/*! The weights entry is empty when there is no weights file. \throw
    IOException when the PDB or the trajectory is missing. */
IMPBFFEXPORT std::vector<std::string> find_reference_rotamer_files(
        const std::string& lib_dir, const std::string& dye_name,
        int cutoff = 30);

//! Draw a rotamer index according to the weights.
/*! \param[in] weights need not be normalised
    \param[in] seed negative draws freely */
IMPBFFEXPORT int sample_rotamer_index(const std::vector<double>& weights,
                                      int seed = -1);

// The two frame writers that were here are in `bin/imp_bff`, with the one
// command that attaches them. They are `IMP::OptimizerState`s, so they could be
// C++ -- but writing RMF needs `IMP.rmf`, and widening this module's dependency
// graph for a program's convenience is the wrong trade: `bin/` is where a
// program lives, and its Python already imports RMF.

IMPBFF_END_NAMESPACE

#endif //IMPBFF_DYESAMPLING_H
