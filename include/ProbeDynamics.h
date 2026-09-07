/**
 *  \file IMP/bff/ProbeDynamics.h
 *  \brief Langevin and Brownian dynamics of an attached label.
 *
 * The label is read as a molecule, not as a cloud: its own force field holds
 * it together, its backbone anchor is held still, and the protein around the
 * labelling site is a wall of soft spheres. What comes back is a trajectory --
 * every frame's coordinates, its potential energy, and (for dynamics) its
 * kinetic energy -- which is what a diffusion coefficient, an orientation
 * correlation or a rotamer library is measured from.
 *
 * \note Nothing here is a dye's. A fluorophore, a spin label and a quencher
 * are one object to this machinery -- a molecule on a linker with an anchor
 * held still -- and the class is named for what it does rather than for the
 * first thing it was pointed at.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_PROBEDYNAMICS_H
#define IMPBFF_PROBEDYNAMICS_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/ProbeLibrary.h>

#include <IMP/Pointer.h>
#include <IMP/ScoringFunction.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/Simulator.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Boltzmann's constant in kcal/mol/K, spelled once.
IMPBFFEXPORT double kb_kcal();

//! Give particles the decorators an integrator needs.
/*!
    `md` wants a velocity, `bd` wants a diffusion coefficient -- and the
    coefficient is Einstein's for the particle's own radius at this
    temperature, so a big atom diffuses more slowly than a small one without
    anybody saying so.

    \param[in] model,particles what to decorate
    \param[in] temperature K
    \param[in] integrator `md` or `bd`
    \param[in] radii one per particle; empty leaves an existing radius alone
               and gives an undecorated particle 1.7 A
*/
IMPBFFEXPORT void prepare_particles(
        IMP::Model* model, const IMP::ParticleIndexes& particles,
        double temperature, const std::string& integrator,
        const std::vector<double>& radii = std::vector<double>());

//! A simulator for \p mobile, decorated as the integrator needs.
/*!
    \param[in] model the model
    \param[in] mobile the particles that move; everything else stays where it
               is, because it is simply not in this list
    \param[in] scoring_function what the dynamics follows
    \param[in] integrator `md` (velocity Verlet with a Langevin thermostat) or
               `bd` (Brownian, overdamped)
    \param[in] temperature K
    \param[in] timestep_fs the step, fs
    \param[in] friction_ps the thermostat's friction, 1/ps; `md` only
    \param[in] seed negative leaves the generator alone
    \throw ValueException for any other \p integrator
*/
IMPBFFEXPORT IMP::atom::Simulator* make_langevin_simulator(
        IMP::Model* model, const IMP::ParticleIndexes& mobile,
        IMP::ScoringFunction* scoring_function,
        const std::string& integrator = "md", double temperature = 300.0,
        double timestep_fs = 2.0, double friction_ps = 10.0, int seed = -1);

//! What a run of #IMP::bff::AttachedProbeDynamics produced.
struct IMPBFFEXPORT LangevinTrajectory {
    //! Flat `n_frames * n_atoms * 3`.
    std::vector<double> coordinates;
    //! One per frame: the time, the potential energy, and the kinetic energy
    //! (NaN under `bd`, which has no velocities).
    std::vector<double> times_fs, potential_energy, kinetic_energy;
    std::vector<std::string> atom_names;
    std::string integrator;
    double temperature, timestep_fs;
    int n_frames, n_atoms;

    LangevinTrajectory()
        : integrator("md"), temperature(300.0), timestep_fs(2.0), n_frames(0),
          n_atoms(0) {}

    void get_coordinates(double** out_view, int* n_out_view) const;
    void get_times_fs(double** out_view, int* n_out_view) const;
    void get_potential_energy(double** out_view, int* n_out_view) const;
    void get_kinetic_energy(double** out_view, int* n_out_view) const;

    IMP_SHOWABLE_INLINE(LangevinTrajectory,
                        out << "LangevinTrajectory(" << n_frames << " frames of "
                            << n_atoms << " atoms)");
};
IMP_VALUES(LangevinTrajectory, LangevinTrajectories);

//! Langevin (`md`) or Brownian (`bd`) dynamics of an attached label.
/*!
    The label's own force field holds it together (#build_probe_restraints), its
    `N`/`CA`/`C`/`O` anchor is held still, and the protein's heavy atoms near
    the site repel it (#IMP::core::SoftSpherePairScore over a bipartite close
    pair container). The labelled residue is not an obstacle to its own label,
    and neither are hydrogens.
*/
class IMPBFFEXPORT AttachedProbeDynamics {
    IMP::Pointer<IMP::Model> model_;
    IMP::atom::Hierarchy protein_, dye_;
    ProbeForceFieldSystem system_;
    IMP::ParticleIndexes label_particles_, mobile_, fixed_, obstacles_;
    std::vector<std::string> site_ids_, atom_names_;
    IMP::Pointer<IMP::ScoringFunction> scoring_function_;
    IMP::Pointer<IMP::atom::Simulator> simulator_;
    std::string integrator_;
    double temperature_, timestep_fs_, friction_ps_;
    IMP::algebra::Vector3D site_ca_;

public:
    //! \param[in] protein_hier the labelled protein (the dye is already placed)
    /*! \param[in] label_hier the dye, as read from \p label_mol2
        \param[in] label_mol2 the dye's MOL2 -- its bonds, angles and torsions
        \param[in] chain,residue the labelled site, whose own atoms are not
                   obstacles
        \param[in] integrator `md` or `bd`
        \param[in] temperature K
        \param[in] timestep_fs the step; negative takes 2 fs for `md` and
                   0.5 fs for `bd`. Overdamped `bd` with the same stiff bonds
                   is stable only below about 2 fs, which is why the two
                   defaults differ.
        \param[in] friction_ps the thermostat's friction, 1/ps
        \param[in] interaction_sphere protein heavy atoms within this of the
                   site's CA are obstacles, A
        \param[in] repulsion_k the label-protein soft-sphere constant
        \param[in] seed negative leaves the generator alone
        \throw ValueException when the MOL2 and the hierarchy disagree about
               the atom count, when the site has no CA, or for an unknown
               integrator */
    AttachedProbeDynamics(IMP::atom::Hierarchy protein_hier,
                       IMP::atom::Hierarchy label_hier,
                       const std::string& label_mol2, const std::string& chain,
                       int residue, const std::string& integrator = "md",
                       double temperature = 300.0, double timestep_fs = -1.0,
                       double friction_ps = 10.0,
                       double interaction_sphere = 25.0,
                       double repulsion_k = 10.0, int seed = -1);

    //! The dye's coordinates as they stand, flat `n_atoms * 3`.
    void get_coordinates(double** out_view, int* n_out_view) const;
    //! The potential energy as it stands.
    double get_energy() const;
    //! Conjugate-gradient relaxation before dynamics; returns the score.
    double minimize(int n_steps = 200);

    //! Integrate, keeping a frame every \p write_every steps.
    /*! \param[in] n_steps how far to integrate
        \param[in] write_every how often to keep a frame */
    LangevinTrajectory run(int n_steps, int write_every = 10);

    //! T from the kinetic energy of the mobile atoms: \f$2E/(3Nk_B)\f$.
    double kinetic_temperature(double kinetic_energy) const;

    //! The label's atoms, those that move, those held, and the wall.
    /*! Particles rather than indexes: these are what a caller decorates and
        reads coordinates from, and `IMP::core::XYZ(p)` should be all it takes. */
    IMP::ParticlesTemp get_probe_particles() const;
    IMP::ParticlesTemp get_mobile() const;
    IMP::ParticlesTemp get_fixed() const;
    IMP::ParticlesTemp get_obstacles() const;
    std::vector<std::string> get_atom_names() const { return atom_names_; }
    IMP::ScoringFunction* get_scoring_function() const {
        return scoring_function_;
    }
    IMP::atom::Simulator* get_simulator() const { return simulator_; }
    ProbeForceFieldSystem get_system() const { return system_; }

    IMP_SHOWABLE_INLINE(AttachedProbeDynamics,
                        out << "AttachedProbeDynamics(" << integrator_ << ", "
                            << mobile_.size() << " mobile atoms)");
};

IMPBFF_END_NAMESPACE

#endif //IMPBFF_PROBEDYNAMICS_H
