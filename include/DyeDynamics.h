/**
 *  \file IMP/bff/DyeDynamics.h
 *  \brief Dye dynamics on a labelled structure, without naming IMP.
 *
 * The dye roads that run on an atomistic model -- placing a dye on a residue
 * and then integrating its motion -- are built on IMP: its hierarchies, its
 * decorators, its integrators (`ProbeAttachment.h`, `ProbeDynamics.h`). That
 * is a fine way to *implement* them and a poor way to *offer* them: a caller
 * who installed one package should not need IMP's Python to run a dye.
 *
 * So this header is the same computation behind file paths and arrays. It is
 * part of the connection layer -- it exists only where IMP is linked -- but
 * nothing in it names an IMP type, which is what lets it be wrapped without
 * IMP's own SWIG interfaces and shipped in a wheel that carries IMP as a
 * private library.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_DYEDYNAMICS_H
#define IMPBFF_DYEDYNAMICS_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/Base.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

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

//! Place a dye on a residue and write the labelled structure.
/*!
    The dye is superimposed on the site's backbone frame and its clashes with
    the residue's own side chain are resolved the way #IMP::bff::attach_probes
    does -- this is that function, reached by file rather than by hierarchy.

    \param[in] protein_pdb the structure to label
    \param[in] dye_pdb the dye, as an atomistic model
    \param[in] chain,residue where to put it
    \param[in] out_pdb where the labelled structure goes; empty writes nothing
    \param[in] strip_site_sidechain remove the labelled residue's side chain
    \return the number of atoms stripped from the site
    \throw ValueException when the site is not in the structure
    \throw IOException when a file cannot be read or written
*/
IMPBFFEXPORT int attach_dye_to_pdb(
        const std::string& protein_pdb, const std::string& dye_pdb,
        const std::string& chain, int residue,
        const std::string& out_pdb = "", bool strip_site_sidechain = true);

//! Langevin (`md`) or Brownian (`bd`) dynamics of a dye on a structure.
/*!
    The whole road in one call: read the structure and the dye, place the dye
    at the site, build the dye's force field from its MOL2, hold the anchor
    still, let the protein's heavy atoms near the site repel it, relax, and
    integrate. What comes back is #IMP::bff::LangevinTrajectory -- the frames,
    the times and the two energies, as arrays.

    \param[in] protein_pdb,dye_pdb the structure and the dye
    \param[in] dye_mol2 the dye's MOL2: its bonds, angles and torsions
    \param[in] chain,residue the labelled site
    \param[in] n_steps how long to integrate
    \param[in] write_every frames are kept every so many steps
    \param[in] minimize_steps conjugate-gradient relaxation first; 0 skips it
    \param[in] integrator `md` (velocity Verlet with a Langevin thermostat)
               or `bd` (Brownian, overdamped)
    \param[in] temperature K
    \param[in] timestep_fs the step, fs; negative picks the integrator's own
    \param[in] friction_ps the thermostat's friction, 1/ps; `md` only
    \param[in] interaction_sphere how far around the site the protein repels, A
    \param[in] repulsion_k the repulsion's force constant
    \param[in] seed negative leaves the generator alone
    \param[in] strip_site_sidechain remove the labelled residue's side chain
    \throw ValueException for an unknown integrator or a site that is not there
    \throw IOException when a file cannot be read
*/
IMPBFFEXPORT LangevinTrajectory run_dye_langevin(
        const std::string& protein_pdb, const std::string& dye_pdb,
        const std::string& dye_mol2, const std::string& chain, int residue,
        int n_steps = 3000, int write_every = 100, int minimize_steps = 200,
        const std::string& integrator = "md", double temperature = 300.0,
        double timestep_fs = -1.0, double friction_ps = 10.0,
        double interaction_sphere = 25.0, double repulsion_k = 10.0,
        int seed = -1, bool strip_site_sidechain = true);

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_DYEDYNAMICS_H
