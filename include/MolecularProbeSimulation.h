/**
 *  \file IMP/bff/MolecularProbeSimulation.h
 *  \brief A dye on a structure, simulated -- without naming IMP.
 *
 * Placing a dye on a residue and integrating its motion is built on IMP: its
 * hierarchies, its decorators, its integrators (`ProbeAttachment.h`,
 * `ProbeDynamics.h`). That is a fine way to *implement* it and a poor way to
 * *offer* it. A caller who installed one package should not need IMP's
 * Python to run a dye, and should not have to write a PDB file to say where
 * the atoms are.
 *
 * So #IMP::bff::MolecularProbeSimulation is the same computation behind
 * #IMP::bff::ProbeSimulation -- coordinates in arrays, `minimize()`,
 * `step()`, `run()` -- and it names no IMP type anywhere in this header.
 * That is what lets it be wrapped without IMP's own SWIG interfaces and
 * shipped in a package that carries IMP as a private library (PRD-139).
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_MOLECULAR_PROBE_SIMULATION_H
#define IMPBFF_MOLECULAR_PROBE_SIMULATION_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/Base.h>
#include <IMP/bff/HierarchyFrame.h>
#include <IMP/bff/ProbeSimulation.h>

#include <memory>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! A dye placed on a structure and moved, as an ordinary simulation object.
/*!
    Construction places the dye and builds the system; after that it is
    #IMP::bff::ProbeSimulation like any other -- `minimize()`, `step()`,
    `run()`, positions in and out.

    The structure and the dye arrive as #IMP::bff::ProteinFrame, which is
    arrays: coordinates, atom names, residue names and indices, chain ids.
    Where they came from is the caller's business -- a PDB read by
    #IMP::bff::load_structure, a frame of a trajectory, a generator. The
    dye's *chemistry* -- its bonds, angles and torsions -- comes from its
    MOL2, because that is a chemistry file and there is nothing to gain by
    re-deriving it from coordinates.

    Everything IMP is kept behind the pointer: this header names none of it,
    and neither does anything a wrapper sees.
*/
class IMPBFFEXPORT MolecularProbeSimulation : public ProbeSimulation {
    struct Impl;
    std::shared_ptr<Impl> impl_;

 public:
    //! Place \p dye on \p residue of \p chain in \p protein and set up.
    /*!
        \param[in] protein the structure to label
        \param[in] dye the dye's atoms, as they sit before placement
        \param[in] dye_mol2 the dye's MOL2: its bonds, angles and torsions
        \param[in] chain,residue the labelled site
        \param[in] parameters a JSON object; empty means the defaults.
                   `integrator` (`md` or `bd`), `temperature` (K),
                   `timestep_fs` (negative lets the integrator choose),
                   `friction_ps` (`md` only), `interaction_sphere` (A, how
                   far around the site the protein pushes back),
                   `repulsion_k`, `seed` (negative leaves the generator
                   alone), `strip_site_sidechain`.
        \throw ValueException when the site is not in the structure, or for
               an integrator that is neither `md` nor `bd`
        \throw IOException when the MOL2 cannot be read
    */
    MolecularProbeSimulation(const ProteinFrame& protein, const ProteinFrame& dye,
                  const std::string& dye_mol2, const std::string& chain,
                  int residue, const std::string& parameters = "");

    //! The dye's atoms -- the ones that move.
    int get_n_atoms() const override;
    void get_positions(double** out_view, int* n_out_view) const override;
    void set_positions(const std::vector<double>& xyz) override;
    double minimize(int n_steps = 200) override;
    void step(int n_steps) override;
    ProbeSimulationTrajectory run(int n_steps, int write_every = 1) override;
    double get_potential_energy() const override;
    bool has_energy() const override { return true; }
    std::string get_type() const override { return "dye-langevin"; }
    std::string get_parameters() const override;
    //! Only before the system is built: the integrator and its constants are
    //! fixed once the dye is placed, so this refuses afterwards.
    /*! \throw ValueException after construction */
    void set_parameters(const std::string& json) override;

    //! The dye's atom names, in the order the positions are in.
    std::vector<std::string> get_atom_names() const;
    //! How many atoms the labelled residue's side chain lost.
    int get_n_stripped() const;
    //! The labelled structure: the protein's atoms, then the dye's.
    ProteinFrame get_labelled_frame() const;
    //! The temperature a kinetic energy implies, \f$2E/(3Nk_B)\f$.
    /*! Takes the energy rather than reading one, so that a caller can ask it
        of any frame of a trajectory, not only of the state as it stands. */
    double kinetic_temperature(double kinetic_energy) const;
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_MOLECULAR_PROBE_SIMULATION_H
