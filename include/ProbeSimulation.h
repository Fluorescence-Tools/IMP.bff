/**
 *  \file IMP/bff/ProbeSimulation.h
 *  \brief What every probe simulation in this module can do.
 *
 * There is more than one way to move a label: a dye on its linker integrated
 * with Langevin dynamics on an atomistic model, and a probe diffusing on the
 * accessible volume's grid. They are different physics on different state,
 * and until now they were also different *shapes* -- one an object with a
 * constructor full of parameters and a `run` that returns frames, the other
 * an object with a constructor full of grids and a `run` that returns a
 * count. A caller who wanted to swap one for the other rewrote their code.
 *
 * #IMP::bff::ProbeSimulation is the shape they share, and it is deliberately
 * OpenMM's: positions in, `minimize()`, `step()`, `run()`, a state to read
 * back out. What a caller writes against it works for either.
 *
 * Coordinates are always `(n_atoms, 3)` in Angstrom and always plain arrays;
 * no structure file appears in this interface. Times are femtoseconds,
 * because that is the unit an integrator's step is naturally written in and
 * one unit is better than two.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_PROBE_SIMULATION_H
#define IMPBFF_PROBE_SIMULATION_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/IMPCompatibility.h>

#include <limits>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! One recorded run: the frames, and what was known at each of them.
/*!
    The same record for every simulation. `potential_energy` and
    `kinetic_energy` are empty where the simulation has no such quantity --
    a grid walk has neither -- rather than filled with a polite zero, so that
    a caller can tell "not measured" from "measured and small".
*/
struct IMPBFFEXPORT ProbeSimulationTrajectory {
    //! Flat `n_frames * n_atoms * 3`, Angstrom.
    std::vector<double> coordinates;
    //! One per frame: the time in femtoseconds.
    std::vector<double> times_fs;
    //! One per frame, or empty where the simulation has no energy.
    std::vector<double> potential_energy, kinetic_energy;
    //! One per atom, where the atoms have names.
    std::vector<std::string> atom_names;
    //! What produced this: the simulation's type, and its integrator where
    //! it has one (`md`, `bd`, `grid-walk`).
    std::string integrator;
    double temperature, timestep_fs;
    int n_frames, n_atoms;

    ProbeSimulationTrajectory()
        : integrator("md"), temperature(300.0), timestep_fs(2.0), n_frames(0),
          n_atoms(0) {}

    void get_coordinates(double** out_view, int* n_out_view) const;
    void get_times_fs(double** out_view, int* n_out_view) const;
    void get_potential_energy(double** out_view, int* n_out_view) const;
    void get_kinetic_energy(double** out_view, int* n_out_view) const;

    IMP_SHOWABLE_INLINE(ProbeSimulationTrajectory,
                        out << "ProbeSimulationTrajectory(" << n_frames
                            << " frames of " << n_atoms << " atoms, "
                            << integrator << ")");
};
IMP_VALUES(ProbeSimulationTrajectory, ProbeSimulationTrajectories);

//! The shape every probe simulation has: set positions, advance, read back.
/*!
    Modelled on OpenMM, because that is the interface people already know:
    the object holds the system and its parameters, `minimize()` relaxes it,
    `step()` advances it, `run()` advances it while recording, and the
    positions are readable and writable at any point.

    **Parameters are JSON**, not a struct of fields, and that is the whole
    reason one interface can carry these at all. A Langevin run wants a
    temperature, a timestep and a friction; a grid walk wants a diffusion
    coefficient, a slow factor and a voxel edge; a lattice walk wants a
    stencil; a learned forward model wants weights and a cutoff. Naming all
    of those in one struct would give every simulation the fields of every
    other, and adding a kind would move everybody's arguments. A JSON object
    says only what this simulation has, `set_parameters` takes a partial one
    and leaves the rest alone, and `get_parameters` always reports what is
    actually in force -- which is also what belongs in a file beside the
    results.

    What is *not* JSON: the state and the system. Coordinates, grids and
    topologies arrive as arrays through the constructor, because they are
    large and numeric and JSON is neither.

    A simulation that has no minimiser returns NaN from `minimize()` rather
    than pretending; the same for `get_potential_energy()`. Ask
    `has_energy()` when it matters.
*/
class IMPBFFEXPORT ProbeSimulation {
 public:
    virtual ~ProbeSimulation() {}

    //! What kind this is: `dye-langevin`, `grid-diffusion`, ...
    /*! The name a caller switches on, and the one `get_parameters()` should
        be read against -- each kind documents its own keys. */
    virtual std::string get_type() const = 0;

    //! Everything in force, as a JSON object.
    virtual std::string get_parameters() const = 0;

    //! Change some of them; keys that are not named keep their values.
    /*! \throw ValueException on JSON that will not parse, on a key this
               simulation does not have, or on a value out of its range --
               a silently ignored parameter is a silently wrong run. */
    virtual void set_parameters(const std::string& json) = 0;

    //! How many atoms (or walkers) the state has.
    virtual int get_n_atoms() const = 0;

    //! The current positions, `(n_atoms, 3)` in Angstrom.
    virtual void get_positions(double** out_view, int* n_out_view) const = 0;

    //! Put the state somewhere; flat `n_atoms * 3`.
    /*! \throw ValueException on a length that is not `3 * get_n_atoms()` */
    virtual void set_positions(const std::vector<double>& xyz) = 0;

    //! Relax the state; returns the score reached, or NaN where there is none.
    virtual double minimize(int n_steps = 200) {
        (void)n_steps;
        return std::numeric_limits<double>::quiet_NaN();
    }

    //! Advance without recording.
    virtual void step(int n_steps) = 0;

    //! Advance, keeping a frame every \p write_every steps.
    virtual ProbeSimulationTrajectory run(int n_steps, int write_every = 1) = 0;

    //! The potential energy as it stands, or NaN where there is none.
    virtual double get_potential_energy() const {
        return std::numeric_limits<double>::quiet_NaN();
    }

    //! Whether this simulation has a potential energy at all.
    virtual bool has_energy() const { return false; }

    IMP_SHOWABLE_INLINE(ProbeSimulation,
                        out << "ProbeSimulation(" << get_type() << ", "
                            << get_n_atoms() << " atoms)");
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_PROBE_SIMULATION_H
