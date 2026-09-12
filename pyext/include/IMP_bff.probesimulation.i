/*
 * The shape every probe simulation shares: `ProbeSimulation` and the record
 * a run produces, `ProbeSimulationTrajectory`.
 *
 * Wrapped before anything that implements or returns one -- the grid walk in
 * `IMP_bff.sampling.i`, the dye in `IMP_bff.molecularprobesimulation.i`. Parameters cross
 * as JSON text, so a simulation that gains a knob gains nothing here.
 */

IMP_SWIG_VALUE(IMP::bff, ProbeSimulationTrajectory, ProbeSimulationTrajectories);

%include "IMP/bff/ProbeSimulation.h"

// The trajectory's shapes: `(n_frames, n_atoms, 3)` for the coordinates and
// one value per frame for the rest.
%attribute_np3(IMP::bff::ProbeSimulationTrajectory, std::vector<double>, coordinates,
               get_coordinates, n_frames, 3);
%attribute_np(IMP::bff::ProbeSimulationTrajectory, std::vector<double>, times_fs,
              get_times_fs);
%attribute_np(IMP::bff::ProbeSimulationTrajectory, std::vector<double>,
              potential_energy, get_potential_energy);
%attribute_np(IMP::bff::ProbeSimulationTrajectory, std::vector<double>,
              kinetic_energy, get_kinetic_energy);
%attribute(IMP::bff::ProbeSimulationTrajectory, int, n_frames, n_frames);
%attribute(IMP::bff::ProbeSimulationTrajectory, int, n_atoms, n_atoms);
%attribute_py(IMP::bff::ProbeSimulationTrajectory, std::vector<std::string>,
              atom_names, atom_names);

