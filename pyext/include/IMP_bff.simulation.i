/*
 * The shape every probe simulation shares: `ProbeSimulation` and the record
 * a run produces, `SimulationTrajectory`.
 *
 * Wrapped before anything that implements or returns one -- the grid walk in
 * `IMP_bff.sampling.i`, the dye in `IMP_bff.dyedynamics.i`. Parameters cross
 * as JSON text, so a simulation that gains a knob gains nothing here.
 */

IMP_SWIG_VALUE(IMP::bff, SimulationTrajectory, SimulationTrajectories);

%include "IMP/bff/Simulation.h"

// The trajectory's shapes: `(n_frames, n_atoms, 3)` for the coordinates and
// one value per frame for the rest.
%attribute_np3(IMP::bff::SimulationTrajectory, std::vector<double>, coordinates,
               get_coordinates, n_frames, 3);
%attribute_np(IMP::bff::SimulationTrajectory, std::vector<double>, times_fs,
              get_times_fs);
%attribute_np(IMP::bff::SimulationTrajectory, std::vector<double>,
              potential_energy, get_potential_energy);
%attribute_np(IMP::bff::SimulationTrajectory, std::vector<double>,
              kinetic_energy, get_kinetic_energy);
%attribute(IMP::bff::SimulationTrajectory, int, n_frames, n_frames);
%attribute(IMP::bff::SimulationTrajectory, int, n_atoms, n_atoms);
%attribute_py(IMP::bff::SimulationTrajectory, std::vector<std::string>,
              atom_names, atom_names);

%pythoncode %{
# LangevinTrajectory was this record's name while the dye run was the only
# thing that produced one. It is the same type.
LangevinTrajectory = SimulationTrajectory
LangevinTrajectories = SimulationTrajectories
%}
