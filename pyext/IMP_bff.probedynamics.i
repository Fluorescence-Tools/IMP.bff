/*
 * The dynamics of an explicit probe -- IMP::atom::Simulator over IMP
 * particles, so connection layer. `make_langevin_simulator` builds the
 * integrator (`md` gets a Langevin thermostat, `bd` gets Einstein's
 * coefficient for each particle's own radius), `AttachedProbeDynamics` sets up
 * the dye's force field, its held anchor and the wall of protein spheres
 * around the site, and `LangevinTrajectory` is what a run returns.
 */

IMP_SWIG_VALUE(IMP::bff, LangevinTrajectory, LangevinTrajectories);

%feature("kwargs") IMP::bff::make_langevin_simulator;
%feature("kwargs") IMP::bff::AttachedProbeDynamics::AttachedProbeDynamics;
%feature("kwargs") IMP::bff::AttachedProbeDynamics::run;
%feature("kwargs") IMP::bff::prepare_particles;

%include "IMP/bff/ProbeDynamics.h"

// The trajectory's shapes: `(n_frames, n_atoms, 3)` for the coordinates and
// one value per frame for the rest.
%attribute_np3(IMP::bff::LangevinTrajectory, std::vector<double>, coordinates,
               get_coordinates, n_frames, 3);
%attribute_np(IMP::bff::LangevinTrajectory, std::vector<double>, times_fs,
              get_times_fs);
%attribute_np(IMP::bff::LangevinTrajectory, std::vector<double>,
              potential_energy, get_potential_energy);
%attribute_np(IMP::bff::LangevinTrajectory, std::vector<double>,
              kinetic_energy, get_kinetic_energy);
%attribute(IMP::bff::LangevinTrajectory, int, n_frames, n_frames);
%attribute(IMP::bff::LangevinTrajectory, int, n_atoms, n_atoms);

%attribute_np2(IMP::bff::AttachedProbeDynamics, std::vector<double>, coordinates,
               get_coordinates, 3);
%attribute_py(IMP::bff::AttachedProbeDynamics, ParticlesTemp, probe_particles,
              get_probe_particles);
%attribute_py(IMP::bff::AttachedProbeDynamics, ParticlesTemp, mobile, get_mobile);
%attribute_py(IMP::bff::AttachedProbeDynamics, ParticlesTemp, fixed, get_fixed);
%attribute_py(IMP::bff::AttachedProbeDynamics, ParticlesTemp, obstacles,
              get_obstacles);
%attribute_py(IMP::bff::AttachedProbeDynamics, std::vector<std::string>, atom_names,
              get_atom_names);
