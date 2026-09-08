/*
 * The dynamics of an explicit probe -- IMP::atom::Simulator over IMP
 * particles, so connection layer. `make_langevin_simulator` builds the
 * integrator (`md` gets a Langevin thermostat, `bd` gets Einstein's
 * coefficient for each particle's own radius), `AttachedProbeDynamics` sets up
 * the dye's force field, its held anchor and the wall of protein spheres
 * around the site. `LangevinTrajectory` and the file doors that return one
 * are `IMP_bff.dyedynamics.i`, which is wrapped before this.
 */

%feature("kwargs") IMP::bff::make_langevin_simulator;
%feature("kwargs") IMP::bff::AttachedProbeDynamics::AttachedProbeDynamics;
%feature("kwargs") IMP::bff::AttachedProbeDynamics::run;
%feature("kwargs") IMP::bff::prepare_particles;

%include "IMP/bff/ProbeDynamics.h"

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
