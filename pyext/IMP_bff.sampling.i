/*
 * Sampling an explicit probe. What each piece is:
 *
 *  - **The dynamics** are `ProbeDynamics.h`: `make_langevin_simulator` builds
 *    the integrator (`md` gets a Langevin thermostat, `bd` gets Einstein's
 *    coefficient for each particle's own radius), `AttachedProbeDynamics` sets up
 *    the dye's force field, its held anchor and the wall of protein spheres
 *    around the site, and `LangevinTrajectory` is what a run returns.
 *  - **The linker sampler** is `LinkerSampling.h`. It has no IMP model at
 *    all: routing each trial through `XYZ` decorators would use the
 *    particles as a scratch buffer for numbers `LinkerGeometry::apply` has
 *    just returned.
 *  - **The planners** are `RRT.h` -- one tree, two metrics. A caller's
 *    collision test crosses back through the `RRTCollision` director.
 *  - **The Markov walk** (`reconstruct_rotamer_trajectory`) is
 *    `markov_state_trajectory` in `Clustering.h`, beside the transition matrix
 *    it draws from.
 *  - **The clustering and Markov API** is `Clustering.h`, which takes the
 *    flat shapes directly.
 *
 * Two conveniences are deliberately absent. `is_collision_sphere` is eight
 * lines of arithmetic a caller writes inside its own `RRTCollision`; and the
 * sampler has no `out_rmf=`, because a trajectory comes back as coordinates
 * and writing it is `write_rmf`'s job, not the engine's.
 */

IMP_SWIG_VALUE(IMP::bff, LangevinTrajectory, LangevinTrajectories);
IMP_SWIG_VALUE(IMP::bff, LinkerSamplingResult, LinkerSamplingResults);
IMP_SWIG_VALUE(IMP::bff, RRTTree, RRTTrees);
IMP_SWIG_OBJECT(IMP::bff, RRTCollision, RRTCollisions);

// A Python caller subclasses `RRTCollision` to say what clashes; the director
// is what lets the C++ growth loop ask it.
%feature("director") IMP::bff::RRTCollision;

%feature("kwargs") IMP::bff::make_langevin_simulator;
%feature("kwargs") IMP::bff::AttachedProbeDynamics::AttachedProbeDynamics;
%feature("kwargs") IMP::bff::AttachedProbeDynamics::run;
%feature("kwargs") IMP::bff::prepare_particles;
%feature("kwargs") IMP::bff::sample_linker;
%feature("kwargs") IMP::bff::generate_linker_rotamers;
%feature("kwargs") IMP::bff::linker_geometry_from_mol2;
%feature("kwargs") IMP::bff::grow_torsion_rrt;
%feature("kwargs") IMP::bff::grow_rigid_body_rrt;
%feature("kwargs") IMP::bff::markov_state_trajectory;

%include "IMP/bff/ProbeDynamics.h"
%include "IMP/bff/LinkerSampling.h"
%include "IMP/bff/RRT.h"

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

%attribute_np3(IMP::bff::LinkerSamplingResult, std::vector<double>,
               coordinates, get_coordinates, n_frames, 3);
%attribute_np(IMP::bff::LinkerSamplingResult, std::vector<double>, energies,
              get_energies);
%attribute(IMP::bff::LinkerSamplingResult, int, n_frames, n_frames);
%attribute(IMP::bff::LinkerSamplingResult, int, n_atoms, n_atoms);
%attribute(IMP::bff::LinkerSamplingResult, int, n_accepted, n_accepted);
%attribute(IMP::bff::LinkerSamplingResult, double, acceptance, get_acceptance);

%attribute_np2(IMP::bff::RRTTree, std::vector<double>, configurations,
               get_configurations, n_dof);
%attribute(IMP::bff::RRTTree, int, n_nodes, n_nodes);
%attribute(IMP::bff::RRTTree, int, n_dof, n_dof);
%attribute(IMP::bff::RRTTree, int, goal_node, goal_node);
// `parents` is a public member and SWIG already exposes it; wrapping it in a
// `%attribute_py` would make a property of a property.
