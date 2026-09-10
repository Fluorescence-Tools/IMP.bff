/*
 * Sampling an explicit probe. What each piece is:
 *
 *  - **The dynamics** (`ProbeDynamics.h`, over IMP particles) are the
 *    connection layer's, in probedynamics.i.
 *  - **The linker sampler** is `Linker.h`, beside the torsion geometry it
 *    applies. It has no IMP model at
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

// A Python caller subclasses `RRTCollision` to say what clashes; the director
// is what lets the C++ growth loop ask it.
IMP_SWIG_VALUE(IMP::bff, LinkerSamplingResult, LinkerSamplingResults);
IMP_SWIG_VALUE(IMP::bff, RRTTree, RRTTrees);
IMP_SWIG_OBJECT(IMP::bff, RRTCollision, RRTCollisions);
%feature("kwargs") IMP::bff::sample_linker;
%feature("kwargs") IMP::bff::generate_linker_rotamers;
%feature("kwargs") IMP::bff::linker_geometry_from_mol2;
%feature("kwargs") IMP::bff::grow_torsion_rrt;
%feature("kwargs") IMP::bff::grow_rigid_body_rrt;
%feature("kwargs") IMP::bff::markov_state_trajectory;

%feature("director") IMP::bff::RRTCollision;

/* Linker.h needs RotamerLibrary (wrapped above, in core.i): the sampler
   returns one. Its geometry half used to be wrapped earlier, on its own. */
%include "IMP/bff/Linker.h"
%include "IMP/bff/RRT.h"

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
