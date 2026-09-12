/*
 * FRET-restrained rigid-body docking.
 *
 *  - The **values** a run is told and reports (`DockingParameters`,
 *    `PairDistance`, `DockingResult`) and the **assembly** it starts from are
 *    `Docking.h`. `create_docking_assembly` resamples each volume, attaches it
 *    to the rigid body of the atom it hangs off, gives it a radius and a mass
 *    and adds the mean-distance and excluded-volume terms -- in plain IMP, so
 *    the scoring path pulls no sampler dependency.
 *  - `score_structures`, `dock_minimize`, `refine_docking` and
 *    `screen_structures` are an API, not a program: an application drives
 *    them with a cancellation callback (#IMP::bff::DockingStop, a director
 *    class a Python caller subclasses) and reads the result back.
 *  - The **Monte-Carlo sampler** and the **repeated-trial driver** are
 *    `imp_bff dock` and `imp_bff dock-errors`. Both write directories of
 *    results, the first drives IMP.pmi (Python-only, and not a dependency of
 *    this module) and the second forks workers. Those are programs.
 *  - The derivative-enabled mean-distance restraint the minimiser uses is
 *    `ProbeAccessibleVolumeMeanDistanceRestraint`, which has the gradient in C++. There is one.
 */

/* A trap worth knowing about, measured 2026-08-31 and **not** fixed here.
 *
 * `pairs` is a `std::vector` *member*, and SWIG returns a member of class type
 * as a pointer into the owning object -- `& ((arg1)->pairs)`. Reading one off a
 * **temporary** therefore reads freed memory, and a freed `std::vector` reports
 * size 0 rather than crashing:
 *
 *     IMP.bff.score_structures(...).pairs[0]   -> IndexError: index out of range
 *     r = IMP.bff.score_structures(...); r.pairs[0]   -> the pair
 *
 * which reads as "this structure has no distances" -- the most misleading shape
 * a bug can take in a diagnostics table. **Bind the result first.**
 *
 * `%naturalvar` is the documented cure and does not work here: IMP already
 * passes `-naturalvar` globally (`tools/build/make_swig_wrapper.py:38`) and
 * `IMP_bff.types.i` says `%naturalvar;` again, and the generated getter is
 * still the pointer form. The reason is ordering -- SWIG needs
 * `std::vector<PairDistance>` to be a *known, copyable* class when it parses
 * `DockingResult::pairs`, and the `%template` below necessarily comes after
 * `%include "IMP/bff/Docking.h"`, which is where `PairDistance` is declared in
 * the first place. Fixing it means splitting the header or forward-declaring
 * into a `%template`, neither of which belongs in a change about FPS scoring.
 */

IMP_SWIG_VALUE(IMP::bff, PairDistance, PairDistances);
IMP_SWIG_VALUE(IMP::bff, ReferenceAtom, ReferenceAtoms);
IMP_SWIG_VALUE(IMP::bff, ReferenceFit, ReferenceFits);
IMP_SWIG_VALUE(IMP::bff, DockingParameters, DockingParametersList);
IMP_SWIG_VALUE(IMP::bff, DockingResult, DockingResults);
IMP_SWIG_VALUE(IMP::bff, ScreenedStructure, ScreenedStructures);
IMP_SWIG_VALUE(IMP::bff, BootstrapParameters, BootstrapParametersList);
IMP_SWIG_VALUE(IMP::bff, BootstrapReplica, BootstrapReplicas);
IMP_SWIG_VALUE(IMP::bff, BootstrapResult, BootstrapResults);
IMP_SWIG_OBJECT(IMP::bff, DockingStop, DockingStops);
IMP_SWIG_OBJECT(IMP::bff, ScoreTrace, ScoreTraces);

// A Python caller subclasses `DockingStop` to cancel a long run; the director
// is what lets C++ call back into that subclass.
%feature("director") IMP::bff::DockingStop;

%feature("kwargs") IMP::bff::DockingResult::DockingResult;
%feature("kwargs") IMP::bff::create_docking_assembly;
%feature("kwargs") IMP::bff::score_assembly;
%feature("kwargs") IMP::bff::score_structures;
%feature("kwargs") IMP::bff::dock_minimize;
%feature("kwargs") IMP::bff::refine_docking;
%feature("kwargs") IMP::bff::screen_structures;
%feature("kwargs") IMP::bff::collect_pair_distances;
%feature("kwargs") IMP::bff::pair_distances_at_positions;
%feature("kwargs") IMP::bff::fit_reference_atoms;
%feature("kwargs") IMP::bff::fit_reference_positions;
%feature("kwargs") IMP::bff::fps_bootstrap;
%feature("kwargs") IMP::bff::sample_distance_perturbations;
%feature("kwargs") IMP::bff::pose_rmsd;
%feature("kwargs") IMP::bff::pose_superposition;

%include "IMP/bff/Docking.h"

%template(PairDistanceList) std::vector<IMP::bff::PairDistance>;
%template(ScreenedStructureList) std::vector<IMP::bff::ScreenedStructure>;
// `std::vector<double>` is **not** templated here: `IMP.saxs` already wraps it
// as `DistBase` and SWIG wraps a type once across a module and its imports
// (`IMP_bff.types.i`). `sample_distance_perturbations` therefore hands back an
// `IMP.saxs.DistBase`, which indexes and iterates like a list.
%template(BootstrapReplicaList) std::vector<IMP::bff::BootstrapReplica>;

// The derived numbers as attributes: a residual is not a field a caller may
// set, it is what the two distances say.
%attribute(IMP::bff::PairDistance, double, residual, get_residual);
%attribute(IMP::bff::PairDistance, double, chi2, get_chi2);
%attribute(IMP::bff::PairDistance, double, efficiency_model,
           get_efficiency_model);
%attribute(IMP::bff::PairDistance, double, efficiency_exp, get_efficiency_exp);
%attribute_py(IMP::bff::ReferenceAtom, Vector3D, coordinates, get_coordinates);

// The assembly's parts, as attributes: a caller reads them, a run sets them.
%attribute_py(IMP::bff::DockingAssembly, Model, model, get_model);
%attribute_py(IMP::bff::DockingAssembly, Hierarchy, root, get_root);
%attribute_py(IMP::bff::DockingAssembly, ProbeNetworkRestraint, network,
              get_network);
%attribute_py(IMP::bff::DockingAssembly, RestraintSet, restraints,
              get_restraints);
%attribute_py(IMP::bff::DockingAssembly, RestraintsScoringFunction,
              scoring_function, get_scoring_function);
%attribute_py(IMP::bff::DockingAssembly, RigidBodies, rigid_bodies,
              get_rigid_bodies);
%attribute_py(IMP::bff::DockingAssembly, VectorInt, body_of_pdb,
              get_body_of_pdb);
%attribute(IMP::bff::DockingAssembly, bool, mean_position, get_mean_position);
%attribute(IMP::bff::DockingAssembly, double, sigma_da, get_sigma_da);
