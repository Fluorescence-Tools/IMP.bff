/*
 * FRET-restrained rigid-body docking.
 *
 * This file was 1,427 lines of `%pythoncode` -- an engine, as the maintainer
 * put it -- and holds no Python at all now. Where it went, and why:
 *
 *  - The **values** a run is told and reports (`DockingParameters`,
 *    `PairDistance`, `DockingResult`) and the **assembly** it starts from are
 *    `Docking.h`. `build_docking_assembly` does what the IMP.pmi wrapper did
 *    -- resample each volume, attach it to the rigid body of the atom it
 *    hangs off, give it a radius and a mass, add the mean-distance and
 *    excluded-volume terms -- in plain IMP, so the scoring path no longer
 *    pulls a sampler's dependency.
 *  - `score`, `dock_minimize`, `refine` and `screen` are C++ under
 *    self-describing names (`score_structures`, `dock_minimize`,
 *    `refine_docking`, `screen_structures`): an application drives them with
 *    a cancellation callback and reads the result back, so they are an API
 *    and not a program. The callback is #IMP::bff::DockingStop, a director
 *    class a Python caller subclasses.
 *  - The **Monte-Carlo sampler** and the **repeated-trial driver** are
 *    `imp_bff dock` and `imp_bff dock-errors`. Both write directories of
 *    results, the first drives IMP.pmi (Python-only, and not a dependency of
 *    this module) and the second forks workers. Those are programs.
 *  - A *third* copy of the derivative-enabled mean-distance restraint was
 *    defined here; `AVMeanDistanceRestraint` has had the gradient in C++
 *    since the AV batch and is what the minimiser uses.
 *
 * Every entry point in this file raised `NameError` on its first line before
 * this move -- four separate dead names, one per path -- because nothing
 * called it and nothing tested it. `test/test_docking_values.py` is the
 * coverage that would have caught them.
 */

IMP_SWIG_VALUE(IMP::bff, PairDistance, PairDistances);
IMP_SWIG_VALUE(IMP::bff, DockingParameters, DockingParametersList);
IMP_SWIG_VALUE(IMP::bff, DockingResult, DockingResults);
IMP_SWIG_VALUE(IMP::bff, ScreenedStructure, ScreenedStructures);
IMP_SWIG_OBJECT(IMP::bff, DockingStop, DockingStops);
IMP_SWIG_OBJECT(IMP::bff, ScoreTrace, ScoreTraces);

// A Python caller subclasses `DockingStop` to cancel a long run; the director
// is what lets C++ call back into that subclass.
%feature("director") IMP::bff::DockingStop;

%feature("kwargs") IMP::bff::DockingResult::DockingResult;
%feature("kwargs") IMP::bff::build_docking_assembly;
%feature("kwargs") IMP::bff::score_assembly;
%feature("kwargs") IMP::bff::score_structures;
%feature("kwargs") IMP::bff::dock_minimize;
%feature("kwargs") IMP::bff::refine_docking;
%feature("kwargs") IMP::bff::screen_structures;
%feature("kwargs") IMP::bff::collect_pair_distances;
%feature("kwargs") IMP::bff::pair_distances_at_positions;

%include "IMP/bff/Docking.h"

%template(PairDistanceList) std::vector<IMP::bff::PairDistance>;
%template(ScreenedStructureList) std::vector<IMP::bff::ScreenedStructure>;

// The derived numbers as attributes: a residual is not a field a caller may
// set, it is what the two distances say.
%attribute(IMP::bff::PairDistance, double, residual, get_residual);
%attribute(IMP::bff::PairDistance, double, chi2, get_chi2);
%attribute(IMP::bff::PairDistance, double, efficiency_model,
           get_efficiency_model);
%attribute(IMP::bff::PairDistance, double, efficiency_exp, get_efficiency_exp);

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
