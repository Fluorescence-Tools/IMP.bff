/*
 * The project document: what binds one FRET docking or screening run together.
 *
 * All C++. The file layout of both doors -- the legacy `OptionsManager.Export`
 * text twin and the extended fps.json -- and every default value live in
 * `FPSProject.h`/`.cpp`, because a default that exists in two languages is a
 * default that drifts, and FPS's shipped per-mode parameters are the settings
 * behind published numbers.
 *
 * `bin/imp_bff_fps project` picks files and prints tables; it decides nothing
 * about what a project *is*.
 *
 * Must be wrapped **after** `IMP_bff.docking.i`: `FPSProject` returns a
 * `DockingParameters` by value, and SWIG needs that class already declared to
 * wrap the return rather than an opaque pointer.
 *
 * The `std::vector` member trap of `IMP_bff.docking.i` applies here too --
 * `FPSProject::structures`, `selected_distances` and `selected_flags` are
 * vector members, so bind the project to a name before reading one off it.
 */

IMP_SWIG_VALUE(IMP::bff, FPSModeParameters, FPSModeParametersList);
IMP_SWIG_VALUE(IMP::bff, FPSConversionParameters, FPSConversionParametersList);
IMP_SWIG_VALUE(IMP::bff, FPSAVGlobalParameters, FPSAVGlobalParametersList);
IMP_SWIG_VALUE(IMP::bff, FPSProject, FPSProjects);

%feature("kwargs") IMP::bff::fps_mode_parameters;
%feature("kwargs") IMP::bff::fps_project_from_json;
%feature("kwargs") IMP::bff::read_fps_project_txt;
%feature("kwargs") IMP::bff::read_fps_project;
%feature("kwargs") IMP::bff::write_fps_project;
%feature("kwargs") IMP::bff::FPSProject::get_docking_parameters;

%include "IMP/bff/FPSProject.h"

/* No `%template` for `std::vector<FPSModeParameters>`: nothing returns one.
   The five blocks are named members, and `get_parameters(mode)` returns one by
   value -- an FPS project has exactly five modes, not a list of them. */
