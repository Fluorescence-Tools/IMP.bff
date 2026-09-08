/*
 * The dye roads by file path: `attach_dye_to_pdb` and `run_dye_langevin`,
 * with `LangevinTrajectory` -- what a run returns.
 *
 * These are the connection layer's, and they exist only where IMP is linked,
 * but nothing in their signatures names an IMP type. That is what lets them
 * be wrapped without IMP's own SWIG interfaces: a module wrapping
 * `IMP::atom::Hierarchy` has to `%import` IMP's kernel interface, and the
 * extension then calls `PyImport_ImportModule("_IMP_kernel")` at load
 * (IMP_kernel.exceptions.i). Behind paths and arrays, none of that follows,
 * so a wheel that carries IMP as a private C++ library can offer the dye
 * dynamics with no IMP in Python at all (PRD-139).
 *
 * Wrapped ahead of `IMP_bff.probedynamics.i`, which needs
 * `LangevinTrajectory` declared before `ProbeDynamics.h` uses it.
 */

IMP_SWIG_VALUE(IMP::bff, LangevinTrajectory, LangevinTrajectories);

%feature("kwargs") IMP::bff::attach_dye_to_pdb;
%feature("kwargs") IMP::bff::run_dye_langevin;

%include "IMP/bff/DyeDynamics.h"

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
%attribute_py(IMP::bff::LangevinTrajectory, std::vector<std::string>, atom_names,
              atom_names);
