/*
 * A dye on a structure, as an ordinary simulation object: coordinates in
 * arrays, `minimize()`, `step()`, `run()`, positions out. `DyeSimulation`
 * is the connection layer's -- it exists only where IMP is linked -- but
 * nothing in its declaration names an IMP type.
 *
 * That is what lets it be wrapped without IMP's own SWIG interfaces: a
 * module wrapping `IMP::atom::Hierarchy` has to `%import` IMP's kernel
 * interface, and the extension then calls
 * `PyImport_ImportModule("_IMP_kernel")` at load
 * (IMP_kernel.exceptions.i). Behind arrays, none of that follows, so a
 * package that carries IMP as a private C++ library can offer the dye
 * dynamics with no IMP in Python at all (PRD-139).
 */

%feature("kwargs") IMP::bff::DyeSimulation::DyeSimulation;
%feature("kwargs") IMP::bff::DyeSimulation::run;

%include "IMP/bff/DyeDynamics.h"

%attribute_np2(IMP::bff::DyeSimulation, std::vector<double>, positions,
               get_positions, 3);
%attribute(IMP::bff::DyeSimulation, int, n_atoms, get_n_atoms);
%attribute(IMP::bff::DyeSimulation, int, n_stripped, get_n_stripped);
%attribute_py(IMP::bff::DyeSimulation, std::vector<std::string>, atom_names,
              get_atom_names);
