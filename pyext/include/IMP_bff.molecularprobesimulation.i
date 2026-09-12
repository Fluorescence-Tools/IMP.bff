/*
 * A dye on a structure, as an ordinary simulation object: coordinates in
 * arrays, `minimize()`, `step()`, `run()`, positions out. `MolecularProbeSimulation`
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

%feature("kwargs") IMP::bff::MolecularProbeSimulation::MolecularProbeSimulation;
%feature("kwargs") IMP::bff::MolecularProbeSimulation::run;

%include "IMP/bff/MolecularProbeSimulation.h"

%attribute_np2(IMP::bff::MolecularProbeSimulation, std::vector<double>, positions,
               get_positions, 3);
%attribute(IMP::bff::MolecularProbeSimulation, int, n_atoms, get_n_atoms);
%attribute(IMP::bff::MolecularProbeSimulation, int, n_stripped, get_n_stripped);
%attribute_py(IMP::bff::MolecularProbeSimulation, std::vector<std::string>, atom_names,
              get_atom_names);
