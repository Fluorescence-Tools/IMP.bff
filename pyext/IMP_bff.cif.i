/*
 * io/cif.py to C++: the FF writer, template CIF, rotamer library IO, and
 * string utilities. The FF reader is already in ForceFieldCIF.h.
 *
 * The C++ owns all of it now -- CifIO.h carries the `ComponentTemplate` and
 * `RotamerLibraryData` values the dict-shaped Python used to return, and the
 * writer/reader pairs under their own names. The one dict-literal convenience
 * that remains is `forcefield_system_from_json` (a nlohmann bridge in
 * DyeForceField.cpp); `as_forcefield_system`, which decides *whether* to call
 * it, lives in the `IMP.bff.io.cif` package module -- it is duck typing on a
 * Python object, which has no C++ spelling.
 */

IMP_SWIG_VALUE(IMP::bff, ComponentTemplate, ComponentTemplates);
IMP_SWIG_VALUE(IMP::bff, RotamerLibraryData, RotamerLibraryDatas);

// The nested Feature/FeatureAtom/Improper structs, flattened to module level
// so the `features` map and its rows marshal instead of sitting opaque.
%feature("flatnested") IMP::bff::ComponentTemplate::FeatureAtom;
%feature("flatnested") IMP::bff::ComponentTemplate::Feature;
%feature("flatnested") IMP::bff::ComponentTemplate::Improper;

%include "IMP/bff/CifIO.h"

%template(FeatureAtomList) std::vector<IMP::bff::ComponentTemplate::FeatureAtom>;
%template(ImproperList) std::vector<IMP::bff::ComponentTemplate::Improper>;
%template(FeatureMap) std::map<std::string, IMP::bff::ComponentTemplate::Feature>;