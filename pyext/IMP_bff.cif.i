/*
 * io/cif.py to C++: the FF writer, template CIF, rotamer library IO, and
 * string utilities. The FF reader is already in ForceFieldCIF.h.
 *
 * There is no Python left here. What was:
 *
 *  - `read_dye_forcefield_cif(path)` was `read_forcefield_cif(str(path))` --
 *    a second name for the C++ reader, spelled in Python so that it could
 *    take a `pathlib.Path`. There is one name for the reader now, and a
 *    caller holding a `Path` says `str()` where it says it.
 *  - `write_dye_forcefield_cif(path, system)` shadowed the C++ writer of the
 *    same name, to `str()` its path and coerce a dict system -- and so had to
 *    call `_IMP_bff.` to avoid recursing into itself. A function whose only
 *    body is a call to the thing it shadows is a hazard, not a convenience.
 *  - `forcefield_system_from_dict` and `as_forcefield_system` turned a dict
 *    into a typed system through `json.dumps`. That is what a *caller* does
 *    with a dict it built; the library takes systems. `bin/imp_bff` and the
 *    three tests that build systems as dict literals coerce their own.
 */

IMP_SWIG_VALUE(IMP::bff, ComponentTemplate, ComponentTemplates);

// The nested Feature/FeatureAtom/Improper structs, flattened to module level
// so the `features` map and its rows marshal instead of sitting opaque.
%feature("flatnested") IMP::bff::ComponentTemplate::FeatureAtom;
%feature("flatnested") IMP::bff::ComponentTemplate::Feature;
%feature("flatnested") IMP::bff::ComponentTemplate::Improper;

%include "IMP/bff/CifIO.h"

%template(FeatureAtomList) std::vector<IMP::bff::ComponentTemplate::FeatureAtom>;
%template(ImproperList) std::vector<IMP::bff::ComponentTemplate::Improper>;
%template(FeatureMap) std::map<std::string, IMP::bff::ComponentTemplate::Feature>;
