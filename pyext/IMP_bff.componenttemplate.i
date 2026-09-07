/*
 * The cgprobe component template: what a component's features, impropers and
 * probe metadata look like in Python.
 *
 * Nothing in this module parses CIF -- that is IMP's vendored ihm reader --
 * and the writer it does carry is one internal class, `internal/CifWriter.h`.
 */

IMP_SWIG_VALUE(IMP::bff, ComponentTemplate, ComponentTemplates);

// The nested Feature/FeatureAtom/Improper structs, flattened to module level
// so the `features` map and its rows marshal instead of sitting opaque.
%feature("flatnested") IMP::bff::ComponentTemplate::FeatureAtom;
%feature("flatnested") IMP::bff::ComponentTemplate::Feature;
%feature("flatnested") IMP::bff::ComponentTemplate::Improper;

%include "IMP/bff/ComponentTemplate.h"

%template(FeatureAtomList) std::vector<IMP::bff::ComponentTemplate::FeatureAtom>;
%template(ImproperList) std::vector<IMP::bff::ComponentTemplate::Improper>;
%template(FeatureMap) std::map<std::string, IMP::bff::ComponentTemplate::Feature>;
