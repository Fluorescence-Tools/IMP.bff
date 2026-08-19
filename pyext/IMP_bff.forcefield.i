/*
 * The coarse-grained system, typed.
 *
 * It was a nested dictionary with fifteen keys, read by six modules, holding
 * lists whose elements were positional -- a bond was
 * `[site_a, site_b, length, type_id]` and the third slot was a length by
 * convention only. Every consumer wrote `system.get("bonds", [])`, because
 * nothing guaranteed the key was there.
 *
 * The templates below are what give Python `system.sites[0].mass` and
 * `system.bond_types["b1"]` instead. They are verbose and they are the whole
 * point: each one is a container whose element type is now stated.
 */

IMP_SWIG_VALUE(IMP::bff, FFComponent, FFComponents);
IMP_SWIG_VALUE(IMP::bff, FFSite, FFSites);
IMP_SWIG_VALUE(IMP::bff, FFBond, FFBonds);
IMP_SWIG_VALUE(IMP::bff, FFAngle, FFAngles);
IMP_SWIG_VALUE(IMP::bff, FFTorsion, FFTorsions);
IMP_SWIG_VALUE(IMP::bff, FFTorsionType, FFTorsionTypes);
IMP_SWIG_VALUE(IMP::bff, FFLJType, FFLJTypes);
IMP_SWIG_VALUE(IMP::bff, FFProbe, FFProbes);
IMP_SWIG_VALUE(IMP::bff, FFNonbonded, FFNonbondeds);
IMP_SWIG_VALUE(IMP::bff, FFSampling, FFSamplings);
IMP_SWIG_VALUE(IMP::bff, DyeForceFieldSystem, DyeForceFieldSystems);

%include "IMP/bff/DyeForceField.h"

// The containers, instantiated *after* the header so their element types are
// known. Each one is what turns a positional list into a typed sequence:
// `system.bonds[0].length` rather than `system["bonds"][0][2]`.
%template(FFComponentMap) std::map<std::string, IMP::bff::FFComponent>;
%template(FFTorsionTypeMap) std::map<std::string, IMP::bff::FFTorsionType>;
%template(FFLJTypeMap) std::map<std::string, IMP::bff::FFLJType>;
%template(MapStringDouble) std::map<std::string, double>;
%template(MapStringVectorString) std::map<std::string, std::vector<std::string> >;
%template(FFSiteVector) std::vector<IMP::bff::FFSite>;
%template(FFBondVector) std::vector<IMP::bff::FFBond>;
%template(FFAngleVector) std::vector<IMP::bff::FFAngle>;
%template(FFTorsionVector) std::vector<IMP::bff::FFTorsion>;
%template(FFProbeVector) std::vector<IMP::bff::FFProbe>;

// The molecular graph's return types.
%template(PairStringString) std::pair<std::string, std::string>;
%template(VectorPairStringString) std::vector<std::pair<std::string, std::string> >;
%template(VectorVectorString) std::vector<std::vector<std::string> >;

// One Python-side method, because one is needed. `get_bonded_neighbors` and
// `find_rings` return containers a caller can use directly -- the map indexes
// and the list iterates -- so they get no wrapper. `exclusions` does: its
// caller tests membership inside an O(n^2) loop over site pairs, and a SWIG
// vector would make each test a linear scan.
%extend IMP::bff::DyeForceFieldSystem {
  %pythoncode %{
    def exclusions(self, include_impropers=True):
        """The 1-2, 1-3 and 1-4 pairs as a set of sorted ``(a, b)`` tuples.

        ``include_impropers`` adds each improper's a-c, a-d and b-d pairs. It
        is the one thing the Python copies disagreed on, and only in principle:
        no builder of a combined system fills ``impropers``. See
        okf/validation/impropers_are_dropped.md.
        """
        return {(a, b) for a, b in self.get_exclusions(include_impropers)}

  %}
}


// Read-only properties, so `system.sites` reads like the dictionary it
// replaces without `.get()` and without a default that can never be needed.
%attribute_py(IMP::bff::DyeForceFieldSystem, std::string, name, get_name, set_name);
%attribute_py(IMP::bff::DyeForceFieldSystem, PyObject*, components, get_components, set_components);
%attribute_py(IMP::bff::DyeForceFieldSystem, PyObject*, sites, get_sites, set_sites);
%attribute_py(IMP::bff::DyeForceFieldSystem, PyObject*, groups, get_groups, set_groups);
%attribute_py(IMP::bff::DyeForceFieldSystem, PyObject*, rb_groups, get_rb_groups, set_rb_groups);
%attribute_py(IMP::bff::DyeForceFieldSystem, PyObject*, md_fixed_groups, get_md_fixed_groups, set_md_fixed_groups);
%attribute_py(IMP::bff::DyeForceFieldSystem, PyObject*, fixed_groups, get_fixed_groups, set_fixed_groups);
%attribute_py(IMP::bff::DyeForceFieldSystem, PyObject*, bond_types, get_bond_types, set_bond_types);
%attribute_py(IMP::bff::DyeForceFieldSystem, PyObject*, angle_types, get_angle_types, set_angle_types);
%attribute_py(IMP::bff::DyeForceFieldSystem, PyObject*, torsion_types, get_torsion_types, set_torsion_types);
%attribute_py(IMP::bff::DyeForceFieldSystem, PyObject*, improper_types, get_improper_types, set_improper_types);
%attribute_py(IMP::bff::DyeForceFieldSystem, PyObject*, lj_types, get_lj_types, set_lj_types);
%attribute_py(IMP::bff::DyeForceFieldSystem, PyObject*, bonds, get_bonds, set_bonds);
%attribute_py(IMP::bff::DyeForceFieldSystem, PyObject*, angles, get_angles, set_angles);
%attribute_py(IMP::bff::DyeForceFieldSystem, PyObject*, dihedrals, get_dihedrals, set_dihedrals);
%attribute_py(IMP::bff::DyeForceFieldSystem, PyObject*, impropers, get_impropers, set_impropers);
%attribute_py(IMP::bff::DyeForceFieldSystem, PyObject*, probes, get_probes, set_probes);
%attribute_py(IMP::bff::DyeForceFieldSystem, PyObject*, nonbonded, get_nonbonded, set_nonbonded);
%attribute_py(IMP::bff::DyeForceFieldSystem, PyObject*, sampling, get_sampling, set_sampling);
