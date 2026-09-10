/*
 * A probe as a species -- a dye, a fluorescent protein, or a spin label --
 * its library file, and the coarse-grained system it is simulated as. One
 * header (`ProbeLibrary.h`, PRD-138), one interface file.
 *
 * `Probe` and `Spectrum` were frozen dataclasses holding numpy arrays, a CIF
 * reader built on `ihm.format`'s Python parser with a hand-rolled mtime cache,
 * and a Förster radius computed with `np.trapezoid`. None of that is arithmetic
 * that needed C++; what needed it is the same thing `LifetimeSpectrum` needed,
 * which is that the *object* be a C++ value the way `atom`'s and `core`'s are.
 * Everything that reads or computes is C++ -- the name-keyed front door the
 * rotamer code uses, the `_cgprobe_metadata` read, and the two flrCIF item
 * maps, which live in ProbeLibrary.cpp as plain data exposed over SWIG.
 *
 * The probe container (`.mmfdb.pto`, the former `ProbeContainer.h`) speaks
 * paths and returns values; its signatures name nothing wrapped later.
 *
 * The coarse-grained system (the former `ProbeForceField.h`) is typed: the
 * templates below are what give Python `system.sites[0].mass` and
 * `system.bond_types["b1"]` -- each one is a container whose element type is
 * stated, so a bond's third slot is `length` by declaration rather than by
 * convention. It is wrapped here, early, because Mol2IO.h, Scoring.h and
 * ForceFieldCIF.h, all wrapped soon after, take the system.
 */
IMP_SWIG_VALUE(IMP::bff, Spectrum, Spectrums);
IMP_SWIG_VALUE(IMP::bff, Probe, Probes);

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
IMP_SWIG_VALUE(IMP::bff, ProbeForceFieldSystem, ProbeForceFieldSystems);

// The molecular graph's return types, declared before the header so the
// `std::set` below is already a known template when the method is read.
%template(PairStringString) std::pair<std::string, std::string>;
%template(VectorPairStringString) std::vector<std::pair<std::string, std::string> >;
%template(VectorVectorString) std::vector<std::vector<std::string> >;

// `std_set.i`+`std_pair.i` (in types.i) turn this into a native Python `set`
// of `(a, b)` tuples -- exactly the surface a caller expects from a set of
// excluded pairs -- with no hand-written typemap.
%template(SetPairStringString)
        std::set<std::pair<std::string, std::string> >;

// Spelled `exclusions` in Python: it names the pairs, not the getter. The
// rename has to come before the header declares the method.
%rename(exclusions) IMP::bff::ProbeForceFieldSystem::get_exclusions;

%include "IMP/bff/ProbeLibrary.h"

%template(ProbeMap) std::map<std::string, IMP::bff::Probe>;

// The two flrCIF item maps are data about a dictionary, now returned from C++.
// They declare themselves in ProbeLibrary.h; over SWIG they surface as Python
// dicts via the MapStringString template. A bff-native category with no
// flrCIF item simply has no key.
//   PROBE_FLRCIF_ITEMS              -> probe_flrcif_items()
//   FORSTER_RADIUS_FLRCIF_ITEMS     -> forster_radius_flrcif_items()

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

// Read-only properties, so `system.sites` reads like the dictionary it
// replaces without `.get()` and without a default that can never be needed.
%attribute_py(IMP::bff::ProbeForceFieldSystem, std::string, name, get_name, set_name);
%attribute_py(IMP::bff::ProbeForceFieldSystem, PyObject*, components, get_components, set_components);
%attribute_py(IMP::bff::ProbeForceFieldSystem, PyObject*, sites, get_sites, set_sites);
%attribute_py(IMP::bff::ProbeForceFieldSystem, PyObject*, groups, get_groups, set_groups);
%attribute_py(IMP::bff::ProbeForceFieldSystem, PyObject*, rb_groups, get_rb_groups, set_rb_groups);
%attribute_py(IMP::bff::ProbeForceFieldSystem, PyObject*, md_fixed_groups, get_md_fixed_groups, set_md_fixed_groups);
%attribute_py(IMP::bff::ProbeForceFieldSystem, PyObject*, fixed_groups, get_fixed_groups, set_fixed_groups);
%attribute_py(IMP::bff::ProbeForceFieldSystem, PyObject*, bond_types, get_bond_types, set_bond_types);
%attribute_py(IMP::bff::ProbeForceFieldSystem, PyObject*, angle_types, get_angle_types, set_angle_types);
%attribute_py(IMP::bff::ProbeForceFieldSystem, PyObject*, torsion_types, get_torsion_types, set_torsion_types);
%attribute_py(IMP::bff::ProbeForceFieldSystem, PyObject*, improper_types, get_improper_types, set_improper_types);
%attribute_py(IMP::bff::ProbeForceFieldSystem, PyObject*, lj_types, get_lj_types, set_lj_types);
%attribute_py(IMP::bff::ProbeForceFieldSystem, PyObject*, bonds, get_bonds, set_bonds);
%attribute_py(IMP::bff::ProbeForceFieldSystem, PyObject*, angles, get_angles, set_angles);
%attribute_py(IMP::bff::ProbeForceFieldSystem, PyObject*, dihedrals, get_dihedrals, set_dihedrals);
%attribute_py(IMP::bff::ProbeForceFieldSystem, PyObject*, impropers, get_impropers, set_impropers);
%attribute_py(IMP::bff::ProbeForceFieldSystem, PyObject*, probes, get_probes, set_probes);
%attribute_py(IMP::bff::ProbeForceFieldSystem, PyObject*, nonbonded, get_nonbonded, set_nonbonded);
%attribute_py(IMP::bff::ProbeForceFieldSystem, PyObject*, sampling, get_sampling, set_sampling);
