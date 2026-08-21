/*
 * Photoinduced electron transfer: which moieties quench, and how hard.
 *
 * The three values -- `Quencher` (identity), `PETParameters` (the pair's rate)
 * and `ResidueQuenching` (what a diffusing dye sees) -- were frozen dataclasses
 * split across `label.py` and `quenching.py`, with the tables in one and the
 * types that read them in the other. They are one header now, because a rate
 * that takes two partners to define cannot live in only one of them.
 *
 * `quench_radius` and `attenuation_length` are **NaN** where the Python had
 * `None`: NaN is what a C++ double carries, and for the radius it means
 * "inherit the model-wide critical distance" exactly as `None` did.
 *
 * The tables come back as SWIG `std::map` proxies, which iterate and support
 * `[]`, `in`, `len()` and `items()` but not `.get()`.
 */

IMP_SWIG_VALUE(IMP::bff, Quencher, Quenchers);
IMP_SWIG_VALUE(IMP::bff, PETParameters, PETParametersList);
IMP_SWIG_VALUE(IMP::bff, ResidueQuenching, ResidueQuenchings);
IMP_SWIG_VALUE(IMP::bff, PETReference, PETReferences);
IMP_SWIG_VALUE(IMP::bff, ResidueSites, ResidueSitesList);

// Two output views from one call: `atomic_quenching_parameters` returns kQ and
// rC together because they are read together and looked up once.
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_kQ, int* n_out_kQ)};
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_rC, int* n_out_rC)};

// The atoms, straight from numpy.
%apply(double* IN_ARRAY2, int DIM1, int DIM2) {
    (double* coords, int n_atoms, int n_dim)
};

// Scalar bool getters map straight to read-only attributes; the matching
// `get_is_typed()` method is suppressed (SWIG %attribute owns the name).
%attribute(IMP::bff::Quencher, bool, is_typed, get_is_typed);
%attribute(IMP::bff::PETParameters, bool, is_transferred, get_is_transferred);

// A build helper over raw pointers; nothing Python-side appends to a result.
%ignore IMP::bff::ResidueSites::add;

%include "IMP/bff/PETQuenching.h"

%template(QuencherMap) std::map<std::string, IMP::bff::Quencher>;
%template(PETParametersMap) std::map<std::string, IMP::bff::PETParameters>;
%template(ResidueQuenchingMap) std::map<std::string, IMP::bff::ResidueQuenching>;
%template(PETReferenceMap) std::map<std::string, IMP::bff::PETReference>;
%template(MapStringVectorString) std::map<std::string, std::vector<std::string> >;
%template(VectorInt2) std::vector<int>;