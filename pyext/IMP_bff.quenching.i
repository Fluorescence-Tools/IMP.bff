/*
 * Collisional quenching of a diffusing dye: the tables, the fields, the decay.
 *
 * `Quenching.h` is four former headers in one, and this is their one interface
 * file. The tables (`Quencher`, `PETParameters`, `ResidueQuenching`) were
 * frozen dataclasses split across `label.py` and `quenching.py`; they are one
 * header now because a rate that takes two partners to define cannot live in
 * only one of them. `quench_radius` and `attenuation_length` are **NaN** when
 * unset: NaN is what a C++ double carries, and for the radius it means
 * "inherit the model-wide critical distance". The tables come back as SWIG
 * `std::map` proxies, which iterate and support `[]`, `in`, `len()` and
 * `items()` but not `.get()`.
 *
 * The fields are named compositions of a kernel further up the same header --
 * `slow_factor_grid` is `stamp_spheres` with the multiplying combine,
 * `quenching_rate_map` is `quenching_map` with the axis built for it -- and
 * the names are the ones the physics uses. The kernels are the SWIG surface:
 * flat in, flat out, with the `(ng, ng, ng)` reshape the caller's.
 * `radial_diffusion_map` takes its profile pre-evaluated on an
 * integer-Angstrom radius grid rather than as a callable of the distance from
 * the anchor.
 *
 * `quenched_donor_photons` takes its occupancy, mobility and rate-map arrays
 * through the typemaps in `IMP_bff.types.i`, which are global, so the header
 * can sit here rather than next to `PhotonSimulation.h` where it used to.
 *
 * The solvent-accessible surface (how buried a quencher is) is wrapped
 * alongside: the fields read it. The FRET-rate trace, which they also read,
 * is in IMP_bff.fret.i since PRD-138 put it in FRET.h with the pair values.
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

%include "IMP/bff/SolventAccessibleSurface.h"
%include "IMP/bff/Quenching.h"

%template(QuencherMap) std::map<std::string, IMP::bff::Quencher>;
%template(PETParametersMap) std::map<std::string, IMP::bff::PETParameters>;
%template(ResidueQuenchingMap) std::map<std::string, IMP::bff::ResidueQuenching>;
%template(PETReferenceMap) std::map<std::string, IMP::bff::PETReference>;
%template(MapStringVectorString) std::map<std::string, std::vector<std::string> >;
%template(VectorInt2) std::vector<int>;
