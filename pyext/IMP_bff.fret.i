/*
 * FRET between two labels: every pair of two point clouds, and the rate along
 * one dye's trajectory. One header (`FRET.h`, PRD-138), one interface file.
 *
 * The pair geometry and its efficiencies are C++ values whose arrays are read
 * back through the get_*() methods as managed numpy views (ARGOUTVIEWM_ARRAY1)
 * and whose scalar statistics are native attributes via SWIG's `%attribute`:
 * `g.get_R()`, `e.get_E()`, `e.static_efficiency`.
 *
 * The rate trace (the former `FRETRateTrace.h`) takes its orientation factor
 * as a double, so it needs nothing of OrientationFactor.h wrapped first; it
 * used to be wrapped early, beside the quenching fields, and is here now
 * because the pair values must come after their value declarations below.
 */

IMP_SWIG_VALUE(IMP::bff, FRETPairGeometry, FRETPairGeometries);
IMP_SWIG_VALUE(IMP::bff, FRETPairEfficiencies, FRETPairEfficienciesList);

// The matrix members are hidden from SWIG's own member wrapping and
// republished below as shaped attributes: a flat `R` and a `(n1, n2)` `R`
// cannot both be called `R`, and the shape is the one a caller wants.
%ignore IMP::bff::FRETPairGeometry::R;
%ignore IMP::bff::FRETPairGeometry::kappa2;
%ignore IMP::bff::FRETPairGeometry::weight;
%ignore IMP::bff::FRETPairEfficiencies::R;
%ignore IMP::bff::FRETPairEfficiencies::kappa2;
%ignore IMP::bff::FRETPairEfficiencies::weight;
%ignore IMP::bff::FRETPairEfficiencies::E;
%ignore IMP::bff::FRETPairEfficiencies::rate_ratio;
%ignore IMP::bff::FRETPairEfficiencies::k_fret;

%include "IMP/bff/FRET.h"

// The FRETPair values' statistics are public member variables, which SWIG
// exposes directly as members; their arrays are read back through the existing
// get_*() methods (which return managed numpy views, one per matrix) and, as
// attributes, in the `(n1, n2)` shape the pair has. Every consumer reshaped
// those by hand from `n1` and `n2` -- six call sites in the ensemble code
// alone -- and a matrix reshaped wrongly is a matrix transposed in silence.
%attribute_np2v(IMP::bff::FRETPairGeometry, std::vector<double>, R, get_R, n2);
%attribute_np2v(IMP::bff::FRETPairGeometry, std::vector<double>, kappa2, get_kappa2, n2);
%attribute_np2v(IMP::bff::FRETPairGeometry, std::vector<double>, weight, get_weight, n2);
%attribute_np2v(IMP::bff::FRETPairEfficiencies, std::vector<double>, R, get_R, n2);
%attribute_np2v(IMP::bff::FRETPairEfficiencies, std::vector<double>, kappa2, get_kappa2, n2);
%attribute_np2v(IMP::bff::FRETPairEfficiencies, std::vector<double>, weight, get_weight, n2);
%attribute_np2v(IMP::bff::FRETPairEfficiencies, std::vector<double>, E, get_E, n2);
%attribute_np2v(IMP::bff::FRETPairEfficiencies, std::vector<double>, rate_ratio, get_rate_ratio, n2);
%attribute_np2v(IMP::bff::FRETPairEfficiencies, std::vector<double>, k_fret, get_k_fret, n2);
