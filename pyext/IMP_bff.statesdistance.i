/*
 * Distances between two labels, whatever represents them.
 *
 * `representation/distance.py` asked all of these of a `States`, and every one
 * of them was a `.ravel()`, a call, and a `.reshape()`. The pair geometry and
 * its efficiencies are C++ values whose arrays are read back through the
 * get_*() methods as managed numpy views (ARGOUTVIEWM_ARRAY1), and whose scalar
 * statistics are native attributes via SWIG's `%attribute`. No dict-indexing
 * sugar remains: a caller writes `g.get_R()`, `e.get_E()`, `e.static_efficiency`.
 */

IMP_SWIG_VALUE(IMP::bff, FRETPairGeometry, FRETPairGeometries);
IMP_SWIG_VALUE(IMP::bff, FRETDistanceConverter, FRETDistanceConverters);
IMP_SWIG_VALUE(IMP::bff, FRETPairEfficiencies, FRETPairEfficienciesList);

// Scalar state becomes native attributes (no Python in the wrapper). The two
// string fields (position_name, simulation_type) stay methods: SWIG's
// `%attribute` cannot compile a std::string getter.
%attribute(IMP::bff::LabelDistribution, double, simulation_grid_resolution,
          get_simulation_grid_resolution);
%attribute(IMP::bff::LabelDistribution, int, n_points, get_n_points);
%attribute(IMP::bff::DyeDistributionNormal, double, width, get_width);
%attribute(IMP::bff::FRETDistanceConverter, double, forster_radius,
          get_forster_radius, set_forster_radius);
%attribute(IMP::bff::FRETDistanceConverter, double, sigma, get_sigma, set_sigma);

// The distributions header is shared; include it here as the include order
// requires (its renames preceded it historically).
%include "IMP_bff.distributions.i"

// The free functions are C++ under their public names -- the grid/cloud inputs
// arrive as flat vectors (the shared `const std::vector<double>&` typemap
// ravels a contiguous array), and `None` for an optional dipole is the empty
// vector. No wrapper dispatches; `_xyz`-style flattening is the typemap's.

// The constructors are already all-flat and all-default; no shadow is needed.
%include "IMP/bff/DistanceCalibration.h"
%include "IMP/bff/FRETPair.h"
%include "IMP/bff/StatesDistance.h"

// The FRETPair values' statistics are public member variables, which SWIG
// exposes directly as members; their arrays are read back through the existing
// get_*() methods (which return managed numpy views, one per matrix).
