/*
 * Distances between two labels, whatever represents them.
 *
 * `LabelDistribution` and the converter are C++ values; their scalar state
 * becomes native attributes. The FRET pair values these distances feed are
 * in IMP_bff.fret.i, wrapped just before this file.
 */

IMP_SWIG_VALUE(IMP::bff, FRETDistanceConverter, FRETDistanceConverters);

// Scalar state becomes native attributes. The two
// string fields (position_name, simulation_type) stay methods: SWIG's
// `%attribute` cannot compile a std::string getter.
%attribute(IMP::bff::LabelDistribution, double, simulation_grid_resolution,
          get_simulation_grid_resolution);
%attribute(IMP::bff::LabelDistribution, int, n_points, get_n_points);
%attribute(IMP::bff::ProbeDistributionNormal, double, width, get_width);
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
%include "IMP/bff/StatesDistance.h"
