// Fluorescence correlation spectroscopy kernels and graph models.
%ignore IMP::bff::fcs_mdf_g_diff(const std::vector<double>&, double, double, double, double, double, double, double, int, double, bool, int, double, double**);
%ignore IMP::bff::fcs_bunching_factor(const std::vector<double>&, double, const std::vector<double>&, const std::vector<double>&, int, const std::vector<double>&, const std::vector<double>&, double**);
%ignore IMP::bff::fcs_saturated_curve_shape(const std::vector<double>&, double, double, const std::vector<double>&, const std::vector<double>&, int, const std::vector<double>&, double, double, double, bool, int, int, double, const std::vector<double>&, double**);
%ignore IMP::bff::fcs_gaussian_g_diff(const std::vector<double>&, double, double, double, double**);
// C++-only: output-by-reference vectors do not wrap usefully, and Python
// callers have numpy's own hermgauss (which this agrees with to machine
// precision -- the parity test says so).
%ignore IMP::bff::hermgauss;
// The two halves fcs_mdf_g_diff was split into so `FCSMdfCurve` can cache
// them across evaluations. C++-only for the same reason: one writes through
// output-by-reference vectors, and the other takes six of them. A Python
// caller wants the whole curve, which is fcs_mdf_g_diff or the node.
%ignore IMP::bff::fcs_mdf_axial_profiles;
%ignore IMP::bff::fcs_mdf_g_raw;
// The MDF diffusion shape as a graph node, so the `"mdf"` FCS mode joins the
// closed-form ones on the graph instead of returning to Python per iteration.
// Both curves are `GraphNode` subclasses, so they need the same shared_ptr holder
// their base has -- declared before the header that defines them.
%apply(double* IN_ARRAY1, int DIM1) {(double* in_axis, int n_axis)};
%shared_ptr(IMP::bff::FCSMdfCurve);
%shared_ptr(IMP::bff::FCSSaturationCurve);
%include "IMP/bff/FCS.h"
