/*
 * The shapes a label distribution takes, as numpy arrays.
 *
 * The kernels publish managed numpy views (ARGOUTVIEWM_ARRAY1) of their
 * density over the distance axis, so a caller gets an ndarray directly --
 * no `np.asarray(...)` wrapper states the shape a second time.
 */

// The axis->density kernels publish their result as a 1-D managed view.
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {
    (double** out_view, int* n_out_view),
    (double** out_pairs, int* n_out_pairs)
};

// The default-argument split of the trailing (double** out_view,
// int* n_out_view) pair grows a middle overload -- the view pointer without
// its dimension -- that the multi-argument typemap cannot bind; dispatched,
// it converts a bare Python object into the pointer and segfaults (the
// module build hit it in histogram_rda; density_to_points' ignores in
// states.i and proberestraints.i are the same story). The no-output and
// full-view forms stay; the middle form is refused at wrap time.
%ignore IMP::bff::normal_distribution(const std::vector<double>&, double, double, bool, double**);
%ignore IMP::bff::generalized_normal_distribution(const std::vector<double>&, double, double, double, bool, double**);
%ignore IMP::bff::distance_between_gaussian(const std::vector<double>&, double, double, bool, double**);
%ignore IMP::bff::gaussian_distance_mixture(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, int, bool, bool, double**);
%ignore IMP::bff::fcs_mdf_g_diff(const std::vector<double>&, double, double, double, double, double, double, double, int, double, bool, int, double, double**);
%ignore IMP::bff::fcs_bunching_factor(const std::vector<double>&, double, const std::vector<double>&, const std::vector<double>&, int, const std::vector<double>&, const std::vector<double>&, double**);
%ignore IMP::bff::fcs_saturated_curve_shape(const std::vector<double>&, double, double, const std::vector<double>&, const std::vector<double>&, int, const std::vector<double>&, double, double, double, bool, int, int, double, const std::vector<double>&, double**);
%ignore IMP::bff::fcs_gaussian_g_diff(const std::vector<double>&, double, double, double, double**);
%ignore IMP::bff::i0_array(const std::vector<double>&, double**);
%ignore IMP::bff::worm_like_chain(const std::vector<double>&, double, double, bool, bool, double**);
%ignore IMP::bff::worm_like_chain_linker(const std::vector<double>&, double, double, double, bool, double**);
%ignore IMP::bff::ising_chain(const std::vector<double>&, int, double, double, double, double, int, double**);
%ignore IMP::bff::saw_nu(const std::vector<double>&, double, double, double, double**);
%include "IMP/bff/SpecialFunctions.h"
%include "IMP/bff/Distributions.h"
%include "IMP/bff/PolymerChain.h"
// C++-only: output-by-reference vectors do not wrap usefully, and Python
// callers have numpy's own hermgauss (which this agrees with to machine
// precision -- the parity test says so).
%ignore IMP::bff::hermgauss;
// The two halves fcs_mdf_g_diff was split into so `FcsMdfCurve` can cache
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
%shared_ptr(IMP::bff::FcsMdfCurve);
%shared_ptr(IMP::bff::FcsSaturationCurve);
%include "IMP/bff/Fcs.h"

// The transfer-polynomial vector evaluator also publishes a 1-D view. It is
// defined here, after the apply above. (This file was once included from
// statesdistance.i; since PRD-138 the distances live in avmodel.i and this
// file is included from swig.i-in directly, at the same position.)
