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
// Both curves are `Node` subclasses, so they need the same shared_ptr holder
// their base has -- declared before the header that defines them.
%apply(double* IN_ARRAY1, int DIM1) {(double* in_axis, int n_axis)};
%shared_ptr(IMP::bff::FcsMdfCurve);
%shared_ptr(IMP::bff::FcsSaturationCurve);
%include "IMP/bff/Fcs.h"

// The transfer-polynomial vector evaluator also publishes a 1-D view. It is
// defined here, after the apply above. (This file was once included from
// statesdistance.i; since PRD-138 the distances live in avmodel.i and this
// file is included from swig.i-in directly, at the same position.)
