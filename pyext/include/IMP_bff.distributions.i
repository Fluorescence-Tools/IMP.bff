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
%ignore IMP::bff::i0_array(const std::vector<double>&, double**);
%ignore IMP::bff::worm_like_chain(const std::vector<double>&, double, double, bool, bool, double**);
%ignore IMP::bff::worm_like_chain_linker(const std::vector<double>&, double, double, double, bool, double**);
%ignore IMP::bff::ising_chain(const std::vector<double>&, int, double, double, double, double, int, double**);
%ignore IMP::bff::saw_nu(const std::vector<double>&, double, double, double, double**);
%include "IMP/bff/SpecialFunctions.h"
%include "IMP/bff/Distributions.h"
%include "IMP/bff/PolymerChain.h"
// The transfer-polynomial vector evaluator also publishes a 1-D view. It is
// defined here, after the apply above. (This file was once included from
// statesdistance.i; since PRD-138 the distances live in avmodel.i and this
// file is included from swig.i-in directly, at the same position.)
