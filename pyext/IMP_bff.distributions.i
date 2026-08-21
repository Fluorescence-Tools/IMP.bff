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

%include "IMP/bff/Distributions.h"
%include "IMP/bff/PolymerChain.h"

// The transfer-polynomial vector evaluator also publishes a 1-D view. It is
// defined here (after the apply above and before statesdistance.i, which
// includes this file and exposes the scalar dispatch functions from
// DistanceCalibration.h / StatesDistance.h).
%include "IMP/bff/DistanceCalibration.h"
