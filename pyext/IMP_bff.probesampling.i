/*
 * Sampling a tethered dye: the walk, the photons, and the equilibrium.
 *
 * `sampling.py` was the shapes these kernels' callers want -- the trajectory is
 * a C++ value, the mobility field is built where the walk is, the photon trace
 * is split into delay and emitted-flag views, and the decay curve accumulates
 * into the caller's histogram. The reshape-vs-flat choice is the caller's.
 */

IMP_SWIG_VALUE(IMP::bff, ProbeDiffusionTrajectory, ProbeDiffusionTrajectories);

// The equilibrium map, published as a flat view (caller reshapes to the cube).
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_view, int* n_out_view)};

%include "IMP/bff/ProbeSampling.h"

// The trajectory's scalar counters and ratio, as read-only attributes (the
// %extend properties they replace are gone; the `get_*()` methods are owned by
// the %attribute and no longer emitted).
%attribute(IMP::bff::ProbeDiffusionTrajectory, int, n_frames, get_n_frames);
%attribute(IMP::bff::ProbeDiffusionTrajectory, int, n_accepted, get_n_accepted);
%attribute(IMP::bff::ProbeDiffusionTrajectory, int, n_rejected, get_n_rejected);
%attribute(IMP::bff::ProbeDiffusionTrajectory, double, acceptance_ratio,
           get_acceptance_ratio);