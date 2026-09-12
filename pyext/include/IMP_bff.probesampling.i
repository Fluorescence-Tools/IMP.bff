/* Probe sampling kernels, trajectory values and library helpers. */
IMP_SWIG_VALUE(IMP::bff, ProbeDiffusionTrajectory, ProbeDiffusionTrajectories);
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
