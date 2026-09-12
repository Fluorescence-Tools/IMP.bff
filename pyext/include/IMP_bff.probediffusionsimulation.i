/* A probe diffusion simulation; fields remain in C++ until requested. */
IMP_SWIG_VALUE(IMP::bff, ProbeDiffusionSimulation, ProbeDiffusionSimulations);
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_view, int* n_out_view)};

// Scalar statistics as native attributes.
%attribute(IMP::bff::ProbeDiffusionSimulation, int, n_frames, get_n_frames);
%attribute(IMP::bff::ProbeDiffusionSimulation, int, n_accepted, get_n_accepted);
%attribute(IMP::bff::ProbeDiffusionSimulation, int, n_rejected, get_n_rejected);
%attribute(IMP::bff::ProbeDiffusionSimulation, int, ng, get_ng);
%attribute(IMP::bff::ProbeDiffusionSimulation, double, dg, get_dg);
%attribute(IMP::bff::ProbeDiffusionSimulation, double, t_step, get_t_step, set_t_step);
%attribute(IMP::bff::ProbeDiffusionSimulation, double, collision_fraction, get_collision_fraction);

%include "IMP/bff/ProbeDiffusionSimulation.h"
