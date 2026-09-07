/*
 * The particle picture of a tethered probe, as a C++ object.
 *
 * The occupancy grid, the mobility field, the rate map and the trajectory stay
 * on the C++ side -- `sample_grid` alone is tens of millions of indexed reads
 * along a trajectory, and a boundary crossing per read is what makes that
 * expensive. The walk happens in `run()`, the grid fields and the rate trace
 * are read back as managed numpy views, and the scalar statistics are native
 * attributes via `%attribute`.
 */

IMP_SWIG_VALUE(IMP::bff, ProbeDiffusionSimulation, ProbeDiffusionSimulations);

// Array-backed fields publish as 1-D managed numpy views.
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {
    (double** out_view, int* n_out_view)
};

// Scalar statistics as native attributes.
%attribute(IMP::bff::ProbeDiffusionSimulation, int, n_frames, get_n_frames);
%attribute(IMP::bff::ProbeDiffusionSimulation, int, n_accepted, get_n_accepted);
%attribute(IMP::bff::ProbeDiffusionSimulation, int, n_rejected, get_n_rejected);
%attribute(IMP::bff::ProbeDiffusionSimulation, int, ng, get_ng);
%attribute(IMP::bff::ProbeDiffusionSimulation, double, dg, get_dg);
%attribute(IMP::bff::ProbeDiffusionSimulation, double, t_step, get_t_step, set_t_step);
%attribute(IMP::bff::ProbeDiffusionSimulation, double, collision_fraction, get_collision_fraction);

%include "IMP/bff/ProbeDiffusion.h"