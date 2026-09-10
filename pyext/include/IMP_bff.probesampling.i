/*
 * Sampling a tethered dye: the walk, the photons, and the equilibrium -- as
 * kernels and as an object. One header (`ProbeSampling.h`, PRD-138), one
 * interface file.
 *
 * The kernels: `sampling.py` was the shapes their callers want -- the
 * trajectory is a C++ value, the mobility field is built where the walk is,
 * the photon trace is split into delay and emitted-flag views, and the decay
 * curve accumulates into the caller's histogram. The reshape-vs-flat choice
 * is the caller's.
 *
 * The object (`ProbeDiffusionSimulation`, the former `ProbeDiffusion.h`): the
 * occupancy grid, the mobility field, the rate map and the trajectory stay on
 * the C++ side -- `sample_grid` alone is tens of millions of indexed reads
 * along a trajectory, and a boundary crossing per read is what makes that
 * expensive. The walk happens in `run()`, the grid fields and the rate trace
 * are read back as managed numpy views, and the scalar statistics are native
 * attributes via `%attribute`.
 */

IMP_SWIG_VALUE(IMP::bff, ProbeDiffusionSimulation, ProbeDiffusionSimulations);
IMP_SWIG_VALUE(IMP::bff, ProbeDiffusionTrajectory, ProbeDiffusionTrajectories);

// Array-backed fields and the equilibrium map publish as 1-D managed numpy
// views (the caller reshapes to the cube).
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_view, int* n_out_view)};

// Scalar statistics as native attributes.
%attribute(IMP::bff::ProbeDiffusionSimulation, int, n_frames, get_n_frames);
%attribute(IMP::bff::ProbeDiffusionSimulation, int, n_accepted, get_n_accepted);
%attribute(IMP::bff::ProbeDiffusionSimulation, int, n_rejected, get_n_rejected);
%attribute(IMP::bff::ProbeDiffusionSimulation, int, ng, get_ng);
%attribute(IMP::bff::ProbeDiffusionSimulation, double, dg, get_dg);
%attribute(IMP::bff::ProbeDiffusionSimulation, double, t_step, get_t_step, set_t_step);
%attribute(IMP::bff::ProbeDiffusionSimulation, double, collision_fraction, get_collision_fraction);

%include "IMP/bff/ProbeSampling.h"

// The trajectory's scalar counters and ratio, as read-only attributes (the
// %extend properties they replace are gone; the `get_*()` methods are owned by
// the %attribute and no longer emitted).
%attribute(IMP::bff::ProbeDiffusionTrajectory, int, n_frames, get_n_frames);
%attribute(IMP::bff::ProbeDiffusionTrajectory, int, n_accepted, get_n_accepted);
%attribute(IMP::bff::ProbeDiffusionTrajectory, int, n_rejected, get_n_rejected);
%attribute(IMP::bff::ProbeDiffusionTrajectory, double, acceptance_ratio,
           get_acceptance_ratio);
