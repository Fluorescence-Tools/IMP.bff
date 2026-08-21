/*
 * The particle picture of a tethered dye, as a C++ object.
 *
 * The Python held the occupancy grid, the mobility field, the rate map and the
 * trajectory as numpy arrays and drove the walk kernel from a
 * `ThreadPoolExecutor`. Every array crossed the boundary on every access, and
 * `sample_grid` -- tens of millions of indexed reads along the trajectory --
 * was 0.18 s of a profiled suite slice.
 *
 * The threading moved with it, to `std::thread`; the grid fields and the rate
 * trace are read back as managed numpy views through `get_*()` methods, the
 * walk happens in `run()`, and the scalar statistics are native attributes via
 * `%attribute`.
 */

IMP_SWIG_VALUE(IMP::bff, DyeDiffusionSimulation, DyeDiffusionSimulations);

// Array-backed fields publish as 1-D managed numpy views.
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {
    (double** out_view, int* n_out_view)
};

// Scalar statistics as native attributes (no Python in the wrapper).
%attribute(IMP::bff::DyeDiffusionSimulation, int, n_frames, get_n_frames);
%attribute(IMP::bff::DyeDiffusionSimulation, int, n_accepted, get_n_accepted);
%attribute(IMP::bff::DyeDiffusionSimulation, int, n_rejected, get_n_rejected);
%attribute(IMP::bff::DyeDiffusionSimulation, int, ng, get_ng);
%attribute(IMP::bff::DyeDiffusionSimulation, double, dg, get_dg);
%attribute(IMP::bff::DyeDiffusionSimulation, double, t_step, get_t_step, set_t_step);
%attribute(IMP::bff::DyeDiffusionSimulation, double, collision_fraction, get_collision_fraction);

%include "IMP/bff/DyeDiffusion.h"