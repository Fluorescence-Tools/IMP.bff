/*
 * The fields a quenched dye lives in: stickiness, PET rate, FRET rate.
 *
 * Four headers under one shim because the Python that used to wrap them was one
 * module: every function here is a *named composition* of a kernel two headers
 * down -- `slow_factor_grid` is `stamp_spheres` with the multiplying combine,
 * `quenching_rate_map` is `quenching_map` with the axis built for it -- and the
 * names are the ones the physics uses. The kernels are the SWIG surface now:
 * flat in, flat out. The `(ng, ng, ng)` reshape and the `_flat()` ravel that
 * used to live here are the caller's.
 *
 * `radial_diffusion_map` was the one function that stayed Python, because it
 * took a *callable* of the distance from the anchor. It is C++ now and takes
 * the profile pre-evaluated on an integer-Angstrom radius grid instead.
 */

%include "IMP/bff/SolventAccessibleSurface.h"
%include "IMP/bff/QuenchingGrid.h"
%include "IMP/bff/QuenchingMap.h"
%include "IMP/bff/FRETRateTrace.h"