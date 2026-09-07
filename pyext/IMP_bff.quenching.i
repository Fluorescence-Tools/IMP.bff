/*
 * The fields a quenched probe lives in: stickiness, PET rate, FRET rate.
 *
 * Four headers under one interface file because every function here is a
 * *named composition* of a kernel two headers down -- `slow_factor_grid` is
 * `stamp_spheres` with the multiplying combine, `quenching_rate_map` is
 * `quenching_map` with the axis built for it -- and the names are the ones the
 * physics uses. The kernels are the SWIG surface: flat in, flat out, with the
 * `(ng, ng, ng)` reshape the caller's. `radial_diffusion_map` takes its
 * profile pre-evaluated on an integer-Angstrom radius grid rather than as a
 * callable of the distance from the anchor.
 */

%include "IMP/bff/SolventAccessibleSurface.h"
%include "IMP/bff/QuenchingGrid.h"
%include "IMP/bff/QuenchingMap.h"
%include "IMP/bff/FRETRateTrace.h"