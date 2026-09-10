/*
 * The field picture of a tethered probe, as a C++ object.
 *
 * The solver takes its four ng^3 grids straight from numpy and folds `dt` into
 * them itself -- `d = D dt/dg^2` and `decay = exp(-k dt)` -- so a caller never
 * carries pre-scaled grids across the boundary. `run()` and `gradient()`
 * return C++ values whose arrays are managed numpy views.
 */

IMP_SWIG_VALUE(IMP::bff, GridDiffusionSolver, GridDiffusionSolvers);
IMP_SWIG_VALUE(IMP::bff, GridDiffusionResult, GridDiffusionResults);
IMP_SWIG_VALUE(IMP::bff, GridDiffusionGradient, GridDiffusionGradients);

// Every grid field and result array publishes as a 1-D managed numpy view:
// get_density(), get_time(), get_fluorescence(), get_d_diffusion(), ... all
// return ndarrays via the ARGOUTVIEWM_ARRAY1 typemap, with no property sugar:
// a caller writes `solver.get_density()` rather than `solver.density`.
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {
    (double** out_view, int* n_out_view)
};

%include "IMP/bff/GridDiffusionSolver.h"