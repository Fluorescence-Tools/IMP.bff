/*
 * The field picture of a tethered dye, as a C++ object.
 *
 * The Python held four ng^3 grids and folded `dt` into two of them on every
 * call -- `d = D dt/dg^2` and `decay = exp(-k dt)` -- then handed all four
 * across the boundary. All of that is C++ now: the solver takes the grids
 * straight from numpy, and `run()`/`gradient()` return the result as C++ value
 * types whose arrays are read back as numpy views through `get_*()` methods.
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