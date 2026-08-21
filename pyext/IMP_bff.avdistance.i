/*
 * The array kernels an accessible volume is measured with.
 *
 * `representation/av.py` carried a `_kernels` section that was six functions
 * doing nothing but `.ravel()` on the way in and `.reshape()` on the way out of
 * these. The shape belongs to the caller's picture of the data, so it is stated
 * once, on the kernel itself -- the C++ functions publish `(n, 2)`/`(n, 4)`/
 * `(ng, ng, ng)` numpy views whose shape is part of their contract.
 *
 * All bindings are the standard numpy.i suites -- no Python wrapper states a
 * shape a second time, and no custom typemap re-implements one.
 */

// A point cloud is an `(n, 4)` array: x, y, z, weight. random_distances reads
// both clouds row-by-row, so both bind through the stock 2-D input suite.
%apply(double* IN_ARRAY2, int DIM1, int DIM2) {(double* p1, int n_p1, int n_p1c)};
%apply(double* IN_ARRAY2, int DIM1, int DIM2) {(double* p2, int n_p2, int n_p2c)};

// The shaped managed views, bound before the header declares the functions.
// `random_distances` and `density_to_points` bind to the ARRAY2 output typemap
// claimed in types.i on `(double** output, int* n_output1, int* n_output2)`.
// The three-view kernels below need bindings types.i does not have: an IN_ARRAY3
// input for the density cube, an int ARRAY3 for the contact label, and two
// unsigned-char ARRAY3s for the (contact, free) mask pair.
%apply(double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(double* density, int nx, int ny, int nz)};
%apply(double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(double* density, int ng, int ng2, int ng3)};
%apply(double* IN_ARRAY2, int DIM1, int DIM2) {(double* rs, int n_rs, int n_rsc)};
%apply(int** ARGOUTVIEWM_ARRAY3, int* DIM1, int* DIM2, int* DIM3) {(int** output_i, int* dim1, int* dim2, int* dim3)};
%apply(unsigned char** ARGOUTVIEWM_ARRAY3, int* DIM1, int* DIM2, int* DIM3) {(unsigned char** contact, int* contact_dim1, int* contact_dim2, int* contact_dim3)};
%apply(unsigned char** ARGOUTVIEWM_ARRAY3, int* DIM1, int* DIM2, int* DIM3) {(unsigned char** free, int* free_dim1, int* free_dim2, int* free_dim3)};

%include "IMP/bff/AVDistance.h"
