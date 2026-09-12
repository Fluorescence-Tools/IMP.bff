/*
 * A label's states, whatever represents them -- the cloud as a C++ value, the
 * array kernels a cloud is measured with, and the distances between two.
 * One header (`States.h`, PRD-138), one interface file, wrapped where `States`
 * always was: the representations (avmodel.i, the rotamer ensemble) and
 * everything that consumes a `States` come after.
 *
 * The kernels (the former `AVDistance.h`): `representation/av.py` carried a
 * `_kernels` section that was six functions doing nothing but `.ravel()` on
 * the way in and `.reshape()` on the way out of these. The shape belongs to
 * the caller's picture of the data, so it is stated once, on the kernel
 * itself -- the C++ functions publish `(n, 2)`/`(n, 4)`/`(ng, ng, ng)` numpy
 * views whose shape is part of their contract. All bindings are the standard
 * numpy.i suites -- no Python wrapper states a shape a second time, and no
 * custom typemap re-implements one. They used to be wrapped early, before
 * ProbeAccessibleVolumeDecorator.h, whose defaults name this header's distance enum; a default argument
 * is filled in by the C++ compiler, not by SWIG, so the early position was
 * never load-bearing.
 *
 * `States`: the array fields are read back as managed numpy views through
 * `get_*()` methods (ARGOUTVIEWM_ARRAY1) and written through `set_*()`
 * methods taking a flat vector. A caller reads the cloud with
 * `states.get_points()` and reshapes to `(n, 4)` themselves, which is where
 * that shape belongs.
 *
 * The distances (the former `StatesDistance.h`): free functions, C++ under
 * their public names -- the grid/cloud inputs arrive as flat vectors and
 * `None` for an optional dipole is the empty vector; no wrapper dispatches.
 * The converter is a C++ value whose scalar state becomes native attributes.
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
// IMP's kernel registers refusal typemaps for bare `int*`/`double*`
// parameters (values_like_*: they fail to compile on purpose -- a pointer
// into Python-owned memory is a bug waiting). The 2-D output trio above
// covers the full form of `density_to_points`; the default-argument split
// also generates a partial overload (output + one dimension), whose lone
// `int* n_output1` matches the refusal typemap. Python callers always take
// the full `(n, 4)` view, so the partial form is refused at wrap time.
%ignore IMP::bff::density_to_points(double*,int,int,int,double,const std::vector<double>&,double,double**,int*);
%ignore IMP::bff::density_to_points(double*,int,int,int,double,const std::vector<double>&,double,double**);
%apply(unsigned char** ARGOUTVIEWM_ARRAY3, int* DIM1, int* DIM2, int* DIM3) {(unsigned char** contact, int* contact_dim1, int* contact_dim2, int* contact_dim3)};
%apply(unsigned char** ARGOUTVIEWM_ARRAY3, int* DIM1, int* DIM2, int* DIM3) {(unsigned char** free, int* free_dim1, int* free_dim2, int* free_dim3)};

IMP_SWIG_VALUE(IMP::bff, States, StatesList);
IMP_SWIG_VALUE(IMP::bff, FRETDistanceConverter, FRETDistanceConverters);

%attribute(IMP::bff::FRETDistanceConverter, double, forster_radius,
          get_forster_radius, set_forster_radius);
%attribute(IMP::bff::FRETDistanceConverter, double, sigma, get_sigma, set_sigma);

// The `*_vector` accessor is the C++ side of the same buffer the numpy view
// publishes. Wrapping it would give Python a second, slower way to ask the
// same question, so it stays out of the binding.
%ignore IMP::bff::States::get_points_vector;

// The default-argument split of the trailing (double** out_view,
// int* n_out_view) pair grows a middle overload -- the view pointer without
// its dimension -- that the multi-argument typemap cannot bind; dispatched,
// it converts a bare Python object into the pointer and segfaults (the
// module build hit it in histogram_rda; density_to_points' ignores in
// states.i and proberestraints.i are the same story). The no-output and
// full-view forms stay; the middle form is refused at wrap time.
%ignore IMP::bff::histogram_rda(const States&, const States&, const std::vector<double>&, int, bool, double**);
%include "IMP/bff/States.h"

// Scalar state is exposed as native attributes through SWIG's `%attribute`,
// which maps a C++ getter (and setter) pair to an attribute with no Python in
// the wrapper. The array fields keep their get_*()/set_*() methods
// (ARGOUTVIEWM_ARRAY1), which is how a view stays zero-copy and what the C++
// callers of `get_points` want; `%attribute` cannot express a void+out-param
// getter at all.
//
// What the array *attributes* below add is only the table shape: the cloud is
// `(n, 4)` and a dipole set is `(n, 3)`, and a caller that reshapes by hand is
// a caller that can reshape wrongly -- this repository has paid for exactly
// that twice (a density transposed by a C-order reshape, a `radius1` written
// into three radii). They were hand-written Python properties on the
// `ProbeRotamerEnsemble` subclass, so `States` itself did not have them and the
// shape lived in one representation but not the other; declaring them here
// puts it on the base class, once, for every representation.
%attribute(IMP::bff::States, int, n_points, get_n_points);
%attribute(IMP::bff::States, bool, has_volume, get_has_volume);
%attribute(IMP::bff::States, bool, has_orientations, get_has_orientations);
// `%attributestring`, not `%attribute`: the getter returns the string by
// value, and `%attribute` would hand out the address of that temporary.
%attributestring(IMP::bff::States, std::string, position_name,
                 get_position_name, set_position_name);
%attribute_np2(IMP::bff::States, std::vector<double>, points, get_points, 4,
               set_points);
%attribute_np2(IMP::bff::States, std::vector<double>, orientations,
               get_orientations, 3, set_orientations);
// `mu` is what the rotamer code calls the dipoles, and it *is* `orientations`
// -- one getter, so a caller cannot set one and read a stale other.
%attribute_np2(IMP::bff::States, std::vector<double>, mu, get_orientations, 3,
               set_orientations);
%attribute_np(IMP::bff::States, std::vector<double>, attachment_point,
              get_attachment_point, set_attachment_point);
%attribute_np(IMP::bff::States, std::vector<double>, mean_position,
              get_mean_position);
// The provenance record. `%attribute_py` hands back the wrapped
// `std::map<std::string, std::string>`, which indexes and iterates like a
// dict without copying into one.
%attribute_py(IMP::bff::States, MapStringString, params, get_params,
              set_params);
