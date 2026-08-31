/*
 * What a label's configuration space is, as C++ values: `States` ->
 * `AccessibleVolume` -> `ACV`, one description of one cloud.
 *
 * The array fields are read back as managed numpy views through `get_*()`
 * methods (ARGOUTVIEWM_ARRAY1) and written through `set_*()` methods taking a
 * flat vector. A caller reads the cloud with `states.get_points()` and
 * reshapes to `(n, 4)` themselves, which is where that shape belongs. The
 * constructors carry everything as flat vectors (via the shared
 * `const std::vector<double>&` typemap, which accepts a contiguous array of
 * any shape or `None` for an absent grid); the only conversion left is `ACV`'s
 * scalar `slow_radius`, which C++ broadcasts over the centres.
 */

IMP_SWIG_VALUE(IMP::bff, States, StatesList);
IMP_SWIG_VALUE(IMP::bff, AccessibleVolume, AccessibleVolumes);
IMP_SWIG_VALUE(IMP::bff, ACV, ACVs);


// The `*_vector` accessors are the C++ side of the same buffers the numpy
// views above publish. Wrapping them would give Python a second, slower way to
// ask the same question, so they stay out of the binding.
%ignore IMP::bff::States::get_points_vector;
%ignore IMP::bff::AccessibleVolume::get_density_vector;
%ignore IMP::bff::AccessibleVolume::get_grid_origin_vector;

%include "IMP/bff/AVModel.h"

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
// `RotamerEnsemble` subclass, so `States` itself did not have them and the
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
%attribute(IMP::bff::AccessibleVolume, int, ng, get_ng);
%attribute(IMP::bff::AccessibleVolume, double, grid_step, get_grid_step);
%attribute(IMP::bff::ACV, double, trapped_fraction, get_trapped_fraction);
