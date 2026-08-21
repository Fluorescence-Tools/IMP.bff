/*
 * What a label's configuration space is, as C++ values.
 *
 * `States` is the surface every representation supplies -- positions, weights,
 * orientations -- and it was a Python dataclass in
 * `representation/distance.py`, with `AccessibleVolume` deriving from it and
 * `BasicAV`/`ACV` sitting beside them in C++ describing the same cloud a second
 * time. Two objects for one thing is how the two came to disagree about the
 * axis order of `density` without anything noticing (PRD-113 stage 3a). There
 * is one now: `States` -> `AccessibleVolume` -> `ACV`, all C++.
 *
 * The array fields are read back as managed numpy views through `get_*()`
 * methods (ARGOUTVIEWM_ARRAY1) and written through `set_*()` methods taking a
 * flat vector -- no property sugar, no `%pythoncode`. A caller reads the cloud
 * with `states.get_points()` and reshapes to `(n, 4)` themselves, which is
 * where that shape belongs. The constructors carry everything as flat vectors
 * (via the shared `const std::vector<double>&` typemap, which accepts a
 * contiguous array of any shape or `None` for an absent grid); the only
 * conversion left is `ACV`'s scalar `slow_radius`, which C++ broadcasts over
 * the centres.
 */

IMP_SWIG_VALUE(IMP::bff, States, StatesList);
IMP_SWIG_VALUE(IMP::bff, AccessibleVolume, AccessibleVolumes);
IMP_SWIG_VALUE(IMP::bff, ACV, ACVs);

%include "IMP/bff/AVModel.h"

// Scalar state is exposed as native attributes through SWIG's `%attribute`,
// which maps a C++ getter (and setter) pair to an attribute with no Python in
// the wrapper. The array fields are read back as managed numpy views through
// the get_*()/set_*() methods (ARGOUTVIEWM_ARRAY1), which is how a view keeps
// its zero-copy property; `%attribute` cannot express a void+out-param getter.
%attribute(IMP::bff::States, int, n_points, get_n_points);
%attribute(IMP::bff::States, bool, has_volume, get_has_volume);
%attribute(IMP::bff::States, bool, has_orientations, get_has_orientations);
%attribute(IMP::bff::AccessibleVolume, int, ng, get_ng);
%attribute(IMP::bff::AccessibleVolume, double, grid_step, get_grid_step);
%attribute(IMP::bff::ACV, double, trapped_fraction, get_trapped_fraction);
