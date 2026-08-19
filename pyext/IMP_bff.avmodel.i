/*
 * The accessible volume as a C++ value, with the Python surface it had.
 *
 * `BasicAV` and `ACV` were 362 lines of Python holding a point cloud and a
 * density grid as numpy arrays, calling C++ for every actual computation. The
 * arithmetic was already across the boundary; the *object* was not, which is
 * what keeps a module Python-carried.
 *
 * `points`, `mean_position` and the density grids are returned as managed numpy
 * views over the object's own buffers, so reading them is free -- and the
 * `.reshape(-1, 4)` / `.reshape(ng, ng, ng)` that callers expect stays here,
 * because a view over a flat buffer is flat and the shape belongs to the
 * caller's picture of the data, not to the storage.
 */

IMP_SWIG_VALUE(IMP::bff, BasicAV, BasicAVs);
IMP_SWIG_VALUE(IMP::bff, ACV, ACVs);

// Both constructors take every argument by keyword and every one is optional,
// which SWIG cannot express for a function with defaults. The shadows below
// also do the two conversions the Python constructors did: `None` for an
// absent array, and a scalar `slow_radius` broadcast over the centres.
%feature("shadow") IMP::bff::BasicAV::BasicAV %{
def __init__(self, points=None, density=None, grid_origin=None,
             grid_step=1.5, position_name=""):
    import numpy as _np
    def _flat(a):
        return (_np.zeros(0) if a is None
                else _np.ascontiguousarray(_np.asarray(a, dtype=_np.float64)).ravel())
    _IMP_bff.BasicAV_swiginit(self, _IMP_bff.new_BasicAV(
        _flat(points), _flat(density), _flat(grid_origin),
        float(grid_step), str(position_name)))
%}

%feature("shadow") IMP::bff::ACV::ACV %{
def __init__(self, points=None, density=None, grid_origin=None, grid_step=1.5,
             slow_centers=None, slow_radius=10.0, trapped_fraction=0.8,
             position_name=""):
    import numpy as _np
    def _flat(a):
        return (_np.zeros(0) if a is None
                else _np.ascontiguousarray(_np.asarray(a, dtype=_np.float64)).ravel())
    centres = _flat(slow_centers)
    radius = _np.atleast_1d(_np.asarray(slow_radius, dtype=_np.float64)).ravel()
    if radius.size == 1 and centres.size:
        radius = _np.full(centres.size // 3, float(radius[0]))
    _IMP_bff.ACV_swiginit(self, _IMP_bff.new_ACV(
        _flat(points), _flat(density), _flat(grid_origin), float(grid_step),
        centres, _np.ascontiguousarray(radius), float(trapped_fraction),
        str(position_name)))
%}

%attribute_py(IMP::bff::BasicAV, int, n_points, get_n_points);
%attribute_py(IMP::bff::BasicAV, double, grid_step, get_grid_step);
%attribute_py(IMP::bff::BasicAV, std::string, position_name,
              get_position_name, set_position_name);
%attribute_py(IMP::bff::ACV, double, trapped_fraction, get_trapped_fraction);

%include "IMP/bff/AVModel.h"

%extend IMP::bff::BasicAV {
    %pythoncode %{
        @property
        def points(self):
            """The cloud as ``(n, 4)`` -- x, y, z, weight."""
            return _IMP_bff.BasicAV_get_points(self).reshape(-1, 4)

        @property
        def density(self):
            """The grid as ``(ng, ng, ng)``, or ``None`` when none is carried.

            ``None`` rather than an empty array because callers test it -- the
            Python this replaces stored ``None`` and every consumer asks
            ``if av.density is not None``."""
            flat = _IMP_bff.BasicAV_get_density(self)
            if flat.size == 0:
                return None
            ng = self.get_ng()
            return flat.reshape(ng, ng, ng)

        @property
        def grid_origin(self):
            o = _IMP_bff.BasicAV_get_grid_origin(self)
            return None if o.size == 0 else o

        @property
        def mean_position(self):
            """Density-weighted mean position, ``(3,)``."""
            return _IMP_bff.BasicAV_get_mean_position(self)

        def pRDA(self, other, axis=None, n_samples=50000):
            """``p(R_DA)`` and the bin centres, as ``(y, x)``.

            The default axis is 0-200 A in 0.5 A steps, which is the range a
            FRET pair can report over; it is stated here rather than in C++ so
            the two are not free to drift apart.
            """
            if axis is None:
                axis = np.arange(0.0, 200.0, 0.5)
            axis = np.ascontiguousarray(np.asarray(axis, dtype=np.float64))
            y = _IMP_bff.BasicAV_pRDA(self, other, axis, int(n_samples))
            return y, 0.5 * (axis[:-1] + axis[1:])
    %}
}

%extend IMP::bff::ACV {
    %pythoncode %{
        @property
        def slow_centers(self):
            c = _IMP_bff.ACV_get_slow_centers(self)
            return None if c.size == 0 else c.reshape(-1, 3)

        @property
        def slow_radius(self):
            return _IMP_bff.ACV_get_slow_radius(self)

        @property
        def contact_density(self):
            flat = _IMP_bff.ACV_get_contact_density(self)
            if flat.size == 0:
                return None
            ng = self.get_ng()
            return flat.reshape(ng, ng, ng)
    %}
}
