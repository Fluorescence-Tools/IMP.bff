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
 * `points`, `orientations`, `mean_position` and the density grids are returned
 * as managed numpy views over the object's own buffers, so reading them is
 * free -- and the `.reshape(-1, 4)` / `.reshape(ng, ng, ng)` that callers
 * expect stays here, because a view over a flat buffer is flat and the shape
 * belongs to the caller's picture of the data, not to the storage.
 */

IMP_SWIG_VALUE(IMP::bff, States, StatesList);
IMP_SWIG_VALUE(IMP::bff, AccessibleVolume, AccessibleVolumes);
IMP_SWIG_VALUE(IMP::bff, ACV, ACVs);

// Every constructor takes every argument by keyword and every one is optional,
// which SWIG cannot express for a function with defaults. The shadows below
// also do the conversions the Python constructors did: `None` for an absent
// array, a scalar `slow_radius` broadcast over the centres, and a `params`
// mapping whose values are stringified, because it is provenance rather than
// numbers to compute with.
%feature("shadow") IMP::bff::States::States %{
def __init__(self, points=None, attachment_point=None, orientations=None,
             position_name="", params=None):
    _IMP_bff.States_swiginit(self, _IMP_bff.new_States(
        _av_flat(points), _av_flat(attachment_point), _av_flat(orientations),
        str(position_name), _av_params(params)))
%}

%feature("shadow") IMP::bff::AccessibleVolume::AccessibleVolume %{
def __init__(self, points=None, density=None, grid_origin=None, grid_step=1.5,
             position_name="", attachment_point=None, orientations=None,
             params=None):
    _IMP_bff.AccessibleVolume_swiginit(self, _IMP_bff.new_AccessibleVolume(
        _av_flat(points), _av_flat(density), _av_flat(grid_origin),
        float(grid_step), str(position_name), _av_flat(attachment_point),
        _av_flat(orientations), _av_params(params)))
%}

%feature("shadow") IMP::bff::ACV::ACV %{
def __init__(self, points=None, density=None, grid_origin=None, grid_step=1.5,
             slow_centers=None, slow_radius=10.0, trapped_fraction=0.8,
             position_name=""):
    centres = _av_flat(slow_centers)
    radius = np.atleast_1d(np.asarray(slow_radius, dtype=np.float64)).ravel()
    if radius.size == 1 and centres.size:
        radius = np.full(centres.size // 3, float(radius[0]))
    _IMP_bff.ACV_swiginit(self, _IMP_bff.new_ACV(
        _av_flat(points), _av_flat(density), _av_flat(grid_origin),
        float(grid_step), centres, np.ascontiguousarray(radius),
        float(trapped_fraction), str(position_name)))
%}

%attribute_py(IMP::bff::States, int, n_points, get_n_points);
%attribute_py(IMP::bff::States, bool, has_volume, get_has_volume);
%attribute_py(IMP::bff::States, bool, has_orientations, get_has_orientations);
%attribute_py(IMP::bff::States, std::string, position_name,
              get_position_name, set_position_name);
%attribute_py(IMP::bff::AccessibleVolume, double, grid_step, get_grid_step);
%attribute_py(IMP::bff::ACV, double, trapped_fraction, get_trapped_fraction);

%include "IMP/bff/AVModel.h"

%template(MapStringString) std::map<std::string, std::string>;

%pythoncode %{
def _av_flat(a):
    """An array argument as a flat float64 buffer; `None` is empty."""
    return (np.zeros(0) if a is None
            else np.ascontiguousarray(np.asarray(a, dtype=np.float64)).ravel())


def _av_params(p):
    """A provenance mapping with every value stringified.

    `params` records *how* a set of states was produced -- backend, linker
    length, library name, temperature. It is written by whatever built the
    states and read back as a record, never computed with, so one type for the
    values is enough and a string is the one that survives a nested value.
    """
    return {} if not p else {str(k): str(v) for k, v in dict(p).items()}
%}

%extend IMP::bff::States {
    %pythoncode %{
        @property
        def points(self):
            """The cloud as ``(n, 4)`` -- x, y, z, weight."""
            return _IMP_bff.States_get_points(self).reshape(-1, 4)

        @points.setter
        def points(self, value):
            self.set_points(_av_flat(value))

        @property
        def positions(self):
            """``(n, 3)`` state coordinates."""
            return self.points[:, :3]

        @property
        def weights(self):
            """``(n,)`` state weights, unnormalised."""
            return self.points[:, 3]

        @property
        def orientations(self):
            """``(n, 3)`` transition dipoles, or ``None`` when unresolved.

            ``None`` rather than an empty array because that is the question
            callers ask: a representation without dipoles forces an isotropic
            kappa^2, and it should have to be *told* that rather than reading a
            zero-length array as an answer.
            """
            mu = _IMP_bff.States_get_orientations(self)
            return None if mu.size == 0 else mu.reshape(-1, 3)

        @orientations.setter
        def orientations(self, value):
            self.set_orientations(_av_flat(value))

        @property
        def attachment_point(self):
            """``(3,)`` where the label is tied, or ``None``."""
            a = _IMP_bff.States_get_attachment_point(self)
            return None if a.size == 0 else a

        @attachment_point.setter
        def attachment_point(self, value):
            self.set_attachment_point(_av_flat(value))

        @property
        def mean_position(self):
            """Weight-averaged position, ``(3,)``.

            Falls back to the attachment point for an empty or zero-weight
            cloud, because that is the one position a label always has.
            """
            return _IMP_bff.States_get_mean_position(self)

        @property
        def params(self):
            """How these states were produced. Provenance, values as strings."""
            return dict(self.get_params())

        @params.setter
        def params(self, value):
            self.set_params(_av_params(value))

        def pRDA(self, other, axis=None, n_samples=50000):
            """``p(R_DA)`` and the bin centres, as ``(y, x)``.

            The default axis is 0-200 A in 0.5 A steps, which is the range a
            FRET pair can report over; it is stated here rather than in C++ so
            the two are not free to drift apart.
            """
            if axis is None:
                axis = np.arange(0.0, 200.0, 0.5)
            axis = np.ascontiguousarray(np.asarray(axis, dtype=np.float64))
            y = _IMP_bff.States_pRDA(self, other, axis, int(n_samples))
            return y, 0.5 * (axis[:-1] + axis[1:])
    %}
}

%extend IMP::bff::AccessibleVolume {
    %pythoncode %{
        @property
        def density(self):
            """The grid as ``(ng, ng, ng)``, or ``None`` when none is carried.

            ``None`` rather than an empty array because callers test it -- the
            Python this replaces stored ``None`` and every consumer asks
            ``if av.density is not None``."""
            flat = _IMP_bff.AccessibleVolume_get_density(self)
            if flat.size == 0:
                return None
            ng = self.get_ng()
            return flat.reshape(ng, ng, ng)

        @property
        def grid_origin(self):
            o = _IMP_bff.AccessibleVolume_get_grid_origin(self)
            return None if o.size == 0 else o

        @property
        def grid_shape(self):
            """``(ng, ng, ng)``. The grid is cubic -- ``PathMapHeader`` sizes it
            from the linker length in all three axes."""
            ng = self.get_ng()
            return (ng, ng, ng)
    %}
}

%extend IMP::bff::ACV {
    %pythoncode %{
        @staticmethod
        def from_accessible_volume(av, slow_centers=None, slow_radius=10.0,
                                   trapped_fraction=0.8):
            """Split a solved volume into contact and free, on its own grid.

            The same two conversions the constructor does -- `None` for an
            absent array and a scalar radius broadcast over the centres --
            because a caller that has one radius for every centre should not
            have to say it once per centre.
            """
            centres = _av_flat(slow_centers)
            radius = np.atleast_1d(
                np.asarray(slow_radius, dtype=np.float64)).ravel()
            if radius.size == 1 and centres.size:
                radius = np.full(centres.size // 3, float(radius[0]))
            return _IMP_bff.ACV_from_accessible_volume(
                av, centres, np.ascontiguousarray(radius),
                float(trapped_fraction))

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
