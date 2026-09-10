/*
 * LabelLib's interface, name for name (PRD-137 step 6d).
 *
 * LabelLib (https://github.com/Fluorescence-Tools/LabelLib) is the small
 * library that computes an accessible volume and the observables read off
 * one. bff computes the same volumes -- and the rest of the fluorescence
 * stack besides -- so a caller should not have to install both, nor rewrite
 * a working script to move over. `IMP.bff.labellib` is therefore LabelLib's
 * API exactly: the same seven functions, the same argument order, the same
 * array conventions ((4, N) or (N, 4) atoms, (3,) source), and a `Grid3D`
 * with LabelLib's `shape`, `originXYZ`, `discStep`, `grid` (flat, x fastest)
 * and `points()` ((4, n): x, y, z, weight).
 *
 *     import LabelLib as ll        ->    from IMP.bff import labellib as ll
 *
 * The numbers are bff's own -- the same lattice search `IMP.bff.get_av` runs
 * -- so they agree with LabelLib to within the discretisation, not bit for
 * bit; `test/test_labellib_interface.py` measures that against the real
 * LabelLib wherever it is installed.
 */
%pythoncode %{
import sys as _ll_sys
import types as _ll_types

_labellib = _ll_types.ModuleType(
    __name__ + ".labellib",
    """LabelLib's API, computed by bff.

Drop-in for `import LabelLib`: same names, same signatures, same array
conventions. See IMP.bff's documentation for what differs (nothing in the
interface; the volumes are bff's lattice search).
""")


def _ll_atoms(atoms_xyzr):
    """LabelLib takes (4, N) and its `_arr` twins take (N, 4) too. bff's
    get_av wants (N, 4), so this settles the orientation once. A (4, 4) is
    read as (4, N), which is LabelLib's own reading."""
    import numpy as np
    a = np.ascontiguousarray(atoms_xyzr, dtype=np.float64)
    if a.ndim != 2:
        raise ValueError("atomsXyzr must be 2-dimensional, (4, N) or (N, 4)")
    if a.shape[0] == 4:
        return np.ascontiguousarray(a.T)
    if a.shape[1] == 4:
        return a
    raise ValueError("atomsXyzr must be (4, N) or (N, 4), got %r" % (a.shape,))


def _ll_vec3(v, name):
    import numpy as np
    a = np.asarray(v, dtype=np.float64).ravel()
    if a.size != 3:
        raise ValueError("%s must have three components" % name)
    return [float(a[0]), float(a[1]), float(a[2])]


class Grid3D(object):
    """LabelLib's Grid3D: a scalar cube with an origin and a spacing.

    `grid` is flat with **x fastest** (LabelLib's order, numpy's "F"), so
    `numpy.array(g.grid).reshape(g.shape, order="F")` indexes as [ix, iy, iz].
    """

    def __init__(self, shape, originXYZ, discStep, grid=None):
        self.shape = [int(s) for s in shape]
        self.originXYZ = [float(o) for o in originXYZ]
        self.discStep = float(discStep)
        n = self.shape[0] * self.shape[1] * self.shape[2]
        self.grid = [0.0] * n if grid is None else list(grid)

    def points(self):
        """The occupied voxels as (4, n): x, y, z, weight."""
        import numpy as np
        g = np.asarray(self.grid, dtype=np.float64).reshape(self.shape, order="F")
        idx = np.argwhere(g > 0.0)
        if idx.size == 0:
            return np.zeros((4, 0), dtype=np.float64)
        xyz = np.asarray(self.originXYZ) + idx * self.discStep
        w = g[idx[:, 0], idx[:, 1], idx[:, 2]]
        return np.vstack([xyz.T, w[None, :]])

    def __repr__(self):
        return "Grid3D(shape=%r, originXYZ=%r, discStep=%r)" % (
            self.shape, self.originXYZ, self.discStep)


def _ll_grid_from_av(av):
    """An IMP.bff AccessibleVolume as a LabelLib Grid3D. bff writes its cube
    in (x, y, z) order (x slowest); LabelLib's is x fastest."""
    import numpy as np
    d = np.asarray(av.get_density(), dtype=np.float64)
    o = list(np.asarray(av.get_grid_origin(), dtype=np.float64).ravel())
    step = av.get_grid_step()
    n = int(round(len(d) ** (1.0 / 3.0)))
    while n * n * n < len(d):
        n += 1
    shape = [n, n, n]
    return Grid3D(shape, o, step, d.reshape(shape).ravel(order="F"))


def dyeDensityAV1(atomsXyzr, sourceXyz, linkerLength, linkerDiameter,
                  dyeRadius, discStep):
    """Accessible volume of a spherical probe (LabelLib's AV1)."""
    atoms = _ll_atoms(atomsXyzr)
    av = get_av(atoms, _ll_vec3(sourceXyz, "sourceXyz"),
                float(linkerLength), float(linkerDiameter),
                float(dyeRadius), 0.0, 0.0, float(discStep))
    return _ll_grid_from_av(av)


def dyeDensityAV3(atomsXyzr, sourceXyz, linkerLength, linkerDiameter,
                  dyeRadii, discStep):
    """Accessible volume of a three-radius probe (LabelLib's AV3)."""
    import numpy as np
    r = np.asarray(dyeRadii, dtype=np.float64).ravel()
    if r.size != 3:
        raise ValueError("dyeRadii must have three components")
    atoms = _ll_atoms(atomsXyzr)
    av = get_av(atoms, _ll_vec3(sourceXyz, "sourceXyz"),
                float(linkerLength), float(linkerDiameter),
                float(r[0]), float(r[1]), float(r[2]), float(discStep))
    return _ll_grid_from_av(av)


def minLinkerLength(atomsXyzr, sourceXyz, linkerLength, linkerDiameter,
                    dyeRadius, discStep):
    """The shortest obstacle-free linker path to every voxel; unreachable
    voxels are negative."""
    import numpy as np
    # the array door, which reads the obstacles and the attachment site the
    # way get_av does, so that this grid and the volume above are two
    # read-outs of one search rather than two searches
    atoms = _ll_atoms(atomsXyzr)
    grid = get_linker_path_lengths(
        atoms, _ll_vec3(sourceXyz, "sourceXyz"),
        float(linkerLength), float(linkerDiameter), float(dyeRadius),
        float(discStep))
    h = grid.get_header()
    shape = [h.get_nx(), h.get_ny(), h.get_nz()]
    d = np.asarray([grid.get_value(i) for i in range(shape[0] * shape[1] * shape[2])],
                   dtype=np.float64)
    origin = list(grid.get_origin())
    return Grid3D(shape, origin, grid.get_spacing(), d.reshape(shape).ravel(order="F"))


def addWeights(grid, xyzRQ):
    """Add a local weight Q to every occupied voxel within R of (x, y, z).

    `xyzRQ` is (5, N) or (N, 5): x, y, z, radius, weight."""
    import numpy as np
    a = np.ascontiguousarray(xyzRQ, dtype=np.float64)
    if a.ndim != 2:
        raise ValueError("xyzRQ must be 2-dimensional, (5, N) or (N, 5)")
    if a.shape[0] == 5:
        a = a.T
    elif a.shape[1] != 5:
        raise ValueError("xyzRQ must be (5, N) or (N, 5), got %r" % (a.shape,))
    g = np.asarray(grid.grid, dtype=np.float64).reshape(grid.shape, order="F")
    o = np.asarray(grid.originXYZ)
    step = grid.discStep
    ii = np.indices(grid.shape).reshape(3, -1).T
    xyz = o + ii * step
    flat = g.ravel(order="F").copy()
    occupied = flat > 0.0
    for x, y, z, r, q in a:
        d2 = ((xyz - np.array([x, y, z])) ** 2).sum(axis=1)
        flat[occupied & (d2 <= r * r)] += q
    return Grid3D(grid.shape, grid.originXYZ, step, flat)


def meanDistance(g1, g2, nsamples=100000):
    """The mean distance between two clouds."""
    return average_distance(list(g1.points().T.ravel()),
                            list(g2.points().T.ravel()), int(nsamples), 0)


def meanEfficiency(g1, g2, R0, nsamples=100000):
    """The mean FRET efficiency between two clouds."""
    r = mean_fret_distance(list(g1.points().T.ravel()),
                           list(g2.points().T.ravel()), float(R0),
                           int(nsamples), 0)
    return fret_efficiency(r, float(R0))


def sampleDistanceDistInv(g1, g2, nsamples=1000000):
    """Distance samples drawn from the two clouds."""
    import numpy as np
    p1 = np.ascontiguousarray(g1.points().T, dtype=np.float64)
    p2 = np.ascontiguousarray(g2.points().T, dtype=np.float64)
    out = np.asarray(random_distances(p1, p2, int(nsamples), 0))
    return out[:, 0] if out.ndim == 2 else out


# LabelLib's `_arr` twins take the same arrays either way round, which the
# functions above already do; the names exist so that a script written
# against either spelling runs unchanged.
dyeDensityAV1_arr = dyeDensityAV1
dyeDensityAV3_arr = dyeDensityAV3
minLinkerLength_arr = minLinkerLength
addWeights_arr = addWeights
meanDistance_arr = meanDistance
meanEfficiency_arr = meanEfficiency
sampleDistanceDistInv_arr = sampleDistanceDistInv

for _n in ("Grid3D", "dyeDensityAV1", "dyeDensityAV3", "minLinkerLength",
           "addWeights", "meanDistance", "meanEfficiency",
           "sampleDistanceDistInv", "dyeDensityAV1_arr", "dyeDensityAV3_arr",
           "minLinkerLength_arr", "addWeights_arr", "meanDistance_arr",
           "meanEfficiency_arr", "sampleDistanceDistInv_arr"):
    setattr(_labellib, _n, locals()[_n])
del _n

labellib = _labellib
_ll_sys.modules[__name__ + ".labellib"] = _labellib
%}
