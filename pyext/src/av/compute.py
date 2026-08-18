"""Standalone AV computation and data structures.

Provides the :func:`compute_av` function that computes an accessible
volume from raw atomic coordinates, without requiring a full IMP or
chisurf setup.
"""

from __future__ import annotations

import threading
from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Backend detection
# ---------------------------------------------------------------------------

_HAS_IMP_BFF = False
try:
    import IMP
    import IMP.algebra
    import IMP.atom
    import IMP.bff
    import IMP.core
    # Verify that the AV decorator is exposed in this build
    if hasattr(IMP.bff, "AV"):
        _HAS_IMP_BFF = True
except ImportError:
    pass

#: Radius around the attachment atom inside which obstacles are ignored, so the
#: linker can leave the atom it is tied to.
#:
#: Below roughly 2 Å the attachment atom's own neighbours block every starting
#: voxel and IMP's path search returns an **empty volume without raising** — the
#: reference test in this repository uses 2.0 and expects zero points at 1.0.
#: 2.1 Å is carbon's CHARMM radius, so it sits just clear of that cliff and
#: matches the atom a linker is normally tied to.
DEFAULT_ALLOWED_SPHERE_RADIUS = 2.1

#: Serialises IMP object *construction*, not IMP computation.
#:
#: IMP builds its decorators through SWIG, which is not safe from several
#: threads at once. A caller running labelling sites in a thread pool — which is
#: the obvious way to use this function — otherwise dies inside
#: ``XYZR_setup_particle``, receiving a particle *name* where it expects a
#: particle, because two threads are part-way through construction in the same
#: interpreter.
#:
#: The lock deliberately does **not** cover ``resample()``. Measured on 148l E36
#: at 0.5 Å: resample is 636 ms of a 662 ms build (96 %) while all the SWIG
#: construction is 22 ms (3 %). Serialising the compute would cost a threaded
#: caller almost all of its parallelism in order to guard the phase that is not
#: the one failing. Each thread resamples its own model.
_IMP_BUILD_LOCK = threading.Lock()


def _path_map_header(path_map):
    """The PathMapHeader, tolerating an incompletely wrapped build.

    IMP 2.24 binds ``PathMap.get_path_map_header`` to a C symbol
    (``DensityMap_get_path_map_header``) that is absent from the compiled
    ``_IMP_bff``, so calling it raises ``AttributeError``. The writable accessor
    is bound to a different symbol and works, and the header is only read here.
    """
    try:
        return path_map.get_path_map_header()
    except AttributeError:
        return path_map.get_path_map_header_writable()


def _path_map_accessible_density(path_map):
    """``(density, edge_length, origin)`` for a path map's accessible volume.

    Two things here are easy to get wrong and silent when wrong:

    **Axis order.** IMP numbers voxels with *x* fastest
    (``i = x + nx*y + nx*ny*z``), so reshaping the flat tile values C-order into
    ``(nx, ny, nz)`` — which makes the *last* axis fastest — transposes the
    volume, exchanging x and z. An accessible volume is globular enough that the
    transpose still overlaps the truth by ~83 % of its voxels, so the point
    count, the bounding box and the total volume all look right while the shape
    is mirrored.

    **Origin.** IMP sizes the map itself, and the attachment atom lands on the
    middle voxel only when the edge length is odd. The corner must be read from
    the header rather than reconstructed from the attachment point.
    """
    header = path_map.get_header()
    nx, ny, nz = header.get_nx(), header.get_ny(), header.get_nz()
    if not (nx == ny == nz):
        raise ValueError(
            f"IMP PathMap grid is not cubic ({nx}, {ny}, {nz}); "
            "this backend assumes equal nx, ny, nz."
        )
    tile_values = path_map.get_tile_values(
        IMP.bff.PM_TILE_ACCESSIBLE_DENSITY,
        (0.0, _path_map_header(path_map).get_max_path_length()),
    )
    density = (
        np.asarray(tile_values, dtype=np.float64)
        .reshape((nz, ny, nx), order="C")
        .transpose(2, 1, 0)
    )
    density = np.ascontiguousarray(np.where(density > 0.0, 1.0, 0.0))
    origin = np.array(
        [header.get_xorigin(), header.get_yorigin(), header.get_zorigin()],
        dtype=np.float64,
    )
    return density, int(nx), origin


def _particle_index(particle):
    """The index of *particle*, whichever spelling this IMP uses.

    ``Particle.get_particle_index()`` exists in some IMP builds and
    ``get_index()`` in others; calling the wrong one is an ``AttributeError``
    at AV-construction time, which is the last place anyone looks.
    """
    for name in ("get_particle_index", "get_index"):
        getter = getattr(particle, name, None)
        if getter is not None:
            return getter()
    raise AttributeError(
        f"{type(particle).__name__} exposes neither get_particle_index() nor "
        "get_index(); cannot address it in the model."
    )


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class AccessibleVolume:
    """Result of an AV computation.

    Attributes
    ----------
    points : (N, 4) float64
        Point cloud ``(x, y, z, weight)``.
    density : (nx, ny, nz) float64
        3-D voxel density.
    grid_origin : (3,) float64
        Coordinate of the first voxel centre.
    grid_step : float
        Voxel spacing (Å).
    grid_shape : tuple[int, int, int]
        Voxel grid dimensions ``(nx, ny, nz)``.
    attachment_point : (3,) float64
        Coordinates of the attachment site.
    position_name : str
        Optional human-readable label.
    params : dict
        Computation parameters used.
    """

    points: np.ndarray
    density: np.ndarray
    grid_origin: np.ndarray
    grid_step: float
    grid_shape: Tuple[int, int, int]
    attachment_point: np.ndarray
    position_name: str = ""
    params: dict = field(default_factory=dict)

    @property
    def n_points(self) -> int:
        """Number of points in the point cloud."""
        return self.points.shape[0] if self.points.ndim == 2 else 0

    @property
    def mean_position(self) -> np.ndarray:
        """Density-weighted mean position."""
        if self.n_points == 0:
            return self.attachment_point.copy()
        w = self.points[:, 3]
        ws = w.sum()
        if ws == 0:
            return self.attachment_point.copy()
        return np.average(self.points[:, :3], axis=0, weights=w)

    @property
    def has_volume(self) -> bool:
        """Whether the AV contains any points."""
        return self.n_points > 0


# ---------------------------------------------------------------------------
# IMP.bff backend
# ---------------------------------------------------------------------------

def _av_imp_bff(
    atoms_xyz: np.ndarray,
    atoms_vdw: np.ndarray,
    source_xyz: np.ndarray,
    linker_length: float = 20.0,
    linker_width: float = 0.5,
    dye_radii: Tuple[float, float, float] = (3.5, 0.0, 0.0),
    grid_resolution: float = 1.5,
    allowed_sphere_radius: float = DEFAULT_ALLOWED_SPHERE_RADIUS,
) -> AccessibleVolume:
    """Compute AV using IMP.bff.

    Parameters
    ----------
    atoms_xyz : (N, 3) float64
        Atomic coordinates.
    atoms_vdw : (N,) float64
        Van der Waals radii.
    source_xyz : (3,) float64
        Attachment site coordinates.
    linker_length, linker_width, dye_radii, grid_resolution
        See :func:`compute_av`.
    allowed_sphere_radius : float
        See :func:`compute_av`.

    Returns
    -------
    AccessibleVolume

    Notes
    -----
    This backend did not work at all before 2026-07-27, and could not have:
    it read its density by handing the *source* particle to an
    ``IMP.em.SampledDensityMap``, which raises ``UsageException`` because that
    particle carries no ``IMP.atom.Mass`` — and would have sampled a Gaussian
    blob around the attachment point rather than the accessible volume even
    with one. It also never passed ``allowed_sphere_radius``, leaving it below
    the value at which IMP's path search can start. Nothing noticed because no
    test exercised this branch.

    It now reads ``PM_TILE_ACCESSIBLE_DENSITY`` off the path map, which is the
    quantity being asked for.
    """
    if not _HAS_IMP_BFF:
        raise ImportError("IMP.bff is required for the IMP.bff AV backend.")

    source = np.asarray(source_xyz, dtype=np.float64)
    coords = np.asarray(atoms_xyz, dtype=np.float64)
    radii = np.asarray(atoms_vdw, dtype=np.float64)

    # Everything that *constructs* SWIG objects is serialised; `resample()` is
    # not. See `_IMP_BUILD_LOCK` — the split is what makes a threaded caller
    # worth having, since resample is ~96 % of the wall clock.
    with _IMP_BUILD_LOCK:
        model = IMP.Model()
        root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model))

        # Every atom is an obstacle, carrying a radius. The attachment site is
        # one of them rather than a separate massless marker: the source has to
        # be a real XYZR particle in the same hierarchy, and
        # `allowed_sphere_radius` handles that it is also an obstacle to itself.
        particles = []
        for i in range(len(coords)):
            ap = IMP.Particle(model)
            xyzr = IMP.core.XYZR.setup_particle(ap)
            xyzr.set_coordinates(IMP.algebra.Vector3D(*coords[i]))
            xyzr.set_radius(float(radii[i]))
            ah = IMP.atom.Hierarchy.setup_particle(ap)
            ah.set_name(f"atom_{i}")
            root.add_child(ah)
            particles.append(ap)

        coincident = np.flatnonzero(np.all(np.isclose(coords, source), axis=1))
        if coincident.size:
            source_particle = particles[int(coincident[0])]
        else:
            # The attachment site is not one of the atoms; add it as a
            # zero-radius particle so it is still a real XYZR to anchor to.
            source_particle = IMP.Particle(model)
            xyzr = IMP.core.XYZR.setup_particle(source_particle)
            xyzr.set_coordinates(IMP.algebra.Vector3D(*source))
            xyzr.set_radius(0.0)
            ah = IMP.atom.Hierarchy.setup_particle(source_particle)
            ah.set_name("source")
            root.add_child(ah)

        # The AV is decorated onto its **own** particle, with the source passed
        # separately. Setting it up on the source particle itself — which this
        # did — leaves the resampled map at the coordinate origin: header origin
        # (0, 0, 0), obstacles nowhere near the search region, and most of the
        # grid reported accessible (259 788 points against 67 889, centred on
        # nothing).
        rad = IMP.algebra.Vector3D(dye_radii[0], dye_radii[1], dye_radii[2])
        av_particle = IMP.Particle(model)
        IMP.bff.AV.do_setup_particle(
            model,
            av_particle,
            source_particle,
            linker_length=linker_length,
            radii=rad,
            linker_width=linker_width,
            # Without this the search starts inside the attachment atom's own
            # neighbourhood and returns an empty volume, reporting nothing.
            allowed_sphere_radius=allowed_sphere_radius,
            contact_volume_thickness=0.0,
            contact_volume_trapped_fraction=-1,
            simulation_grid_resolution=grid_resolution,
        )
        av_decorator = IMP.bff.AV(model, av_particle)

    # The C++ path search, on this thread's own model. Unlocked on purpose.
    av_decorator.resample()

    with _IMP_BUILD_LOCK:
        path_map = av_decorator.get_map()
        density, nx, origin = _path_map_accessible_density(path_map)
        # The point cloud comes from IMP directly rather than from the grid
        # above. That is deliberate: it makes the cloud and the density
        # **independent** readings of the same volume, so a test comparing them
        # can catch a mis-indexed grid. Deriving the points from the density
        # instead would make any indexing error self-consistent, and therefore
        # invisible.
        xyz_density = np.asarray(path_map.get_xyz_density(), dtype=np.float64)

    shape = (nx, nx, nx)
    voxel_size = float(grid_resolution)

    if xyz_density.size:
        pts = np.ascontiguousarray(
            xyz_density.reshape((-1, xyz_density.shape[-1]))[:, :4]
        )
    else:
        pts = np.zeros((0, 4), dtype=np.float64)

    return AccessibleVolume(
        points=pts,
        density=density,
        grid_origin=origin,
        grid_step=voxel_size,
        grid_shape=shape,
        attachment_point=source_xyz,
        params={
            "backend": "imp_bff",
            "linker_length": linker_length,
            "linker_width": linker_width,
            "dye_radii": list(dye_radii),
            "grid_resolution": grid_resolution,
            "allowed_sphere_radius": allowed_sphere_radius,
        },
    )


def compute_av(
    atoms_xyz: np.ndarray,
    atoms_vdw: np.ndarray,
    source_xyz: np.ndarray,
    linker_length: float = 20.0,
    linker_width: float = 0.5,
    dye_radii: Tuple[float, float, float] = (3.5, 0.0, 0.0),
    grid_resolution: float = 1.5,
    allowed_sphere_radius: float = DEFAULT_ALLOWED_SPHERE_RADIUS,
) -> AccessibleVolume:
    """Compute an accessible volume from raw atomic coordinates.

    The **array** front door: no PDB file and no fps position definition, which
    is what ``restraints`` needs. :func:`IMP.bff.compute_av` is the other one,
    taking a structure plus an fps position; both drive the same
    ``IMP.bff.AV`` / ``PathMap`` core.

    A ``backend`` argument used to select between this and LabelLib. It is gone
    (PRD-112 stage 1): **IMP.bff's AV is the only backend**, decided 2026-08-11
    and applied to ChiSurf by PRD-97 and to QuEst by PRD-109. This module was
    the last place still shipping the alternative.

    Parameters
    ----------
    atoms_xyz : (N, 3) float64
        Cartesian coordinates of all atoms in the structure.
    atoms_vdw : (N,) float64
        Van der Waals radii of all atoms (Å).
    source_xyz : (3,) float64
        Coordinates of the attachment site (Å).
    linker_length : float
        Length of the dye linker (Å).  Default 20.0.
    linker_width : float
        Width of the linker (Å).  Default 0.5.
    dye_radii : tuple of float
        Dye radii ``(r1, r2, r3)`` (Å).  Default ``(3.5, 0.0, 0.0)``
        corresponds to the AV1 (single-sphere) model.
    grid_resolution : float
        Voxel grid spacing (Å).  Default 1.5.
    allowed_sphere_radius : float
        Radius around the attachment atom inside which obstacles are ignored,
        so the linker can leave the atom it is tied to (Å). Default 2.1.

    Returns
    -------
    AccessibleVolume
        The computed dye distribution.

    Raises
    ------
    ImportError
        If this build does not expose IMP.bff's AV decorator.

    Examples
    --------
    >>> import numpy as np
    >>> from IMP.bff.av import compute_av
    >>> # Toy system: 2 atoms + attachment site
    >>> atoms_xyz = np.array([[0.,0.,0.],[3.,0.,0.]])
    >>> atoms_vdw = np.array([1.5, 1.5])
    >>> source_xyz = np.array([0.,0.,0.])
    >>> try:
    ...     av = compute_av(atoms_xyz, atoms_vdw, source_xyz)
    ...     print(av.n_points)
    ... except ImportError:
    ...     print("IMP.bff AV not available in this build")
    """
    if not _HAS_IMP_BFF:
        raise ImportError(
            "IMP.bff's AV decorator is not available in this build; it is the "
            "only accessible-volume backend."
        )
    return _av_imp_bff(
        atoms_xyz, atoms_vdw, source_xyz,
        linker_length, linker_width, dye_radii, grid_resolution,
        allowed_sphere_radius,
    )
