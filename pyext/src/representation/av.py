"""The accessible volume: where a tethered dye can reach.

The implicit representation — a region the dye may occupy, uniform on that
region, with a mobility field over it — against the rotamer library's discrete
states and the coarse-grained dye's explicit coordinates. All three are
representations, which is why they are siblings here rather than three
top-level packages.

Two front doors onto the same C++ core, and **both are needed**:

* :func:`~IMP.bff.representation.av.compute.compute_av` takes **arrays** —
  coordinates and radii — and knows nothing about files. This is the entry point
  for a caller that already has a structure in memory.
* :mod:`~IMP.bff.representation.av.structure` takes a **structure and an fps
  position**, which is what a labelling file describes.

They grew independently (imp-tricks and ChiSurf), and each carries fixes the
other lacks: the array door has the AV-on-its-own-particle handling and the
build/resample lock split, the structure door has the ``disc_step`` trap and
source clearance. PRD-113 stage 3 set out to merge them and did not; what it did
do is put them side by side, which is the precondition.

:class:`~IMP.bff.representation.av.basic.BasicAV` is the point cloud with its
distance methods, and :class:`~IMP.bff.representation.av.acv.ACV` splits it into
contact and free volumes for a dye that sticks.
"""

from __future__ import annotations

from dataclasses import field
from functools import lru_cache
from typing import Dict, Optional, Tuple
import json
import logging
import os
import threading

import numpy as np

import IMP.bff.io.structure as _structure
from IMP.bff.label import default_strip_mask, strip_pdb_lines
# (was: import _kernels) -- now in this module
from IMP.bff.representation.pathmap import resample_av
from IMP.bff.representation.distance import AccessibleVolume
import IMP
import IMP.algebra
import IMP.atom
import IMP.bff
import IMP.core
import IMP.em

__all__ = [
    'ACV',
    'AccessibleVolume',
    'BasicAV',
    'compute_av',
]

# --------------------------------------------------------------------------
# _kernels
# --------------------------------------------------------------------------
r"""Distance and density kernels for accessible volumes, in C++.

These were numba (``@njit(cache=True, nogil=True)``) until PRD-113 stage 3;
numba is a prototyping tool in this package, not a runtime dependency, so the
numerics are now :mod:`IMP.bff`'s compiled ``AVDistance`` and this module is the
array adapter that reshapes for it. The signatures are unchanged.

The kernels take **point arrays**, not the ``AV`` decorators the ``av_distance``
family takes, so they serve a rotamer library, a coarse-grained ensemble or an
MD trajectory as readily as an accessible volume.

One behavioural note that no rewrite could avoid: :func:`random_distances` draws
from a different generator than numba's, so *individual samples differ*. Every
quantity built on it -- :func:`average_distance`, :func:`mean_fret_distance`,
``pRDA`` -- is a Monte-Carlo estimator and agrees with the numba version only to
the sampling error, about :math:`1/\sqrt{n}`. Tests on these must compare
distributions, never recorded numbers.
"""

def random_distances(
    p1: np.ndarray,
    p2: np.ndarray,
    n_samples: int,
    seed: int = 0,
) -> np.ndarray:
    """Draw random distance-weight pairs from two AV point clouds.

    Each sample picks one random point from each AV and computes the
    Euclidean distance and the product of their density weights. The two
    clouds are drawn independently, so the pairs sample the joint distribution.

    Parameters
    ----------
    p1 : (n1, 4) float64
        Points of the first AV (x, y, z, weight).
    p2 : (n2, 4) float64
        Points of the second AV (x, y, z, weight).
    n_samples : int
        Number of random samples to draw.
    seed : int
        Seed for reproducibility. Reproducible run to run; **not** the same
        stream as the numba version this replaces.

    Returns
    -------
    (n_samples, 2) float64
        Column 0: Euclidean distances (Å).
        Column 1: weight products.
    """
    out = IMP.bff.random_distances(
        np.ascontiguousarray(p1, dtype=np.float64).ravel(),
        np.ascontiguousarray(p2, dtype=np.float64).ravel(),
        int(n_samples), int(seed),
    )
    return np.asarray(out, dtype=np.float64).reshape(int(n_samples), 2)


def density2points(
    nx: int,
    ny: int,
    nz: int,
    dg: float,
    density: np.ndarray,
    r0: np.ndarray,
    threshold: float = 0.0,
) -> tuple[int, np.ndarray]:
    """Convert a 3-D density grid to a point cloud.

    Parameters
    ----------
    nx, ny, nz : int
        Grid dimensions.
    dg : float
        Grid spacing (Å).
    density : (nx, ny, nz) float64
        Density values.
    r0 : (3,) float64
        Grid origin (coordinates of the first voxel centre).
    threshold : float
        Minimum density for a voxel to be included.

    Returns
    -------
    n : int
        Number of retained points.
    points : (n, 4) float64
        Point cloud ``(x, y, z, weight)``. Exactly *n* rows -- the numba
        version returned a grid-sized buffer with unused trailing entries, and
        callers sliced it to *n*, which is still correct.
    """
    flat = IMP.bff.density_to_points(
        np.ascontiguousarray(density, dtype=np.float64).ravel(),
        int(nx), int(ny), int(nz), float(dg),
        np.asarray(r0, dtype=np.float64).ravel(), float(threshold),
    )
    points = np.asarray(flat, dtype=np.float64).reshape(-1, 4)
    return points.shape[0], points


def weighted_mean(points: np.ndarray, n: int) -> np.ndarray:
    """Weighted mean position of a point cloud.

    Parameters
    ----------
    points : (n_points, 4) float64
        Point cloud ``(x, y, z, weight)``.
    n : int
        Number of valid points.

    Returns
    -------
    (3,) float64
        Weighted mean coordinate; the origin for an empty cloud.
    """
    if points is None or n == 0:
        return np.zeros(3, dtype=np.float64)
    m = IMP.bff.points_weighted_mean(
        np.ascontiguousarray(points[:n], dtype=np.float64).ravel())
    return np.asarray(m, dtype=np.float64)


def average_distance(points1: np.ndarray, n1: int,
                     points2: np.ndarray, n2: int,
                     n_samples: int = 50000) -> float:
    """Weighted mean inter-point distance between two AVs.

    Parameters
    ----------
    points1 : (n_points, 4) float64
    n1 : int
        Number of valid rows in *points1*.
    points2 : (n_points, 4) float64
    n2 : int
        Number of valid rows in *points2*.
    n_samples : int
        Number of random pairs to sample.

    Returns
    -------
    float
        Mean distance :math:`\\langle R_{DA}\\rangle` (Å), to the sampling error.
    """
    return IMP.bff.average_distance(
        np.ascontiguousarray(points1[:n1], dtype=np.float64).ravel(),
        np.ascontiguousarray(points2[:n2], dtype=np.float64).ravel(),
        int(n_samples), 0)


def mean_fret_distance(points1: np.ndarray, n1: int,
                       points2: np.ndarray, n2: int,
                       forster_radius: float = 52.0,
                       n_samples: int = 50000) -> float:
    """FRET-averaged distance :math:`R_E` between two AVs.

    The *efficiency* is averaged and converted back, which is not the same as
    averaging the distance: :math:`1/r^6` weights close pairs far more heavily,
    so this is always the shorter of the two.

    Parameters
    ----------
    points1, points2 : AV point clouds.
    n1, n2 : Number of valid rows.
    forster_radius : float
        Förster radius :math:`R_0` (Å).
    n_samples : int

    Returns
    -------
    float
        :math:`R_E` (Å); 0 if the mean efficiency saturates at 1, infinity if
        it reaches 0.
    """
    return IMP.bff.mean_fret_distance(
        np.ascontiguousarray(points1[:n1], dtype=np.float64).ravel(),
        np.ascontiguousarray(points2[:n2], dtype=np.float64).ravel(),
        float(forster_radius), int(n_samples), 0)


def split_av_acv(
    density: np.ndarray,
    dg: float,
    radius: np.ndarray,
    rs: np.ndarray,
    r0: np.ndarray,
) -> tuple[int, int, np.ndarray, np.ndarray]:
    """Split an AV density grid into contact (slow) and non-contact parts.

    Grid points within *radius* of any slow center *rs* are assigned to
    the contact volume; all others to the non-contact volume.

    Parameters
    ----------
    density : (ng, ng, ng) float64
        Density grid.
    dg : float
        Grid spacing (Å).
    radius : (n,) float64
        Slow radius per center. A length-1 array broadcasts to all centres.
    rs : (n, 3) float64
        Slow-center coordinates (Å).
    r0 : (3,) float64
        **The grid anchor**: the position of voxel ``(ng - 1) // 2`` on each
        axis, not the corner and not the attachment atom. Voxel ``i`` sits at
        ``r0 + (i - (ng - 1) // 2) * dg``, and the inverse used here is
        ``floor((p - r0) / dg) + (ng - 1) // 2``.

        This said "grid origin (attachment-site coordinates)", which is
        ambiguous between three different points, and the code used the float
        corner ``(ng - 1) / 2`` with ``int()`` truncation. Both were wrong, and
        both were fixed on 2026-07-28 after the identical pair of defects was
        found and measured in QuEst:

        * the float corner differs from the integer offset on every **even**
          ``ng`` -- the normal case, 80 and 86 on the reference sites -- which
          puts the stamped spheres half a voxel from the density they mask;
        * ``int()`` truncates **toward zero**, so a centre on the negative side
          of ``r0`` rounds *up* while every other index map rounds down: one
          voxel per axis for half the grid, whatever the parity.

        Registration check to repeat if this is ever changed: every occupied
        voxel of ``density``, converted to Angstrom and back through this map,
        must return to itself.

    Returns
    -------
    n_contact : int
        Number of contact-voxel candidates.
    n_non_contact : int
        Number of non-contact voxel candidates.
    contact_mask : (ng, ng, ng) uint8
        Binary mask for contact voxels (1 where contact).
    non_contact_mask : (ng, ng, ng) uint8
        Binary mask for non-contact voxels (1 where non-contact).

        ``uint8``, not ``float64``: these are binary masks, and at ``ng = 92``
        a float64 pair costs 12.5 MB per labelling site against 1.6 MB. A
        residue scan builds one per site.
    """
    ng = int(density.shape[0])
    rs = np.ascontiguousarray(rs, dtype=np.float64).reshape(-1, 3)
    radius = np.asarray(radius, dtype=np.float64).ravel()
    if radius.size != rs.shape[0]:
        radius = np.full(rs.shape[0], radius[0], dtype=np.float64)

    label = np.asarray(IMP.bff.split_contact_volume(
        np.ascontiguousarray(density, dtype=np.float64).ravel(),
        ng, float(dg), radius, rs.ravel(),
        np.asarray(r0, dtype=np.float64).ravel(),
    ), dtype=np.int32).reshape(ng, ng, ng)

    contact = (label == IMP.bff.AV_VOXEL_CONTACT).astype(np.uint8)
    non_contact = (label == IMP.bff.AV_VOXEL_FREE).astype(np.uint8)
    return int(contact.sum()), int(non_contact.sum()), contact, non_contact


# --------------------------------------------------------------------------
# basic
# --------------------------------------------------------------------------
"""BasicAV — accessible volume for a single dye labelling site.

BasicAV stores the accessible volume as a point cloud and provides
distance calculations (R_DA, R_mp, R_E, pRDA) and I/O methods.
"""

BACKENDS_AVAILABLE = False
try:
    import IMP
    import IMP.bff
    BACKENDS_AVAILABLE = True
except ImportError:
    pass


# ``BasicAV`` and ``ACV`` are **C++**.
#
# They were 362 lines here holding a point cloud and a density grid as numpy
# arrays and calling C++ for every actual computation -- the arithmetic was
# already across the boundary, the *object* was not, and that is what keeps a
# module Python-carried rather than C++-carried.
#
# The surface is unchanged: keyword construction, ``points`` as ``(n, 4)``,
# ``density`` as ``(ng, ng, ng)`` or ``None``, ``mean_position``, ``n_points``,
# ``dRmp`` / ``dRDA`` / ``dRDAE`` / ``pRDA``, and ``ACV.from_basic_av``. Gated
# against the Python it replaces: points and all three distances bit-identical,
# the ACV density to 1e-14, ``pRDA`` to 2e-15 (a different accumulation order in
# the histogram).
#
# One thing did *not* survive, deliberately. The Python kept ``_padded_points``
# alongside ``n_points`` because the numba kernel it once called returned a
# grid-sized buffer with unused trailing rows. ``density_to_points`` has
# returned exactly the points it found since it started publishing a numpy
# view, so the two could no longer disagree and the count is simply the cloud's
# length.
#
# See ``include/IMP/bff/AVModel.h`` and ``pyext/IMP_bff.avmodel.i``.
from IMP.bff import ACV, BasicAV  # noqa: F401


# --------------------------------------------------------------------------
# compute
# --------------------------------------------------------------------------
"""Standalone AV computation and data structures.

Provides the :func:`compute_av` function that computes an accessible
volume from raw atomic coordinates, without requiring a full IMP or
chisurf setup.
"""

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


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

def _av_from_arrays(
    atoms_xyz: np.ndarray,
    atoms_vdw: np.ndarray,
    source_xyz: np.ndarray,
    linker_length: float = 20.0,
    linker_width: float = 0.5,
    dye_radii: Tuple[float, float, float] = (3.5, 0.0, 0.0),
    grid_resolution: float = 1.5,
    allowed_sphere_radius: float = DEFAULT_ALLOWED_SPHERE_RADIUS,
    search_stencil: int = 0,
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

    reading = resample_av(
        model, source_particle,
        linker_length=linker_length, linker_width=linker_width,
        radii=dye_radii, disc_step=grid_resolution,
        allowed_sphere_radius=allowed_sphere_radius,
        search_stencil=search_stencil,
    )
    nx = reading.shape[0]
    origin = reading.origin
    xyz_density = reading.points_xyzw
    # This door's conventions: float64, and the point weights are whatever IMP
    # reported. The structure door in `IMP.bff.representation.av.structure`
    # keeps the raw float32 values and forces the weights to one.
    #
    # The density used to be **binarised** here, `np.where(d > 0, 1, 0)`. For
    # AV1 that is a no-op -- the carve already emits exactly 0 or 1 -- but it
    # destroys AV3, whose whole content is the grading: the density is the
    # fraction of the three probe radii that fit, {1/3, 2/3, 1}, matching
    # LabelLib's `excludeConcentricSpheres`. Binarising it returns the AV1
    # volume of the *smallest* radius and silently discards the other two, which
    # is the same defect as `get_radii()` dropping `radius2`: an AV3 that looks
    # like it worked. The raw values are kept, so AV1 is unchanged and AV3 means
    # something.
    density = np.ascontiguousarray(np.asarray(reading.density, dtype=float))

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
            "search_stencil": search_stencil or 74,
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
    search_stencil: int = 0,
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
    search_stencil : int
        Dijkstra neighbour stencil. 0 (default) keeps the AV's own default of
        **74**, the LabelLib reference metric. 26 is the speed option: ~1.9×
        faster for ~20 % less volume. 30 is the historical variant.

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
    >>> from IMP.bff.representation.av import compute_av
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
    return _av_from_arrays(
        atoms_xyz, atoms_vdw, source_xyz,
        linker_length, linker_width, dye_radii, grid_resolution,
        allowed_sphere_radius, search_stencil,
    )


# --------------------------------------------------------------------------
# structure
# --------------------------------------------------------------------------
"""Accessible-volume computation on IMP.bff's own AV machinery.

Self-contained on ``IMP.bff.AV``/``PathMap`` — there is no alternative backend
and no optional dependency. (The pre-move ChiSurf implementation carried a
LabelLib backend and selection logic; that was dropped deliberately when the
module moved here: nothing in the IMP stack depends on LabelLib.)
"""

logger = logging.getLogger(__name__)



# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------




def _active_backend_name() -> str:
    """Return the active backend name for diagnostics."""
    return "imp-bff"


# ---------------------------------------------------------------------------
# AV computation
# ---------------------------------------------------------------------------

def _av_from_structure(
    pdb_path: str,
    source_info: Dict,
    linker_length: float,
    linker_width: float,
    radii: Tuple[float, float, float],
    disc_step: float,
) -> AccessibleVolume:
    """Compute an AV with ``IMP.bff.AV``."""
    chain = source_info.get("chain_identifier", "")
    resseq = source_info.get("residue_seq_number", 0)
    aname = source_info.get("atom_name", "CA")

    # The FPS strip: positions are calibrated for a structure whose
    # attachment residue does not wall in its own dye, so the side chain
    # goes (minus the attachment atom) before anything is measured. A
    # declared ``strip_mask`` outside the fps dialect raises here -- loud,
    # because the alternative is computing against obstacles the document
    # said to remove.
    strip_mask = str(source_info.get("strip_mask") or "").strip()
    pdb_path = _stripped_pdb_for(pdb_path, chain, resseq, aname, strip_mask)

    model = IMP.Model()
    hierarchy = IMP.atom.read_pdb(pdb_path, model, IMP.atom.NonWaterPDBSelector())

    sel = IMP.atom.Selection(hierarchy)
    if chain:
        sel.set_chain_id(chain)
    sel.set_residue_index(int(resseq))
    sel.set_atom_type(IMP.atom.AtomType(str(aname)))
    particles = sel.get_selected_particles()
    if not particles:
        raise ValueError(f"Attachment site {chain}:{resseq}:{aname} not found")
    attachment_particle = particles[0]

    # Source clearance. The path search inflates obstacles by half the linker
    # width, so the free sphere around the attachment atom has to clear that
    # inflation (plus a grid step of slack) or the source tile is walled in
    # and the AV comes back empty. The strip above already removes the
    # attachment residue's side chain, which is what lets FPS-calibrated
    # small clearances (``allowed_sphere_radius: 1``) compute a real cloud.
    # A position that declares allowed_sphere_radius (schema field) keeps
    # its own value -- honoured exactly, never escalated: an empty AV at
    # the declared parameters is the answer, not a signal to retry at
    # invented ones.
    default_clearance = max(1.5, 0.5 * linker_width + 0.5 * disc_step)
    allowed_sphere_radius = float(
        source_info.get("allowed_sphere_radius", default_clearance))

    reading = resample_av(
        model, attachment_particle,
        linker_length=linker_length, linker_width=linker_width, radii=radii,
        disc_step=disc_step, allowed_sphere_radius=allowed_sphere_radius,
        contact_volume_thickness=float(
            source_info.get("contact_volume_thickness", 0.0)),
        contact_volume_trapped_fraction=float(
            source_info.get("contact_volume_trapped_fraction", -1)),
    )
    nx, ny, nz = reading.shape
    # This door's conventions, deliberately kept as they were: float32 density
    # values (not binarised) and uniform point weights. The array door in
    # `IMP.bff.representation.av.compute` chooses differently on both counts; reconciling them
    # is a behaviour change and belongs to a later stage, not to this move.
    density = np.ascontiguousarray(reading.density, dtype=np.float32)
    points = (
        np.column_stack([reading.points_xyzw[:, :3],
                         np.ones(reading.points_xyzw.shape[0])])
        if reading.points_xyzw.size else np.zeros((0, 4), dtype=np.float64))
    origin = reading.origin
    step = reading.step
    att_xyz = reading.source_xyz

    return AccessibleVolume(
        points=points,
        density=density,
        grid_origin=origin,
        grid_step=float(step),
        grid_shape=(nx, ny, nz),
        attachment_point=att_xyz,
    )


def compute_av_from_structure(
    atoms: np.ndarray,
    source_xyz: np.ndarray,
    linker_length: float,
    linker_width: float,
    radii: Tuple[float, float, float],
    disc_step: float = 1.5,
    pdb_path: Optional[str] = None,
    source_info: Optional[Dict] = None,
) -> AccessibleVolume:
    """Compute an accessible volume from a structure and an fps position.

    The **structure** front door, and the other of the pair: :func:`compute_av`
    takes raw arrays, this one takes a PDB plus a position definition in the
    fps dialect. They were both called ``compute_av`` while they sat in
    different modules of ``representation/av/``; merging those modules made one
    shadow the other silently, and since the package re-exported the array one,
    the shadow would have swapped the public function for its sibling.

    Parameters
    ----------
    atoms : (N, 4) float64
        Columns: x, y, z, vdw_radius. Kept for API compatibility — the AV is
        computed from ``pdb_path``, which carries atom identity.
    source_xyz : (3,) float64
        Attachment point coordinates (informational; the attachment atom is
        resolved from ``source_info``).
    pdb_path : str
        Structure the AV is computed on. Required.
    source_info : dict
        Position definition (fps.json position fields). Required.
    """
    if pdb_path is None or source_info is None:
        raise ValueError(
            "compute_av_from_structure requires pdb_path and source_info")
    # `source_info` is a position definition in the fps dialect, and every other
    # AV parameter in it is honoured (`allowed_sphere_radius`, `strip_mask`,
    # `contact_volume_*`). `simulation_grid_resolution` is the exception: it is
    # *written* from `disc_step` below rather than read, so a caller who states
    # it here and leaves `disc_step` at its default silently gets 1.5 A. That is
    # the fps field a reader is most likely to trust, so disagreement is an
    # error rather than a preference -- a resolution is the one AV parameter
    # whose being wrong is invisible in the result.
    declared = source_info.get("simulation_grid_resolution")
    if declared is not None and abs(float(declared) - float(disc_step)) > 1e-9:
        raise ValueError(
            f"source_info declares simulation_grid_resolution={float(declared)} A "
            f"but disc_step={float(disc_step)} A was passed. The AV is built at "
            "disc_step; pass the resolution there, or drop it from source_info.")
    return _av_from_structure(
        pdb_path,
        source_info,
        linker_length,
        linker_width,
        radii,
        disc_step,
    )


def compute_avs_for_structure(
    atoms: np.ndarray,
    positions: Dict,
    pdb_path: str | list[str] | None = None,
    disc_step: Optional[float] = None,
) -> Dict[str, AccessibleVolume]:
    """Compute AVs for all positions in an fps.json ``Positions`` dict.

    Parameters
    ----------
    atoms : (N, 4) float64
        xyzr from :func:`load_structure_with_vdw`. Used as fallback if pdb_path is not given.
    positions : dict
        fps.json Positions section.
    """
    if isinstance(pdb_path, (list, tuple)):
        pdb_paths = list(pdb_path)
    elif isinstance(pdb_path, str) and "," in pdb_path:
        pdb_paths = [p.strip() for p in pdb_path.split(",")]
    elif isinstance(pdb_path, str):
        pdb_paths = [pdb_path]
    else:
        pdb_paths = []

    avs: Dict[str, AccessibleVolume] = {}
    for pname, pdef in positions.items():
        bi = int(pdef.get("body_id", 0))
        curr_pdb = pdb_paths[bi] if bi < len(pdb_paths) else (pdb_paths[0] if pdb_paths else None)

        if curr_pdb is not None:
            curr_atoms = load_structure_with_vdw(curr_pdb)
        else:
            curr_atoms = atoms

        ll = float(pdef.get("linker_length", 20.0))
        lw = float(pdef.get("linker_width", 1.0))
        r1 = float(pdef.get("radius1", 3.5))
        r2 = float(pdef.get("radius2", 0.0))
        r3 = float(pdef.get("radius3", 0.0))
        ds = float(disc_step or pdef.get("simulation_grid_resolution", 1.5))

        chain = pdef.get("chain_identifier", "")
        resseq = pdef.get("residue_seq_number", 0)
        aname = pdef.get("atom_name", "CA")

        source_xyz = _find_attachment_point(curr_atoms, chain, resseq, aname, pdb_path=curr_pdb)
        if source_xyz is None:
            avs[pname] = AccessibleVolume(
                points=np.zeros((0, 4), dtype=np.float64),
                density=np.zeros((1, 1, 1), dtype=np.float32),
                grid_origin=np.zeros(3),
                grid_step=ds,
                grid_shape=(1, 1, 1),
                attachment_point=np.zeros(3),
                position_name=pname,
            )
            continue

        av = compute_av_from_structure(
            atoms=curr_atoms,
            source_xyz=source_xyz,
            linker_length=ll,
            linker_width=lw,
            radii=(r1, r2, r3),
            disc_step=ds,
            pdb_path=curr_pdb,
            source_info=pdef,
        )
        av.position_name = pname
        av.params = pdef
        avs[pname] = av
    return avs


# ---------------------------------------------------------------------------
# Helper: vdW radii
# ---------------------------------------------------------------------------

# From FPS data/vdW.txt (selected common elements)
VDW_RADII = {
    1: 1.20, 2: 1.40, 3: 1.82, 4: 1.53, 5: 1.92, 6: 1.70, 7: 1.55,
    8: 1.52, 9: 1.47, 12: 1.73, 14: 2.10, 15: 1.80, 16: 1.80,
    17: 1.75, 19: 2.27, 20: 1.97, 26: 1.56, 30: 1.39,
}
_DEFAULT_VDW = 1.70
_ELEMENT_NUMBERS = {
    "H": 1,
    "HE": 2,
    "LI": 3,
    "BE": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "MG": 12,
    "SI": 14,
    "P": 15,
    "S": 16,
    "CL": 17,
    "K": 19,
    "CA": 20,
    "FE": 26,
    "ZN": 30,
}


def _element_symbol_from_pdb_line(line: str) -> str:
    """Return an element symbol parsed from a PDB ATOM/HETATM line.

    Columns 77-78 carry the element in a modern PDB and are used verbatim when
    present. Older files (the shipped FPS screening structures are 66-character
    records) leave them empty, so the element is read from the atom-name field
    instead. There the element is *right-justified in columns 13-14*: a blank or
    numeric column 13 means a one-letter element, so ``" CA "`` is an
    α-carbon while ``"CA  "`` is calcium. A two-letter reading is additionally
    rejected when columns 15-16 contain a digit, which is how a four-character
    hydrogen name such as ``"HE21"`` is written — helium would otherwise win.

    Parameters
    ----------
    line : str
        PDB ATOM or HETATM record.

    Returns
    -------
    str
        Uppercase element symbol, or an empty string when it cannot be parsed.
    """
    symbol = line[76:78].strip().upper() if len(line) >= 78 else ""
    if symbol:
        return symbol

    name_field = line[12:16].ljust(4).upper()
    candidate = name_field[:2].strip()
    if (
        len(candidate) == 2
        and candidate in _ELEMENT_NUMBERS
        and not any(ch.isdigit() for ch in name_field[2:])
    ):
        return candidate
    letters = "".join(ch for ch in name_field if ch.isalpha())
    if letters[:1] not in _ELEMENT_NUMBERS and letters[:2] in _ELEMENT_NUMBERS:
        # Left-padding a two-letter element (" ZN ") breaks the column rule, but
        # here the strict reading is not an element at all, so take the pair.
        return letters[:2]
    return letters[:1]


def _pdb_cache_token(pdb_path: str) -> tuple[str, int, int]:
    """Return a cache token that changes when a PDB file changes.

    Parameters
    ----------
    pdb_path : str
        Path to a PDB file.

    Returns
    -------
    tuple
        Absolute path, modification time in ns, and file size.
    """
    path = os.path.abspath(pdb_path)
    stat = os.stat(path)
    return path, stat.st_mtime_ns, stat.st_size


@lru_cache(maxsize=32)
def _load_pdb_records_cached(
    pdb_path: str,
    mtime_ns: int,
    size: int,
) -> tuple[tuple[str, int, str, float, float, float, float], ...]:
    """Load ATOM/HETATM records from a PDB file.

    Parameters
    ----------
    pdb_path : str
        Absolute path to a PDB file.
    mtime_ns : int
        File modification timestamp used as part of the cache key.
    size : int
        File size used as part of the cache key.

    Returns
    -------
    tuple
        Records containing chain, residue number, atom name, xyz, and vdW radius.
    """
    del mtime_ns, size
    rows = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                resseq = int(line[22:26].strip())
            except ValueError:
                continue
            chain = line[21].strip()
            atom_name = line[12:16].strip()
            element = _element_symbol_from_pdb_line(line)
            atomic_number = _ELEMENT_NUMBERS.get(element, 0)
            rows.append((
                chain,
                resseq,
                atom_name,
                xyz[0],
                xyz[1],
                xyz[2],
                VDW_RADII.get(atomic_number, _DEFAULT_VDW),
            ))
    if not rows:
        raise ValueError(f"No ATOM/HETATM coordinates found in '{pdb_path}'")
    return tuple(rows)


def _cached_pdb_records(pdb_path: str) -> tuple[tuple[str, int, str, float, float, float, float], ...]:
    """Return cached PDB records for a path.

    Parameters
    ----------
    pdb_path : str
        Path to a PDB file.

    Returns
    -------
    tuple
        Cached ATOM/HETATM records.
    """
    return _load_pdb_records_cached(*_pdb_cache_token(pdb_path))


def _load_pdb_xyzr_direct(pdb_path: str) -> np.ndarray:
    """Load PDB ATOM/HETATM coordinates and vdW radii without IMP.

    Parameters
    ----------
    pdb_path : str
        Path to a PDB file.

    Returns
    -------
    numpy.ndarray
        ``(N, 4)`` array with ``x, y, z, vdw_radius`` columns.
    """
    records = _cached_pdb_records(pdb_path)
    return np.asarray(
        [(x, y, z, radius) for _, _, _, x, y, z, radius in records],
        dtype=np.float64,
    )


def load_structure_with_vdw(pdb_path: str) -> np.ndarray:
    """Load a PDB and return (N, 4) array: x, y, z, vdw_radius.

    Parses PDB records directly to avoid IMP/CHARMM warnings for unsupported
    HETATM residues, then falls back to IMP.atom if direct parsing fails.
    """
    try:
        return _load_pdb_xyzr_direct(pdb_path)
    except Exception:
        pass

    coords, particles, _model, _hier = _structure.load_structure_with_particles(pdb_path)
    vdw = np.full(coords.shape[0], _DEFAULT_VDW, dtype=np.float64)
    for i, p in enumerate(particles):
        try:
            at = IMP.atom.Atom(p)
            elem = at.get_element()
            vdw[i] = VDW_RADII.get(elem, _DEFAULT_VDW)
        except Exception:
            pass
    return np.column_stack([coords, vdw])


def _find_attachment_point(
    atoms: np.ndarray,
    chain: str,
    resseq: int,
    atom_name: str,
    pdb_path: Optional[str] = None,
) -> Optional[np.ndarray]:
    """Find the coordinates of an attachment atom.

    If pdb_path is provided, the atom is resolved by identity — chain, residue
    number and atom name — and a miss stays a miss: an unresolvable site
    returns ``None`` rather than a positional guess, because a dye attached to
    an unrelated atom yields a plausible and entirely wrong accessible volume.
    Without a PDB file the atoms array carries no identity at all, so the
    residue sequence number is used as a proxy index into it.

    Parameters
    ----------
    atoms : (N, 4) ndarray
        The atoms array (x, y, z, vdw_radius).
    chain : str
        The chain identifier.
    resseq : int
        The residue sequence number.
    atom_name : str
        The attachment atom name (e.g., 'CA', 'CB').
    pdb_path : str, optional
        Path to the PDB file for exact matching.

    Returns
    -------
    ndarray or None
        The (3,) coordinates of the attachment atom, or None if not found.
    """
    if pdb_path and os.path.exists(pdb_path):
        try:
            records = _cached_pdb_records(pdb_path)
        except Exception:
            logger.warning("Could not read atom records from %s", pdb_path, exc_info=True)
            return None
        for line_chain, line_resseq, line_atom_name, x, y, z, _ in records:
            if line_resseq == resseq and line_atom_name == atom_name:
                if not chain or line_chain == chain:
                    return np.array([x, y, z], dtype=np.float64)
        logger.warning(
            "Attachment atom '%s:%s:%s' does not exist in %s",
            chain, resseq, atom_name, pdb_path,
        )
        return None

    return atoms[resseq - 1, :3] if resseq > 0 and resseq <= atoms.shape[0] else None


# ---------------------------------------------------------------------------
# The FPS strip: obstacles removed around the attachment site
# ---------------------------------------------------------------------------

#: The strip grammar and the AV default live in :mod:`IMP.bff.label`
#: (PRD-106): the AV build keeps the backbone plus the attachment atom
#: (``default_strip_mask``), imported above.

#: Stripped-PDB cache, keyed like the record cache plus the site and mask.
#: One file per distinct (structure, site, mask) per process.
_STRIPPED_PDB_CACHE: Dict[tuple, str] = {}


def _stripped_pdb_for(
    pdb_path: str,
    chain: str,
    resseq: int,
    atom_name: str,
    strip_mask: Optional[str] = None,
) -> str:
    """A copy of *pdb_path* with the atoms selected by *strip_mask* removed.

    The strip is the FPS convention: fps.json positions are calibrated for a
    structure whose attachment residue does not wall in its own dye, so the
    default removes the attachment residue's side chain minus the attachment
    atom (the backbone stays). A declared mask in the fps dialect is
    honoured as given; one outside the dialect raises rather than being
    ignored.

    The attachment atom itself is always kept, whatever the mask selects --
    the attachment is resolved by ``(chain, residue, atom name)`` from this
    file.

    Parameters
    ----------
    pdb_path : str
        Structure to strip.
    chain : str
        Attachment chain identifier (empty matches any chain).
    resseq : int
        Attachment residue sequence number.
    atom_name : str
        Attachment atom name.
    strip_mask : str, optional
        Declared fps ``strip_mask``; empty means the default strip.

    Returns
    -------
    str
        Path of the stripped copy (the original path on I/O failure -- a
        cloud computed against the unstripped structure beats no cloud).
    """
    mask = (strip_mask or "").strip()
    key = (*_pdb_cache_token(pdb_path), chain, int(resseq), atom_name, mask)
    cached = _STRIPPED_PDB_CACHE.get(key)
    if cached is not None and os.path.exists(cached):
        return cached
    try:
        with open(pdb_path) as source:
            lines = source.readlines()
    except OSError:
        return str(pdb_path)

    strip = mask or default_strip_mask(chain, int(resseq), atom_name)
    kept = strip_pdb_lines(lines, strip, keep_attachment=(chain, int(resseq), atom_name))

    try:
        import tempfile

        with tempfile.NamedTemporaryFile(
            "w", suffix=".pdb", prefix="bff_av_", delete=False
        ) as handle:
            handle.writelines(kept)
            result = handle.name
    except OSError:
        return str(pdb_path)
    _STRIPPED_PDB_CACHE[key] = result
    return result


# --------------------------------------------------------------------------
# Reporting AV positions
#
# Was in `tools.py`. It reads AVs and prints their mean positions, which is
# this module's subject rather than a path helper's.
# --------------------------------------------------------------------------
def display_mean_av_positions(
        used_avs: typing.List[IMP.bff.AV],
        dye_radius: float = 3.5
):
    """Decorate the accessible volume particles
    with a mass and radii, so that they appear
    in the rmf file"""
    for dye in used_avs:
        # Add radius and mass to dye
        p_dye = dye.get_particle()
        if not IMP.core.XYZR.get_is_setup(p_dye):
            p_dye = IMP.core.XYZR.setup_particle(p_dye)
            p_dye.set_radius(dye_radius)
        else:
            p_dye = IMP.core.XYZR(p_dye)
        if not IMP.atom.Mass.get_is_setup(p_dye):
            IMP.atom.Mass.setup_particle(p_dye, 100)

        # add dye to atom hier
        p_att = dye.get_source()
        h_att = IMP.atom.Hierarchy(p_att)
        IMP.atom.Hierarchy.setup_particle(p_dye)
        h_dye = IMP.atom.Hierarchy(p_dye)
        h_att.add_child(h_dye)
        # create bond between dye and source
        if not IMP.atom.Bonded.get_is_setup(p_dye):
            IMP.atom.Bonded.setup_particle(p_dye)
        if not IMP.atom.Bonded.get_is_setup(p_att):
            IMP.atom.Bonded.setup_particle(p_att)
        if not IMP.atom.get_bond(IMP.atom.Bonded(p_dye), IMP.atom.Bonded(p_att)):
            IMP.atom.create_bond(
                IMP.atom.Bonded(p_dye), IMP.atom.Bonded(p_att), 1)

