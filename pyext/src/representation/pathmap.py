"""The one place ``IMP.bff.AV`` is driven and its path map is read.

Two builders grew independently — ``IMP.bff.av.compute`` from imp-tricks (raw
arrays in) and ``IMP.bff.fret.av`` from ChiSurf (a structure plus an fps
position). Both set up the same C++ decorator, call the same ``resample()`` and
read the same ``PathMap``, and both had to learn the same lessons; each learned
some of them.

The array path had found that IMP numbers voxels with *x* fastest, so a C-order
reshape into ``(nx, ny, nz)`` returns the volume **transposed**. The structure
path had not, and shipped a mirrored density for as long as it existed — 71 % of
its own point cloud landed on an occupied voxel instead of 100 %, and the mean
donor lifetime on T4L A132 came out 4.9 % wrong (PRD-113 stage 3a).

The structure path had found that the AV must be decorated onto its **own**
particle with the source passed separately, and that the build must be
serialised while ``resample()`` need not be. The array path had found the same
independently, and documents it.

Sharing this function is what stops that from happening a third time. What it
deliberately does **not** share is the conventions the two doors apply to the
result — see :class:`PathMapReading`.
"""

from __future__ import annotations

import threading
from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import numpy as np

import IMP
import IMP.algebra
import IMP.atom
import IMP.core
import IMP.bff

__all__ = ["PathMapReading", "resample_av", "IMP_BUILD_LOCK"]

#: Everything that *constructs* SWIG objects is serialised; ``resample()`` is
#: not. That split is what makes a threaded caller worth having, since resample
#: is ~96 % of the wall clock.
IMP_BUILD_LOCK = threading.RLock()


@dataclass(frozen=True)
class PathMapReading:
    """The raw readings of one resampled path map.

    Deliberately *raw*: the two front doors disagree about what to do with these
    and this stage does not adjudicate.

    * ``density`` carries the accessible-density values. The array door
      binarises them to 0/1 and returns float64; the structure door keeps the
      values as float32.
    * ``points_xyzw`` is IMP's ``get_xyz_density()`` output, ``(N, 4)``. The
      array door keeps the fourth column as a weight; the structure door
      overwrites it with ones.
    * ``source_xyz`` is ``AV.get_source_coordinates()``. The array door instead
      returns the attachment coordinate it was handed.

    Those differences are **behaviour**, and PRD-113 stage 3 forbids behaviour
    changes — the gate is that each door reproduces its own output byte for
    byte. Choosing between them is a deliberate decision for a later stage, and
    it is recorded here rather than papered over.
    """

    density: np.ndarray
    points_xyzw: np.ndarray
    origin: np.ndarray
    step: float
    shape: Tuple[int, int, int]
    source_xyz: np.ndarray


def _accessible_density(path_map, nx: int, ny: int, nz: int) -> np.ndarray:
    """Tile values as an ``(nx, ny, nz)`` array, in the *coordinate* axis order.

    IMP orders the flat values ``i = x + nx*y + nx*ny*z`` — *x* fastest — so a
    C-order reshape into ``(nx, ny, nz)``, which makes the last axis fastest,
    returns the volume transposed. A mirrored volume keeps the right voxel
    count, the right bounding box and the right total volume, so only a
    voxel-by-voxel comparison against the point cloud catches it; see
    ``test/fret/test_av_grid_registration.py``.
    """
    values = path_map.get_tile_values(
        IMP.bff.PM_TILE_ACCESSIBLE_DENSITY,
        (0.0, path_map.get_path_map_header().get_max_path_length()),
    )
    return (
        np.asarray(values, dtype=np.float64)
        .reshape((nz, ny, nx), order="C")
        .transpose(2, 1, 0)
    )


def resample_av(
    model: "IMP.Model",
    source_particle,
    *,
    linker_length: float,
    linker_width: float,
    radii: Sequence[float],
    disc_step: float,
    allowed_sphere_radius: float,
    contact_volume_thickness: float = 0.0,
    contact_volume_trapped_fraction: float = -1.0,
) -> PathMapReading:
    """Decorate a fresh particle as an AV, resample it, and read the map.

    :param model: holds the obstacle hierarchy and the attachment particle.
    :param source_particle: the attachment atom, a real ``XYZR`` in that model.
    :param allowed_sphere_radius: obstacles inside this radius of the attachment
        are ignored, so the linker can leave the atom it is tied to. Without it
        the search starts inside the attachment atom's own neighbourhood and
        returns an empty volume, reporting nothing.

    The AV is decorated onto its **own** particle with the source passed
    separately. Setting it up on the source particle leaves the resampled map at
    the coordinate origin — header origin (0, 0, 0), obstacles nowhere near the
    search region, and most of the grid reported accessible.
    """
    r1, r2, r3 = (list(radii) + [0.0, 0.0, 0.0])[:3]
    with IMP_BUILD_LOCK:
        av_particle = IMP.Particle(model)
        IMP.bff.AV.do_setup_particle(
            model,
            av_particle,
            source_particle,
            linker_length=float(linker_length),
            linker_width=float(linker_width),
            radii=IMP.algebra.Vector3D(float(r1), float(r2), float(r3)),
            allowed_sphere_radius=float(allowed_sphere_radius),
            contact_volume_thickness=float(contact_volume_thickness),
            contact_volume_trapped_fraction=float(contact_volume_trapped_fraction),
            simulation_grid_resolution=float(disc_step),
        )
        av = IMP.bff.AV(model, av_particle)

    # The C++ path search, on this thread's own model. Unlocked on purpose.
    av.resample()

    with IMP_BUILD_LOCK:
        path_map = av.get_map()
        header = path_map.get_header()
        nx, ny, nz = header.get_nx(), header.get_ny(), header.get_nz()
        density = _accessible_density(path_map, nx, ny, nz)
        # The point cloud comes from IMP directly rather than from the grid
        # above. That is deliberate: it makes the cloud and the density
        # **independent** readings of the same volume, so a test comparing them
        # can catch a mis-indexed grid. Deriving the points from the density
        # would make any indexing error self-consistent, and invisible -- which
        # is exactly how the transposed density survived in the structure door.
        xyz_density = np.asarray(path_map.get_xyz_density(), dtype=np.float64)
        if xyz_density.size:
            points = np.ascontiguousarray(
                xyz_density.reshape((-1, xyz_density.shape[-1]))[:, :4])
        else:
            points = np.zeros((0, 4), dtype=np.float64)
        origin = np.array(
            [header.get_origin(0), header.get_origin(1), header.get_origin(2)],
            dtype=np.float64)
        source_xyz = np.array(av.get_source_coordinates(), dtype=np.float64)

    return PathMapReading(
        density=density,
        points_xyzw=points,
        origin=origin,
        step=float(header.get_spacing()),
        shape=(nx, ny, nz),
        source_xyz=source_xyz,
    )
