"""Distance and density kernels for accessible volumes, in C++.

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
the sampling error, about :math:`1/\\sqrt{n}`. Tests on these must compare
distributions, never recorded numbers.
"""

from __future__ import annotations

import numpy as np

import IMP.bff


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
