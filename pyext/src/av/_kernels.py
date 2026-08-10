"""Numba-accelerated kernels for accessible volume operations.

These kernels provide the low-level computations used by AV classes.
They are pure Numba (no IMP or LabelLib dependency) and can be used
standalone.

Every kernel here is ``@_njit(cache=True, nogil=True)``, and both flags are
load-bearing rather than decoration:

``nogil``
    ``density2points`` visits every voxel of the grid and is the single most
    expensive step of a labelling-site computation. A caller scanning many
    sites in a thread pool -- which is what QuEst's web backend does -- would
    otherwise serialise on the GIL precisely where the work is.
``cache``
    Without it every one of these compiles afresh in each new process, which a
    short CLI run pays in full and never amortises.

A caution that goes with ``cache=True``: a warm cache is only re-read, never
re-validated against the compiler configuration, so a green suite proves
nothing about compilation after these are edited. Clear ``*.nbi``/``*.nbc``
under ``__pycache__`` before trusting a run.
"""

from __future__ import annotations

import math
import numpy as np
from IMP.bff._jit import njit as _njit, jit as _jit


@_njit(cache=True, nogil=True)
def random_distances(
    p1: np.ndarray,
    p2: np.ndarray,
    n_samples: int,
    seed: int = 0,
) -> np.ndarray:
    """Draw random distance-weight pairs from two AV point clouds.

    Each sample picks one random point from each AV and computes the
    Euclidean distance and the product of their density weights.

    Parameters
    ----------
    p1 : (n1, 4) float64
        Points of the first AV (x, y, z, weight).
    p2 : (n2, 4) float64
        Points of the second AV (x, y, z, weight).
    n_samples : int
        Number of random samples to draw.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    (n_samples, 2) float64
        Column 0: Euclidean distances (Å).
        Column 1: weight products.
    """
    np.random.seed(seed)
    n1 = p1.shape[0]
    n2 = p2.shape[0]
    result = np.empty((n_samples, 2), dtype=np.float64)
    for i in range(n_samples):
        i1 = np.random.randint(0, n1) if n1 > 0 else 0
        i2 = np.random.randint(0, n2) if n2 > 0 else 0
        dx = p1[i1, 0] - p2[i2, 0]
        dy = p1[i1, 1] - p2[i2, 1]
        dz = p1[i1, 2] - p2[i2, 2]
        result[i, 0] = math.sqrt(dx * dx + dy * dy + dz * dz)
        result[i, 1] = p1[i1, 3] * p2[i2, 3]
    return result


@_njit(cache=True, nogil=True)
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
    points : (n_points, 4) float64
        Point cloud ``(x, y, z, weight)``.  The array may be larger than
        *n*; unused trailing entries should be discarded.
    """
    max_points = nx * ny * nz
    points = np.zeros((max_points, 4), dtype=np.float64)
    n = 0
    for ix in range(nx):
        x = dg * ix + r0[0]
        for iy in range(ny):
            y = dg * iy + r0[1]
            for iz in range(nz):
                val = density[ix, iy, iz]
                if val > threshold:
                    points[n, 0] = x
                    points[n, 1] = y
                    points[n, 2] = dg * iz + r0[2]
                    points[n, 3] = val
                    n += 1
    return n, points


@_njit(cache=True, nogil=True)
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
        Weighted mean coordinate.
    """
    if n == 0:
        return np.zeros(3, dtype=np.float64)
    w_sum = 0.0
    mean = np.zeros(3, dtype=np.float64)
    for i in range(n):
        w = points[i, 3]
        w_sum += w
        mean[0] += points[i, 0] * w
        mean[1] += points[i, 1] * w
        mean[2] += points[i, 2] * w
    if w_sum > 0.0:
        mean /= w_sum
    return mean


@_njit(cache=True, nogil=True)
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
        Mean distance <R_DA> (Å).
    """
    d = random_distances(points1[:n1], points2[:n2], n_samples)
    w_sum = 0.0
    rda = 0.0
    for i in range(n_samples):
        w_sum += d[i, 1]
        rda += d[i, 0] * d[i, 1]
    if w_sum > 0.0:
        return rda / w_sum
    return 0.0


@_njit(cache=True, nogil=True)
def mean_fret_distance(points1: np.ndarray, n1: int,
                       points2: np.ndarray, n2: int,
                       forster_radius: float = 52.0,
                       n_samples: int = 50000) -> float:
    """FRET-averaged distance R_E between two AVs.

    Parameters
    ----------
    points1, points2 : AV point clouds.
    n1, n2 : Number of valid rows.
    forster_radius : float
        Förster radius R0 (Å).
    n_samples : int

    Returns
    -------
    float
        FRET-averaged distance R_E (Å).
    """
    d = random_distances(points1[:n1], points2[:n2], n_samples)
    w_sum = 0.0
    mean_e = 0.0
    for i in range(n_samples):
        r = d[i, 0]
        w = d[i, 1]
        e = 1.0 / (1.0 + (r / forster_radius) ** 6.0)
        w_sum += w
        mean_e += e * w
    if w_sum > 0.0:
        mean_e /= w_sum
    if mean_e <= 0.0 or mean_e >= 1.0:
        return 0.0 if mean_e >= 1.0 else math.inf
    return forster_radius * ((1.0 / mean_e) - 1.0) ** (1.0 / 6.0)


@_njit(cache=True, nogil=True)
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
        Slow radius per center (broadcast scalar to all if needed).
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
    ng = density.shape[0]
    n_radii = rs.shape[0]

    if len(radius) != n_radii:
        rad = np.zeros(n_radii, dtype=np.float64) + radius[0]
    else:
        rad = radius

    contact = np.zeros((ng, ng, ng), dtype=np.uint8)
    non_contact = np.zeros((ng, ng, ng), dtype=np.uint8)
    n_contact = 0
    n_non = 0

    # Integer offset, and `floor` below rather than `int`. See `r0` above.
    half = (ng - 1) // 2

    for ix in range(ng):
        for iy in range(ng):
            for iz in range(ng):
                if density[ix, iy, iz] <= 0.0:
                    continue

                overlapped = 0
                for isa in range(n_radii):
                    ix0 = int(np.floor((rs[isa, 0] - r0[0]) / dg)) + half
                    iy0 = int(np.floor((rs[isa, 1] - r0[1]) / dg)) + half
                    iz0 = int(np.floor((rs[isa, 2] - r0[2]) / dg)) + half
                    r_idx = int(rad[isa] / dg)
                    dx = ix - ix0
                    dy = iy - iy0
                    dz = iz - iz0
                    if dx * dx + dy * dy + dz * dz < r_idx * r_idx:
                        overlapped = 1
                        break

                if overlapped > 0:
                    contact[ix, iy, iz] = 1
                    n_contact += 1
                else:
                    non_contact[ix, iy, iz] = 1
                    n_non += 1

    return n_contact, n_non, contact, non_contact
