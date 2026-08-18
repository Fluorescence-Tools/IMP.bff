"""Solvent-accessible surface area (Shrake-Rupley) for selected atoms.

Used by the PET model to weight a quencher by how exposed it is: a tryptophan
buried in the core cannot be reached by a dye however close its centre comes.
Cold -- once per simulation -- so it stays a straightforward numba loop.

Moved here from QuEst (``quest/core/av.py``) by PRD-109.
"""

from __future__ import annotations

import numpy as np

import IMP.bff


__all__ = ["sphere_points", "solvent_accessible_surface"]

#: Sphere-point count QuEst used, kept as the default so numbers do not move.
DEFAULT_SPHERE_POINTS = 590


def sphere_points(n: int) -> np.ndarray:
    """*n* roughly equidistant points on the unit sphere (golden spiral).

    **C++** (:file:`include/IMP/bff/SolventAccessibleSurface.h`). The Python it
    replaces built the spiral in ``float32``; this is ``float64`` and agrees
    with an independent float64 construction to 1e-15, where the old one was
    off by up to 2.4e-6. The port is the more accurate of the two.
    """
    return np.asarray(IMP.bff.sphere_points(int(n)),
                      dtype=np.float64).reshape(-1, 3)


def _asa(r, vdw, probe_atom_indices, points, probe, radius):
    """Shrake-Rupley area of the selected atoms. **C++.**

    The neighbour cutoff is ``(2*radius + probe)**2``, which is the tight
    criterion here: a sample sits ``radius`` from its own centre and is occluded
    by a neighbour within ``radius + probe`` of it. QuEst's kernel tested a
    *squared* distance against an *unsquared* sum, and against ``vdw`` rather
    than the sampling radius -- a 2.45 A cutoff where 6.0 A was meant, so nearly
    every occluding neighbour was missed. Fixed on the move (PRD-109); nothing
    downstream changed, because nothing consumed the result.
    """
    return np.asarray(
        IMP.bff.solvent_accessible_surface_area(
            np.ascontiguousarray(r, dtype=np.float64).ravel(),
            np.ascontiguousarray(vdw, dtype=np.float64).ravel(),
            [int(i) for i in np.asarray(probe_atom_indices).ravel()],
            np.ascontiguousarray(points, dtype=np.float64).ravel(),
            float(probe), float(radius)),
        dtype=np.float64)


def solvent_accessible_surface(
    xyz,
    vdw,
    probe_atom_indices,
    points=None,
    probe: float = 1.0,
    radius: float = 2.5,
) -> np.ndarray:
    """Accessible surface area (A^2) of each atom in *probe_atom_indices*.

    :param xyz: ``(n, 3)`` atom coordinates.
    :param vdw: ``(n,)`` van der Waals radii.
    :param probe_atom_indices: indices of the atoms to measure.
    :param points: unit-sphere sample points; :func:`sphere_points` by default.
    :param probe: probe radius.
    :param radius: radius at which the sphere points are placed.
    """
    xyz = np.asarray(xyz, dtype=np.float32)
    vdw = np.asarray(vdw, dtype=np.float32)
    probe_atom_indices = np.asarray(probe_atom_indices, dtype=np.uint32)
    if points is None or np.asarray(points).ndim != 2:
        points = sphere_points(DEFAULT_SPHERE_POINTS)
    points = np.asarray(points, dtype=np.float32)
    return _asa(xyz, vdw, probe_atom_indices, points, float(probe), float(radius))
