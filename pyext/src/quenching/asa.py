"""Solvent-accessible surface area (Shrake-Rupley) for selected atoms.

Used by the PET model to weight a quencher by how exposed it is: a tryptophan
buried in the core cannot be reached by a dye however close its centre comes.
Cold -- once per simulation -- so it stays a straightforward numba loop.

Moved here from QuEst (``quest/core/av.py``) by PRD-109.
"""

from __future__ import annotations

import numpy as np

from .._jit import njit

__all__ = ["sphere_points", "solvent_accessible_surface"]

#: Sphere-point count QuEst used, kept as the default so numbers do not move.
DEFAULT_SPHERE_POINTS = 590


def sphere_points(n: int) -> np.ndarray:
    """*n* roughly equidistant points on the unit sphere (golden-spiral)."""
    inc = np.pi * (3.0 - np.sqrt(5.0))
    offset = 2.0 / float(n)
    k = np.arange(n, dtype=np.float32)
    y = k * offset - 1.0 + (offset / 2.0)
    rd = np.sqrt(1.0 - y * y)
    phi = k * inc
    x = np.cos(phi) * rd
    z = np.sin(phi) * rd
    return np.stack([x, y, z], axis=1).astype(np.float32)


@njit(cache=True)
def _asa(r, vdw, probe_atom_indices, points, probe, radius):
    n_sphere = points.shape[0]
    n_probe_atoms = probe_atom_indices.shape[0]
    n_atoms = r.shape[0]
    if n_probe_atoms == 0 or n_sphere == 0 or n_atoms == 0:
        return np.zeros(n_probe_atoms, dtype=np.float32)

    neighbor_indices = np.zeros(n_atoms, dtype=np.uint32)
    probe_atom_asa = np.zeros(n_probe_atoms, dtype=np.float32)
    c = 4.0 * np.pi / float(n_sphere)
    # A sample point sits `radius` from its own centre and is occluded by a
    # neighbour within `radius + probe` of it, so a neighbour can matter only
    # within `2 * radius + probe`. This is the tight criterion for *this*
    # kernel, which places every sample at a uniform `radius` and tests
    # occlusion at a uniform `radius + probe`.
    #
    # QuEst's kernel tested `dist2 < 2 * (vdw[i] + probe)` -- a **squared**
    # distance against an unsquared sum, and against `vdw` rather than the
    # radius the samples actually sit at. At the defaults that is a cutoff of
    # sqrt(6) = 2.45 A instead of 6.0 A, so almost every occluding neighbour
    # was missed and the area came back far too large. Fixed in the move
    # (PRD-109); nothing downstream moved, because nothing consumed the result
    # -- QuEst computed it on every quencher setup and its `quencher_asa`
    # accessor had no `return` statement.
    neighbor_cutoff2 = (2.0 * radius + probe) ** 2
    for i in range(n_probe_atoms):
        probe_atom_index = probe_atom_indices[i]
        n_neighbor = 0
        aX, aY, aZ = r[probe_atom_index]
        for atom_index in range(n_atoms):
            if probe_atom_index != atom_index:
                bX, bY, bZ = r[atom_index]
                dist2 = (aX - bX) ** 2 + (aY - bY) ** 2 + (aZ - bZ) ** 2
                if dist2 < neighbor_cutoff2:
                    neighbor_indices[n_neighbor] = atom_index
                    n_neighbor += 1

        n_accessible_point = 0
        for j in range(n_sphere):
            is_accessible = True
            aX = points[j, 0] * radius + r[probe_atom_index, 0]
            aY = points[j, 1] * radius + r[probe_atom_index, 1]
            aZ = points[j, 2] * radius + r[probe_atom_index, 2]
            for k in range(n_neighbor):
                bX, bY, bZ = r[neighbor_indices[k]]
                dist2 = (aX - bX) ** 2 + (aY - bY) ** 2 + (aZ - bZ) ** 2
                if dist2 < (radius + probe) * (radius + probe):
                    is_accessible = False
                    break
            if is_accessible:
                n_accessible_point += 1
        probe_atom_asa[i] = c * n_accessible_point * vdw[probe_atom_index] ** 2
    return probe_atom_asa


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
