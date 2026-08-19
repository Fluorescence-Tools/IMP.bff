"""Per-voxel stickiness and quenching-rate maps stamped onto an AV grid.

Both kernels stamp spheres of influence -- one per quencher or slow centre --
onto the occupied voxels of an accessible-volume density grid, so the Brownian
walk can read a local slow-down factor and a local quenching rate by indexing
rather than by testing every centre every frame.

Moved here from QuEst (``quest/core/av.py``) by PRD-109.

**One indexing convention, everywhere** -- see :func:`grid_center_index`. Both
registration defects this code once had were *upstream's* originally (float
corner, truncate-toward-zero); the fixes travel with the code.
"""

from __future__ import annotations

import numpy as np

import IMP.bff


__all__ = [
    "grid_center_index",
    "slow_factor_grid",
    "quenching_rate_grid",
    "av_contact_mask",
]


def grid_center_index(ng: int) -> int:
    """The voxel index that the grid anchor ``x0`` sits on.

    **One convention, everywhere.** A cubic grid of edge *ng* anchored at ``x0``
    places voxel ``i`` at ``x0 + (i - grid_center_index(ng)) * dg``, and that is
    the *definition*, not one of two defensible readings.

    It follows from how the grid is built. On the ``IMP.bff`` path
    ``x0 = grid_origin + ((ng - 1) // 2) * dg``, so voxel ``i`` sits at
    ``grid_origin + i * dg`` by construction, and the integer offset is the one
    that inverts that. The float corner ``(ng - 1) / 2`` differs whenever *ng*
    is **even** -- the normal case, not an edge case -- and integer voxel
    indexing cannot express it.

    QuEst had four call sites disagreeing about this until 2026-07-28: the grid
    stamps used the integer offset and the trajectory kernels the float one, so
    on an even *ng* the walk read the quenching map half a voxel from where the
    quenchers had been placed.

    Numba cannot call this from an ``njit`` kernel, so the kernels below spell
    ``(ng - 1) // 2`` inline. If you add a call site, use this function; if you
    cannot, quote this docstring.
    """
    return (int(ng) - 1) // 2


def _center_grid_indices(rs, r0, dg, ng, radius):
    """Voxel coordinates and integer radii of the sphere centres. **C++.**

    Uses ``floor``, not truncation. ``int()`` truncates *toward zero*, so a
    centre on the negative side of ``r0`` would round up while every other map
    here rounds down -- the walk's occupancy test and the trajectory sampler
    both floor. The two disagreed by one voxel per axis for every centre with a
    negative offset, which is half the grid.
    """
    ix0, iy0, iz0, ridx = IMP.bff.center_grid_indices(
        np.ascontiguousarray(rs, dtype=np.float64).ravel(),
        np.ascontiguousarray(r0, dtype=np.float64).ravel(),
        float(dg), int(ng),
        np.ascontiguousarray(radius, dtype=np.float64).ravel())
    return (np.asarray(ix0, dtype=np.int64), np.asarray(iy0, dtype=np.int64),
            np.asarray(iz0, dtype=np.int64), np.asarray(ridx, dtype=np.int64))


def _stamp(density, ng, radius, rs, r0, dg, values, combine):
    """One stamping kernel for both combines. **C++.**

    Stickiness *multiplies* (identity 1) and rates *add* (identity 0) -- the
    only difference between what used to be two near-identical Python kernels.
    Rates add because parallel channels do.
    """
    out = IMP.bff.stamp_spheres(
        np.ascontiguousarray(density, dtype=np.float64).ravel(), int(ng),
        np.ascontiguousarray(radius, dtype=np.float64).ravel(),
        np.ascontiguousarray(rs, dtype=np.float64).ravel(),
        np.ascontiguousarray(r0, dtype=np.float64).ravel(), float(dg),
        np.ascontiguousarray(values, dtype=np.float64).ravel(), int(combine))
    return np.asarray(out, dtype=np.float64).reshape(int(ng), int(ng), int(ng))


def _slow_factor_grid(density, ng, slow_radius, rs, r0, dg, slow_fact):
    """Stickiness field: factors multiply."""
    return _stamp(density, ng, slow_radius, rs, r0, dg, slow_fact,
                  IMP.bff.GRID_COMBINE_MULTIPLY)


def _additive_factor_grid(density, ng, radius, rs, r0, dg, values):
    """Rate field: parallel channels add."""
    return _stamp(density, ng, radius, rs, r0, dg, values,
                  IMP.bff.GRID_COMBINE_ADD)


def slow_factor_grid(density, ng, dg, slow_radius, rs, r0, slow_fact):
    """Per-voxel diffusion scaling from overlapping sticky spheres.

    Stickiness **multiplies** where spheres overlap. Voxels outside the AV keep
    1.0, so the factor is only meaningful where the walk can go.

    :param density: ``(ng, ng, ng)`` binary occupancy of the accessible volume.
    :param dg: voxel edge in Angstrom.
    :param slow_radius: ``(n,)`` radius of each sticky sphere.
    :param rs: ``(n, 3)`` sphere centres.
    :param r0: the grid anchor (see :func:`grid_center_index`).
    :param slow_fact: ``(n,)`` factor in [0, 1] per centre.
    """
    return _slow_factor_grid(
        np.asarray(density, dtype=np.uint8),
        int(ng),
        np.asarray(slow_radius, dtype=np.float64),
        np.asarray(rs, dtype=np.float64),
        np.asarray(r0, dtype=np.float64),
        float(dg),
        np.asarray(slow_fact, dtype=np.float64),
    )


def quenching_rate_grid(density, ng, dg, radius, rs, r0, values):
    """Per-voxel quenching rate (1/ns) from overlapping quencher spheres.

    Rates **add up** where contact spheres overlap, which is the physical
    composition rule for independent PET channels. Same geometry and the same
    indexing convention as :func:`slow_factor_grid`; only the accumulator
    differs.
    """
    return _additive_factor_grid(
        np.asarray(density, dtype=np.uint8),
        int(ng),
        np.asarray(radius, dtype=np.float64),
        np.asarray(rs, dtype=np.float64),
        np.asarray(r0, dtype=np.float64),
        float(dg),
        np.asarray(values, dtype=np.float64),
    )


def av_contact_mask(density, ng, dg, slow_radius, rs, r0):
    """The contact ("slow") part of an AV density grid, as a binary mask.

    A thin adapter over :func:`IMP.bff.representation.av._kernels.split_av_acv`, which returns
    both parts plus their counts; for a binary density the contact part is what
    the PET model wants.
    """
    from IMP.bff.representation.av import split_av_acv

    density = np.asarray(density, dtype=np.uint8)
    ng = int(ng)
    _, _, contact, _ = split_av_acv(
        density.reshape(ng, ng, ng).astype(np.float64),
        float(dg),
        np.asarray(slow_radius, dtype=np.float64),
        np.asarray(rs, dtype=np.float64),
        np.asarray(r0, dtype=np.float64),
    )
    return np.ascontiguousarray(contact.reshape(density.shape), dtype=np.uint8)
