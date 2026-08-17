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

from .._jit import njit, prange

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


@njit(cache=True)
def _center_grid_indices(rs, r0, dg, ng, radius):
    """Voxel coordinates and integer radii of the sphere centres."""
    n = rs.shape[0]
    ix0 = np.empty(n, dtype=np.int64)
    iy0 = np.empty(n, dtype=np.int64)
    iz0 = np.empty(n, dtype=np.int64)
    radius_idx = np.empty(n, dtype=np.int64)
    offset = (ng - 1) // 2
    for i in range(n):
        # `floor`, not `int`. `int()` truncates **toward zero**, so for a centre
        # on the negative side of `r0` it rounds *up* while every other map here
        # rounds down -- the walk's own occupancy test is `int(pos)` with
        # `pos >= 0`, which is floor, and so is the trajectory sampler. The two
        # disagreed by one voxel per axis for every quencher with a negative
        # offset, which is half the grid.
        ix0[i] = int(np.floor((rs[i, 0] - r0[0]) / dg)) + offset
        iy0[i] = int(np.floor((rs[i, 1] - r0[1]) / dg)) + offset
        iz0[i] = int(np.floor((rs[i, 2] - r0[2]) / dg)) + offset
        radius_idx[i] = int(radius[i] / dg)
    return ix0, iy0, iz0, radius_idx


# The two kernels below visit only the voxels inside each centre's bounding box
# instead of testing every voxel against every centre. The outer loop stays over
# x-slabs so the work is parallel and free of write races.
@njit(cache=True, parallel=True)
def _slow_factor_grid(density, ng, slow_radius, rs, r0, dg, slow_fact):
    factors = np.ones(density.shape, dtype=np.float64)
    n_slow_center = rs.shape[0]
    ix0s, iy0s, iz0s, radius_idxs = _center_grid_indices(rs, r0, dg, ng, slow_radius)

    for ix in prange(ng):
        slab = np.ones((ng, ng), dtype=np.float64)
        touched = False
        for isa in range(n_slow_center):
            r_idx = radius_idxs[isa]
            dx = ix - ix0s[isa]
            remaining = r_idx * r_idx - dx * dx
            if remaining <= 0:
                continue
            factor = slow_fact[isa]
            iy0 = iy0s[isa]
            iz0 = iz0s[isa]
            y_lo = max(0, iy0 - r_idx)
            y_hi = min(ng - 1, iy0 + r_idx)
            for iy in range(y_lo, y_hi + 1):
                dy = iy - iy0
                span2 = remaining - dy * dy
                if span2 <= 0:
                    continue
                span = int(np.sqrt(span2))
                # The membership test is strict (d^2 < r^2), so drop the
                # boundary voxel when span2 is a perfect square.
                if span * span >= span2:
                    span -= 1
                if span < 0:
                    continue
                z_lo = max(0, iz0 - span)
                z_hi = min(ng - 1, iz0 + span)
                for iz in range(z_lo, z_hi + 1):
                    slab[iy, iz] *= factor
                    touched = True
        if touched:
            for iy in range(ng):
                for iz in range(ng):
                    if density[ix, iy, iz] != 0:
                        factors[ix, iy, iz] = slab[iy, iz]
    return factors


@njit(cache=True, parallel=True)
def _additive_factor_grid(density, ng, radius, rs, r0, dg, values):
    factors = np.zeros(density.shape, dtype=np.float64)
    n_center = rs.shape[0]
    ix0s, iy0s, iz0s, radius_idxs = _center_grid_indices(rs, r0, dg, ng, radius)

    for ix in prange(ng):
        slab = np.zeros((ng, ng), dtype=np.float64)
        touched = False
        for isa in range(n_center):
            value = values[isa]
            if value == 0.0:
                continue
            r_idx = radius_idxs[isa]
            dx = ix - ix0s[isa]
            remaining = r_idx * r_idx - dx * dx
            if remaining <= 0:
                continue
            iy0 = iy0s[isa]
            iz0 = iz0s[isa]
            y_lo = max(0, iy0 - r_idx)
            y_hi = min(ng - 1, iy0 + r_idx)
            for iy in range(y_lo, y_hi + 1):
                dy = iy - iy0
                span2 = remaining - dy * dy
                if span2 <= 0:
                    continue
                span = int(np.sqrt(span2))
                if span * span >= span2:
                    span -= 1
                if span < 0:
                    continue
                z_lo = max(0, iz0 - span)
                z_hi = min(ng - 1, iz0 + span)
                for iz in range(z_lo, z_hi + 1):
                    slab[iy, iz] += value
                    touched = True
        if touched:
            for iy in range(ng):
                for iz in range(ng):
                    if density[ix, iy, iz] != 0:
                        factors[ix, iy, iz] = slab[iy, iz]
    return factors


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

    A thin adapter over :func:`IMP.bff.av._kernels.split_av_acv`, which returns
    both parts plus their counts; for a binary density the contact part is what
    the PET model wants.
    """
    from ..av._kernels import split_av_acv

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
