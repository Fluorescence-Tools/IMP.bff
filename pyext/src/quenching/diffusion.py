"""Brownian dynamics of a dye sphere inside its accessible volume.

The dye is a sphere diffusing on the AV's occupancy grid: a Gaussian step is
proposed each frame, accepted if it lands on an occupied voxel and rejected
(position held) otherwise, with the step width scaled down inside sticky
regions. This is the coarse counterpart of
:class:`IMP.bff.LangevinDyeSampler`, which integrates an explicit all-atom dye
under a force field -- a different model at a different scale.

Moved here from QuEst (``quest/core/av.py``) by PRD-109.

Two stickiness forms are supported and share one entry point: a **scalar**
factor applied wherever a separate binary "slow" grid is occupied, and a
**per-voxel** factor grid from :func:`IMP.bff.slow_factor_grid`. The step width
scales with ``sqrt(factor)`` in both, because the factor scales the diffusion
coefficient and the step width goes as ``sqrt(D)``.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np

from .._jit import njit

__all__ = ["DyeDiffusionTrajectory", "simulate_dye_diffusion"]


class DyeDiffusionTrajectory(NamedTuple):
    """A Brownian trajectory of the dye centre, in Angstrom.

    ``xyz`` is relative to the grid anchor: add the attachment point to place it
    in the structure's frame.
    """

    xyz: np.ndarray
    accepted: np.ndarray
    n_accepted: int
    n_rejected: int

    @property
    def n_frames(self) -> int:
        return int(self.xyz.shape[0])

    @property
    def acceptance_ratio(self) -> float:
        total = self.n_accepted + self.n_rejected
        return float(self.n_accepted) / total if total else 0.0


# Sequential on purpose: the walk is a strict chain of dependent steps, so
# `parallel=True` buys no speed-up and splits the random draws across threads,
# which makes `random_seed` non-reproducible.
#
# `nogil` is what makes running *whole* trajectories concurrently worthwhile:
# without it a caller's thread pool serialises on the GIL. Numba's RNG state is
# thread-local, so each worker still reproduces its own `random_seed`.
@njit(cache=True, nogil=True)
def _simulate_scalar(d, ds, dg, t_max, t_step, D, slow_fact, random_seed):
    if random_seed >= 0:
        np.random.seed(random_seed)
    n_samples = int(t_max / t_step)
    ng = d.shape[0]
    sigma = np.sqrt(2 * D * 3 * t_step) / dg
    pos = np.zeros((n_samples, 3), dtype=np.float64)
    accepted = np.zeros(n_samples, dtype=np.uint8)
    r = np.random.normal(0.0, sigma, size=(n_samples, 3))
    slow_fact_sqrt = np.sqrt(slow_fact)

    found = False
    for _ in range(1000):
        x_idx = np.random.randint(0, ng)
        y_idx = np.random.randint(0, ng)
        z_idx = np.random.randint(0, ng)
        if d[x_idx, y_idx, z_idx] > 0:
            pos[0, 0] = x_idx
            pos[0, 1] = y_idx
            pos[0, 2] = z_idx
            found = True
            break
    if not found:
        return pos, accepted, 0, 0

    n_accepted = 1
    n_rejected = 0
    i_accepted = 0
    accepted[0] = 1
    for i in range(1, n_samples):
        x_idx = int(pos[i_accepted, 0])
        y_idx = int(pos[i_accepted, 1])
        z_idx = int(pos[i_accepted, 2])
        if ds[x_idx, y_idx, z_idx] > 0:
            r[i, 0] *= slow_fact_sqrt
            r[i, 1] *= slow_fact_sqrt
            r[i, 2] *= slow_fact_sqrt

        pos[i, 0] = pos[i_accepted, 0] + r[i, 0]
        pos[i, 1] = pos[i_accepted, 1] + r[i, 1]
        pos[i, 2] = pos[i_accepted, 2] + r[i, 2]

        x_idx = int(pos[i, 0])
        y_idx = int(pos[i, 1])
        z_idx = int(pos[i, 2])
        if 0 <= x_idx < ng and 0 <= y_idx < ng and 0 <= z_idx < ng and d[x_idx, y_idx, z_idx] > 0:
            i_accepted = i
            n_accepted += 1
            accepted[i] = 1
        else:
            pos[i, 0] = pos[i_accepted, 0]
            pos[i, 1] = pos[i_accepted, 1]
            pos[i, 2] = pos[i_accepted, 2]
            n_rejected += 1
            accepted[i] = 0
    # Integer offset, matching `grids._center_grid_indices` -- see
    # `grids.grid_center_index`. Numba cannot call that helper from here.
    return (pos - (ng - 1) // 2) * dg, accepted, n_accepted, n_rejected


# Sequential on purpose, and `nogil` for concurrency: see `_simulate_scalar`.
@njit(cache=True, nogil=True)
def _simulate_grid(d, slow_factor, dg, t_max, t_step, D, random_seed):
    if random_seed >= 0:
        np.random.seed(random_seed)
    n_samples = int(t_max / t_step)
    ng = d.shape[0]
    sigma = np.sqrt(2 * D * 3 * t_step) / dg
    pos = np.zeros((n_samples, 3), dtype=np.float64)
    accepted = np.zeros(n_samples, dtype=np.uint8)
    r = np.random.normal(0.0, sigma, size=(n_samples, 3))

    found = False
    for _ in range(1000):
        x_idx = np.random.randint(0, ng)
        y_idx = np.random.randint(0, ng)
        z_idx = np.random.randint(0, ng)
        if d[x_idx, y_idx, z_idx] > 0:
            pos[0, 0] = x_idx
            pos[0, 1] = y_idx
            pos[0, 2] = z_idx
            found = True
            break
    if not found:
        return pos, accepted, 0, 0

    n_accepted = 1
    n_rejected = 0
    i_accepted = 0
    accepted[0] = 1
    for i in range(1, n_samples):
        x_idx = int(pos[i_accepted, 0])
        y_idx = int(pos[i_accepted, 1])
        z_idx = int(pos[i_accepted, 2])
        factor = slow_factor[x_idx, y_idx, z_idx]
        if factor < 1.0:
            factor_sqrt = np.sqrt(factor)
            r[i, 0] *= factor_sqrt
            r[i, 1] *= factor_sqrt
            r[i, 2] *= factor_sqrt

        pos[i, 0] = pos[i_accepted, 0] + r[i, 0]
        pos[i, 1] = pos[i_accepted, 1] + r[i, 1]
        pos[i, 2] = pos[i_accepted, 2] + r[i, 2]

        x_idx = int(pos[i, 0])
        y_idx = int(pos[i, 1])
        z_idx = int(pos[i, 2])
        if 0 <= x_idx < ng and 0 <= y_idx < ng and 0 <= z_idx < ng and d[x_idx, y_idx, z_idx] > 0:
            i_accepted = i
            n_accepted += 1
            accepted[i] = 1
        else:
            pos[i, 0] = pos[i_accepted, 0]
            pos[i, 1] = pos[i_accepted, 1]
            pos[i, 2] = pos[i_accepted, 2]
            n_rejected += 1
            accepted[i] = 0
    # Integer offset, matching `grids._center_grid_indices`.
    return (pos - (ng - 1) // 2) * dg, accepted, n_accepted, n_rejected


def simulate_dye_diffusion(
    density,
    slow_density=None,
    dg: float = 0.5,
    t_max: float = 10000.0,
    t_step: float = 0.002,
    D: float = 40.0,
    slow_fact=1.0,
    random_seed=None,
) -> DyeDiffusionTrajectory:
    """Run one Brownian trajectory of the dye centre in its accessible volume.

    :param density: ``(ng, ng, ng)`` binary occupancy of the accessible volume.
    :param slow_density: binary "slow" (contact) grid, used with a **scalar**
        *slow_fact*. Ignored when *slow_fact* is a 3-D grid.
    :param dg: voxel edge in Angstrom.
    :param t_max: total simulated time in ns.
    :param t_step: time step in ns.
    :param D: diffusion coefficient in A^2/ns.
    :param slow_fact: either a scalar factor applied inside *slow_density*, or a
        ``(ng, ng, ng)`` per-voxel factor grid from
        :func:`IMP.bff.slow_factor_grid`.
    :param random_seed: seed for a reproducible walk; ``None`` draws freely.

    The returned coordinates are relative to the grid anchor.
    """
    seed = -1 if random_seed is None else int(random_seed)
    slow = np.asarray(slow_fact)
    if slow.ndim == 3:
        result = _simulate_grid(
            np.asarray(density, dtype=np.uint8),
            np.asarray(slow_fact, dtype=np.float64),
            float(dg),
            float(t_max),
            float(t_step),
            float(D),
            seed,
        )
    else:
        if slow_density is None:
            slow_density = np.zeros_like(np.asarray(density, dtype=np.uint8))
        result = _simulate_scalar(
            np.asarray(density, dtype=np.uint8),
            np.asarray(slow_density, dtype=np.uint8),
            float(dg),
            float(t_max),
            float(t_step),
            float(D),
            float(slow_fact),
            seed,
        )
    return DyeDiffusionTrajectory(*result)
