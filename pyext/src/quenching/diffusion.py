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

import IMP.bff

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
    occupancy = np.ascontiguousarray(np.asarray(density, dtype=np.int32))
    ng = int(occupancy.shape[0])
    slow = np.asarray(slow_fact, dtype=np.float64)

    # Mobility is a *field*: one scaling per voxel. The scalar-plus-mask form
    # is just one way to build it, and building it here collapses what used to
    # be two near-identical kernels into one.
    if slow.ndim == 3:
        mobility = np.ascontiguousarray(slow, dtype=np.float64)
    elif float(slow) == 1.0 or slow_density is None:
        mobility = np.empty(0, dtype=np.float64)   # uniform medium
    else:
        mask = np.asarray(slow_density, dtype=bool)
        mobility = np.where(mask, float(slow), 1.0)

    accepted_out = IMP.bff.VectorInt()
    counts = IMP.bff.VectorInt()
    flat = IMP.bff.brownian_walk_in_volume(
        occupancy.ravel(), mobility.ravel(), ng, float(dg), float(t_max),
        float(t_step), float(D), seed, accepted_out, counts)

    xyz = np.asarray(flat, dtype=np.float64).reshape(-1, 3)
    accepted = np.asarray(accepted_out, dtype=np.uint8)
    n_acc, n_rej = (int(counts[0]), int(counts[1])) if len(counts) == 2 else (0, 0)
    if xyz.size == 0:
        # No accessible starting voxel: an empty trajectory of the right length.
        n_steps = int(t_max / t_step)
        xyz = np.zeros((n_steps, 3), dtype=np.float64)
        accepted = np.zeros(n_steps, dtype=np.uint8)
    return DyeDiffusionTrajectory(xyz, accepted, n_acc, n_rej)
