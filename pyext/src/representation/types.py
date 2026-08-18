"""The accessible volume: a region a dye can reach, as one of several representations."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np

from .states import States

__all__ = ["AccessibleVolume"]


@dataclass(kw_only=True)
class AccessibleVolume(States):
    """States enumerated as a voxel grid, plus the grid itself.

    One definition. There were two identical ones -- ``IMP.bff.av.compute`` and
    ``IMP.bff.fret.av``, same seven fields in the same order -- which is how they
    came to disagree about the axis order of ``density`` without anything
    noticing (PRD-113 stage 3a).

    :param density: ``(nx, ny, nz)`` accessible density, in the **coordinate**
        axis order. IMP orders its flat tile values with *x* fastest; a C-order
        reshape into ``(nx, ny, nz)`` returns the volume transposed, and a
        mirrored volume has the right voxel count, bounding box and total volume,
        so only a voxel-by-voxel comparison against ``points`` catches it.
    :param grid_origin: ``(3,)`` coordinate of the first voxel centre.
    :param grid_step: voxel edge in Angstrom.
    :param grid_shape: ``(nx, ny, nz)``.
    """

    density: np.ndarray
    grid_origin: np.ndarray
    grid_step: float
    grid_shape: Tuple[int, int, int]
