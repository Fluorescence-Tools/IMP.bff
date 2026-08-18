"""States: positions, weights and orientations -- what every representation supplies.

An accessible-volume grid point, a rotamer, a coarse-grained conformer and an MD
frame are the same kind of thing: *a state the dye can occupy, with a weight, a
position and possibly an orientation*. What differs is only how the states were
generated. Everything downstream -- distances, kappa^2, the interaction terms --
consumes states, so it is written once and works for all of them.

This is the abstraction that was missing. ``RotamerEnsemble`` inherited from the
concrete :class:`~IMP.bff.representation.AccessibleVolume` instead, which is why
distance code happened to work for rotamers: by inheritance, not by design. The
cost showed up in the fields a rotamer library then had to carry and could not
fill -- ``density=zeros((0, 0, 0))``, ``grid_step=0.0``, ``grid_shape=(0, 0, 0)``
-- placeholders for a grid that does not exist. No consumer ever read them:
every one of them used ``points``, ``mean_position``, ``n_points`` or
``has_volume``, which is exactly this surface.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional

import numpy as np

__all__ = ["States"]


@dataclass(kw_only=True)
class States:
    """A weighted set of states of one label.

    :param points: ``(N, 4)`` -- ``x, y, z, weight``. One row per state.
    :param attachment_point: ``(3,)`` where the label is tied to the structure.
    :param orientations: ``(N, 3)`` transition dipoles, when the representation
        resolves them. ``None`` for a positional-only model such as an AV, which
        is why kappa^2 from an AV needs an isotropic assumption and kappa^2 from
        a rotamer library does not.
    :param position_name: human-readable label for the site.
    :param params: how these states were produced -- representation parameters,
        not dye or site properties.
    """

    points: np.ndarray
    attachment_point: np.ndarray
    orientations: Optional[np.ndarray] = None
    position_name: str = ""
    params: Dict = field(default_factory=dict)

    @property
    def positions(self) -> np.ndarray:
        """``(N, 3)`` state coordinates."""
        return self.points[:, :3]

    @property
    def weights(self) -> np.ndarray:
        """``(N,)`` state weights, unnormalised."""
        return self.points[:, 3]

    @property
    def n_points(self) -> int:
        return self.points.shape[0] if self.points.ndim == 2 else 0

    @property
    def has_volume(self) -> bool:
        return self.n_points > 0

    @property
    def has_orientations(self) -> bool:
        return self.orientations is not None and len(self.orientations) > 0

    @property
    def mean_position(self) -> np.ndarray:
        """Weight-averaged position, falling back to the attachment point."""
        if self.n_points == 0:
            return np.asarray(self.attachment_point).copy()
        w = self.points[:, 3]
        if w.sum() == 0:
            return np.asarray(self.attachment_point).copy()
        return np.average(self.points[:, :3], axis=0, weights=w)
