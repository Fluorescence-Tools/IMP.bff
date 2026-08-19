"""ACV — Accessible Contact Volume with trapped dye fraction.

The ACV extends BasicAV by splitting the density into a **contact**
(slow-diffusion) and a **non-contact** (free-diffusion) part.  A
``trapped_fraction`` parameter controls how much of the density
resides in the contact volume near specified slow centres.

Examples
--------
>>> import numpy as np
>>> from IMP.bff.representation.av import BasicAV, ACV
>>> pts = np.random.randn(500, 4).astype(np.float64)
>>> pts[:, 3] = np.abs(pts[:, 3]) + 0.01  # positive weights
>>> av = BasicAV(points=pts)
>>> acv = ACV.from_basic_av(av, trapped_fraction=0.5, slow_radius=10.0)
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from IMP.bff.representation.av.basic import BasicAV
from IMP.bff.representation.av import _kernels


class ACV(BasicAV):
    """Accessible Contact Volume with trapped dye fraction.

    Parameters
    ----------
    points : (n, 4) ndarray, optional
        Pre-computed point cloud.
    density : (ng, ng, ng) ndarray, optional
        Density grid (will be split into contact / non-contact).
    grid_origin : (3,) ndarray, optional
        Coordinates of the first voxel centre.
    grid_step : float
        Voxel spacing (Å).
    slow_centers : (n, 3) ndarray, optional
        Coordinates of slow-diffusion centres.
    slow_radius : float or (n,) ndarray
        Radius around each centre (Å).  Default 10.0.
    trapped_fraction : float
        Fraction of dye trapped in the contact volume.  Default 0.8.
    position_name : str
        Human-readable label.
    """

    def __init__(
        self,
        points: Optional[np.ndarray] = None,
        density: Optional[np.ndarray] = None,
        grid_origin: Optional[np.ndarray] = None,
        grid_step: float = 1.5,
        slow_centers: Optional[np.ndarray] = None,
        slow_radius: float | np.ndarray = 10.0,
        trapped_fraction: float = 0.8,
        position_name: str = "",
    ):
        super().__init__(
            points=points,
            density=density,
            grid_origin=grid_origin,
            grid_step=grid_step,
            position_name=position_name,
        )
        self._slow_centers = slow_centers
        self._slow_radius = np.atleast_1d(
            np.asarray(slow_radius, dtype=np.float64)
        )
        self._trapped_fraction = float(trapped_fraction)
        self._contact_density = None

        if density is not None and slow_centers is not None:
            self._update_contact_density()

    # ------------------------------------------------------------------
    # Factory
    # ------------------------------------------------------------------
    @classmethod
    def from_basic_av(
        cls,
        av: BasicAV,
        slow_centers: np.ndarray,
        slow_radius: float | np.ndarray = 10.0,
        trapped_fraction: float = 0.8,
    ) -> "ACV":
        """Create an ACV from an existing BasicAV.

        Parameters
        ----------
        av : BasicAV
            Source AV (needs ``density`` and ``grid_origin``).
        slow_centers : (n, 3) ndarray
            Slow-centre coordinates.
        slow_radius : float
            Radius around centres (Å).
        trapped_fraction : float
            Dye fraction trapped in contact volume.

        Returns
        -------
        ACV
        """
        return cls(
            density=av.density,
            grid_origin=av._grid_origin,
            grid_step=av.grid_step,
            slow_centers=slow_centers,
            slow_radius=slow_radius,
            trapped_fraction=trapped_fraction,
            position_name=av.position_name,
        )

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------
    @property
    def slow_centers(self) -> Optional[np.ndarray]:
        """Coordinates of slow-diffusion centres."""
        return self._slow_centers

    @slow_centers.setter
    def slow_centers(self, value: Optional[np.ndarray]):
        self._slow_centers = value
        if value is not None and self._density is not None:
            self._update_contact_density()

    @property
    def slow_radius(self) -> np.ndarray:
        """Radius (Å) around each slow centre."""
        return self._slow_radius

    @slow_radius.setter
    def slow_radius(self, value: float | np.ndarray):
        self._slow_radius = np.atleast_1d(np.asarray(value, dtype=np.float64))
        if self._slow_centers is not None and self._density is not None:
            self._update_contact_density()

    @property
    def trapped_fraction(self) -> float:
        """Fraction of dye trapped in the contact volume."""
        return self._trapped_fraction

    @trapped_fraction.setter
    def trapped_fraction(self, value: float):
        self._trapped_fraction = float(value)
        if self._slow_centers is not None and self._density is not None:
            self._update_contact_density()

    @property
    def contact_density(self) -> Optional[np.ndarray]:
        """Portion of the density that belongs to the contact volume."""
        return self._contact_density

    # ------------------------------------------------------------------
    # Contact-volume logic
    # ------------------------------------------------------------------
    def set_slow_centers_from_atoms(
        self, atoms_xyz: np.ndarray, atom_name: str = "CB"
    ):
        """Set slow centres from an atom coordinate array.

        This convenience method mimics the chisurf behaviour where
        ``slow_centers='CB'`` uses all ``CB`` atoms as slow centres.

        Parameters
        ----------
        atoms_xyz : (N, 3) or (N, 4+) ndarray
            Atom coordinates.  If 4+ columns, the 4th column is
            expected to contain vdW radii (currently unused for
            targeting).
        atom_name : str
            Unused here; the caller is expected to pre-filter.
        """
        self._slow_centers = np.asarray(atoms_xyz[:, :3], dtype=np.float64)
        if self._density is not None:
            self._update_contact_density()

    def _update_contact_density(self):
        """Recompute the split between contact and non-contact density."""
        if self._density is None or self._grid_origin is None:
            return
        if self._slow_centers is None or len(self._slow_centers) == 0:
            return

        n_idx, n_non_idx, contact_mask, non_contact_mask = _kernels.split_av_acv(
            self._density,
            self.grid_step,
            self._slow_radius,
            self._slow_centers,
            self._grid_origin,
        )

        cd = self._density * contact_mask
        nd = self._density * non_contact_mask

        if n_idx > 0:
            cd *= self._trapped_fraction / n_idx
        if n_non_idx > 0:
            nd *= (1.0 - self._trapped_fraction) / n_non_idx

        self._contact_density = cd
        new_density = cd + nd
        new_density /= new_density.sum()
        self._density = new_density
        self.update_points()

    def update(self):
        """Recompute the contact split and the point cloud."""
        self._update_contact_density()
        BasicAV.update_points(self)

    # ------------------------------------------------------------------
    # String representation
    # ------------------------------------------------------------------
    def __repr__(self) -> str:
        return (
            f"ACV(name={self.position_name!r}, n_points={self.n_points}, "
            f"trapped={self._trapped_fraction:.2f})"
        )
