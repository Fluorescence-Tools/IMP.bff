"""BasicAV — accessible volume for a single dye labelling site.

BasicAV stores the accessible volume as a point cloud and provides
distance calculations (R_DA, R_mp, R_E, pRDA) and I/O methods.
"""

from __future__ import annotations

import json
import os
from typing import Optional

import numpy as np

from IMP.bff.av import _kernels

BACKENDS_AVAILABLE = False
try:
    import IMP
    import IMP.bff
    BACKENDS_AVAILABLE = True
except ImportError:
    pass


class BasicAV:
    """Accessible volume for a single dye labelling site.

    Parameters
    ----------
    points : (n, 4) ndarray, optional
        Pre-computed point cloud (x, y, z, weight).  If given, the AV
        is defined directly from these points without backend computation.
    density : (nx, ny, nz) ndarray, optional
        3-D density grid (will be converted to points on demand).
    grid_origin : (3,) ndarray, optional
        Coordinates of the first voxel centre when *density* is provided.
    grid_step : float, optional
        Voxel spacing (Å).
    position_name : str, optional
        Human-readable name for this AV (e.g. ``"donor_72"``).

    Attributes
    ----------
    points : (n, 4) ndarray
        Point cloud ``(x, y, z, weight)``.
    n_points : int
        Number of points in the cloud.
    density : (nx, ny, nz) ndarray or None
        Voxel density grid (if available).
    grid_origin, grid_step : see above.
    mean_position : (3,) ndarray
        Density-weighted mean position of the AV.

    Examples
    --------
    >>> import numpy as np
    >>> from IMP.bff.av import BasicAV
    >>> # Create a simple spherical AV
    >>> n = 500
    >>> pts = np.random.randn(n, 4).astype(np.float64)
    >>> pts[:, 3] = 1.0 / n  # uniform weights
    >>> av = BasicAV(points=pts, position_name="test")
    >>> av.n_points == n
    True
    >>> av.mean_position.shape
    (3,)
    """

    def __init__(
        self,
        points: Optional[np.ndarray] = None,
        density: Optional[np.ndarray] = None,
        grid_origin: Optional[np.ndarray] = None,
        grid_step: float = 1.5,
        position_name: str = "",
    ):
        self._points = None
        self._padded_points = None  # may be larger than n_points
        self.n_points = 0
        self._density = density
        self._grid_origin = grid_origin
        self.grid_step = float(grid_step)
        self.position_name = position_name
        self._min_points = 150

        if points is not None:
            self.set_points(points)

        if density is not None and points is None:
            self.update_points()

    # ------------------------------------------------------------------
    # Point cloud
    # ------------------------------------------------------------------
    def set_points(self, points: np.ndarray) -> None:
        """Set the point cloud directly.

        Parameters
        ----------
        points : (n, 4) ndarray
            Columns: ``x, y, z, weight``.
        """
        self._padded_points = np.asarray(points, dtype=np.float64)
        self.n_points = len(points)
        self._points = None

    def update_points(self) -> None:
        """Convert the internal density grid to a point cloud."""
        if self._density is None or self._grid_origin is None:
            raise ValueError("No density grid available to convert.")
        ng = self._density.shape[0]
        n, pts = _kernels.density2points(
            ng, ng, ng, self.grid_step, self._density, self._grid_origin
        )
        self._padded_points = pts
        self.n_points = n

    @property
    def points(self) -> np.ndarray:
        """The valid portion of the point cloud as a view.

        Returns
        -------
        (n_points, 4) ndarray
        """
        if self._padded_points is None or self.n_points == 0:
            return np.empty((0, 4), dtype=np.float64)
        return self._padded_points[: self.n_points]

    @property
    def density(self) -> Optional[np.ndarray]:
        """Density grid (may be None)."""
        return self._density

    # ------------------------------------------------------------------
    # Mean position
    # ------------------------------------------------------------------
    @property
    def mean_position(self) -> np.ndarray:
        """Weighted mean position ``(x, y, z)``."""
        return _kernels.weighted_mean(self._padded_points, self.n_points)

    # ------------------------------------------------------------------
    # Distance calculations
    # ------------------------------------------------------------------
    def dRmp(self, other: "BasicAV") -> float:
        r"""Distance between mean positions (Rmp).

        Parameters
        ----------
        other : BasicAV

        Returns
        -------
        float
            :math:`R_{\\mathrm{mp}}` in Å.
        """
        mp1 = self.mean_position
        mp2 = other.mean_position
        return float(
            np.sqrt(np.sum((mp1 - mp2) ** 2))
        )

    def dRDA(self, other: "BasicAV", n_samples: int = 50000) -> float:
        r"""Mean inter-dye distance :math:`\\langle R_{DA} \\rangle`.

        Sampled via random weighted pairing of points from both AVs.

        Parameters
        ----------
        other : BasicAV
        n_samples : int
            Number of random pairs to sample.

        Returns
        -------
        float
            :math:`\\langle R_{DA} \\rangle` in Å.
        """
        return _kernels.average_distance(
            self._padded_points, self.n_points,
            other._padded_points, other.n_points,
            n_samples,
        )

    def dRDAE(self, other: "BasicAV", forster_radius: float, n_samples: int = 50000) -> float:
        r"""FRET-averaged distance :math:`R_E`.

        Parameters
        ----------
        other : BasicAV
        forster_radius : float
            Förster radius :math:`R_0` (Å).
        n_samples : int

        Returns
        -------
        float
            :math:`R_E` in Å.
        """
        return _kernels.mean_fret_distance(
            self._padded_points, self.n_points,
            other._padded_points, other.n_points,
            forster_radius, n_samples,
        )

    def pRDA(
        self,
        other: "BasicAV",
        axis: Optional[np.ndarray] = None,
        n_samples: int = 50000,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""Distance distribution :math:`p(R_{DA})`.

        Parameters
        ----------
        other : BasicAV
        axis : (n_bins + 1,) ndarray, optional
            Bin edges for the histogram.  Default:
            ``np.arange(0, 200, 0.5)``.
        n_samples : int

        Returns
        -------
        y : ndarray
            Probability density.
        x : ndarray
            Bin centres.
        """
        if axis is None:
            axis = np.arange(0.0, 200.0, 0.5)
        d = _kernels.random_distances(
            self._padded_points[: self.n_points],
            other._padded_points[: other.n_points],
            n_samples,
        )
        y, _ = np.histogram(d[:, 0], bins=axis, weights=d[:, 1])
        if np.sum(y) > 0:
            y = y.astype(np.float64) / np.sum(y)
        x = 0.5 * (axis[:-1] + axis[1:])
        return y, x

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    def save_xyz(self, filename: str) -> None:
        """Save the point cloud as an XYZ file.

        Parameters
        ----------
        filename : str
            Output path.
        """
        pts = self.points
        np.savetxt(
            filename,
            np.column_stack([pts[:, :3], pts[:, 3] / max(pts[:, 3].max(), 1e-12) * 50]),
            fmt="%.6f",
            header=f"{self.n_points} atoms",
        )

    # ------------------------------------------------------------------
    # Pretty printing
    # ------------------------------------------------------------------
    def __repr__(self) -> str:
        return (
            f"BasicAV(name={self.position_name!r}, "
            f"n_points={self.n_points})"
        )
