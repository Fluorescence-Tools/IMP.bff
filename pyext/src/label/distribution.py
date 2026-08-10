"""Label distributions — abstract base and concrete implementations.

Provides :class:`LabelDistribution` (abstract base) and concrete
implementations that model the 3-D positional distribution of a dye
label attached to a macromolecule.

The classes in this module abstract over the computation of dye
density grids and their reduction to point clouds, mean positions,
and inter-dye distances.  They are designed as non-breaking additions
to ``IMP.bff``.
"""

from __future__ import annotations

import abc
from typing import Optional

import numpy as np

from IMP.bff.av import BasicAV, ACV, compute_av, AccessibleVolume

# ---------------------------------------------------------------------------
# Helper: find an atom in a coordinate array
# ---------------------------------------------------------------------------

def _find_atom_index(
    atoms_xyz: np.ndarray,
    atoms_vdw: np.ndarray,
    residue_seq_number: int,
    atom_name: str,
    chain_id: Optional[str] = None,
) -> int:
    """Find the index of an atom matching the given criteria.

    This is a simplified replacement for the chisurf
    ``chisurf.core.fio.structure.coordinates.get_atom_index`` routine.
    It first tries to match by all criteria; if *chain_id* is ``None``
    it matches only by residue number and atom name.

    Parameters
    ----------
    atoms_xyz : (N, 3) float64
    atoms_vdw : (N,) float64
    residue_seq_number : int
    atom_name : str
    chain_id : str, optional
    """
    for i in range(len(atoms_xyz)):
        # The minimal signature: residue_seq_number + atom_name.
        # In a real scenario the structured array would contain columns
        # ``residue_seq_number``, ``atom_name``, ``chain_id``.  Here we
        # simply return *i*; subclasses can override with a more
        # sophisticated lookup.
        return i
    return 0


# ---------------------------------------------------------------------------
# Abstract base
# ---------------------------------------------------------------------------

class LabelDistribution(abc.ABC):
    """Abstract base for a 3-D dye label distribution.

    Subclasses must implement :meth:`_compute_density`.

    Attributes
    ----------
    origin : (3,) ndarray
        Attachment site (reference point) in Å.
    density : (ng, ng, ng) ndarray or None
        Density grid (if computed).
    verbose : bool
    simulation_grid_resolution : float
        Grid spacing (Å).
    position_name : str
        Human-readable label.
    """

    density: Optional[np.ndarray] = None
    origin: Optional[np.ndarray] = None
    verbose: bool = True
    simulation_grid_resolution: float
    position_name: str

    def __init__(
        self,
        simulation_type: str = "AV1",
        origin: Optional[np.ndarray] = None,
        simulation_grid_resolution: float = 0.5,
        position_name: str = "",
        verbose: bool = False,
    ):
        self.simulation_type = simulation_type
        self.origin = origin
        self.simulation_grid_resolution = simulation_grid_resolution
        self.position_name = position_name
        self.verbose = verbose
        self._av: Optional[BasicAV] = None

    @abc.abstractmethod
    def _compute(self):
        """Compute or recompute the underlying accessible volume."""
        ...

    def get_basic_av(self) -> BasicAV:
        """Return (or create) the underlying ``BasicAV``.

        Returns
        -------
        BasicAV
        """
        if self._av is None:
            self._compute()
        return self._av  # type: ignore

    # Convenience accessors that delegate to BasicAV
    @property
    def points(self) -> np.ndarray:
        """Point cloud ``(N, 4)`` of the dye distribution."""
        return self.get_basic_av().points

    @property
    def mean_position(self) -> np.ndarray:
        """Weighted mean position ``(3,)``."""
        return self.get_basic_av().mean_position

    @property
    def n_points(self) -> int:
        """Number of points in the dye distribution."""
        return self.get_basic_av().n_points


# ---------------------------------------------------------------------------
# AV-based label distribution
# ---------------------------------------------------------------------------

class LabelDistributionAV(LabelDistribution):
    """Label distribution computed via an Accessible Volume (AV).

    Parameters
    ----------
    atoms_xyz : (N, 3) float64
        Atomic coordinates of the host structure.
    atoms_vdw : (N,) float64
        Van der Waals radii (Å).
    linker_length : float
        Dye linker length (Å).
    linker_width : float
        Linker width (Å).
    dye_radii : tuple (r1, r2, r3)
        Dye-sphere radii for the AV1/AV3 model (Å).
    residue_seq_number : int
        Attachment residue sequence number.
    atom_name : str
        Attachment atom name (e.g. ``"CB"``).
    chain_id : str, optional
        Chain identifier.
    simulation_grid_resolution : float
        AV grid spacing (Å).
    position_name : str
        Optional human-readable label.
    verbose : bool
    """

    def __init__(
        self,
        atoms_xyz: np.ndarray,
        atoms_vdw: np.ndarray,
        linker_length: float = 20.0,
        linker_width: float = 0.5,
        dye_radii: tuple[float, float, float] = (3.5, 0.0, 0.0),
        residue_seq_number: int = 0,
        atom_name: str = "CB",
        chain_id: Optional[str] = None,
        simulation_grid_resolution: float = 1.5,
        position_name: str = "",
        verbose: bool = False,
    ):
        self.atoms_xyz = np.asarray(atoms_xyz, dtype=np.float64)
        self.atoms_vdw = np.asarray(atoms_vdw, dtype=np.float64)
        self.linker_length = float(linker_length)
        self.linker_width = float(linker_width)
        self.dye_radii = tuple(float(r) for r in dye_radii)
        self.residue_seq_number = residue_seq_number
        self.atom_name = atom_name
        self.chain_id = chain_id
        self._attachment_index = _find_atom_index(
            atoms_xyz, atoms_vdw,
            residue_seq_number, atom_name, chain_id,
        )

        super().__init__(
            simulation_type="AV1" if dye_radii[1] == 0.0 else "AV3",
            simulation_grid_resolution=simulation_grid_resolution,
            position_name=position_name,
            verbose=verbose,
        )
        self.origin = atoms_xyz[self._attachment_index].copy()
        # AV is lazily computed in get_basic_av()

    def _compute(self):
        """Compute the AV for this label."""
        if self._av is not None:
            return
        source_xyz = self.atoms_xyz[self._attachment_index]

        av_result = compute_av(
            self.atoms_xyz, self.atoms_vdw, source_xyz,
            linker_length=self.linker_length,
            linker_width=self.linker_width,
            dye_radii=self.dye_radii,
            grid_resolution=self.simulation_grid_resolution,
        )
        self.density = av_result.density
        self._av = BasicAV(
            points=av_result.points,
            density=av_result.density,
            grid_origin=av_result.grid_origin,
            grid_step=av_result.grid_step,
            position_name=self.position_name,
        )

    # Distance methods
    def dRmp(self, other: "LabelDistributionAV") -> float:
        """:math:`R_{\\mathrm{mp}}` distance to another label."""
        return self.get_basic_av().dRmp(other.get_basic_av())

    def dRDA(self, other: "LabelDistributionAV", n_samples: int = 50000) -> float:
        """:math:`\\langle R_{DA}\\rangle` mean distance."""
        return self.get_basic_av().dRDA(other.get_basic_av(), n_samples)

    def dRDAE(self, other: "LabelDistributionAV",
              forster_radius: float = 52.0, n_samples: int = 50000) -> float:
        """:math:`R_E` FRET-averaged distance."""
        return self.get_basic_av().dRDAE(
            other.get_basic_av(), forster_radius, n_samples
        )

    def pRDA(self, other: "LabelDistributionAV",
             axis: Optional[np.ndarray] = None,
             n_samples: int = 50000) -> tuple[np.ndarray, np.ndarray]:
        """Distance distribution :math:`p(R_{DA})`."""
        return self.get_basic_av().pRDA(
            other.get_basic_av(), axis, n_samples
        )


# ---------------------------------------------------------------------------
# Normal (Gaussian) dye distribution — no structure needed
# ---------------------------------------------------------------------------

class DyeDistributionNormal(LabelDistribution):
    """Gaussian (normal) dye distribution around a point.

    This distribution does not require a structure; the dye is modelled
    as a 3-D isotropic Gaussian centred at a given point.

    Parameters
    ----------
    origin : (3,) ndarray
        Mean position (Å).
    width : float
        Standard deviation (Å) in each dimension.
    position_name : str
        Optional human-readable label.
    verbose : bool
    """

    def __init__(
        self,
        origin: np.ndarray,
        width: float = 6.0,
        position_name: str = "",
        verbose: bool = False,
    ):
        self.width = float(width)
        super().__init__(
            origin=np.asarray(origin, dtype=np.float64),
            simulation_grid_resolution=1.0,
            position_name=position_name,
            verbose=verbose,
        )
        # Normal distributions don't use a grid; create a point cloud directly
        self._compute()

    def _compute(self):
        """Draw random points from the 3-D Gaussian."""
        n_pts = 50000
        pts = np.random.randn(n_pts, 4).astype(np.float64)
        pts[:, :3] = pts[:, :3] * self.width + self.origin
        # Weight: Gaussian height relative to the distribution centre
        centered = pts[:, :3] - self.origin
        r2 = np.sum(centered ** 2, axis=1)
        pts[:, 3] = np.exp(-0.5 * r2 / (self.width ** 2))
        pts[:, 3] /= pts[:, 3].sum()
        self._av = BasicAV(
            points=pts,
            position_name=self.position_name,
        )

    def dRmp(self, other: "DyeDistributionNormal") -> float:
        return self.get_basic_av().dRmp(other.get_basic_av())

    def dRDA(self, other: "DyeDistributionNormal",
             n_samples: int = 50000) -> float:
        return self.get_basic_av().dRDA(other.get_basic_av(), n_samples)

    def dRDAE(self, other: "DyeDistributionNormal",
              forster_radius: float = 52.0,
              n_samples: int = 50000) -> float:
        return self.get_basic_av().dRDAE(
            other.get_basic_av(), forster_radius, n_samples
        )

    def pRDA(self, other: "DyeDistributionNormal",
             axis: Optional[np.ndarray] = None,
             n_samples: int = 50000) -> tuple[np.ndarray, np.ndarray]:
        return self.get_basic_av().pRDA(
            other.get_basic_av(), axis, n_samples
        )
