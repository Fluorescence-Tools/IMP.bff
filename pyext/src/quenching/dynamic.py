"""A dye whose accessible volume carries mobility, quenching and FRET fields.

:class:`DynamicAccessibleVolume` wraps an :class:`IMP.bff.AccessibleVolume` and
adds the three fields the grid solver needs, plus the two things they buy:

* an **equilibrium occupancy** in place of the flat AV density -- with a
  position-dependent diffusion coefficient the dye dwells where it moves
  slowly, so the stationary distribution is not uniform, and every average
  taken over it (⟨R_DA⟩ included) shifts;
* a **donor decay** integrated on the grid, quenching and FRET together, with
  no photon sampling and so no shot noise.

Moved here from ChiSurf (``chisurf/core/structure/av/__init__.py::DynamicAV``)
by PRD-109.

This is the field formulation. :class:`IMP.bff.DyeDiffusionSimulation` is the
particle one: it walks the dye and reads rates along the trajectory. They answer
the same question and disagree in an informative way -- the walk resolves the
dye's history (so it can produce a correlation function), the field resolves the
whole distribution at once (so it is deterministic and cheap to converge).
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from . import maps
from .solver import (
    GridDiffusionSolver,
    diffusion_stability_limit,
    equilibrium_occupancy,
)

__all__ = ["DynamicAccessibleVolume"]


class DynamicAccessibleVolume:
    """An accessible volume with diffusion, quenching and FRET rate fields.

    :param av: the :class:`IMP.bff.AccessibleVolume` to decorate.
    :param atoms: structured array of the obstacle atoms, with ``res_name``,
        ``atom_name`` and ``coord`` fields.
    :param tau0: unquenched donor lifetime in ns.
    :param dye_radius: dye sphere radius in Angstrom.
    :param free_diffusion: unhindered diffusion coefficient in A^2/ns.
    :param contact_distance: dye-to-atom distance counted as contact, for the
        mobility field.
    :param slow_factor: mobility multiplier applied per contacting atom.
    """

    def __init__(
        self,
        av,
        atoms,
        *,
        tau0: float = 4.0,
        dye_radius: float = 3.5,
        free_diffusion: float = 8.0,
        contact_distance: float = 6.5,
        slow_factor: float = 0.985,
    ):
        self.av = av
        self.atoms = atoms
        self.tau0 = float(tau0)
        self.dye_radius = float(dye_radius)
        self.free_diffusion = float(free_diffusion)
        self.contact_distance = float(contact_distance)
        self.slow_factor = float(slow_factor)

        self._diffusion_map: Optional[np.ndarray] = None
        self._quenching_rate_map: Optional[np.ndarray] = None
        self._fret_rate_map: Optional[np.ndarray] = None
        self._occupancy: Optional[np.ndarray] = None

    # -- grid geometry -------------------------------------------------------

    @property
    def density(self) -> np.ndarray:
        """Binary occupancy of the accessible volume."""
        return np.ascontiguousarray(self.av.density, dtype=np.float64)

    @property
    def dg(self) -> float:
        return float(self.av.grid_step)

    @property
    def x0(self) -> np.ndarray:
        """The grid anchor -- the attachment point."""
        return np.asarray(self.av.attachment_point, dtype=np.float64)

    @property
    def bounds(self) -> np.ndarray:
        """Where the dye may be: the domain mask for the solver."""
        return (self.density > 0).astype(np.float64)

    # -- the three fields ----------------------------------------------------

    @property
    def diffusion_map(self) -> np.ndarray:
        if self._diffusion_map is None:
            self.update_diffusion_map()
        return self._diffusion_map

    @property
    def quenching_rate_map(self) -> np.ndarray:
        if self._quenching_rate_map is None:
            raise ValueError(
                "No quenching map yet -- call update_quenching_map(quencher) "
                "with the per-atom PET parameters."
            )
        return self._quenching_rate_map

    @property
    def fret_rate_map(self) -> Optional[np.ndarray]:
        return self._fret_rate_map

    @property
    def rate_map(self) -> np.ndarray:
        """Total decay rate per voxel: quenching, plus FRET if an acceptor is set."""
        total = self.quenching_rate_map
        if self._fret_rate_map is not None:
            total = total + self._fret_rate_map
        return total

    def update_diffusion_map(self, radial_profile=None) -> np.ndarray:
        """Build the mobility field: free diffusion, slowed by nearby atoms."""
        self._diffusion_map = maps.diffusion_coefficient_map(
            self.density,
            self.x0,
            self.dg,
            np.ascontiguousarray(self.atoms["coord"], dtype=np.float64),
            free_diffusion=self.free_diffusion,
            min_distance=self.contact_distance,
            slow_factor=self.slow_factor,
            radial_profile=radial_profile,
        )
        return self._diffusion_map

    def update_quenching_map(self, quencher, rC: Optional[float] = None) -> np.ndarray:
        """Build the PET field from per-atom ``(kQ, rC)`` parameters.

        :param quencher: ``{res_name: {atom_name: (kQ, rC)}}``.
        :param rC: overrides every per-atom characteristic distance with one
            electron-transfer length, which is how ChiSurf drove it.
        """
        kQ, rC_atoms = maps.atomic_quenching_parameters(self.atoms, quencher)
        if rC is not None:
            rC_atoms = np.where(kQ > 0.0, float(rC), 0.0)
        self._quenching_rate_map = maps.quenching_rate_map(
            self.density,
            self.x0,
            self.dg,
            np.ascontiguousarray(self.atoms["coord"], dtype=np.float64),
            kQ,
            rC_atoms,
            tau0=self.tau0,
            dye_radius=self.dye_radius,
        )
        return self._quenching_rate_map

    def update_fret_map(
        self,
        acceptor: "DynamicAccessibleVolume",
        forster_radius: float = 52.0,
        acceptor_step: int = 2,
    ) -> np.ndarray:
        """Build the FRET field against an acceptor's accessible volume.

        The donor's radiative rate is ``1/tau0``. ChiSurf read it from the
        ``foerster_radius`` keyword by a copy-paste slip
        (``kf = kwargs.get('foerster_radius', 1./tau0)``), so passing a Forster
        radius also set the radiative rate to it; taken from ``tau0`` here.
        """
        self._fret_rate_map = maps.fret_rate_map(
            self.density,
            acceptor.density,
            self.x0,
            acceptor.x0,
            self.dg,
            acceptor.dg,
            forster_radius,
            1.0 / self.tau0,
            acceptor_step=acceptor_step,
        )
        return self._fret_rate_map

    # -- propagation ---------------------------------------------------------

    def _solver(self, density, rate_map, t_step: Optional[float]) -> GridDiffusionSolver:
        d_map = self.diffusion_map
        if t_step is None:
            # Half the explicit limit: stable with room for the map to change.
            t_step = 0.5 * diffusion_stability_limit(float(d_map.max()), self.dg)
        return GridDiffusionSolver(
            d_map, self.bounds, density, rate_map, t_step=t_step, dg=self.dg
        )

    @property
    def occupancy(self) -> np.ndarray:
        """The equilibrium distribution of the dye, normalised.

        Not the AV density: that weights every accessible voxel equally, and a
        dye that moves slowly near the surface spends more time there.
        """
        if self._occupancy is None:
            self.update_occupancy()
        return self._occupancy

    def update_occupancy(self, t_step: Optional[float] = None, **kwargs) -> np.ndarray:
        """The equilibrium occupancy, in closed form.

        ``p ∝ 1/D`` — see :func:`IMP.bff.quenching.solver.equilibrium_occupancy`.
        Propagating to it instead is possible but slow and, on a real site where
        the compounding slow factor makes ``D`` span orders of magnitude, may not
        converge at all: on T4L site 132 it was still drifting after 40 000
        iterations. Pass ``iterate=True`` to do it the long way anyway.
        """
        if kwargs.pop("iterate", False):
            self._occupancy = self._solver(
                self.bounds, None, t_step).equilibrium(**kwargs)
        else:
            self._occupancy = equilibrium_occupancy(
                self.diffusion_map, self.bounds)
        return self._occupancy

    def donor_decay(
        self, t_max: float = 50.0, t_step: Optional[float] = None, n_out: int = 10
    ):
        """Integrate the donor decay on the grid, from the equilibrium start.

        :returns: a :class:`IMP.bff.quenching.solver.GridDiffusionResult` --
            time axis in ns, surviving excited-state fraction, final density.
        """
        start = self.occupancy
        solver = self._solver(start, self.rate_map, t_step)
        n_steps = max(1, int(t_max / solver.t_step))
        return solver.run(n_steps, n_out=n_out)
