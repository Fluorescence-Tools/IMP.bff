"""The particle model, end to end: structure and site in, donor decay out.

:class:`DyeDiffusionSimulation` walks a dye sphere in its accessible volume and
reports the quenching rate it experiences frame by frame.
:class:`QuenchedDonorDecay` puts that together with the accessible volume, the
PET chemistry and the photon Monte-Carlo, and is the object a caller drives.

Moved here from QuEst (``quest/core/dye_diffusion.py``:
``SimulateDiffusion``, ``DonorDecay``) by PRD-109. What did **not** come with
them is QuEst's application layer -- the project document, the four surfaces,
scanning, run artefacts, file saving and the verbose printing. These take arrays
and an :class:`IMP.bff.AccessibleVolume`, and return numbers.

The field counterpart is :class:`IMP.bff.DynamicAccessibleVolume`.
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from typing import Optional, Sequence

import numpy as np

from . import sites
from .diffusion import simulate_dye_diffusion
from .fret_trace import fret_rate_pair_trace, fret_rate_trace
from .grids import grid_center_index, quenching_rate_grid, slow_factor_grid
from .pet import DEFAULT_DYE_RADIUS, normalize_amino_acid_quenching
from .photon import simulate_photon_trace

__all__ = ["DyeDiffusionSimulation", "QuenchedDonorDecay", "MAX_PARALLEL_TRAJECTORIES"]

#: Independent trajectories are concatenated before photon sampling. More than
#: this buys little: the walk is already ergodic over a run of ~10 us.
MAX_PARALLEL_TRAJECTORIES = 8


def _trajectory_seeds(random_seed, n_trajectories):
    """One seed per trajectory, derived from the base seed or drawn freshly."""
    if random_seed is None:
        if n_trajectories <= 1:
            return [None]
        return [
            int(seed)
            for seed in np.random.SeedSequence().generate_state(
                n_trajectories, dtype=np.uint32
            )
        ]
    base = int(random_seed)
    max_seed = (2 ** 31) - 1
    # A large prime stride, so trajectories from one base seed do not share a
    # low-order pattern the way `base + i` would.
    return [int((base + i * 104729) % max_seed) for i in range(n_trajectories)]


def _resolve_parallel(parallel_trajectories):
    requested = -1 if parallel_trajectories is None else int(parallel_trajectories)
    if requested == -1:
        return max(1, min(os.cpu_count() or 1, MAX_PARALLEL_TRAJECTORIES))
    return max(1, min(requested, MAX_PARALLEL_TRAJECTORIES))


class DyeDiffusionSimulation:
    """A dye's Brownian walk in its accessible volume, and the rate it sees.

    :param density: binary AV occupancy, ``(ng, ng, ng)``.
    :param dg: voxel edge in Angstrom.
    :param x0: the grid anchor -- the attachment point.
    :param slow_density: binary contact grid, used with a scalar *slow_fact*.
    :param slow_factor_map: per-voxel stickiness, used instead of the pair
        above when given.
    :param quenching_rate_map: per-voxel quenching rate (1/ns), sampled along
        the trajectory to give :attr:`k_quench`.
    """

    def __init__(
        self,
        density,
        dg: float,
        x0,
        slow_density=None,
        slow_factor_map=None,
        quenching_rate_map=None,
    ):
        self.density = np.ascontiguousarray(density, dtype=np.uint8)
        self.dg = float(dg)
        self.x0 = np.asarray(x0, dtype=np.float64)
        self.slow_density = (
            None if slow_density is None
            else np.ascontiguousarray(slow_density, dtype=np.uint8)
        )
        self.slow_factor_map = slow_factor_map
        self.quenching_rate_map = quenching_rate_map

        self.trajectory: Optional[np.ndarray] = None
        self.t_step: Optional[float] = None
        self.n_accepted = 0
        self.n_rejected = 0

    def run(
        self,
        D: float = 40.0,
        slow_fact: float = 0.01,
        t_step: float = 0.002,
        t_max: float = 10000.0,
        n_trajectories: int = -1,
        random_seed=None,
    ) -> Optional[np.ndarray]:
        """Simulate the walk and return the trajectory in the structure's frame.

        *n_trajectories* independent walks are run and **concatenated**, not
        averaged: the photon Monte-Carlo downstream draws a random start frame
        per photon, so a longer concatenated record is exactly what it wants,
        and the walks parallelise because the kernel releases the GIL.

        :returns: ``(n_frames, 3)`` positions in Angstrom, or ``None`` if no
            walk found a starting point (an empty accessible volume).
        """
        self.t_step = float(t_step)
        n_trajectories = _resolve_parallel(n_trajectories)
        seeds = _trajectory_seeds(random_seed, n_trajectories)

        slow = self.slow_factor_map if self.slow_factor_map is not None else slow_fact

        def one(seed):
            return simulate_dye_diffusion(
                self.density,
                self.slow_density,
                self.dg,
                t_max=t_max,
                t_step=t_step,
                D=D,
                slow_fact=slow,
                random_seed=seed,
            )

        if n_trajectories == 1:
            results = [one(seeds[0])]
        else:
            with ThreadPoolExecutor(max_workers=n_trajectories) as pool:
                results = list(pool.map(one, seeds))

        successful = [r for r in results if r.n_accepted > 0]
        if not successful:
            self.trajectory = None
            return None

        parts = [r.xyz for r in successful]
        trajectory = np.concatenate(parts, axis=0) if len(parts) > 1 else parts[0]
        self.trajectory = trajectory + self.x0
        self.n_accepted = int(sum(r.n_accepted for r in successful))
        self.n_rejected = int(sum(r.n_rejected for r in successful))
        return self.trajectory

    @property
    def n_frames(self) -> int:
        return 0 if self.trajectory is None else int(self.trajectory.shape[0])

    @property
    def mean_position(self) -> np.ndarray:
        if self.trajectory is None:
            raise ValueError("Run the simulation first.")
        return self.trajectory.mean(axis=0)

    def sample_grid(self, grid) -> np.ndarray:
        """Read a per-voxel field along the trajectory.

        The index map here **must** be the one the grids were stamped with, or
        the walk reads rates from beside where the quenchers were placed. Two
        details carry that, and both were bugs once:

        * the centre offset is the integer ``(ng - 1) // 2``
          (:func:`IMP.bff.grid_center_index`), not the float ``(ng - 1) / 2``,
          which differ on every even edge length -- and even is the normal case;
        * the conversion is ``floor``, not ``trunc``: ``trunc`` maps ``[-1, 0)``
          to 0, so a position up to one voxel *below* the grid would be treated
          as inside it and read voxel 0's rate.
        """
        if self.trajectory is None:
            raise ValueError("Run the simulation first.")
        values = np.asarray(grid, dtype=np.float32)
        centre = np.array(
            [grid_center_index(n) for n in values.shape], dtype=np.float64
        )
        indices = np.floor(
            (self.trajectory - self.x0) / self.dg + centre
        ).astype(np.int64)
        inside = np.ones(indices.shape[0], dtype=bool)
        for axis, size in enumerate(values.shape):
            inside &= (indices[:, axis] >= 0) & (indices[:, axis] < size)
        out = np.zeros(indices.shape[0], dtype=np.float32)
        if np.any(inside):
            out[inside] = values[tuple(indices[inside].T)]
        return out

    @property
    def k_quench(self) -> np.ndarray:
        """The quenching rate (1/ns) the dye experiences, frame by frame."""
        if self.trajectory is None:
            raise ValueError("Run the simulation first.")
        if self.quenching_rate_map is None:
            return np.zeros(self.n_frames, dtype=np.float32)
        return self.sample_grid(self.quenching_rate_map)

    @property
    def quenched(self) -> np.ndarray:
        """Which frames the dye spent in contact with a quencher."""
        return self.k_quench > 0.0

    @property
    def collision_fraction(self) -> float:
        quenched = self.quenched
        return float(quenched.sum()) / quenched.shape[0] if quenched.size else 0.0


class QuenchedDonorDecay:
    """A labelled site's donor decay, from the structure and the PET chemistry.

    :param av: the donor's :class:`IMP.bff.AccessibleVolume`.
    :param atoms: structured obstacle atoms (``chain``, ``res_id``,
        ``res_name``, ``atom_name``, ``coord``).
    :param tau0: unquenched donor lifetime in ns.
    :param quenching_table: per-residue interaction table; the defaults of
        :func:`IMP.bff.amino_acid_quenching_defaults` when omitted.
    :param critical_distance: contact radius inherited by residues whose table
        entry leaves ``quench_radius`` unset.
    :param slow_radius: radius of each residue's sticky sphere.
    :param n_photons: excitation events for the photon Monte-Carlo.
    """

    def __init__(
        self,
        av,
        atoms,
        *,
        tau0: float = 4.0,
        quenching_table=None,
        critical_distance: float = 7.0,
        slow_radius: float = 10.0,
        dye_radius: float = DEFAULT_DYE_RADIUS,
        diffusion_coefficient: float = 40.0,
        slow_fact: float = 0.05,
        t_step: float = 0.004,
        t_max: float = 10000.0,
        n_photons: int = 100000,
        n_trajectories: int = -1,
        random_seed=None,
    ):
        self.av = av
        self.atoms = atoms
        self.tau0 = float(tau0)
        self.quenching_table = normalize_amino_acid_quenching(quenching_table)
        self.critical_distance = float(critical_distance)
        self.slow_radius = float(slow_radius)
        self.dye_radius = float(dye_radius)
        self.diffusion_coefficient = float(diffusion_coefficient)
        self.slow_fact = min(1.0, max(0.0, float(slow_fact)))
        self.t_step = float(t_step)
        self.t_max = float(t_max)
        self.n_photons = int(n_photons)
        self.n_trajectories = n_trajectories
        self.random_seed = random_seed

        self._sites: Optional[sites.ResidueSites] = None
        self._quenching_rate_map = None
        self._slow_factor_map = None
        self._diffusion: Optional[DyeDiffusionSimulation] = None
        self._photon_trace = None

    # -- grid geometry -------------------------------------------------------

    @property
    def density(self) -> np.ndarray:
        return np.ascontiguousarray(self.av.density, dtype=np.uint8)

    @property
    def dg(self) -> float:
        return float(self.av.grid_step)

    @property
    def x0(self) -> np.ndarray:
        return np.asarray(self.av.attachment_point, dtype=np.float64)

    @property
    def sites(self) -> sites.ResidueSites:
        """The slow and quench centres of every residue in the structure."""
        if self._sites is None:
            self._sites = sites.residue_sites(self.atoms, self.quenching_table)
        return self._sites

    # -- the two grids -------------------------------------------------------

    @property
    def quenching_rate_map(self) -> np.ndarray:
        if self._quenching_rate_map is None:
            self.update_grids()
        return self._quenching_rate_map

    @property
    def slow_factor_map(self) -> np.ndarray:
        if self._slow_factor_map is None:
            self.update_grids()
        return self._slow_factor_map

    def update_grids(self):
        """Stamp the quenching-rate and stickiness fields onto the AV grid."""
        residue_sites = self.sites
        names = residue_sites.residue_names
        density = self.density
        ng = density.shape[0]

        self._quenching_rate_map = quenching_rate_grid(
            density, ng, self.dg,
            sites.quench_radii_for_residues(
                names, self.quenching_table, self.critical_distance
            ),
            residue_sites.quench_centers,
            self.x0,
            sites.quenching_rates_for_residues(names, self.quenching_table),
        )
        self._slow_factor_map = slow_factor_grid(
            density, ng, self.dg,
            np.full(len(names), self.slow_radius, dtype=np.float64),
            residue_sites.slow_centers,
            self.x0,
            sites.slow_factors_for_residues(names, self.quenching_table),
        )
        return self._quenching_rate_map, self._slow_factor_map

    # -- the walk ------------------------------------------------------------

    @property
    def diffusion(self) -> DyeDiffusionSimulation:
        if self._diffusion is None:
            self.simulate_diffusion()
        return self._diffusion

    def simulate_diffusion(self) -> bool:
        """Run the Brownian walk. Returns whether a trajectory was produced."""
        simulation = DyeDiffusionSimulation(
            self.density,
            self.dg,
            self.x0,
            slow_factor_map=self.slow_factor_map,
            quenching_rate_map=self.quenching_rate_map,
        )
        simulation.run(
            D=self.diffusion_coefficient,
            slow_fact=self.slow_fact,
            t_step=self.t_step,
            t_max=self.t_max,
            n_trajectories=self.n_trajectories,
            random_seed=self.random_seed,
        )
        self._diffusion = simulation
        self._photon_trace = None
        return simulation.trajectory is not None

    @property
    def k_quench(self) -> np.ndarray:
        return self.diffusion.k_quench

    # -- photons and the decay ----------------------------------------------

    def _photon_seed(self):
        """Distinct from the walk's, so the photon draws are not its replay."""
        if self.random_seed is None:
            return None
        return int((int(self.random_seed) + 7919) % ((2 ** 31) - 1))

    @property
    def photon_trace(self):
        if self._photon_trace is None:
            self.simulate_photons()
        return self._photon_trace

    def simulate_photons(self, k_quench=None):
        """Race photons against the quenching rate along the trajectory."""
        walk = self.diffusion
        if k_quench is None:
            k_quench = walk.k_quench
        self._photon_trace = simulate_photon_trace(
            self.n_photons,
            k_quench,
            t_step=walk.t_step,
            tau0=self.tau0,
            random_seed=self._photon_seed(),
        )
        return self._photon_trace

    @property
    def quantum_yield(self) -> float:
        """Emitted photons over excitation events."""
        if self.n_photons <= 0:
            return 0.0
        _delays, emitted = self.photon_trace
        return float(emitted.sum()) / self.n_photons

    @property
    def fluorescence_lifetime(self) -> float:
        """The species-averaged lifetime of the emitted photons."""
        delays, emitted = self.photon_trace
        kept = delays[emitted == 1]
        return float(kept.mean()) if kept.size else 0.0

    def lifetime_spectrum(self, n_species: int = 128):
        """The decay as ``(amplitude, rate)`` pairs -- the neutral output.

        Prefer this to :meth:`decay_histogram`. A spectrum carries no bin width
        and no time range, so whatever owns the instrument can convolve, bin and
        add noise on its own terms; a histogram has already chosen all three.
        See :mod:`IMP.bff.observables` for the contract.

        Built from the *per-frame* total rate along the trajectory: each frame
        is a state the dye occupies, weighted equally because the walk visits
        them in proportion to their occupancy. That makes this the **static
        approximation** -- exact only if the dye held each position for a whole
        excited-state lifetime, which is precisely what a diffusion simulation
        exists to deny. The returned spectrum is marked ``exact=False`` for that
        reason, and :meth:`decay_histogram` (which races each photon against the
        moving rate) is the one that resolves the averaging.

        :param n_species: coarse-grain to at most this many species. The
            trajectory has one state per frame -- often 10^5 -- and almost none
            are distinguishable.
        """
        from IMP.bff.observables import lifetime_spectrum_from_rates

        k = np.asarray(self.k_quench, dtype=np.float64).ravel()
        if k.size == 0:
            raise ValueError("no quenching rate trace; run the walk first")
        intrinsic = 1.0 / self.tau0 if self.tau0 > 0.0 else 0.0
        spectrum = lifetime_spectrum_from_rates(k + intrinsic, exact=False)
        return spectrum.coarse_grain(n_species) if n_species else spectrum

    def decay_histogram(self, n_bins: int = 4096, time_range=(0.0, 50.0)):
        """Histogram of the photons that were actually emitted.

        A binned curve is **not** the neutral output -- it has chosen a bin
        width and a time range, both of which are instrument settings. Use
        :meth:`lifetime_spectrum` unless the binning is the point. This is kept
        because the photon trace resolves the *time-dependence* of the quenching
        rate, which a spectrum cannot represent: it is the answer the static
        approximation is an approximation to.

        **Only emitted photons.** A quenched excitation comes back with
        ``dt = 0``, so histogramming the whole trace piles every non-emitted
        photon into the first bin -- a spike of photons that never existed, and
        a curve whose total is the excitation count rather than ``QY * N``.
        """
        delays, emitted = self.photon_trace
        counts, edges = np.histogram(
            delays[emitted == 1], range=time_range, bins=n_bins
        )
        return edges, counts

    # -- FRET ----------------------------------------------------------------

    def fret_rate_trace(
        self,
        acceptor,
        forster_radius: float = 52.0,
        kappa2=None,
        r_min: float = 7.0,
    ) -> np.ndarray:
        """Per-frame FRET rate against an acceptor.

        *acceptor* may be a :class:`QuenchedDonorDecay` whose walk has been run
        -- then both dyes are resolved in time and the frames are paired -- or
        anything with a ``points`` cloud, which is the fast-acceptor limit.
        Pairing requires the same time base; a mismatch raises rather than
        truncating.
        """
        donor_trajectory = self.diffusion.trajectory
        if donor_trajectory is None:
            raise ValueError("The donor walk produced no trajectory.")

        if isinstance(acceptor, QuenchedDonorDecay):
            acceptor_trajectory = acceptor.diffusion.trajectory
            if acceptor_trajectory is None:
                raise ValueError("The acceptor walk produced no trajectory.")
            n = min(donor_trajectory.shape[0], acceptor_trajectory.shape[0])
            return fret_rate_pair_trace(
                donor_trajectory[:n], acceptor_trajectory[:n],
                forster_radius, self.tau0, r_min=r_min, kappa2=kappa2,
            )

        points = np.asarray(getattr(acceptor, "points", acceptor))
        return fret_rate_trace(
            donor_trajectory, points[:, :3], forster_radius, self.tau0,
            r_min=r_min, kappa2=kappa2,
        )

    def fret_efficiency(
        self, acceptor, forster_radius: float = 52.0, kappa2=None
    ) -> float:
        """``1 - QY_DA / QY_D``, both from the photon Monte-Carlo.

        Taking the ratio of two simulated quantum yields rather than an analytic
        formula keeps the quenching in: the donor is quenched by PET in both
        terms, so what is left is the FRET.
        """
        donor_yield = self.quantum_yield
        if donor_yield <= 0.0:
            return 0.0
        k_total = self.k_quench + self.fret_rate_trace(
            acceptor, forster_radius, kappa2
        ).astype(np.float32)
        _delays, emitted = simulate_photon_trace(
            self.n_photons, k_total, t_step=self.diffusion.t_step,
            tau0=self.tau0, random_seed=self._photon_seed(),
        )
        paired_yield = float(emitted.sum()) / self.n_photons
        return 1.0 - paired_yield / donor_yield
