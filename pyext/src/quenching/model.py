"""Rate and mobility maps on an accessible-volume grid.

The *field* formulation of the quenched-dye model: instead of walking a dye and
reading rates along its trajectory (:mod:`IMP.bff.sampling`), assign
every voxel a diffusion coefficient and a decay rate, and propagate the excited
state on the grid (:mod:`IMP.bff.sampling`). Deterministic, no shot
noise, and the natural home for a *distance-dependent* quenching law.

Moved here from ChiSurf (``chisurf/core/structure/av/functions.py``) by PRD-109.

Two quenching laws, deliberately both kept:

:func:`quenching_rate_map`
    ``k(r) = 1/tau0 + sum_a kQ_a * exp(-(|r - r_a| - r_dye) / rC_a)`` -- a
    Dexter-like exponential in the dye-to-atom distance, per **atom**.
:func:`IMP.bff.quenching_rate_grid`
    A hard contact sphere per quencher: full rate inside ``quench_radius``,
    zero outside, rates adding where spheres overlap. This is QuEst's law.

They are different models, not two spellings of one. The exponential has no
sharp contact radius and reaches further at low amplitude; the step model is
what QuEst's published parameters were calibrated against.
"""

from __future__ import annotations

from .fret_trace import fret_rate_pair_trace, fret_rate_trace
from .grids import grid_center_index, quenching_rate_grid, slow_factor_grid
from .pet import DEFAULT_DYE_RADIUS, QUENCHER_ATOMS, normalize_amino_acid_quenching
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from typing import NamedTuple, Optional, Sequence
import os

import numpy as np

from IMP.bff.sampling import simulate_dye_diffusion, simulate_photon_trace
import IMP.bff

__all__ = [
    'DyeDiffusionSimulation',
    'MAX_PARALLEL_TRAJECTORIES',
    'QuenchedDonorDecay',
    'ResidueSites',
    'atomic_quenching_parameters',
    'diffusion_coefficient_map',
    'fret_rate_map',
    'grid_axis',
    'quench_radii_for_residues',
    'quenching_rate_map',
    'quenching_rates_for_residues',
    'radial_diffusion_map',
    'residue_sites',
    'slow_diffusion_near_atoms',
    'slow_factors_for_residues',
]

# --------------------------------------------------------------------------
# maps
# --------------------------------------------------------------------------
"""Rate and mobility maps on an accessible-volume grid.

The *field* formulation of the quenched-dye model: instead of walking a dye and
reading rates along its trajectory (:mod:`IMP.bff.sampling`), assign
every voxel a diffusion coefficient and a decay rate, and propagate the excited
state on the grid (:mod:`IMP.bff.sampling`). Deterministic, no shot
noise, and the natural home for a *distance-dependent* quenching law.

Moved here from ChiSurf (``chisurf/core/structure/av/functions.py``) by PRD-109.

Two quenching laws, deliberately both kept:

:func:`quenching_rate_map`
    ``k(r) = 1/tau0 + sum_a kQ_a * exp(-(|r - r_a| - r_dye) / rC_a)`` -- a
    Dexter-like exponential in the dye-to-atom distance, per **atom**.
:func:`IMP.bff.quenching_rate_grid`
    A hard contact sphere per quencher: full rate inside ``quench_radius``,
    zero outside, rates adding where spheres overlap. This is QuEst's law.

They are different models, not two spellings of one. The exponential has no
sharp contact radius and reaches further at low amplitude; the step model is
what QuEst's published parameters were calibrated against.
"""

def grid_axis(ng: int, dg: float) -> np.ndarray:
    """Voxel-centre offsets from the grid anchor, in Angstrom.

    ``(i - (ng - 1) // 2) * dg`` -- the same integer centre offset as
    :func:`IMP.bff.grid_center_index`, so a map built here indexes the way the
    Brownian walk samples it.
    """
    return (np.arange(int(ng), dtype=np.float64) - (int(ng) - 1) // 2) * float(dg)


def atomic_quenching_parameters(atoms, quencher):
    """Per-atom quenching rate ``kQ`` and characteristic distance ``rC``.

    :param atoms: structured array with ``res_name`` and ``atom_name`` fields.
    :param quencher: ``{res_name: {atom_name: (kQ, rC)}}``. Atoms not named in
        it do not quench.
    :returns: ``(kQ, rC)``, each ``(n_atoms,)`` float64, zero where absent.
    """
    n = len(atoms)
    kQ = np.zeros(n, dtype=np.float64)
    rC = np.zeros(n, dtype=np.float64)
    res_names = atoms["res_name"]
    atom_names = atoms["atom_name"]
    for i in range(n):
        by_atom = quencher.get(str(res_names[i]))
        if not by_atom:
            continue
        entry = by_atom.get(str(atom_names[i]))
        if entry is None:
            continue
        kQ[i], rC[i] = float(entry[0]), float(entry[1])
    return kQ, rC


# Every kernel below iterates the FULL grid, `range(ng)`.
#
# ChiSurf's originals ran `for ix in range(-npm, npm)` with `npm = (ng - 1) // 2`
# -- one short on every axis, so the outermost slab was never written. In
# `assign_diffusion_to_grid_1` the output was `np.empty_like(...)`, so that slab
# was **uninitialised memory** read back as a diffusion coefficient, and that
# variant is the one `DynamicAV` used. In the rate maps the output was zeroed,
# so the outer face silently reported no quenching and no FRET. For an even
# `ng` -- which is what the IMP AV path produces -- it missed two slabs, not
# one. Fixed on the move (PRD-109).


def _slow_near_atoms(d_map, density, axis, r0, atoms_xyz, min_distance_sq, factor):
    """Slow the mobility wherever the dye contacts an atom. **C++.**

    The factor is applied *once per contacting atom*, so it compounds: a voxel
    touching 45 atoms at 0.985 ends at 0.5, and at 0.9 at 9e-3. Only values
    within a whisker of 1.0 are meaningful -- a property of the model, not a
    taste.
    """
    ng = int(density.shape[0])
    return np.asarray(IMP.bff.slow_near_atoms(
        np.ascontiguousarray(d_map, dtype=np.float64).ravel(),
        np.ascontiguousarray(density, dtype=np.float64).ravel(),
        np.ascontiguousarray(axis, dtype=np.float64).ravel(),
        np.ascontiguousarray(r0, dtype=np.float64).ravel(),
        np.ascontiguousarray(atoms_xyz, dtype=np.float64).ravel(),
        float(min_distance_sq), float(factor)),
        dtype=np.float64).reshape(ng, ng, ng)


def slow_diffusion_near_atoms(
    d_map, density, r0, dg, atoms_xyz, min_distance, slow_factor
) -> np.ndarray:
    """Scale a diffusion map down once per atom in contact with each voxel.

    The more atoms crowd a voxel the slower the dye moves there; the factor is
    applied once per contacting atom, so it compounds. Voxels outside the
    accessible volume are set to zero -- the dye cannot be there, and a non-zero
    coefficient would let the solver leak population into the protein.

    J. Chem. Phys. 126, 044707 (2007), eq. (3).

    :param d_map: base diffusion coefficient per voxel (A^2/ns).
    :param density: binary occupancy of the accessible volume.
    :param r0: grid anchor (the attachment point).
    :param dg: voxel edge in Angstrom.
    :param atoms_xyz: ``(n, 3)`` atom coordinates.
    :param min_distance: contact distance, roughly dye radius + 2 * vdW.
    :param slow_factor: factor applied per contacting atom, in [0, 1].
    """
    density = np.ascontiguousarray(density, dtype=np.float64)
    ng = density.shape[0]
    return _slow_near_atoms(
        np.ascontiguousarray(d_map, dtype=np.float64),
        density,
        grid_axis(ng, dg),
        np.asarray(r0, dtype=np.float64),
        np.ascontiguousarray(atoms_xyz, dtype=np.float64),
        float(min_distance) ** 2,
        float(slow_factor),
    )


def radial_diffusion_map(density, dg, f) -> np.ndarray:
    """A diffusion map from a radial profile ``f(distance from the anchor)``.

    Used to give the dye a mobility that depends on how far the linker is
    stretched, before :func:`slow_diffusion_near_atoms` adds the local crowding.
    Vectorised, so *f* must accept an array.
    """
    density = np.asarray(density)
    ng = density.shape[0]
    axis = grid_axis(ng, dg)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    return np.asarray(f(np.sqrt(x * x + y * y + z * z)), dtype=np.float64)


def diffusion_coefficient_map(
    density, r0, dg, atoms_xyz, *, free_diffusion, min_distance, slow_factor,
    radial_profile=None,
) -> np.ndarray:
    """The mobility field: a base coefficient, slowed by nearby atoms.

    :param free_diffusion: the unhindered dye diffusion coefficient (A^2/ns),
        used where *radial_profile* is not given.
    :param radial_profile: optional callable of the distance from the anchor,
        replacing the constant base coefficient.
    """
    density = np.asarray(density)
    if radial_profile is None:
        base = np.full(density.shape, float(free_diffusion), dtype=np.float64)
    else:
        base = radial_diffusion_map(density, dg, radial_profile)
    return slow_diffusion_near_atoms(
        base, density, r0, dg, atoms_xyz, min_distance, slow_factor
    )


def _quenching_map(density, axis, r0, atoms_xyz, kQ, rC, dye_radius, inv_tau0):
    """PET rate field. **C++.**

    Voxels outside the accessible volume stay at **zero**, not at ``1/tau0``:
    the dye cannot be there, so it has no decay rate there.
    """
    ng = int(density.shape[0])
    return np.asarray(IMP.bff.quenching_map(
        np.ascontiguousarray(density, dtype=np.float64).ravel(),
        np.ascontiguousarray(axis, dtype=np.float64).ravel(),
        np.ascontiguousarray(r0, dtype=np.float64).ravel(),
        np.ascontiguousarray(atoms_xyz, dtype=np.float64).ravel(),
        np.ascontiguousarray(kQ, dtype=np.float64).ravel(),
        np.ascontiguousarray(rC, dtype=np.float64).ravel(),
        float(dye_radius), float(inv_tau0)),
        dtype=np.float64).reshape(ng, ng, ng)


def quenching_rate_map(
    density, r0, dg, atoms_xyz, kQ, rC, *, tau0, dye_radius
) -> np.ndarray:
    """Total decay rate per voxel: intrinsic plus exponential PET.

    ``k(r) = 1/tau0 + sum_a kQ_a * exp(-(|r - r_a| - r_dye) / rC_a)``

    The distance is measured from the dye **surface** (the atom distance less
    the dye radius), and the sum runs over every atom with a non-zero ``kQ`` and
    ``rC``. Unlike the contact-sphere law in
    :func:`IMP.bff.quenching_rate_grid` there is no cut-off: a distant atom
    contributes exponentially little rather than nothing.

    Voxels outside the accessible volume stay at zero -- **not** at ``1/tau0``.
    The solver reads this as the rate field on a domain masked by the same
    occupancy, so an unreachable voxel carries no rate at all.
    """
    density = np.ascontiguousarray(density, dtype=np.float64)
    ng = density.shape[0]
    return _quenching_map(
        density,
        grid_axis(ng, dg),
        np.asarray(r0, dtype=np.float64),
        np.ascontiguousarray(atoms_xyz, dtype=np.float64),
        np.ascontiguousarray(kQ, dtype=np.float64),
        np.ascontiguousarray(rC, dtype=np.float64),
        float(dye_radius),
        1.0 / float(tau0) if tau0 > 0.0 else 0.0,
    )


def _fret_map(density_d, density_a, axis_d, axis_a, r0_d, r0_a, r0_6, kf, step):
    """FRET rate field, donor volume against acceptor volume. **C++.**

    **The harmonic mean**: the acceptor-weighted mean transfer *time* is
    accumulated and inverted, which is the static limit -- the acceptor does not
    move within the donor's excited-state lifetime. ``fret_rate_trace``
    accumulates the arithmetic mean of *rates* along a trajectory, which is the
    fast-exchange limit. Different physics, not variants.
    """
    ng = int(density_d.shape[0])
    return np.asarray(IMP.bff.fret_map(
        np.ascontiguousarray(density_d, dtype=np.float64).ravel(),
        np.ascontiguousarray(density_a, dtype=np.float64).ravel(),
        np.ascontiguousarray(axis_d, dtype=np.float64).ravel(),
        np.ascontiguousarray(axis_a, dtype=np.float64).ravel(),
        np.ascontiguousarray(r0_d, dtype=np.float64).ravel(),
        np.ascontiguousarray(r0_a, dtype=np.float64).ravel(),
        float(r0_6), float(kf), int(step)),
        dtype=np.float64).reshape(ng, ng, ng)


def fret_rate_map(
    density_donor,
    density_acceptor,
    r0_donor,
    r0_acceptor,
    dg_donor,
    dg_acceptor,
    forster_radius: float,
    kf: float,
    acceptor_step: int = 2,
) -> np.ndarray:
    """An effective FRET rate for every donor voxel, from the acceptor cloud.

    A fixed donor position does not have *a* FRET rate: it has a distribution of
    them, one per accessible acceptor position, and the FRET-induced donor decay
    is ``sum_i x_i exp(-k_i t)``. This approximates that sum by a single
    exponential whose rate is the reciprocal of the **mean transfer time**::

        k_eff = 1 / <1 / k_i>_A

    which is the *harmonic* mean of the rates, not the arithmetic one. The two
    differ, and the choice matters: the harmonic mean is dominated by the slow
    (distant) acceptor positions, the arithmetic mean by the fast (close) ones.
    :func:`IMP.bff.fret_rate_trace` takes the arithmetic mean instead, because
    it models a fast-exchanging acceptor where the *rates* average.

    :param kf: the donor's radiative rate ``1/tau0``.
    :param acceptor_step: stride over the acceptor grid; the cost is the product
        of the two grids' sizes.
    """
    density_donor = np.ascontiguousarray(density_donor, dtype=np.float64)
    density_acceptor = np.ascontiguousarray(density_acceptor, dtype=np.float64)
    if float(kf) <= 0.0:
        raise ValueError("kf (the donor's radiative rate) must be positive.")
    return _fret_map(
        density_donor,
        density_acceptor,
        grid_axis(density_donor.shape[0], dg_donor),
        grid_axis(density_acceptor.shape[0], dg_acceptor),
        np.asarray(r0_donor, dtype=np.float64),
        np.asarray(r0_acceptor, dtype=np.float64),
        float(forster_radius) ** 6,
        float(kf),
        max(1, int(acceptor_step)),
    )


# --------------------------------------------------------------------------
# sites
# --------------------------------------------------------------------------
"""Per-residue geometry for the PET model: where a residue sticks and quenches.

A residue enters the quenching model at two different points, and they are not
the same point in space:

*slow centre*
    The coarse side-chain position -- CB, else CA -- that drives the unspecific
    stickiness. Stickiness is a bulk property of the side chain, so its centroid
    is the right handle.
*quench centre*
    The centroid of the residue's PET-active atoms. Electron transfer happens at
    the redox-active moiety, and for tryptophan the indole ring sits ~3.3 A from
    CB -- which matters when contact distances are a few Angstrom.

Moved here from QuEst (``quest/core/dye_diffusion.py``) by PRD-109.
"""

_FALLBACK = {
    "slow_factor": 1.0,
    "kQ": 0.0,
    "quench_radius": None,
    "quench_atoms": ("CB",),
}


class ResidueSites(NamedTuple):
    """One slow centre and one quench centre per residue, plus its type."""

    slow_centers: np.ndarray
    quench_centers: np.ndarray
    residue_names: list

    def __len__(self) -> int:
        return len(self.residue_names)


def _name(value) -> str:
    if isinstance(value, bytes):
        return value.decode("ascii", errors="ignore").strip().upper()
    return str(value).strip().upper()


def residue_sites(atoms, quenching_table=None) -> ResidueSites:
    """Group *atoms* by residue and locate its slow and quench centres.

    :param atoms: structured array with ``chain``, ``res_id``, ``res_name``,
        ``atom_name`` and ``coord`` fields.
    :param quenching_table: the per-residue interaction table, whose
        ``quench_atoms`` decide which atoms define each quench centre. The
        defaults from :data:`IMP.bff.QUENCHER_ATOMS` are used when omitted.

    Residues are keyed by ``(chain, res_id, res_name)``. **Keying on ``res_id``
    alone is wrong** and was a real defect in QuEst: residue numbers restart per
    chain, so in a homodimer every number occurs twice and two residues' atoms
    were folded into one centre.
    """
    table = normalize_amino_acid_quenching(quenching_table)

    def wanted_atoms(residue_name):
        params = table.get(residue_name)
        if params is None:
            return frozenset(QUENCHER_ATOMS.get(residue_name, ("CB",)))
        return frozenset(params.get("quench_atoms") or ("CB",))

    by_residue = OrderedDict()
    for index, atom in enumerate(atoms):
        key = (_name(atom["chain"]), int(atom["res_id"]), _name(atom["res_name"]))
        by_residue.setdefault(key, []).append(index)

    slow_centers = []
    quench_centers = []
    residue_names = []
    for (_chain, _res_id, res_name), indices in by_residue.items():
        block = atoms[np.asarray(indices, dtype=np.int64)]
        atom_names = [_name(n) for n in block["atom_name"]]
        coords = np.asarray(block["coord"], dtype=np.float64)
        residue_name = _name(res_name)

        if "CB" in atom_names:
            selected = atom_names.index("CB")
        elif "CA" in atom_names:
            selected = atom_names.index("CA")
        else:
            selected = 0
        slow_center = coords[selected]

        matched = [i for i, n in enumerate(atom_names) if n in wanted_atoms(residue_name)]
        quench_center = coords[matched].mean(axis=0) if matched else slow_center

        slow_centers.append(slow_center)
        quench_centers.append(quench_center)
        residue_names.append(residue_name)

    return ResidueSites(
        np.asarray(slow_centers, dtype=np.float64).reshape(-1, 3),
        np.asarray(quench_centers, dtype=np.float64).reshape(-1, 3),
        residue_names,
    )


def _lookup(table, residue_name):
    return table.get(_name(residue_name), _FALLBACK)


def slow_factors_for_residues(residue_names: Sequence[str], table) -> np.ndarray:
    """The stickiness factor of each residue, in ``residue_names`` order."""
    return np.asarray(
        [_lookup(table, r)["slow_factor"] for r in residue_names], dtype=np.float64
    )


def quenching_rates_for_residues(residue_names: Sequence[str], table) -> np.ndarray:
    """The quenching rate (1/ns) of each residue."""
    return np.asarray(
        [_lookup(table, r)["kQ"] for r in residue_names], dtype=np.float64
    )


def quench_radii_for_residues(
    residue_names: Sequence[str], table, critical_distance: float = 0.0
) -> np.ndarray:
    """The contact radius of each residue, inheriting *critical_distance*.

    A ``quench_radius`` of ``None`` in the table means "use the model-wide
    critical distance", which is how a project sets one radius for everything
    and overrides it per residue type where it matters.
    """
    default = float(critical_distance or 0.0)
    radii = []
    for residue in residue_names:
        value = _lookup(table, residue).get("quench_radius")
        radii.append(default if value is None else float(value))
    return np.asarray(radii, dtype=np.float64)


# --------------------------------------------------------------------------
# model
# --------------------------------------------------------------------------
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


# ``DyeDiffusionSimulation`` is **C++**.
#
# It held the occupancy grid, the mobility field, the rate map and the
# trajectory as numpy arrays and drove the walk kernel from a
# ``ThreadPoolExecutor``. The kernel was already across the boundary; the state
# was not, and ``sample_grid`` -- tens of millions of indexed reads along the
# trajectory -- was the single largest numpy cost left in the package.
#
# Gated bit-for-bit against the Python before it was deleted, for one walk and
# for four concatenated: trajectory, accepted/rejected counts, ``mean_position``,
# ``k_quench``, ``quenched`` and ``sample_grid`` all identical. That was only
# possible because the C++ derives its per-walk seeds with the same
# ``(base + i * 104729) % (2**31 - 1)`` the Python used; a different derivation
# would have forced a distributional comparison instead.
#
# ``sample_grid`` over 400 000 frames: 14.1 ms -> 1.22 ms.
#
# See ``include/IMP/bff/DyeDiffusion.h`` and ``pyext/IMP_bff.dyediffusion.i``.
from IMP.bff import DyeDiffusionSimulation  # noqa: F401


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

        self._sites: Optional[ResidueSites] = None
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
    def sites(self) -> ResidueSites:
        """The slow and quench centres of every residue in the structure."""
        if self._sites is None:
            self._sites = residue_sites(self.atoms, self.quenching_table)
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
            quench_radii_for_residues(
                names, self.quenching_table, self.critical_distance
            ),
            residue_sites.quench_centers,
            self.x0,
            quenching_rates_for_residues(names, self.quenching_table),
        )
        self._slow_factor_map = slow_factor_grid(
            density, ng, self.dg,
            np.full(len(names), self.slow_radius, dtype=np.float64),
            residue_sites.slow_centers,
            self.x0,
            slow_factors_for_residues(names, self.quenching_table),
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

    def photons_fused(self):
        """The walk, the rate along it, and the photon race -- in one call.

        Identical to running :meth:`simulate_diffusion` and then
        :meth:`simulate_photons`, and asserted to be so photon-for-photon in
        ``test/quenching/test_fused_decay.py``. What it skips is the
        **trajectory**: the three-call path hands a ``(n_frames, 3)`` array to
        Python, reads a rate map along it and hands the trace back, and at the
        default ``t_max`` that is 20 million doubles crossing the boundary to
        produce a few thousand photons. Nothing downstream of the decay wants
        the coordinates.

        Use the three calls when the trajectory *is* the point -- a correlation
        function, a visualisation, an inspection of the rate trace. Use this one
        inside a fitting loop.

        :returns: ``(delay_times, emitted)``, the same pair
            :attr:`photon_trace` gives.
        """
        # NOT `self.diffusion` -- that property runs the split walk on first
        # access, so reading it here would do the whole simulation twice and
        # make the fused path the slower one. Measured: 4.9 s against 1.3 s at
        # a million steps, until this was an explicit construction.
        if self._diffusion is None:
            self._diffusion = DyeDiffusionSimulation(
                self.density, self.dg, self.x0,
                slow_factor_map=self.slow_factor_map,
                quenching_rate_map=self.quenching_rate_map,
            )
        walk = self._diffusion
        D, t_step, t_max = self.diffusion_coefficient, self.t_step, self.t_max
        walk.t_step = float(t_step)
        seeds = _trajectory_seeds(self.random_seed,
                                  _resolve_parallel(self.n_trajectories))
        slow = (walk.slow_factor_map if walk.slow_factor_map is not None
                else self.slow_fact)
        mobility = (np.ascontiguousarray(np.asarray(slow, dtype=np.float64)).ravel()
                    if np.ndim(slow) == 3 else
                    (np.empty(0) if float(slow) == 1.0 or walk.slow_density is None
                     else np.where(np.asarray(walk.slow_density, dtype=bool),
                                   float(slow), 1.0).ravel()))
        rate_map = self.quenching_rate_map
        ng = int(np.asarray(walk.density).shape[0])

        stats = IMP.bff.VectorDouble()
        flat = IMP.bff.quenched_donor_photons(
            np.ascontiguousarray(np.asarray(walk.density, dtype=np.int32)).ravel(),
            mobility,
            np.ascontiguousarray(np.asarray(rate_map, dtype=np.float64)).ravel(),
            ng, float(walk.dg), float(t_max), float(t_step), float(D),
            IMP.bff.VectorInt([-1 if s is None else int(s) for s in seeds]),
            float(self.tau0), int(self.n_photons),
            -1 if self._photon_seed() is None else int(self._photon_seed()),
            stats)

        packed = np.asarray(flat, dtype=np.float64).reshape(-1, 2)
        self._photon_trace = (np.ascontiguousarray(packed[:, 0]),
                              packed[:, 1].astype(np.uint8))
        # `list()`, not indexing: this module sets SWIG's `kwargs` feature, and
        # the generated `__getitem__` will not take a positional index under it.
        summary = list(stats)
        if len(summary) == 5:
            # The fused kernel never materialises a trajectory, but it does
            # report the step counts, and the walk record is where callers look
            # for them.
            walk.set_step_counts(int(summary[1]), int(summary[2]))
            self._fused_stats = {
                "n_frames": int(summary[0]),
                "mean_k_quench": float(summary[3]),
                "collision_fraction": float(summary[4]),
            }
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
