"""PET quenching of a dye diffusing in its accessible volume.

The model QuEst was built on, moved here by PRD-109 so that ``IMP.bff`` owns it:
a dye sphere diffuses on an accessible-volume grid, is slowed near sticky
residues, and is quenched by photo-induced electron transfer whenever it comes
within contact distance of a redox-active side chain. The observable is the
donor's fluorescence decay, and comparing it with a measured lifetime is what
calibrates the accessible *contact* volume.

Three layers, each usable on its own:

``pet``, ``asa``
    The chemistry: which residues quench, how fast, where on the residue, and
    how buried it is.
``grids``, ``maps``, ``fret_trace``
    The fields: stamp the quenching rate, the mobility and the FRET rate onto
    the accessible-volume grid, and read them along a trajectory.
``model``
    :class:`DyeDiffusionSimulation` and :class:`QuenchedDonorDecay`, which put
    the layers together for one labelling site.

The **integrators moved out** in PRD-113 stage 7. The Brownian walk, the
Smoluchowski field solver and the excited-state Monte Carlo are now
:mod:`IMP.bff.sampling`: none of them is specific to PET quenching -- a Brownian
walk in a volume is a Brownian walk in a volume -- and filing a general
integrator under the first physics that used it is how it comes to look like a
detail of one model.

Not to be confused with :class:`IMP.bff.LangevinDyeSampler` (PRD-108), which
integrates an **explicit all-atom dye** under a force field. That is a different
model at a different scale; this one is a sphere on a grid, cheap enough to scan
every position in a protein.

Was a six-module subpackage. The six averaged 320 lines and only one of them
(``model``) imported any of the others, so the split bought no separation --
it was the file fragmentation the consolidation exists to remove. One module
of the same size as ``io/cif.py`` or ``scoring.py``, sectioned by what each
part was.
"""

from __future__ import annotations
import numpy as np
import IMP.bff
from IMP.bff import grid_center_index
from collections import OrderedDict
from IMP.bff.photophysics import kappa2_isotropic
from typing import Optional
from IMP.bff.sampling import GridDiffusionSolver, diffusion_stability_limit, equilibrium_occupancy
from typing import NamedTuple, Optional, Sequence
import os
from IMP.bff.sampling import simulate_photon_trace
from IMP.bff import DyeDiffusionSimulation

__all__ = [
    'DEFAULT_DYE_RADIUS',
    'DyeDiffusionSimulation',
    'DynamicAccessibleVolume',
    'MAX_KAPPA2',
    'MAX_PARALLEL_TRAJECTORIES',
    'PET_QUENCHING_REFERENCE',
    'QUENCHER_ATOMS',
    'QuenchedDonorDecay',
    'ResidueSites',
    'STANDARD_AMINO_ACID_RESIDUES',
    'amino_acid_quenching_defaults',
    'atomic_quenching_parameters',
    'av_contact_mask',
    'diffusion_coefficient_map',
    'fret_rate_map',
    'fret_rate_pair_trace',
    'fret_rate_trace',
    'grid_axis',
    'grid_center_index',
    'normalize_amino_acid_quenching',
    'quench_radii_for_residues',
    'quencher_atom_indices',
    'quencher_centers',
    'quenching_rate_grid',
    'quenching_rate_map',
    'quenching_rate_per_frame',
    'quenching_rates_for_residues',
    'radial_diffusion_map',
    'residue_sites',
    'slow_diffusion_near_atoms',
    'slow_factor_grid',
    'slow_factors_for_residues',
    'solvent_accessible_surface',
    'sphere_points',
]


# --------------------------------------------------------------------------
# grids -- Per-voxel stickiness and quenching-rate maps stamped onto an AV grid.
# --------------------------------------------------------------------------

def _center_grid_indices(rs, r0, dg, ng, radius):
    """Voxel coordinates and integer radii of the sphere centres. **C++.**

    Uses ``floor``, not truncation. ``int()`` truncates *toward zero*, so a
    centre on the negative side of ``r0`` would round up while every other map
    here rounds down -- the walk's occupancy test and the trajectory sampler
    both floor. The two disagreed by one voxel per axis for every centre with a
    negative offset, which is half the grid.
    """
    ix0, iy0, iz0, ridx = IMP.bff.center_grid_indices(
        np.ascontiguousarray(rs, dtype=np.float64).ravel(),
        np.ascontiguousarray(r0, dtype=np.float64).ravel(),
        float(dg), int(ng),
        np.ascontiguousarray(radius, dtype=np.float64).ravel())
    return (np.asarray(ix0, dtype=np.int64), np.asarray(iy0, dtype=np.int64),
            np.asarray(iz0, dtype=np.int64), np.asarray(ridx, dtype=np.int64))
def _stamp(density, ng, radius, rs, r0, dg, values, combine):
    """One stamping kernel for both combines. **C++.**

    Stickiness *multiplies* (identity 1) and rates *add* (identity 0) -- the
    only difference between what used to be two near-identical Python kernels.
    Rates add because parallel channels do.
    """
    out = IMP.bff.stamp_spheres(
        np.ascontiguousarray(density, dtype=np.float64).ravel(), int(ng),
        np.ascontiguousarray(radius, dtype=np.float64).ravel(),
        np.ascontiguousarray(rs, dtype=np.float64).ravel(),
        np.ascontiguousarray(r0, dtype=np.float64).ravel(), float(dg),
        np.ascontiguousarray(values, dtype=np.float64).ravel(), int(combine))
    return np.asarray(out, dtype=np.float64).reshape(int(ng), int(ng), int(ng))
def _slow_factor_grid(density, ng, slow_radius, rs, r0, dg, slow_fact):
    """Stickiness field: factors multiply."""
    return _stamp(density, ng, slow_radius, rs, r0, dg, slow_fact,
                  IMP.bff.GRID_COMBINE_MULTIPLY)
def _additive_factor_grid(density, ng, radius, rs, r0, dg, values):
    """Rate field: parallel channels add."""
    return _stamp(density, ng, radius, rs, r0, dg, values,
                  IMP.bff.GRID_COMBINE_ADD)
def slow_factor_grid(density, ng, dg, slow_radius, rs, r0, slow_fact):
    """Per-voxel diffusion scaling from overlapping sticky spheres.

    Stickiness **multiplies** where spheres overlap. Voxels outside the AV keep
    1.0, so the factor is only meaningful where the walk can go.

    :param density: ``(ng, ng, ng)`` binary occupancy of the accessible volume.
    :param dg: voxel edge in Angstrom.
    :param slow_radius: ``(n,)`` radius of each sticky sphere.
    :param rs: ``(n, 3)`` sphere centres.
    :param r0: the grid anchor (see :func:`grid_center_index`).
    :param slow_fact: ``(n,)`` factor in [0, 1] per centre.
    """
    return _slow_factor_grid(
        np.asarray(density, dtype=np.uint8),
        int(ng),
        np.asarray(slow_radius, dtype=np.float64),
        np.asarray(rs, dtype=np.float64),
        np.asarray(r0, dtype=np.float64),
        float(dg),
        np.asarray(slow_fact, dtype=np.float64),
    )
def quenching_rate_grid(density, ng, dg, radius, rs, r0, values):
    """Per-voxel quenching rate (1/ns) from overlapping quencher spheres.

    Rates **add up** where contact spheres overlap, which is the physical
    composition rule for independent PET channels. Same geometry and the same
    indexing convention as :func:`slow_factor_grid`; only the accumulator
    differs.
    """
    return _additive_factor_grid(
        np.asarray(density, dtype=np.uint8),
        int(ng),
        np.asarray(radius, dtype=np.float64),
        np.asarray(rs, dtype=np.float64),
        np.asarray(r0, dtype=np.float64),
        float(dg),
        np.asarray(values, dtype=np.float64),
    )
def av_contact_mask(density, ng, dg, slow_radius, rs, r0):
    """The contact ("slow") part of an AV density grid, as a binary mask.

    A thin adapter over :func:`IMP.bff.representation.av._kernels.split_av_acv`, which returns
    both parts plus their counts; for a binary density the contact part is what
    the PET model wants.
    """
    from IMP.bff.representation.av import split_av_acv

    density = np.asarray(density, dtype=np.uint8)
    ng = int(ng)
    _, _, contact, _ = split_av_acv(
        density.reshape(ng, ng, ng).astype(np.float64),
        float(dg),
        np.asarray(slow_radius, dtype=np.float64),
        np.asarray(rs, dtype=np.float64),
        np.asarray(r0, dtype=np.float64),
    )
    return np.ascontiguousarray(contact.reshape(density.shape), dtype=np.uint8)


# --------------------------------------------------------------------------
# asa -- Solvent-accessible surface area (Shrake-Rupley) for selected atoms.
# --------------------------------------------------------------------------

DEFAULT_SPHERE_POINTS = 590
def sphere_points(n: int) -> np.ndarray:
    """*n* roughly equidistant points on the unit sphere (golden spiral).

    **C++** (:file:`include/IMP/bff/SolventAccessibleSurface.h`). The Python it
    replaces built the spiral in ``float32``; this is ``float64`` and agrees
    with an independent float64 construction to 1e-15, where the old one was
    off by up to 2.4e-6. The port is the more accurate of the two.
    """
    return np.asarray(IMP.bff.sphere_points(int(n)),
                      dtype=np.float64).reshape(-1, 3)
def _asa(r, vdw, probe_atom_indices, points, probe, radius):
    """Shrake-Rupley area of the selected atoms. **C++.**

    The neighbour cutoff is ``(2*radius + probe)**2``, which is the tight
    criterion here: a sample sits ``radius`` from its own centre and is occluded
    by a neighbour within ``radius + probe`` of it. QuEst's kernel tested a
    *squared* distance against an *unsquared* sum, and against ``vdw`` rather
    than the sampling radius -- a 2.45 A cutoff where 6.0 A was meant, so nearly
    every occluding neighbour was missed. Fixed on the move (PRD-109); nothing
    downstream changed, because nothing consumed the result.
    """
    return np.asarray(
        IMP.bff.solvent_accessible_surface_area(
            np.ascontiguousarray(r, dtype=np.float64).ravel(),
            np.ascontiguousarray(vdw, dtype=np.float64).ravel(),
            [int(i) for i in np.asarray(probe_atom_indices).ravel()],
            np.ascontiguousarray(points, dtype=np.float64).ravel(),
            float(probe), float(radius)),
        dtype=np.float64)
def solvent_accessible_surface(
    xyz,
    vdw,
    probe_atom_indices,
    points=None,
    probe: float = 1.0,
    radius: float = 2.5,
) -> np.ndarray:
    """Accessible surface area (A^2) of each atom in *probe_atom_indices*.

    :param xyz: ``(n, 3)`` atom coordinates.
    :param vdw: ``(n,)`` van der Waals radii.
    :param probe_atom_indices: indices of the atoms to measure.
    :param points: unit-sphere sample points; :func:`sphere_points` by default.
    :param probe: probe radius.
    :param radius: radius at which the sphere points are placed.
    """
    xyz = np.asarray(xyz, dtype=np.float32)
    vdw = np.asarray(vdw, dtype=np.float32)
    probe_atom_indices = np.asarray(probe_atom_indices, dtype=np.uint32)
    if points is None or np.asarray(points).ndim != 2:
        points = sphere_points(DEFAULT_SPHERE_POINTS)
    points = np.asarray(points, dtype=np.float32)
    return _asa(xyz, vdw, probe_atom_indices, points, float(probe), float(radius))


# --------------------------------------------------------------------------
# pet -- Photo-induced electron transfer (PET) quenching of a tethered dye.
# --------------------------------------------------------------------------

STANDARD_AMINO_ACID_RESIDUES = (
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
)
QUENCHER_ATOMS = OrderedDict((
    ("ALA", ("CB",)),
    ("ARG", ("CZ", "NE", "NH1", "NH2")),
    ("ASN", ("CG", "OD1", "ND2")),
    ("ASP", ("CG", "OD1", "OD2")),
    ("CYS", ("SG",)),
    ("GLN", ("CD", "OE1", "NE2")),
    ("GLU", ("CD", "OE1", "OE2")),
    ("GLY", ("CA",)),
    ("HIS", ("CG", "ND1", "CD2", "CE1", "NE2")),
    ("ILE", ("CB",)),
    ("LEU", ("CB",)),
    ("LYS", ("NZ",)),
    ("MET", ("SD",)),
    ("PHE", ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")),
    ("PRO", ("N", "CB", "CG", "CD")),
    ("SER", ("OG",)),
    ("THR", ("OG1",)),
    ("TRP", ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2", "NE1", "CG", "CD1")),
    ("TYR", ("CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH")),
    ("VAL", ("CB",)),
))
DEFAULT_DYE_RADIUS = 3.5
PET_QUENCHING_REFERENCE = OrderedDict((
    ("TRP", {"kQ": 3.5, "contact_distance": 5.0}),
    ("TYR", {"kQ": 2.0, "contact_distance": 5.0}),
    ("MET", {"kQ": 1.67, "contact_distance": 3.5}),
    ("HIS", {"kQ": 1.0, "contact_distance": 4.7}),
    ("CYS", {"kQ": 0.8, "contact_distance": 3.5}),
    ("PRO", {"kQ": 2.0, "contact_distance": 4.0}),
))
_DEFAULT_TABLE = OrderedDict(
    (
        res,
        {
            "slow_factor": 1.0,
            "kQ": 0.0,
            # ``None`` means "inherit the model-wide critical distance".
            "quench_radius": None,
            "quench_atoms": list(QUENCHER_ATOMS[res]),
        },
    )
    for res in STANDARD_AMINO_ACID_RESIDUES
)
def _residue_name(value) -> str:
    if isinstance(value, bytes):
        return value.decode("ascii", errors="ignore").strip().upper()
    return str(value).strip().upper()
def _clamp_slow_factor(value) -> float:
    return min(1.0, max(0.0, float(value)))
def _quench_radius(value):
    """A positive quench radius, or ``None`` to inherit the global one."""
    if value is None:
        return None
    try:
        radius = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(radius) or radius <= 0.0:
        return None
    return radius
def _quench_atoms(value, default):
    """A de-duplicated list of atom names, falling back to *default*."""
    if value is None:
        return list(default)
    if isinstance(value, (str, bytes)):
        value = [value]
    atoms = []
    for atom in value:
        name = _residue_name(atom)
        if name and name not in atoms:
            atoms.append(name)
    return atoms or list(default)
def normalize_amino_acid_quenching(table=None):
    """The full per-residue interaction table, with defaults filled in.

    Accepts a partial table keyed by residue name, and the legacy shape where
    a bare number meant the slow factor. Unknown residue names are kept, so a
    non-standard residue can be given a rate.
    """
    normalized = OrderedDict(
        (res, dict(params, quench_atoms=list(params["quench_atoms"])))
        for res, params in _DEFAULT_TABLE.items()
    )
    if not table:
        return normalized
    for residue, params in table.items():
        name = _residue_name(residue)
        defaults = normalized.get(name) or {
            "slow_factor": 1.0,
            "kQ": 0.0,
            "quench_radius": None,
            "quench_atoms": list(QUENCHER_ATOMS.get(name, ("CB",))),
        }
        if isinstance(params, dict):
            slow_factor = params.get("slow_factor", defaults["slow_factor"])
            kQ = params.get("kQ", defaults["kQ"])
            quench_radius = params.get("quench_radius", defaults["quench_radius"])
            quench_atoms = params.get("quench_atoms", defaults["quench_atoms"])
        else:
            # Legacy format: a bare number meant the slow factor.
            slow_factor = params
            kQ = defaults["kQ"]
            quench_radius = defaults["quench_radius"]
            quench_atoms = defaults["quench_atoms"]
        normalized[name] = {
            "slow_factor": _clamp_slow_factor(slow_factor),
            "kQ": max(0.0, float(kQ)),
            "quench_radius": _quench_radius(quench_radius),
            "quench_atoms": _quench_atoms(
                quench_atoms, QUENCHER_ATOMS.get(name, ("CB",))
            ),
        }
    return normalized
def amino_acid_quenching_defaults(
    kQ_scale: float = 1.0,
    slow_factor: float = 1.0,
    dye_radius: float = DEFAULT_DYE_RADIUS,
):
    """A full interaction table built from :data:`PET_QUENCHING_REFERENCE`.

    :param kQ_scale:
        Dye-specific multiplier applied to every reference ``kQ``. Dyes that are
        harder to reduce or oxidise use a value below one.
    :param slow_factor:
        Diffusion scaling applied near every residue (unspecific stickiness).
    :param dye_radius:
        Radius of the dye sphere in Angstrom. The trajectory tracks the dye
        *centre*, so the reference surface contact distances are offset by this
        radius to give centre-to-centre quench radii.
    """
    radius = max(0.0, float(dye_radius))
    table = normalize_amino_acid_quenching()
    for residue, params in table.items():
        params["slow_factor"] = _clamp_slow_factor(slow_factor)
        reference = PET_QUENCHING_REFERENCE.get(residue)
        if reference is None:
            continue
        params["kQ"] = max(0.0, float(reference["kQ"]) * float(kQ_scale))
        params["quench_radius"] = _quench_radius(
            radius + float(reference["contact_distance"])
        )
    return table
def quencher_atom_indices(atoms, selection):
    """Atom indices per residue type, for a ``{residue: [atom names]}`` selection.

    *atoms* is a structured array with ``res_name`` and ``atom_name`` fields --
    the shape QuEst's structure reader and :mod:`IMP.bff.representation.av` both produce.
    Returns an ``OrderedDict`` keyed the same way as *selection*, with a
    ``uint32`` index array per residue type (possibly empty).
    """
    res_name = atoms["res_name"]
    atom_name = atoms["atom_name"]
    indices = OrderedDict()
    for residue in selection:
        found = [
            np.where((res_name == residue) & (atom_name == name))[0]
            for name in selection[residue]
        ]
        if found:
            indices[residue] = np.array(np.hstack(found), dtype=np.uint32)
        else:
            indices[residue] = np.array([], dtype=np.uint32)
    return indices
def quencher_centers(atoms, selection):
    """Quenching-centre coordinates per residue type.

    One row per selected atom, grouped by residue type -- the caller decides
    whether to average them into a per-residue centroid or to treat each atom as
    its own centre. Returns an ``OrderedDict`` of ``(n, 3)`` arrays.
    """
    indices = quencher_atom_indices(atoms, selection)
    coord = atoms["coord"]
    return OrderedDict(
        (residue, coord[indices[residue]]) for residue in selection
    )
def _rate_per_frame(collided, k_quench):
    """Total quenching rate per frame. **C++.**

    Sums the rate constants of the quenchers in contact at each frame -- rates
    add, because the channels are parallel.
    """
    collided = np.ascontiguousarray(collided)
    n_frames = int(collided.shape[0])
    return np.asarray(
        IMP.bff.quenching_rate_per_frame(
            [int(v) for v in collided.astype(np.int32).ravel()],
            n_frames,
            np.ascontiguousarray(k_quench, dtype=np.float64).ravel()),
        dtype=np.float64)
def quenching_rate_per_frame(collided, k_quench) -> np.ndarray:
    """Sum the rates of the quenching atoms the dye touched, frame by frame.

    The **per-atom** route to a quenching trace, as against sampling a stamped
    rate grid along the trajectory (:meth:`IMP.bff.DyeDiffusionSimulation.k_quench`).
    Use it when the contact flags are what you have -- from a distance
    calculation against explicit atoms rather than from a voxel map.

    Frames are independent, so each row reduces on its own thread and the
    ``(n_frames, n_atoms)`` product is never materialised -- for a long
    trajectory against a whole protein's quenching atoms that would be the
    largest array in the calculation.

    Moved here from ChiSurf (``chisurf/core/structure/av/dynamic.py``) by
    PRD-109, where it was a plain Python double loop whose docstring already
    claimed the threading it did not have.

    :param collided: ``(n_frames, n_atoms)`` flags, non-zero where the dye was
        within the critical distance of that atom in that frame.
    :param k_quench: ``(n_atoms,)`` quenching rate per atom, in 1/ns.
    :returns: ``(n_frames,)`` total quenching rate.
    """
    collided = np.ascontiguousarray(collided)
    k_quench = np.ascontiguousarray(k_quench, dtype=np.float64)
    if collided.ndim != 2:
        raise ValueError("`collided` must be (n_frames, n_atoms).")
    if collided.shape[1] != k_quench.shape[0]:
        raise ValueError(
            f"`collided` has {collided.shape[1]} atoms but `k_quench` has "
            f"{k_quench.shape[0]}."
        )
    return _rate_per_frame(collided.astype(np.uint8), k_quench)


# --------------------------------------------------------------------------
# fret_trace -- Time-resolved FRET rate along a dye trajectory.
# --------------------------------------------------------------------------

MAX_KAPPA2 = 4.0
def _rate_trace(trajectory, acceptor_points, R0, tau0, r_min2, kappa2_scale):
    """FRET rate per frame against a static acceptor cloud. **C++.**

    The **arithmetic mean of rates** over the cloud -- the fast-exchange limit,
    where the acceptor re-randomises within the donor's excited-state lifetime.
    ``fret_rate_map`` accumulates the mean transfer *time* and inverts it, which
    is the static limit. Different physics, not two spellings.
    """
    # Already a numpy array: the kernel hands back a view over its own buffer.
    # A trajectory runs to millions of frames, and a returned std::vector would
    # cost ~35-40 ns each to build and walk back.
    return IMP.bff.fret_rate_trace_kernel(
        np.ascontiguousarray(trajectory, dtype=np.float64).ravel(),
        np.ascontiguousarray(acceptor_points, dtype=np.float64).ravel(),
        float(R0), float(tau0), float(r_min2), float(kappa2_scale))
def _rate_pair_trace(donor, acceptor, R0, tau0, r_min2, kappa2_scale):
    """FRET rate per frame from two trajectories, paired frame by frame. **C++.**"""
    return IMP.bff.fret_rate_pair_trace_kernel(
        np.ascontiguousarray(donor, dtype=np.float64).ravel(),
        np.ascontiguousarray(acceptor, dtype=np.float64).ravel(),
        float(R0), float(tau0), float(r_min2), float(kappa2_scale))
def _kappa2_scale(kappa2) -> float:
    """Turn an orientation factor into a multiplier on ``R0**6``.

    A published Forster radius is quoted *at* the isotropic 2/3, so the transfer
    rate carries the ratio ``kappa2 / (2/3)`` -- equivalently the ``1.5 *
    kappa2`` of the literature -- and passing 2/3 leaves the rate as it was.

    **This model does not compute kappa2, and cannot**: its dye is a
    structureless sphere diffusing in a volume, with no transition dipole to
    orient. What this does is make the *assumption* explicit and adjustable
    instead of silently baked into R0. Modelling it needs a rotational degree of
    freedom per dye and a correlation time to go with it -- a change to the
    physics, not a call to a library. :func:`IMP.bff.kappa2_from_dipoles` does
    compute it, but it takes the two dipole vectors, which is exactly what a
    point in a volume has not got. Use :class:`IMP.bff.RotamerEnsemble` when
    dipoles matter.
    """
    value = float(kappa2_isotropic() if kappa2 is None else kappa2)
    if not np.isfinite(value) or value < 0.0 or value > MAX_KAPPA2:
        raise ValueError(
            f"kappa2 must be in [0, {MAX_KAPPA2:g}] (0 = perpendicular dipoles, "
            f"4 = collinear); got {value!r}."
        )
    return value / kappa2_isotropic()
def fret_rate_trace(
    trajectory,
    acceptor_points,
    R0: float,
    tau0: float,
    r_min: float = 7.0,
    max_acceptor_points: int = 512,
    kappa2=None,
) -> np.ndarray:
    """Per-frame FRET rate (1/ns) along a donor trajectory.

    For every donor position the transfer rate is averaged over the acceptor's
    accessible volume::

        k_FRET(t) = (1/tau0) * < (R0 / |r_D(t) - r_A|)^6 >_A * (kappa2 / (2/3))

    This resolves the **donor** in time -- the point of simulating its diffusion
    at all -- while treating the acceptor as sampling its own volume quickly
    compared with the donor's excited-state lifetime. That is the *fast-acceptor
    limit*, and it is not valid if the acceptor is immobilised; use
    :func:`fret_rate_pair_trace` then.

    :param trajectory: ``(n_frames, 3)`` donor positions in Angstrom.
    :param acceptor_points: ``(n, 3)`` acceptor AV cloud in Angstrom.
    :param R0: Forster radius in Angstrom.
    :param tau0: unquenched donor lifetime in ns.
    :param r_min: closest approach of the two dye centres, in Angstrom. Both
        clouds are centre positions and can overlap in space, where
        ``(R0/r)^6`` would diverge; two dye spheres cannot interpenetrate, so
        the distance is floored at the sum of their radii.
    :param max_acceptor_points: cap on the acceptor points averaged over, taken
        by an even stride. The average converges long before the full cloud, and
        the cost is ``n_frames * n_acceptor``.
    :param kappa2: orientation factor; the isotropic 2/3 by default.
    """
    trajectory = np.ascontiguousarray(trajectory, dtype=np.float64).reshape((-1, 3))
    points = np.ascontiguousarray(acceptor_points, dtype=np.float64).reshape((-1, 3))
    if points.shape[0] == 0:
        raise ValueError("Acceptor accessible volume is empty; cannot compute FRET rates.")
    if tau0 <= 0.0:
        raise ValueError("tau0 must be positive to derive a FRET rate.")
    if points.shape[0] > max_acceptor_points:
        stride = int(np.ceil(points.shape[0] / float(max_acceptor_points)))
        points = np.ascontiguousarray(points[::stride])
    r_min = max(float(r_min), 1e-6)
    return _rate_trace(
        trajectory, points, float(R0), float(tau0), r_min * r_min, _kappa2_scale(kappa2)
    )
def fret_rate_pair_trace(
    donor_trajectory,
    acceptor_trajectory,
    R0: float,
    tau0: float,
    r_min: float = 7.0,
    kappa2=None,
) -> np.ndarray:
    """Per-frame FRET rate (1/ns) from *two* trajectories.

    Both dyes are resolved in time::

        k_FRET(t) = (1/tau0) * (R0 / |r_D(t) - r_A(t)|)^6 * (kappa2 / (2/3))

    This is the general case. :func:`fret_rate_trace` averages the acceptor over
    its accessible volume instead, which is the fast-acceptor limit.

    The two do not agree, and the reason is **occupancy**, not the averaging
    order. A cloud is weighted uniformly -- every accessible voxel counts the
    same -- while a trajectory is weighted by where the dye actually spends its
    time, and a sticky dye dwells near the protein surface rather than out in
    the free volume. Which way that moves the efficiency depends on the
    geometry; it is not a fixed bias. Measured on 148l E15->E90 with a sticky
    acceptor (D = 2.5 A^2/ns, slow radius 11 A) it raised E from 0.73 to 0.78,
    because the sticky region there lies toward the donor.

    The two walks are independent, so pairing frame *i* with frame *i* samples
    the joint distribution correctly **provided both were simulated on the same
    time base** -- same ``t_step`` and the same number of frames. A mismatch
    raises rather than truncating: a trajectory is a *concatenation of one walk
    per excitation*, so cutting it at an arbitrary frame splits a walk and pairs
    the tail of one excitation against the head of another.

    **Unequal lengths are the normal case, not an error case.** The acceptor sits
    at a different site with different residues around it and is therefore
    quenched differently, so the two dyes do not need the same number of
    excitations to collect the same number of photons. Matching them is the
    *caller's* job and the lever is the **number of walks**, not ``t_max`` or
    ``t_step`` -- those are usually already shared and changing them will not
    help.

    :param donor_trajectory: ``(n_frames, 3)`` donor centres in Angstrom.
    :param acceptor_trajectory: ``(n_frames, 3)`` acceptor centres, frame-aligned.
    :param R0: Forster radius in Angstrom.
    :param tau0: unquenched donor lifetime in ns.
    :param r_min: closest approach of the two dye centres, in Angstrom.
    :param kappa2: orientation factor; the isotropic 2/3 by default.
    """
    donor = np.ascontiguousarray(donor_trajectory, dtype=np.float64).reshape((-1, 3))
    acceptor = np.ascontiguousarray(acceptor_trajectory, dtype=np.float64).reshape((-1, 3))
    if donor.shape[0] == 0 or acceptor.shape[0] == 0:
        raise ValueError("Both trajectories must have frames to pair.")
    if tau0 <= 0.0:
        raise ValueError("tau0 must be positive to derive a FRET rate.")
    if donor.shape[0] != acceptor.shape[0]:
        raise ValueError(
            "Donor and acceptor trajectories must have the same number of "
            f"frames to pair ({donor.shape[0]} vs {acceptor.shape[0]}). This is "
            "expected whenever the two dyes are quenched differently -- they sit "
            "at different sites -- so they need different numbers of excitations "
            "for the same photon count. Simulate the same number of *walks* for "
            "both dyes; t_max and t_step are usually already shared and are not "
            "the lever. Truncating here is refused because a trajectory "
            "concatenates one walk per excitation, so an arbitrary cut splits a "
            "walk and pairs one excitation's tail against another's head."
        )
    r_min = max(float(r_min), 1e-6)
    return _rate_pair_trace(
        donor, acceptor, float(R0), float(tau0), r_min * r_min, _kappa2_scale(kappa2)
    )


# --------------------------------------------------------------------------
# dynamic -- A dye whose accessible volume carries mobility, quenching and FRET fields.
# --------------------------------------------------------------------------

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
        flux_form: str = "smoluchowski",
    ):
        self.av = av
        self.atoms = atoms
        self.tau0 = float(tau0)
        self.dye_radius = float(dye_radius)
        self.free_diffusion = float(free_diffusion)
        self.contact_distance = float(contact_distance)
        self.slow_factor = float(slow_factor)
        if flux_form not in ("smoluchowski", "ito"):
            raise ValueError(
                f"flux_form must be 'smoluchowski' or 'ito', not {flux_form!r}")
        #: Where the dye sits at equilibrium. See
        #: :func:`IMP.bff.sampling.smoluchowski.equilibrium_occupancy`.
        self.flux_form = flux_form

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
        self._diffusion_map = diffusion_coefficient_map(
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
        kQ, rC_atoms = atomic_quenching_parameters(self.atoms, quencher)
        if rC is not None:
            rC_atoms = np.where(kQ > 0.0, float(rC), 0.0)
        self._quenching_rate_map = quenching_rate_map(
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
        self._fret_rate_map = fret_rate_map(
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
            d_map, self.bounds, density, rate_map, t_step=t_step, dg=self.dg,
            flux_form=self.flux_form,
        )

    @property
    def occupancy(self) -> np.ndarray:
        """The equilibrium distribution of the dye, normalised.

        Under the default ``flux_form="smoluchowski"`` this **is** the AV
        density, uniform over the accessible voxels, and does not depend on the
        mobility -- equilibrium is thermodynamics, mobility is kinetics, and a
        dye slowed by friction with no attraction is still found everywhere it
        can reach. Under ``"ito"`` it is ``p ∝ 1/D``: see
        :func:`IMP.bff.sampling.smoluchowski.equilibrium_occupancy` for why that is
        the inherited behaviour rather than the physics.
        """
        if self._occupancy is None:
            self.update_occupancy()
        return self._occupancy

    def update_occupancy(self, t_step: Optional[float] = None, **kwargs) -> np.ndarray:
        """The equilibrium occupancy, in closed form.

        Uniform on the accessible domain under the default flux form,
        ``p ∝ 1/D`` under ``"ito"`` — see
        :func:`IMP.bff.sampling.smoluchowski.equilibrium_occupancy`. Propagating to it
        instead is possible but slow and, in the ``"ito"`` form on a real site
        where the compounding slow factor makes ``D`` span orders of magnitude,
        may not converge at all: on T4L site 132 it was still drifting after
        40 000 iterations. Pass ``iterate=True`` to do it the long way anyway.
        """
        if kwargs.pop("iterate", False):
            self._occupancy = self._solver(
                self.bounds, None, t_step).equilibrium(**kwargs)
        else:
            self._occupancy = equilibrium_occupancy(
                self.diffusion_map, self.bounds, self.flux_form)
        return self._occupancy

    def donor_decay(
        self, t_max: float = 50.0, t_step: Optional[float] = None, n_out: int = 10
    ):
        """Integrate the donor decay on the grid, from the equilibrium start.

        :returns: a :class:`IMP.bff.sampling.smoluchowski.GridDiffusionResult` --
            time axis in ns, surviving excited-state fraction, final density.
        """
        start = self.occupancy
        solver = self._solver(start, self.rate_map, t_step)
        n_steps = max(1, int(t_max / solver.t_step))
        return solver.run(n_steps, n_out=n_out)


# --------------------------------------------------------------------------
# model -- Rate and mobility maps on an accessible-volume grid.
# --------------------------------------------------------------------------

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
