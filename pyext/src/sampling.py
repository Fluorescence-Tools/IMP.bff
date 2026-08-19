"""Stage 3 — how configurations are drawn, enumerated, or propagated.

A representation says *where* the dye can be. Sampling says how you get at it,
and the method follows from the representation:

======================  ==============================================
representation          how it is sampled
======================  ==============================================
accessible volume       pathfinding on the grid, then a walk in it
rotamer library         screening the listed conformers
coarse-grained dye      molecular dynamics under the force field
======================  ==============================================

* :mod:`~IMP.bff.sampling` — one rejection-sampled trajectory at a
  time, which resolves the dye's *history*.
* :mod:`~IMP.bff.sampling` — the same dynamics as a density,
  ``dp/dt = div(D grad p) - k p``. A Langevin trajectory and its Fokker-Planck
  density are **one model** sampled or integrated, not two rival ones; they must
  agree on ``D``, and ``test/dynamics`` asserts they do.
* :mod:`~IMP.bff.sampling` — kinetic Monte Carlo over an
  excitation's fate, which is a different variable being sampled: not where the
  dye is, but whether the photon escapes.
* :mod:`~IMP.bff.sampling` — loading a library and drawing from
  it by weight.

This package was ``dynamics/`` until the layout was reconciled with the four
stages. "Dynamics" named the integrators and left library screening and
pathfinding without a home; sampling covers all three representations, which is
the point of having a stage rather than a package per technique.

The explicit dye's MD samplers are still in ``cgdye/sampling`` — that package is
off the domain layout, and moving its samplers here without the model they
sample would split it in half.
"""

from __future__ import annotations

from pathlib import Path
from typing import NamedTuple, Optional, Tuple
import random

import numpy as np

import IMP
import IMP.atom
import IMP.bff
import IMP.core

__all__ = [
    'DyeDiffusionTrajectory',
    'GridDiffusionGradient',
    'GridDiffusionResult',
    'GridDiffusionSolver',
    'apply_rotamer_coordinates',
    'diffusion_stability_limit',
    'equilibrium_occupancy',
    'load_rotamer_library_dcd',
    'sample_rotamer_index',
    'simulate_dye_diffusion',
    'simulate_photon_trace',
    'simulate_quenched_decay',
]

# --------------------------------------------------------------------------
# brownian
# --------------------------------------------------------------------------
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
    :param D: diffusion coefficient in A^2/ns, in the standard convention --
        the per-Cartesian-component step variance is ``2 D dt``, so
        ``<dx^2> = 2 D t``, the same D that :class:`GridDiffusionSolver` takes.
    :param slow_fact: either a scalar factor applied inside *slow_density*, or a
        ``(ng, ng, ng)`` per-voxel factor grid from
        :func:`IMP.bff.slow_factor_grid`.
    :param random_seed: seed for a reproducible walk; ``None`` draws freely.

    The returned coordinates are relative to the grid anchor.

    .. note::
       **The default is a free-solution value, and it is entered correctly.**
       ``D = 40 A^2/ns`` is 400 um^2/s (1 A^2/ns = 10 um^2/s), which is free
       Alexa488 in water -- and this parameter is documented, back to QuEst, as
       the dye's diffusion *in solution, when not interacting with the surface*.
       Stickiness is applied separately, through the mobility field.

       That matters because the step width was wrong until 2026-08-18:
       ``sqrt(2*D*3*dt)`` per Cartesian component, the total three-dimensional
       MSD used as one component's width, so the walk diffused at 3D. The
       tempting inference -- that ``D`` had been calibrated around it and must
       move by 3x in compensation -- is **false**. The factor of three is a
       transcription error, and its source is still in the docstring it came
       from: QuEst's 2019 ``simulate_traj_point`` cites a Berkeley teaching page
       whose ``k = sqrt(D * dimensions * tau)`` is the *magnitude* of a step
       that is then given a random direction, not a per-component sigma. No
       calibration of ``D`` exists in any repository in this stack; every
       document that mentions these parameters calls them uncalibrated
       transferable starting values.

       So the width was corrected and ``D`` was left alone. See
       ``okf/validation/particle_vs_field_diffusion.md``.
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

    counts = IMP.bff.VectorInt()
    # The grids go through numpy's own buffer: the SWIG typemap takes
    # (pointer, length) from the array, so nothing is copied or converted. The
    # dtypes below are therefore load-bearing -- int32 and float64 are what the
    # typemaps accept, and anything else would force a conversion right back.
    flat = IMP.bff.brownian_walk_in_volume(
        occupancy.ravel(), mobility.ravel(), ng, float(dg), float(t_max),
        float(t_step), float(D), seed, counts)

    # Four columns: x, y, z, accepted. Already a numpy array -- the kernel
    # hands back a view over its own buffer, so there is nothing to convert. A
    # returned std::vector would cost ~35-40 ns per element, which at a million
    # steps is 140 ms; an out-parameter would cost ~340 ns.
    packed = flat.reshape(-1, 4)
    # Strided views, not copies. `ascontiguousarray` here duplicated three
    # million doubles on a million-step walk -- 54 ms of a 200 ms call -- and
    # every consumer either does arithmetic that produces a fresh contiguous
    # array anyway (`trajectory + x0`) or passes through a numpy typemap that
    # handles a stride.
    xyz = packed[:, :3]
    accepted = packed[:, 3].astype(np.uint8)
    n_acc, n_rej = (int(counts[0]), int(counts[1])) if len(counts) == 2 else (0, 0)
    if xyz.size == 0:
        # No accessible starting voxel: an empty trajectory of the right length.
        n_steps = int(t_max / t_step)
        xyz = np.zeros((n_steps, 3), dtype=np.float64)
        accepted = np.zeros(n_steps, dtype=np.uint8)
    return DyeDiffusionTrajectory(xyz, accepted, n_acc, n_rej)


# --------------------------------------------------------------------------
# excited_state
# --------------------------------------------------------------------------
"""Photon and decay-curve Monte-Carlo against a per-frame quenching rate.

Two entry points, both driven by a quenching rate sampled along a trajectory:
:func:`simulate_photon_trace` samples individual photons (giving a photon
histogram with the counting statistics of a real measurement), and
:func:`simulate_quenched_decay` integrates the excited-state population directly
(giving a smooth curve for the same model, much faster, no shot noise).

Moved here from QuEst (``quest/core/photon.py``) by PRD-109.
"""

def simulate_photon_trace(
    n_ph: int,
    k_quench,
    t_step: float = 0.01,
    tau0: float = 0.25,
    random_seed: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Monte-Carlo photon trace against a per-frame quenching rate.

    :param n_ph: number of excitation events to simulate.
    :param k_quench: ``(n_frames,)`` quenching rate in 1/ns along the trajectory.
    :param t_step: trajectory time step in ns.
    :param tau0: intrinsic (unquenched) lifetime in ns.
    :param random_seed: makes the trace reproducible, at the cost of running
        single-threaded (see :func:`_trace_seeded`).
    :returns: ``(delay_times, emitted)`` -- the sampled delay time of each event
        and a ``uint8`` flag saying whether a photon was emitted or the dye was
        quenched first. Delay time is 0 for quenched events.
    """
    flat = IMP.bff.photon_trace(
        int(n_ph), np.ascontiguousarray(np.asarray(k_quench, dtype=np.float64).ravel()),
        float(t_step), float(tau0),
        -1 if random_seed is None else int(random_seed))
    # Two columns: delay, emitted. The flag comes back inside the returned array
    # because a SWIG out-parameter is a wrapper numpy walks one element at a
    # time -- 13.7 ms of a 152 ms trace, measured.
    packed = np.asarray(flat, dtype=np.float64).reshape(-1, 2)
    return (np.ascontiguousarray(packed[:, 0]), packed[:, 1].astype(np.uint8))


def simulate_quenched_decay(
    n_curves: int,
    decay: np.ndarray,
    dt_tac: float,
    k_quench,
    t_step: float,
    tau0: float,
    random_seed: Optional[int] = None,
) -> None:
    """Fill *decay* in place with a trajectory-driven fluorescence decay curve.

    For each of *n_curves* virtual trajectories, a random starting frame is
    picked along *k_quench* and the excited-state population is propagated in
    steps of *t_step* under a total rate ``1/tau0 + k_quench[frame]``, with the
    emitted intensity accumulated into TAC bins of width *dt_tac*.

    The result is the smooth decay of the same model
    :func:`simulate_photon_trace` samples photon-by-photon -- no shot noise, and
    far cheaper when only the curve is wanted.

    :param n_curves: number of virtual decay curves to average over.
    :param decay: 1-D array filled in place.
    :param dt_tac: TAC bin width in ns.
    :param k_quench: ``(n_frames,)`` quenching rate in 1/ns per trajectory frame.
    :param t_step: trajectory time step in ns.
    :param tau0: intrinsic excited-state lifetime in ns.
    :param random_seed: makes the curve reproducible. The numba this replaced
        had no seed at all -- the starting frames were drawn from the global
        state -- so the function could not be pinned in a test.
    """
    decay += np.asarray(IMP.bff.quenched_decay(
        int(n_curves), int(np.asarray(decay).shape[0]), float(dt_tac),
        np.ascontiguousarray(np.asarray(k_quench, dtype=np.float64).ravel()),
        float(t_step), float(tau0),
        -1 if random_seed is None else int(random_seed)), dtype=np.float64)


# --------------------------------------------------------------------------
# rotamer_library
# --------------------------------------------------------------------------
"""Loading a rotamer library, and drawing a state from it.

The mechanics of the discrete-states representation: read the library off disk
(DCD frames plus weights), put one of its conformers onto a hierarchy, and draw
an index according to the weights.

Was ``cgdye/sampling/rotamer.py`` until the 2026-08-18 cleanup. It sat with the
force-field samplers because it uses the same file formats, but a rotamer
library is not something you *sample from a potential* -- it is a list of states
with weights already attached, which is a representation.
"""

def load_rotamer_library_dcd(pdb_path, dcd_path, weights_path=None, max_frames=None):
    """Load a reference rotamer library from PDB+DCD (+ optional weights).

    Atom names come from the PDB through IMP.atom and coordinates from the DCD
    through the in-tree reader, so this needs nothing beyond IMP and numpy.
    """
    from IMP.bff.io.structure import read_dcd

    model = IMP.Model()
    hierarchy = IMP.atom.read_pdb(str(pdb_path), model, IMP.atom.AllPDBSelector())
    atom_names = [
        IMP.atom.Atom(leaf).get_atom_type().get_string().strip()
        for leaf in IMP.atom.get_leaves(hierarchy)
    ]

    coords = read_dcd(dcd_path, max_frames=max_frames)

    if coords.shape[0] == 0:
        raise ValueError(f"No frames found in DCD: {dcd_path}")

    if weights_path is not None and Path(weights_path).exists():
        w = []
        with open(weights_path) as fh:
            for line in fh:
                txt = line.strip()
                if not txt:
                    continue
                w.append(float(txt))
        weights = np.asarray(w, dtype=float)
        if max_frames is not None:
            weights = weights[:max_frames]
        if len(weights) != coords.shape[0]:
            n = min(len(weights), coords.shape[0])
            coords = coords[:n]
            weights = weights[:n]
    else:
        weights = np.ones(coords.shape[0], dtype=float)

    total = float(weights.sum())
    if total <= 0:
        weights = np.ones_like(weights)
        total = float(weights.sum())
    weights = weights / total

    return {
        "coords": coords,
        "weights": weights,
        "atom_names": atom_names,
    }


def find_reference_rotamer_files(lib_dir, dye_name, cutoff=30):
    """Return (pdb, dcd, weights) paths for a reference dye+linker name."""
    lib_dir = Path(lib_dir)
    pdb = lib_dir / f"{dye_name}.pdb"
    dcd = lib_dir / f"{dye_name}_cutoff{cutoff}.dcd"
    weights = lib_dir / f"{dye_name}_cutoff{cutoff}_weights.txt"
    if not pdb.exists() or not dcd.exists():
        raise FileNotFoundError(f"Missing required reference files for {dye_name}")
    return pdb, dcd, weights if weights.exists() else None


def apply_rotamer_coordinates(dye_hier, coords):
    """Apply one rotamer coordinate set to an IMP dye hierarchy in-place."""
    atoms = list(IMP.atom.get_by_type(dye_hier, IMP.atom.ATOM_TYPE))
    arr = np.asarray(coords, dtype=float)
    if arr.shape[0] != len(atoms):
        raise ValueError(
            f"Atom count mismatch: coords={arr.shape[0]} hierarchy={len(atoms)}"
        )
    for a, c in zip(atoms, arr):
        IMP.core.XYZ(a).set_coordinates(
            IMP.algebra.Vector3D(float(c[0]), float(c[1]), float(c[2]))
        )


def sample_rotamer_index(weights, rng=None):
    """Sample a rotamer index from normalized weights."""
    if rng is None:
        rng = random
    idx = list(range(len(weights)))
    return rng.choices(idx, weights=weights, k=1)[0]


# --------------------------------------------------------------------------
# smoluchowski
# --------------------------------------------------------------------------
"""Explicit propagation of the dye's excited state on an AV grid.

The field counterpart of the Brownian walk: instead of sampling trajectories and
histogramming photons, integrate

    dp/dt = div(D grad p) - k p

on the accessible-volume grid, with ``D`` the mobility field
(:func:`IMP.bff.diffusion_coefficient_map`) and ``k`` the decay-rate field
(:func:`IMP.bff.quenching_rate_map`, plus a FRET map if there is an acceptor).
Deterministic and free of shot noise, at the cost of a stability limit on the
time step.

Two things it gives that the walk does not: the **equilibrium** occupancy of the
volume (propagate with ``k = 0`` until the distribution stops moving -- the
sticky regions fill up, which a uniform AV density does not represent), and the
donor decay as a smooth curve in one pass.

Moved here from ChiSurf (``chisurf/core/structure/av/functions.py``,
``DiffusionIterator``) by PRD-109. ChiSurf carried an OpenCL backend beside the
Python one; it is not brought across -- ``IMP.bff`` ships through conda-forge as
part of IMP and cannot take a ``pyopencl`` dependency, and the kernel here is
jitted.
"""

#: Python name -> the C++ enum.
_FLUX = {"smoluchowski": 0, "ito": 1}



class GridDiffusionGradient(NamedTuple):
    """What :meth:`GridDiffusionSolver.gradient` returns -- ``dL/dD``,
    ``dL/dk`` and ``dL/dp0`` per voxel, each ``(ng, ng, ng)``, in the units
    of the solver's inputs (A^2/ns, 1/ns, and the unnormalised initial
    density)."""

    d_diffusion: np.ndarray
    d_rate: np.ndarray
    d_density: np.ndarray


class GridDiffusionResult(NamedTuple):
    """Time axis, surviving excited-state fraction, and the final density."""

    time: np.ndarray
    fluorescence: np.ndarray
    density: np.ndarray


def diffusion_stability_limit(d_max: float, dg: float, k_max: float = 0.0) -> float:
    """Largest stable time step for the explicit scheme.

    The updated voxel keeps the coefficient ``1 - 6 D dt/dg^2 - k dt``, so both
    terms constrain the step:

        dt <= 1 / (6 D / dg^2 + k_max)

    With ``k_max = 0`` this is the familiar ``dt <= dg^2 / (6 D)``.

    **Two thresholds, and they differ by a factor of two.** The amplification of
    the most oscillatory mode is ``|1 - 6 D dt/dg^2 - k dt|``, so the scheme
    *diverges* only once that sum exceeds **2**. This function returns the
    stronger bound, where the sum reaches **1** and the coefficient first goes
    negative: between the two the answer stays bounded but the density can go
    negative, which for a probability density is not an acceptable answer either.
    The old ``dg^2/(6 D)`` was the same positivity bound for pure diffusion.

    **The rate term used to be left out, and it dominates.** Measured on T4L
    site 19 at 2.5 A with ``t_step = 0.5 / (6*25/dg^2)`` (a 25 A^2/ns ceiling on
    ``D``): diffusion contributes 0.16 to that coefficient and quenching 2.01 --
    twelve times more -- and the sum of 2.17 is past the divergence threshold,
    so the decay reached 7e36 while this function called the step safe. Worse,
    the divergence need not look like one: site 124, at 2.10, returned a smooth,
    finite, entirely plausible decay that was 2.6 % wrong at 25 ns.

    Since :class:`GridDiffusionSolver` now integrates the rate term
    exponentially it is unconditionally stable in ``k`` and only the diffusion
    term binds. The criterion is kept honest here for callers sizing a step for
    the plain scheme, and because a caller who passes ``k_max`` deserves an
    answer that accounts for it.
    """
    d_max = float(d_max)
    k_max = float(k_max)
    denominator = 6.0 * d_max / float(dg) ** 2 + k_max
    if denominator <= 0.0:
        return float("inf")
    return 1.0 / denominator


def equilibrium_occupancy(
    diffusion_map, bounds, flux_form: str = "smoluchowski"
) -> np.ndarray:
    """The stationary occupancy of the volume, in closed form.

    Which closed form depends on how the flux is discretised, and the two
    disagree completely:

    ``"smoluchowski"`` (default, and the physics)
        ``p_eq ∝ 1`` on the accessible domain -- **independent of D**.
    ``"ito"`` (the inherited ChiSurf behaviour)
        ``p_eq ∝ 1/D``, so the dye piles up wherever it moves slowly.

    **The default is Smoluchowski because equilibrium is thermodynamics and
    mobility is kinetics.** A dye slowed by friction near the protein surface,
    with no attractive interaction, must still be found uniformly across its
    accessible volume at equilibrium -- it merely takes longer to get around.
    Letting a friction field set the distribution asserts a potential that was
    never specified.

    This is also the form the field's own canonical treatment uses: the
    **Haas-Steinberg** equation for diffusion-modulated FRET,

        ∂N(r,t)/∂t = -[1/tau_D + k_T(r)] N + D d/dr [ p(r) d/dr ( N / p(r) ) ]

    is written with exactly this structure so that its stationary state is the
    *given* distance distribution ``p(r)`` for any ``D``. ``p(r)`` comes from the
    chain statistics -- here, from the accessible volume -- and ``D`` is a
    separate kinetic parameter fitted against it. QuEst's notebook
    ``04_diffusion_modulated_fret.ipynb`` integrates that equation, but with a
    constant ``D`` and a uniform ``p(r)``, where both discretisations coincide
    and neither is tested.

    A genuinely *sticky* dye -- one with an attractive interaction, not merely a
    slower one -- is represented by a non-uniform ``p_eq``, which is a separate
    physical input this signature is shaped to accept later. It is not the same
    thing as a mobility field, and conflating them is what the inherited
    ``"ito"`` behaviour did.

    The ``"ito"`` branch below documents what it is:

    :func:`_step` propagates the flux as ``d[i]*p[i] - d[j]*p[j]``, which
    discretises

        ∂p/∂t = ∇²(D p)                    (the Itô / divergence form)

    **not** ``∇·(D ∇p)``. The stationary state of ``∇²(Dp)`` is ``D p = const``,
    i.e. ``p_eq ∝ 1/D`` -- so the dye accumulates where it moves slowly. It is
    not a small effect: on T4L site 132 the ratio of peak to mean occupancy is
    ~75. Verified against the iterative solver to a maximum relative deviation
    of **1.1e-13** (``test_quenching_field.py``).

    Either way, use this rather than :meth:`GridDiffusionSolver.equilibrium`,
    which spends tens of thousands of iterations converging to it -- and in the
    ``"ito"`` form on a real site, with ``D`` varying by orders of magnitude
    through the compounding slow factor, may not converge in any reasonable
    number at all.
    """
    diffusion_map = np.ascontiguousarray(diffusion_map, dtype=np.float64)
    mask = np.asarray(bounds) > 0
    occupancy = np.zeros_like(diffusion_map)
    if flux_form == "smoluchowski":
        usable = mask
    elif flux_form == "ito":
        usable = mask & (diffusion_map > 0.0)
    else:
        raise ValueError(
            f"flux_form must be 'smoluchowski' or 'ito', not {flux_form!r}")
    if not usable.any():
        return occupancy
    if flux_form == "smoluchowski":
        occupancy[usable] = 1.0
    else:
        occupancy[usable] = 1.0 / diffusion_map[usable]
    total = occupancy.sum()
    return occupancy / total if total else occupancy


def _step_smoluchowski(nxt, cur, d, k, bounds):
    """One explicit Euler step of the **Smoluchowski** form. **C++.**

    The flux between two voxels is ``D_ij (p_i - p_j)`` with ``D_ij`` the
    interface mobility, rather than ``D_i p_i - D_j p_j``. The difference is the
    whole physics: this flux vanishes when ``p`` is uniform *whatever* ``D``
    does in space, so a spatially varying mobility changes how fast the dye
    redistributes and **not where it ends up**.

    *k* carries ``exp(-k dt)``, applied as a factor: that is the exact solution
    of ``dp/dt = -k p`` over the step, so the rate contributes no stability
    constraint. ``1 - k dt`` goes negative and diverges once ``k dt > 1``.
    """
    ng = int(cur.shape[0])
    out = np.asarray(IMP.bff.diffusion_step(
        np.ascontiguousarray(cur, dtype=np.float64).ravel(),
        np.ascontiguousarray(d, dtype=np.float64).ravel(),
        np.ascontiguousarray(k, dtype=np.float64).ravel(),
        np.ascontiguousarray(bounds, dtype=np.float64).ravel(),
        ng, IMP.bff.FLUX_SMOLUCHOWSKI), dtype=np.float64).reshape(ng, ng, ng)
    nxt[...] = out
    return nxt


def _step(nxt, cur, d, k, bounds):
    """One explicit Euler step of the inherited **Ito** form. **C++.**

    Flux ``D_i p_i - D_j p_j``, whose stationary state is ``p ~ 1/D`` -- a
    friction field acting as an attractive potential. Kept so the old behaviour
    can be reproduced and compared; see :func:`equilibrium_occupancy`.

    *k* carries ``exp(-k dt)``, applied as a factor: that is the exact solution
    of ``dp/dt = -k p`` over the step, so the rate contributes no stability
    constraint. ``1 - k dt`` goes negative and diverges once ``k dt > 1``.
    """
    ng = int(cur.shape[0])
    out = np.asarray(IMP.bff.diffusion_step(
        np.ascontiguousarray(cur, dtype=np.float64).ravel(),
        np.ascontiguousarray(d, dtype=np.float64).ravel(),
        np.ascontiguousarray(k, dtype=np.float64).ravel(),
        np.ascontiguousarray(bounds, dtype=np.float64).ravel(),
        ng, IMP.bff.FLUX_ITO), dtype=np.float64).reshape(ng, ng, ng)
    nxt[...] = out
    return nxt


class GridDiffusionSolver:
    """Propagate an excited-state density on a masked 3-D grid.

    :param diffusion_map: ``D`` per voxel, in A^2/ns.
    :param bounds: non-zero where the dye may be; the domain mask.
    :param density: initial population; normalised on the first step.
    :param rate_map: ``k`` per voxel, in 1/ns. Zero (the default) propagates
        pure diffusion, which is how the equilibrium occupancy is obtained.
    :param t_step: integration step in ns. Must satisfy
        :func:`diffusion_stability_limit`; a larger one raises rather than
        returning a diverged field.
    :param dg: voxel edge in Angstrom.
    """

    def __init__(
        self,
        diffusion_map,
        bounds,
        density,
        rate_map=None,
        t_step: float = 1.0,
        dg: float = 1.0,
        check_stability: bool = True,
        flux_form: str = "smoluchowski",
    ):
        self.diffusion_map = np.ascontiguousarray(diffusion_map, dtype=np.float64)
        self.bounds = np.ascontiguousarray(
            np.asarray(bounds) > 0, dtype=np.float64
        )
        self.density = np.ascontiguousarray(density, dtype=np.float64)
        if rate_map is None:
            rate_map = np.zeros_like(self.diffusion_map)
        self.rate_map = np.ascontiguousarray(rate_map, dtype=np.float64)
        self.dg = float(dg)
        self.n_iterations = 0
        self.check_stability = bool(check_stability)
        self.t_step = float(t_step)
        if flux_form not in ("smoluchowski", "ito"):
            raise ValueError(
                f"flux_form must be 'smoluchowski' or 'ito', not {flux_form!r}")
        #: How the flux between voxels is discretised, which decides where the
        #: dye sits at equilibrium. ``"smoluchowski"`` keeps the equilibrium
        #: independent of ``D`` (the Haas-Steinberg structure); ``"ito"`` is the
        #: inherited ChiSurf behaviour, whose equilibrium is ``p ∝ 1/D``. They
        #: coincide when ``D`` is uniform.
        self.flux_form = flux_form
        self._kernel = _step if flux_form == "ito" else _step_smoluchowski

    def _validate(self):
        # The stencil cannot reach the outer shell, so a domain touching it
        # would lose population there with no warning.
        edge = self.bounds.copy()
        edge[1:-1, 1:-1, 1:-1] = 0.0
        if edge.any():
            raise ValueError(
                "The domain reaches the outer shell of the grid, where the "
                "7-point stencil cannot be evaluated and population is "
                "discarded. Enlarge the grid so the accessible volume is "
                "surrounded by at least one empty voxel."
            )
        if not self.check_stability:
            return
        limit = diffusion_stability_limit(float(self.diffusion_map.max()), self.dg)
        if self.t_step > limit:
            raise ValueError(
                f"t_step {self.t_step:g} ns exceeds the explicit-scheme "
                f"stability limit dg^2 / (6 D_max) = {limit:g} ns for "
                f"dg = {self.dg:g} A and D_max = {self.diffusion_map.max():g} "
                "A^2/ns. The scheme would diverge rather than lose accuracy."
            )

    def run(self, n_steps: int, n_out: int = 10) -> GridDiffusionResult:
        """Integrate *n_steps* steps, reporting every *n_out*.

        :returns: the time axis, the surviving excited-state fraction
            ``sum(p) / sum(p_0)`` at each reported step -- the donor decay -- and
            the final density.
        """
        self._validate()
        n_steps = int(n_steps)
        n_out = max(1, int(n_out))

        cur = self.density * self.bounds
        total = cur.sum()
        if total > 0.0:
            cur = cur / total
        nxt = np.zeros_like(cur)

        # Fold dt into the coefficients so the inner loop is pure arithmetic.
        d = self.diffusion_map * (self.t_step / self.dg ** 2)
        # Exact over the step, and unconditionally stable: see `_step`.
        k = np.exp(-self.rate_map * self.t_step)

        ng = int(self.bounds.shape[0])
        n_reports = n_steps // n_out + 1
        time = np.arange(n_reports, dtype=np.float64) * self.t_step * n_out
        fluorescence = np.zeros(n_reports, dtype=np.float64)

        # **The loop runs in C++.** Calling the step kernel from a Python loop
        # crosses the binding once per step and allocates a grid each time; a
        # solve is 1000-10000 steps, and doing that ran 58x slower than keeping
        # the loop inside (1.33 ms/step against 0.023).
        fluo = IMP.bff.VectorDouble()
        flat = IMP.bff.diffusion_propagate(
            np.ascontiguousarray(cur, dtype=np.float64).ravel(),
            np.ascontiguousarray(d, dtype=np.float64).ravel(),
            np.ascontiguousarray(k, dtype=np.float64).ravel(),
            np.ascontiguousarray(self.bounds, dtype=np.float64).ravel(),
            int(ng), _FLUX[self.flux_form], int(n_steps), int(n_out), fluo)
        cur = np.asarray(flat, dtype=np.float64).reshape(ng, ng, ng)
        fluorescence = np.asarray(list(fluo), dtype=np.float64)[:n_reports]
        self.n_iterations += n_steps

        self.density = cur
        return GridDiffusionResult(time, fluorescence, cur)

    def gradient(self, dL_dF, n_steps: int, n_out: int = 10,
                 density=None) -> GridDiffusionGradient:
        """Gradient of a loss on the decay with respect to every voxel of
        ``D``, ``k`` and the initial density -- the adjoint of :meth:`run`.

        For ``L = sum_k dL_dF[k] * fluorescence[k]`` (any loss, once its
        sensitivity to each reported point is known -- a Poisson likelihood's
        residuals, a chi-square's ``2 (model - data) / sigma^2``), returns
        ``dL/dD``, ``dL/dk`` and ``dL/dp0`` for the run that :meth:`run` with
        the same ``n_steps``, ``n_out`` and starting density performs. All
        voxels at once, for about four times the cost of one :meth:`run`
        (measured 4.4x on a 41^3 grid: one forward with checkpoints, one
        forward re-run between checkpoints, and a memory-bound reverse sweep
        of the transposed stencil that reads five fields per neighbour where
        the forward reads three; see ``diffusion_propagate_adjoint``).
        Finite differences need one run *per parameter*.

        The gradient is of the *discrete* scheme actually run, so it agrees
        with a central difference of :meth:`run` to roundoff; that is also why
        the time step must not depend on the parameters being differentiated
        (choose it from the range of ``D`` you fit over, not the current
        value).

        :param dL_dF: sensitivity of the loss to each reported population,
            length ``n_steps // n_out + 1`` (as :meth:`run` reports).
        :param density: the starting density of the run being differentiated.
            Defaults to the solver's current density -- note that :meth:`run`
            replaces that with the final one, so for a completed run pass the
            density it started from.
        :returns: :class:`GridDiffusionGradient` -- ``d_density`` includes
            the normalisation :meth:`run` applies (``p0 = p b / sum(p b)``).
        """
        self._validate()
        n_steps = int(n_steps)
        n_out = max(1, int(n_out))
        n_reports = n_steps // n_out + 1
        dL_dF = np.ascontiguousarray(dL_dF, dtype=np.float64).ravel()
        if dL_dF.shape[0] != n_reports:
            raise ValueError(
                f"dL_dF has {dL_dF.shape[0]} entries; run(n_steps={n_steps}, "
                f"n_out={n_out}) reports {n_reports} points")

        p = self.density if density is None else np.asarray(density, dtype=np.float64)
        p = np.ascontiguousarray(p, dtype=np.float64)
        cur = p * self.bounds
        total = cur.sum()
        if total > 0.0:
            cur = cur / total

        dt = self.t_step
        d = self.diffusion_map * (dt / self.dg ** 2)
        decay = np.exp(-self.rate_map * dt)
        ng = int(self.bounds.shape[0])

        flat = IMP.bff.diffusion_propagate_adjoint(
            np.ascontiguousarray(cur).ravel(),
            np.ascontiguousarray(d).ravel(),
            np.ascontiguousarray(decay).ravel(),
            np.ascontiguousarray(self.bounds).ravel(),
            int(ng), _FLUX[self.flux_form], int(n_steps), int(n_out), dL_dF)
        g = np.asarray(flat, dtype=np.float64).reshape(3, ng, ng, ng)
        g_d, g_decay, g_cur = g[0], g[1], g[2]

        # Chain rule through the folding run() applies:
        #   d = D dt/dg^2                 -> dL/dD = dL/dd * dt/dg^2
        #   decay = exp(-k dt)            -> dL/dk = dL/ddecay * (-dt) * decay
        #   cur = p b / sum(p b)          -> dL/dp = b (g - <g, cur>) / total
        d_diffusion = g_d * (dt / self.dg ** 2)
        d_rate = g_decay * (-dt) * decay
        if total > 0.0:
            d_density = self.bounds * (g_cur - float(np.sum(g_cur * cur))) / total
        else:
            d_density = np.zeros_like(g_cur)
        return GridDiffusionGradient(d_diffusion, d_rate, d_density)

    def equilibrium(
        self, n_steps: int = 20000, tolerance: float = 1e-8, n_check: int = 100
    ) -> np.ndarray:
        """Propagate with no decay until the occupancy stops moving.

        **Prefer :func:`equilibrium_occupancy`**, which is the same answer in
        closed form (``p ∝ 1/D``, agreeing to 1.1e-13). This iterates toward it,
        which on a real site — where ``D`` spans orders of magnitude through the
        compounding slow factor — can fail to converge in any practical number
        of steps: on T4L site 132 the distribution was still drifting after
        40 000 iterations. Kept as the independent check that the closed form is
        right, and for the case where someone changes the flux discretisation
        and the closed form no longer holds.

        :returns: the normalised equilibrium density.
        """
        self._validate()
        cur = self.density * self.bounds
        total = cur.sum()
        if total > 0.0:
            cur = cur / total
        nxt = np.zeros_like(cur)
        zero_rate = np.ones_like(cur)   # exp(-0 * dt): no decay
        d = self.diffusion_map * (self.t_step / self.dg ** 2)

        # Also in C++, and for the same reason as `run`: this loop is the
        # single hottest thing in the package -- tens of thousands of steps,
        # and it was crossing the binding on every one of them.
        ng = int(self.bounds.shape[0])
        d_flat = np.ascontiguousarray(d, dtype=np.float64).ravel()
        one = np.ones(cur.size, dtype=np.float64)      # exp(-0 * dt): no decay
        b_flat = np.ascontiguousarray(self.bounds, dtype=np.float64).ravel()
        previous = cur.copy()
        remaining = int(n_steps)
        chunk = max(1, int(n_check))
        while remaining > 0:
            take = min(chunk, remaining)
            fluo = IMP.bff.VectorDouble()
            flat = IMP.bff.diffusion_propagate(
                np.ascontiguousarray(cur, dtype=np.float64).ravel(),
                d_flat, one, b_flat, ng, _FLUX[self.flux_form],
                take, take, fluo)
            cur = np.asarray(flat, dtype=np.float64).reshape(ng, ng, ng)
            self.n_iterations += take
            remaining -= take
            drift = float(np.abs(cur - previous).sum())
            if drift < tolerance:
                break
            previous = cur.copy()

        total = cur.sum()
        if total > 0.0:
            cur = cur / total
        self.density = cur
        return cur
