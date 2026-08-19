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
    from IMP.bff.io.structure import read_trajectory

    model = IMP.Model()
    hierarchy = IMP.atom.read_pdb(str(pdb_path), model, IMP.atom.AllPDBSelector())
    atom_names = [
        IMP.atom.Atom(leaf).get_atom_type().get_string().strip()
        for leaf in IMP.atom.get_leaves(hierarchy)
    ]

    coords = read_trajectory(dcd_path, n_atoms=len(atom_names),
                             max_frames=max_frames)
    if coords.shape[0] == 0:
        raise ValueError(f"No frames found in {dcd_path}")

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
    dcd = lib_dir / f"{dye_name}_cutoff{cutoff}.bcif"
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



# ``GridDiffusionResult``, ``GridDiffusionGradient`` and
# ``diffusion_stability_limit`` are **C++**, and re-exported here.
#
# The two records moved because the methods that build them did. The stability
# limit moved because it is the rule the solver *enforces*, and a check written
# on one side of the boundary against a scheme run on the other is how the two
# drift -- it already had, once: the rate term was missing from it while the
# scheme it guards has always carried ``k dt`` in the same coefficient, and on
# T4L site 19 that let a decay reach 7e36 while the check called the step safe.
from IMP.bff import (  # noqa: F401
    GridDiffusionGradient, GridDiffusionResult, diffusion_stability_limit,
)

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


# ``GridDiffusionSolver`` is **C++**, and the two per-step kernels went with
# it: ``_step`` and ``_step_smoluchowski`` existed only to be bound as
# ``self._kernel``, and nothing outside the class ever called either.
#
# The solver held four ``ng^3`` grids and folded ``dt`` into two of them on
# every call -- ``d = D dt/dg^2`` and ``decay = exp(-k dt)`` -- then handed all
# four across the boundary. At ``ng = 41`` that is 276 000 doubles per call,
# and ``equilibrium()`` makes up to 200 of them.
#
# The move is not a speed change: the C++ loop already dominated, and ``run``
# measures the same on either side. It is a *shape* change -- the solver is now
# an object IMP owns rather than a Python class holding numpy arrays.
#
# Gated against the Python before it was deleted, on a 25^3 grid with a
# spatially varying ``D`` and ``k``: the reported decay to 6e-16, the final
# density to 1e-19, ``equilibrium`` to 6e-19, and the three adjoint gradients
# to 2e-14 relative -- the last on ``d_density``, whose normalisation term is a
# dot product and so sums in a different order.
#
# See ``include/IMP/bff/GridDiffusionSolver.h`` and
# ``pyext/IMP_bff.griddiffusion.i``.
from IMP.bff import GridDiffusionSolver  # noqa: F401


# --------------------------------------------------------------------------
# Optimizer states that write frames
#
# `IMP.OptimizerState` subclasses that `bin/imp_bff flexfit` attaches to write
# a trajectory as it samples. Sampling machinery, so they sit with the rest of
# it now that the command that uses them lives in `bin/`.
# --------------------------------------------------------------------------


class WriteRMFFrame(IMP.OptimizerState):

    def __init__(
            self,
            filename,
            root_hier: IMP.atom.Hierarchy,
            restraints: typing.List[IMP.Restraint],
            name: str = "WriteRMFFrame"
    ):
        model = root_hier.get_model()
        super().__init__(model, name)
        self.restraints = restraints
        fileName, fileExtension = os.path.splitext(filename)
        fn = pathlib.Path(fileName + ".0.rmf3")
        if pathlib.Path(fn).exists():
            for i in range(100):
                fn = pathlib.Path(fileName + ".{:d}.rmf3".format(i))
                if not fn.exists():
                    break
        self._rmf_filename = str(fn)
        self._rh = RMF.create_rmf_file(str(fn))
        IMP.rmf.add_hierarchies(self._rh, [root_hier])
        IMP.rmf.add_restraints(self._rh, restraints)
        IMP.rmf.save_frame(self._rh)

    def do_update(self, arg0):
        #print(*[r.evaluate(False) for r in self.restraints], sep="\t")
        IMP.rmf.save_frame(self._rh)



class WritePDBFrame(IMP.OptimizerState):

    def __init__(
            self,
            filename,
            root_hier: IMP.atom.Hierarchy,
            restraints: typing.List[IMP.Restraint],
            restraint_filename: str = None,
            name: str = "WriteDCDFrame",
            output_objects: typing.List = None,
            multi_state: bool = False
    ):
        model = root_hier.get_model()
        super().__init__(model, name)

        import os
        self._pdb_basename = filename
        if restraint_filename is None:
            restraint_filename = os.path.splitext(filename)[0] + ".rst.txt"
        self._restraint_filename = restraint_filename
        self.restraints = restraints
        self._hier = root_hier
        self.frame = 0
        self.output_objects = output_objects
        self.multi_state = multi_state

    def do_update(self, arg0):
        with open(self._restraint_filename, "a+") as fp:
            fp.write("%s\t" % self.frame)
            fp.write("\t".join(["{:.3f}".format(r.evaluate(False)) for r in self.restraints]))
            fp.write("\t")
            if isinstance(self.output_objects, list):
                for obj in self.output_objects:
                    fp.write(str(obj) + "\t")
            fp.write("\n")
        lead = os.path.splitext(self._pdb_basename)[0]
        if not self.multi_state:
            hiers = [self._hier]
        else:
            hiers = self._hier.get_children()
        for i, hier in enumerate(hiers):
            out_fn = lead + "_state_" + str(i) + "_" + "{:04d}".format(self.frame) + ".pdb"
            IMP.atom.write_pdb(hier, out=out_fn)
        self.frame += 1

