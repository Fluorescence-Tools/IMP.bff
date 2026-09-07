"""PRD-110 stage 0: is the quenching model's parameter vector identifiable?

    <arm64>/bin/python benchmark/quenching_identifiability.py [--resolution 2.0]

**No network, on purpose.** PRD-110 proposes a physics-informed surrogate whose
whole justification is that calibrating the PET parameters against a measured
decay costs hundreds of forward solves. Before building anything, this asks the
question that decides whether the surrogate is worth having:

    can theta be recovered from a fluorescence decay at all?

If two parameters act on ``F(t)`` in nearly the same way, the fit is
under-determined, and a faster forward model fits an under-determined problem
faster. That is not progress. ``kQ`` and the stickiness factor are the pair to
worry about: both work by keeping the dye near a quencher, one by raising the
rate there and the other by making it dwell.

What it measures, on a real site:

1. **Sensitivity.** ``J = dF/dtheta`` by central differences, in units of the
   noise, so the columns are comparable.
2. **Identifiability.** The eigenvalues of the Fisher information ``J^T J`` and
   the condition number. A tiny eigenvalue means a *direction* in parameter
   space the decay cannot see, and its eigenvector says which combination.
3. **The baseline the surrogate has to beat.** A real least-squares fit to a
   synthetic decay at known parameters: evaluations, wall-clock, and whether
   the known values come back.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

import IMP
import IMP.atom
import IMP.bff
import IMP.core
import IMP.bff as maps
from IMP.bff import (
    GridDiffusionSolver,
    diffusion_stability_limit,
    equilibrium_occupancy,
)

#: Reference site. T4 lysozyme, the structure the AV suite already uses.
STRUCTURE = "structure/T4L/3GUN.pdb"
SITE = dict(chain="A", residue=132, atom="CB")

#: The parameters under test, with the ranges a calibration would search.
#: `tau0` is excluded on purpose: it is measured independently on free dye, and
#: including a parameter that is known from elsewhere would flatter the result.
#: **`slow_factor` is only meaningful within a whisker of 1.0**, and that is a
#: property of the model, not a taste. It is applied *once per contacting atom*
#: and therefore compounds: on T4L site 132 at a 6.5 A contact distance a voxel
#: sees up to 45 atoms, so
#:
#:     0.985^45 = 0.5      0.95^45 = 0.10      0.9^45 = 9e-3      0.6^45 = 1e-10
#:
#: At 0.6 the mobility spans ten orders of magnitude across the volume, and
#: because the equilibrium is `p ∝ 1/D` the dye ends up entirely inside a
#: frozen contact shell — at which point the decay becomes *exactly* independent
#: of `free_diffusion` (measured: 9.5e-9 over an 80x range of D). ChiSurf's
#: default of 0.985 is the sane end of that, and the range below reflects it.
PARAMETERS = (
    # name              start   low    high   unit
    ("free_diffusion",    8.0,   0.5,  40.0, "A^2/ns"),
    ("slow_factor",     0.985,   0.90,  1.0, "-"),
    ("contact_distance",  6.5,   3.0,  12.0, "A"),
    ("kQ_scale",          1.0,   0.05, 10.0, "-"),
    ("rC",                1.5,   0.3,   5.0, "A"),
)

#: Photons in a decent TCSPC decay; sets the noise the sensitivity is measured
#: against. Poisson at the peak, so sigma/F ~ 1/sqrt(N_peak).
N_PEAK_COUNTS = 10_000

TAU0 = 4.0
T_MAX = 25.0
N_TIME = 64


def load_atoms(pdb_path):
    """The obstacle atoms, in the field layout the quenching model reads."""
    model = IMP.Model()
    hierarchy = IMP.atom.read_pdb(
        pdb_path, model, IMP.atom.NonWaterNonHydrogenPDBSelector())
    leaves = IMP.atom.get_leaves(hierarchy)
    dtype = [("chain", "U4"), ("res_id", "i8"), ("res_name", "U4"),
             ("atom_name", "U4"), ("coord", "f8", 3), ("radius", "f8")]
    rows = []
    for leaf in leaves:
        atom = IMP.atom.Atom(leaf)
        residue = IMP.atom.get_residue(atom)
        chain = IMP.atom.get_chain(residue)
        rows.append((
            chain.get_id(), residue.get_index(),
            residue.get_residue_type().get_string(),
            atom.get_atom_type().get_string().strip(),
            IMP.core.XYZ(leaf).get_coordinates(),
            IMP.core.XYZR(leaf).get_radius(),
        ))
    return np.array(rows, dtype=dtype)


def quencher_table(kQ_scale: float, rC: float):
    """`{comp_id: PETParameters}`, the way `atomic_quenching_parameters` wants.

    Built from the same reference chemistry the contact-sphere model uses, so
    the two formulations are calibrated against comparable numbers.
    """
    return IMP.bff.reference_pet_parameters(
        rate_scale=kQ_scale, attenuation_length=rC)


class Site:
    """One labelling site, ready to be evaluated at many parameter vectors."""

    def __init__(self, pdb_path, resolution: float, site=None,
                 linker_length: float = 20.0,
                 flux_form: str = "smoluchowski"):
        site = SITE if site is None else site
        self.site = site
        self.atoms = load_atoms(pdb_path)
        # `disc_step` is the grid resolution, and it is the *only* way in.
        # `get_av` overwrites `source_info["simulation_grid_resolution"]`
        # from this argument (`fret/av.py:167`) rather than reading it, so
        # passing the resolution in `source_info` is silently ignored -- which
        # is how every run before 2026-08-18 was made at the 1.5 A default while
        # reporting whatever `--resolution` said. Asserted below rather than
        # trusted.
        self.av = IMP.bff.get_av(
            np.zeros((1, 4)), np.zeros(3), linker_length, 0.5, (3.5, 3.5, 3.5),
            disc_step=resolution,
            pdb_path=pdb_path,
            source_info={
                "chain_identifier": site["chain"],
                "residue_seq_number": site["residue"],
                "atom_name": site["atom"],
                "simulation_type": "AV1",
                "linker_length": linker_length, "linker_width": 0.5,
                "radius1": 3.5,
                "allowed_sphere_radius": 2.1,
            },
        )
        self.density = np.ascontiguousarray(self.av.density, dtype=np.float64)
        self.dg = float(self.av.grid_step)
        if abs(self.dg - float(resolution)) > 1e-9:
            raise RuntimeError(
                f"asked for a {resolution} A grid and got {self.dg} A -- the "
                "resolution argument is not reaching the AV builder.")
        self.x0 = np.asarray(self.av.attachment_point, dtype=np.float64)
        self.bounds = (self.density > 0).astype(np.float64)
        self.xyz = np.ascontiguousarray(self.atoms["coord"], dtype=np.float64)

        # The stencil cannot reach the outer shell; an AV that touches it would
        # lose population silently, so the solver raises. Pad instead.
        edge = self.bounds.copy()
        edge[1:-1, 1:-1, 1:-1] = 0.0
        if edge.any():
            self._pad()

        self.linker_length = float(linker_length)
        self.flux_form = flux_form
        self.time = np.linspace(0.0, T_MAX, N_TIME)
        self.n_evaluations = 0

        #: **One time step for every theta.** Deriving it from `D_max` — the
        #: obvious thing, since the stability limit is `dg^2/(6 D_max)` — makes
        #: the integrator's discretisation error a function of the parameter
        #: being differentiated, and that error lands in the derivative.
        #: Measured: `|dF/dD|` drifts 0.0273 -> 0.0173 across finite-difference
        #: steps of 2-20 % with an adaptive step (58 %), against 0.01452 ->
        #: 0.01479 (1.9 %) with this one. The adaptive step *inflated* the
        #: sensitivity to D by about a factor of two.
        self.t_step = 0.5 * diffusion_stability_limit(
            max(p[3] for p in PARAMETERS if p[0] == "free_diffusion"), self.dg)

    def _pad(self):
        """Grow the grid by one voxel a side so the AV clears the boundary."""
        for name in ("density", "bounds"):
            grid = getattr(self, name)
            padded = np.zeros(tuple(n + 2 for n in grid.shape), dtype=np.float64)
            padded[1:-1, 1:-1, 1:-1] = grid
            setattr(self, name, np.ascontiguousarray(padded))

    def decay(self, theta) -> np.ndarray:
        """The donor decay `F(t)` at a parameter vector, normalised to F(0)=1.

        Two solves: the equilibrium occupancy (which depends on the mobility
        field, so it moves with theta) and the decay itself.
        """
        free_diffusion, slow_factor, contact_distance, kQ_scale, rC = theta
        self.n_evaluations += 1

        d_map = maps.diffusion_coefficient_map(
            self.density, self.x0, self.dg, self.xyz,
            free_diffusion=free_diffusion,
            min_distance=contact_distance,
            slow_factor=slow_factor,
        )
        kQ, rC_atoms = maps.atomic_quenching_parameters(
            self.atoms["res_name"], self.atoms["atom_name"],
            quencher_table(kQ_scale, rC))
        rate = maps.quenching_rate_map(
            self.density, self.x0, self.dg, self.xyz, kQ, rC_atoms,
            tau0=TAU0, probe_radius=3.5,
        )

        # Closed form, `p ∝ 1/D`. Iterating for it is not merely slower: on
        # this site the compounding slow factor makes `D` span orders of
        # magnitude and the iteration was still drifting after 40 000 steps, so
        # every gradient measured through it was a gradient of an unconverged
        # field.
        equilibrium = equilibrium_occupancy(d_map, self.bounds, self.flux_form)

        solver = GridDiffusionSolver(
            d_map, self.bounds, equilibrium, rate,
            t_step=self.t_step, dg=self.dg, flux_form=self.flux_form)
        result = solver.run(max(1, int(T_MAX / self.t_step)), n_out=64)

        curve = np.interp(self.time, result.time, result.fluorescence)
        return curve / curve[0] if curve[0] > 0 else curve


def jacobian(site, theta, relative_step=0.05):
    """`dF/dtheta` by central differences, one column per parameter."""
    theta = np.asarray(theta, dtype=np.float64)
    columns = []
    for i, (name, _s, low, high, _u) in enumerate(PARAMETERS):
        step = relative_step * max(abs(theta[i]), 1e-3)
        step = min(step, 0.25 * (high - low))
        up, down = theta.copy(), theta.copy()
        up[i] = min(theta[i] + step, high)
        down[i] = max(theta[i] - step, low)
        f_up, f_down = site.decay(up), site.decay(down)
        columns.append((f_up - f_down) / (up[i] - down[i]))
        print(f"    d/d{name:<18} |grad| = {np.linalg.norm(columns[-1]):.4g}",
              flush=True)
    return np.column_stack(columns)


def analyse(site, theta, sigma):
    """Fisher information of the decay with respect to theta."""
    print("  Jacobian (central differences):", flush=True)
    J = jacobian(site, theta) / sigma[:, None]

    # Scale to relative parameter changes so the units do not decide the answer.
    scale = np.array([max(abs(t), 1e-6) for t in theta])
    J_rel = J * scale[None, :]

    fisher = J_rel.T @ J_rel
    eigenvalues, eigenvectors = np.linalg.eigh(fisher)
    order = np.argsort(eigenvalues)[::-1]
    return J_rel, eigenvalues[order], eigenvectors[:, order]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolution", type=float, default=2.0,
                        help="AV grid resolution in Angstrom (default 2.0).")
    parser.add_argument("--fit", action="store_true",
                        help="also run the least-squares baseline (slow).")
    parser.add_argument("--out", type=Path, default=None,
                        help="write the numbers as JSON.")
    args = parser.parse_args()

    pdb_path = IMP.bff.get_example_path(STRUCTURE)
    print(f"site: {STRUCTURE} {SITE['chain']}{SITE['residue']}.{SITE['atom']}, "
          f"resolution {args.resolution} A", flush=True)

    site = Site(pdb_path, args.resolution)
    print(f"  grid {site.bounds.shape}, {int(site.bounds.sum())} accessible voxels",
          flush=True)

    theta0 = np.array([p[1] for p in PARAMETERS], dtype=np.float64)
    names = [p[0] for p in PARAMETERS]

    start = time.perf_counter()
    reference = site.decay(theta0)
    one_evaluation = time.perf_counter() - start
    print(f"  one forward evaluation: {one_evaluation:.2f} s "
          f"(equilibrium + decay)", flush=True)

    # Poisson counting noise on a decay peaking at N_PEAK_COUNTS.
    counts = reference * N_PEAK_COUNTS
    sigma = np.sqrt(np.maximum(counts, 1.0)) / N_PEAK_COUNTS

    J_rel, eigenvalues, eigenvectors = analyse(site, theta0, sigma)

    print("\n  Fisher information (relative parameter changes):", flush=True)
    print(f"    eigenvalues: {np.array2string(eigenvalues, precision=3)}")
    condition = eigenvalues[0] / max(eigenvalues[-1], 1e-300)
    print(f"    condition number: {condition:.4g}")

    print("\n  eigen-directions, strongest first "
          "(what the decay measures, and what it cannot see):")
    for rank in range(len(names)):
        vector = eigenvectors[:, rank]
        terms = " ".join(
            f"{w:+.2f}*{n}" for n, w in
            sorted(zip(names, vector), key=lambda kv: -abs(kv[1]))
            if abs(w) > 0.15
        )
        print(f"    lambda={eigenvalues[rank]:10.3g}  {terms}")
    weakest = eigenvectors[:, -1]

    print("\n  per-parameter sensitivity (sigma of a 1-parameter fit, relative):")
    try:
        covariance = np.linalg.inv(J_rel.T @ J_rel)
        errors = np.sqrt(np.diag(covariance))
        for name, error in zip(names, errors):
            verdict = "identifiable" if error < 0.5 else "POORLY DETERMINED"
            print(f"    {name:<20} +/- {error * 100:7.1f} %   {verdict}")
    except np.linalg.LinAlgError:
        errors = None
        print("    Fisher matrix is singular: theta is NOT identifiable.")

    payload = {
        "structure": STRUCTURE, "site": SITE,
        "resolution": args.resolution,
        "grid_shape": list(site.bounds.shape),
        "accessible_voxels": int(site.bounds.sum()),
        "seconds_per_evaluation": one_evaluation,
        "parameters": names,
        "theta0": theta0.tolist(),
        "fisher_eigenvalues": eigenvalues.tolist(),
        "condition_number": float(condition),
        "weakest_direction": dict(zip(names, weakest.tolist())),
        "relative_errors": None if errors is None else errors.tolist(),
    }

    if args.fit:
        payload["fit"] = run_fit(site, theta0, reference, sigma)

    if args.out:
        args.out.write_text(json.dumps(payload, indent=2) + "\n")
        print(f"\nwrote {args.out}")
    return 0


def run_fit(site, theta_true, reference, sigma):
    """The baseline PRD-110 stage 3 has to beat: a fit through the solver."""
    from scipy.optimize import least_squares

    rng = np.random.default_rng(0)
    observed = reference + rng.normal(0.0, sigma)

    start_theta = np.array([
        theta_true[i] * (1.6 if i % 2 == 0 else 0.6)
        for i in range(len(theta_true))
    ])
    bounds = (np.array([p[2] for p in PARAMETERS]),
              np.array([p[3] for p in PARAMETERS]))
    start_theta = np.clip(start_theta, bounds[0], bounds[1])

    print("\n  least-squares baseline (through the solver):", flush=True)
    site.n_evaluations = 0
    t0 = time.perf_counter()
    result = least_squares(
        lambda theta: (site.decay(theta) - observed) / sigma,
        start_theta, bounds=bounds, diff_step=0.05, xtol=1e-6,
    )
    elapsed = time.perf_counter() - t0

    print(f"    evaluations: {site.n_evaluations}, wall clock: {elapsed:.1f} s")
    recovery = {}
    for i, (name, *_rest) in enumerate(PARAMETERS):
        true, found = theta_true[i], result.x[i]
        error = (found - true) / true * 100.0
        recovery[name] = {"true": true, "found": found, "percent": error}
        print(f"    {name:<20} true {true:8.3f}  found {found:8.3f} "
              f"({error:+7.1f} %)")
    return {
        "evaluations": site.n_evaluations,
        "seconds": elapsed,
        "cost": float(result.cost),
        "recovery": recovery,
    }


if __name__ == "__main__":
    raise SystemExit(main())
