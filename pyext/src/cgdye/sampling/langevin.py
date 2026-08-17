"""Langevin / Brownian dynamics of an explicit dye on a labelled site (PRD-108).

Real stochastic dynamics -- forces from the dye force field (bonds, angles,
dihedrals, repulsive LJ lower bounds; :func:`topology.combined.dye_forcefield_system`),
a soft-sphere repulsion against the protein's heavy atoms around the site,
and either

* ``integrator="md"``: velocity-Verlet molecular dynamics with a Langevin
  thermostat (``IMP.atom.MolecularDynamics`` + ``LangevinThermostatOptimizerState``,
  time step in fs, friction in 1/ps), or
* ``integrator="bd"``: overdamped Brownian dynamics (``IMP.atom.BrownianDynamics``,
  per-atom Stokes–Einstein diffusion coefficients from the vdW radius).

The dye's backbone anchor atoms (``N CA C O`` of the dye's own residue, which
``attach_dyes`` placed on the protein backbone) stay fixed; everything else
moves. Both integrators sample the Boltzmann distribution of the same energy
function -- the thermodynamic pins in ``test/cgdye/test_langevin_sampler.py``
(harmonic-oscillator temperature, free-particle MSD, a torsion's Boltzmann
histogram, md ⇄ bd agreement) are what "real Langevin" means here. The
collision-gated random walk (``sample-dof-walk``) is a different, cheaper
thing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable, Optional, Sequence

import numpy as np

import IMP
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.core

KB_KCAL = 0.0019872041  # kcal/(mol K)
_VDW_RADIUS = {"H": 1.1, "C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8, "P": 1.8, "F": 1.47, "CL": 1.75}


def _element(name: str) -> str:
    n = str(name).strip().upper()
    for e in ("CL", "BR", "S", "P", "N", "O", "H", "C", "F"):
        if n.startswith(e):
            return e
    return "C"


# ---------------------------------------------------------------------------
# integrator drivers (also used by the unit tests on synthetic systems)
# ---------------------------------------------------------------------------

def _prepare_particles(particles, temperature: float, integrator: str, radii=None):
    """XYZR / Mass / LinearVelocity or Diffusion decorators as the integrator needs."""
    for i, p in enumerate(particles):
        if not IMP.core.XYZR.get_is_setup(p):
            IMP.core.XYZR.setup_particle(p, radii[i] if radii is not None else 1.7)
        elif radii is not None:
            IMP.core.XYZR(p).set_radius(float(radii[i]))
        if not IMP.atom.Mass.get_is_setup(p):
            IMP.atom.Mass.setup_particle(p, 12.011)
        if integrator == "md":
            if not IMP.atom.LinearVelocity.get_is_setup(p):
                IMP.atom.LinearVelocity.setup_particle(p, IMP.algebra.Vector3D(0, 0, 0))
        else:
            r = IMP.core.XYZR(p).get_radius()
            d_coef = IMP.atom.get_einstein_diffusion_coefficient(r, temperature)  # A^2/fs
            if not IMP.atom.Diffusion.get_is_setup(p):
                IMP.atom.Diffusion.setup_particle(p, d_coef)
            else:
                IMP.atom.Diffusion(p).set_diffusion_coefficient(d_coef)


def make_simulator(
    model,
    mobile_particles,
    scoring_function,
    *,
    integrator: str = "md",
    temperature: float = 300.0,
    timestep_fs: float = 2.0,
    friction_ps: float = 10.0,
    seed: Optional[int] = None,
):
    """An ``IMP.atom.Simulator`` for ``mobile_particles`` (decorated as needed).

    ``md``: velocity Verlet + Langevin thermostat (``friction_ps`` in 1/ps);
    ``bd``: Brownian dynamics (overdamped; ``timestep_fs`` is the BD step).
    Fixed particles are simply not in ``mobile_particles`` (their
    coordinates must not be optimized).
    """
    if seed is not None:
        IMP.random_number_generator.seed(int(seed))
    for p in mobile_particles:
        IMP.core.XYZ(p).set_coordinates_are_optimized(True)
    _prepare_particles(mobile_particles, temperature, integrator)
    if integrator == "md":
        sim = IMP.atom.MolecularDynamics(model)
        sim.set_particles(mobile_particles)
        sim.set_scoring_function(scoring_function)
        sim.set_maximum_time_step(float(timestep_fs))
        sim.set_temperature(float(temperature))
        thermostat = IMP.atom.LangevinThermostatOptimizerState(model, mobile_particles, float(temperature), float(friction_ps))
        thermostat.set_period(1)
        sim.add_optimizer_state(thermostat)
        sim.assign_velocities(float(temperature))
        return sim
    if integrator == "bd":
        sim = IMP.atom.BrownianDynamics(model)
        sim.set_particles(mobile_particles)
        sim.set_scoring_function(scoring_function)
        sim.set_maximum_time_step(float(timestep_fs))
        sim.set_temperature(float(temperature))
        return sim
    raise ValueError(f"integrator must be 'md' or 'bd', got {integrator!r}")


# ---------------------------------------------------------------------------
# the dye sampler
# ---------------------------------------------------------------------------

@dataclass
class LangevinTrajectory:
    """What :meth:`LangevinDyeSampler.run` returns."""

    coordinates: np.ndarray            #: (n_frames, n_dye_atoms, 3)
    times_fs: np.ndarray               #: (n_frames,)
    kinetic_energy: np.ndarray         #: (n_frames,) md only (nan for bd)
    potential_energy: np.ndarray       #: (n_frames,)
    atom_names: tuple = ()
    integrator: str = "md"
    temperature: float = 300.0
    timestep_fs: float = 2.0

    @property
    def n_frames(self) -> int:
        return int(self.coordinates.shape[0])


class LangevinDyeSampler:
    """Langevin (``md``) or Brownian (``bd``) dynamics of an attached explicit dye.

    Parameters
    ----------
    protein_hier : IMP.atom.Hierarchy
        The protein (already labelled: ``attach_dyes`` has placed ``dye_hier``).
    dye_hier : IMP.atom.Hierarchy
        The dye read from ``dye_mol2`` (``IMP.atom.read_mol2``), attached at the site.
    dye_mol2 : str
        The dye's MOL2 (topology: bonds, angles, dihedrals, LJ types).
    chain, residue : str, int
        The labelled site (its atoms are not obstacles).
    integrator : "md" | "bd"
    temperature : float (K)
    timestep_fs : float
        2 fs (md) is safe with the 2000 kcal/mol/Å² bonds. Overdamped bd with
        the same stiff bonds is stable only for dt < k_B T / (D k) ≈ 2 fs
        (D ≈ 1.5e-4 Å²/fs for a 1.7 Å atom at 300 K); the default is 0.5 fs.
    friction_ps : float
        Langevin friction (1/ps), md only.
    interaction_sphere : float (Å)
        Protein heavy atoms within this distance of the site CA are obstacles.
    repulsion_k : float
        Soft-sphere spring (kcal/mol/Å²) of the dye–protein and dye–dye clash term.
    seed : int, optional
    """

    def __init__(
        self,
        protein_hier,
        dye_hier,
        dye_mol2: str,
        chain: str,
        residue: int,
        *,
        integrator: str = "md",
        temperature: float = 300.0,
        timestep_fs: Optional[float] = None,
        friction_ps: float = 10.0,
        interaction_sphere: float = 25.0,
        repulsion_k: float = 10.0,
        seed: Optional[int] = None,
    ):
        from IMP.bff.cgdye.sim.dye_restraints import build_dye_restraints
        from IMP.bff.cgdye.topology.combined import dye_forcefield_system

        if integrator not in ("md", "bd"):
            raise ValueError("integrator must be 'md' or 'bd'")
        self.integrator = integrator
        self.temperature = float(temperature)
        self.timestep_fs = float(timestep_fs) if timestep_fs is not None else (2.0 if integrator == "md" else 0.5)
        self.friction_ps = float(friction_ps)
        self.model = protein_hier.get_model()
        self.protein = protein_hier
        self.dye = dye_hier
        self.chain, self.residue = chain, int(residue)

        # dye particles in MOL2 serial order (the order dye_forcefield_system uses)
        dye_atoms = list(IMP.atom.get_by_type(dye_hier, IMP.atom.ATOM_TYPE))
        dye_atoms.sort(key=lambda a: IMP.atom.Atom(a).get_input_index())
        self.dye_particles = [a.get_particle() for a in dye_atoms]
        self.system = dye_forcefield_system(dye_mol2, "dye")
        sites = self.system["sites"]
        if len(sites) != len(self.dye_particles):
            raise ValueError(f"MOL2 has {len(sites)} atoms, the dye hierarchy {len(self.dye_particles)}")
        self.site_particles = {s["id"]: p for s, p in zip(sites, self.dye_particles)}
        self.atom_names = tuple(s["atom_name"] for s in sites)
        anchor_ids = set(self.system["groups"]["dye_anchor"])
        self.mobile = [p for s, p in zip(sites, self.dye_particles) if s["id"] not in anchor_ids]
        self.fixed = [p for s, p in zip(sites, self.dye_particles) if s["id"] in anchor_ids]

        radii = [_VDW_RADIUS.get(_element(s["atom_name"]), 1.7) for s in sites]
        masses = [float(s.get("mass", 12.011)) for s in sites]
        for p, r, m in zip(self.dye_particles, radii, masses):
            if not IMP.core.XYZR.get_is_setup(p):
                IMP.core.XYZR.setup_particle(p, r)
            else:
                IMP.core.XYZR(p).set_radius(r)
            if not IMP.atom.Mass.get_is_setup(p):
                IMP.atom.Mass.setup_particle(p, m)
            else:
                IMP.atom.Mass(p).set_mass(m)
        for p in self.fixed:
            IMP.core.XYZ(p).set_coordinates_are_optimized(False)

        # protein obstacles: heavy atoms within the interaction sphere, not the labelled residue
        ca = None
        obstacles = []
        for a in IMP.atom.get_by_type(protein_hier, IMP.atom.ATOM_TYPE):
            atom = IMP.atom.Atom(a)
            name = atom.get_atom_type().get_string().replace("HET:", "").strip()
            res_p = atom.get_parent()
            res = IMP.atom.Residue(res_p) if IMP.atom.Residue.get_is_setup(res_p) else None
            ch = IMP.atom.Chain(res_p.get_parent()).get_id() if res is not None and IMP.atom.Chain.get_is_setup(res_p.get_parent()) else ""
            same_site = res is not None and res.get_index() == self.residue and (not chain or ch == chain)
            if same_site and name.upper() == "CA":
                ca = IMP.core.XYZ(a).get_coordinates()
            if same_site or name.upper().startswith("H"):
                continue
            obstacles.append(a.get_particle())
        if ca is None:
            raise ValueError(f"site {chain}:{residue} has no CA")
        self.site_ca = np.array(ca)
        obstacles = [p for p in obstacles
                     if (IMP.core.XYZ(p).get_coordinates() - ca).get_magnitude() <= interaction_sphere]
        for p in obstacles:
            if not IMP.core.XYZR.get_is_setup(p):
                IMP.core.XYZR.setup_particle(p, 1.7)
            IMP.core.XYZ(p).set_coordinates_are_optimized(False)
        self.obstacles = obstacles

        # restraints: bonded + intra-dye LJ lower bounds, plus dye-protein soft spheres
        self.restraints = list(build_dye_restraints(self.model, self.system, self.site_particles))
        if obstacles:
            dye_lc = IMP.container.ListSingletonContainer(self.model, [p.get_index() for p in self.mobile])
            prot_lc = IMP.container.ListSingletonContainer(self.model, [p.get_index() for p in obstacles])
            cbpc = IMP.container.CloseBipartitePairContainer(dye_lc, prot_lc, 3.0, 1.0)
            self.restraints.append(IMP.container.PairsRestraint(IMP.core.SoftSpherePairScore(float(repulsion_k)), cbpc, "dye-protein"))
        self.scoring_function = IMP.core.RestraintsScoringFunction(self.restraints)
        self.simulator = make_simulator(
            self.model, self.mobile, self.scoring_function, integrator=integrator,
            temperature=self.temperature, timestep_fs=self.timestep_fs, friction_ps=self.friction_ps, seed=seed)

    # -- running -----------------------------------------------------------
    def coordinates(self) -> np.ndarray:
        return np.array([IMP.core.XYZ(p).get_coordinates() for p in self.dye_particles], dtype=np.float64)

    def energy(self) -> float:
        return float(self.scoring_function.evaluate(False))

    def minimize(self, n_steps: int = 200) -> float:
        """Conjugate-gradient relaxation of the mobile atoms before dynamics."""
        cg = IMP.core.ConjugateGradients(self.model)
        cg.set_scoring_function(self.scoring_function)
        return float(cg.optimize(int(n_steps)))

    def run(self, n_steps: int, write_every: int = 10, out_rmf: Optional[str] = None) -> LangevinTrajectory:
        """Integrate ``n_steps`` steps, keeping a frame every ``write_every`` steps."""
        n_frames = max(1, int(n_steps) // max(1, int(write_every)))
        coords = np.empty((n_frames, len(self.dye_particles), 3))
        ke = np.full(n_frames, np.nan)
        pe = np.empty(n_frames)
        rmf_fh = None
        if out_rmf is not None:
            import IMP.rmf
            import RMF
            root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.model, "root"))
            root.add_child(self.protein)
            root.add_child(self.dye)
            rmf_fh = RMF.create_rmf_file(out_rmf)
            IMP.rmf.add_hierarchies(rmf_fh, [root])
        for frame in range(n_frames):
            self.simulator.optimize(int(write_every))
            coords[frame] = self.coordinates()
            pe[frame] = self.energy()
            if self.integrator == "md":
                ke[frame] = float(self.simulator.get_kinetic_energy())
            if rmf_fh is not None:
                IMP.rmf.save_frame(rmf_fh, str(frame))
        del rmf_fh
        return LangevinTrajectory(
            coordinates=coords, times_fs=(np.arange(1, n_frames + 1) * write_every * self.timestep_fs),
            kinetic_energy=ke, potential_energy=pe, atom_names=self.atom_names,
            integrator=self.integrator, temperature=self.temperature, timestep_fs=self.timestep_fs)

    def kinetic_temperature(self, kinetic_energy: float) -> float:
        """T from the kinetic energy of the mobile atoms (md): 2 KE / (3 N k_B)."""
        n = len(self.mobile)
        return float(2.0 * kinetic_energy / (3.0 * n * KB_KCAL))


__all__ = ["LangevinDyeSampler", "LangevinTrajectory", "make_simulator", "KB_KCAL"]
