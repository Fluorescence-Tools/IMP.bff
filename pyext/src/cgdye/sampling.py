"""Samplers for explicit dyes: linker Metropolis, rotamer libraries, RRT, kinetics.

Sub-modules are imported on demand; import them explicitly
(``from IMP.bff.cgdye.sampling import library_gen``).
"""

from __future__ import annotations

from .topology import parse_dye_mol2, build_graph
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json
import logging
import math
import random

import numpy as np

from IMP.bff.scoring import DyeInternalEnergyEvaluator, boltzmann_weights, dye_internal_system, rotamer_cluster_weights, rotamer_mean_field_weights
import IMP
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.core
import IMP.bff

__all__ = [
    'KB_KCAL',
    'LangevinDyeSampler',
    'LangevinTrajectory',
    'make_langevin_simulator',
]

# --------------------------------------------------------------------------
# clustering
# --------------------------------------------------------------------------
"""Conformational clustering for rotamer library generation."""

def rmsd_no_align(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """RMSD between two conformers, atom order as given, no superposition.

    :func:`IMP.bff.rmsd_no_align`. The conformers already share a frame -- this
    is the metric a rotamer library is clustered under.
    """
    return IMP.bff.rmsd_no_align(
        np.ascontiguousarray(coords1, dtype=np.float64).ravel(),
        np.ascontiguousarray(coords2, dtype=np.float64).ravel())


def cluster_frames_leader(coords: np.ndarray, threshold: float) -> list[int]:
    """Greedy leader clustering; returns the representative frame indices.

    :func:`IMP.bff.cluster_frames_leader`. *coords* is ``(n_frames, n_atoms,
    3)``. The sweep order is part of the answer -- frame 0 leads, and each
    later frame joins the first leader within *threshold* -- so this is not a
    k-medoids that would find better centres. It is the algorithm FRETpredict's
    libraries are built with, which is why it is reproduced rather than
    improved.

    Was Python, and quadratic in frames with a numpy call per pair: 6.6 s for
    1500 frames of 60 atoms against 118 ms here.
    """
    coords = np.ascontiguousarray(coords, dtype=np.float64)
    if coords.ndim != 3:
        raise ValueError("coords must be (n_frames, n_atoms, 3)")
    if coords.shape[0] == 0:
        return []
    return [int(i) for i in IMP.bff.cluster_frames_leader(coords, float(threshold))]


def assign_frames_to_clusters(coords: np.ndarray, centers: list[int]) -> np.ndarray:
    """The nearest leader for each frame, as indices *into centers*.

    :func:`IMP.bff.assign_frames_to_clusters`.
    """
    coords = np.ascontiguousarray(coords, dtype=np.float64)
    return np.asarray(
        IMP.bff.assign_frames_to_clusters(
            coords, np.ascontiguousarray(centers, dtype=np.int32)),
        dtype=int)


# --------------------------------------------------------------------------
# kinetic
# --------------------------------------------------------------------------
"""Utilities for kinetic trajectory reconstruction from rotamer libraries."""

def rotamer_transition_matrix(transition_counts: list[list[int]]) -> np.ndarray:
    """Convert raw transition counts to probabilities.
    
    Returns matrix P where P[i,j] is the probability of jumping from state i to j.
    Rows sum to 1.
    """
    counts = np.array(transition_counts, dtype=float)
    row_sums = counts.sum(axis=1)
    
    # Handle states with no outgoing transitions (shouldn't happen in long walks)
    # by making them sink states (jump to self)
    for i in range(len(row_sums)):
        if row_sums[i] == 0:
            counts[i, i] = 1.0
            row_sums[i] = 1.0
            
    return counts / row_sums[:, np.newaxis]


def rotamer_correlation_times(transition_counts: list[list[int]], timestep: float) -> np.ndarray:
    """Calculate relaxation times from the transition matrix.
    
    Args:
        transition_counts: Matrix of transition counts.
        timestep: Time between frames in the sampling walk.
        
    Returns:
        Array of relaxation times.
    """
    p = rotamer_transition_matrix(transition_counts)
    vals, _ = np.linalg.eig(p.T)
    # Sort eigenvalues by magnitude
    vals = np.sort(np.abs(vals))[::-1]
    
    # Relaxation times: t_i = -timestep / ln(lambda_i)
    # lambda_1 should be 1.0 (stationary state)
    times = []
    for v in vals[1:]:
        if v > 1e-10 and v < 0.99999999:
            times.append(-timestep / np.log(v))
        else:
            times.append(0.0)
    return np.array(times)


def rotamer_rotational_correlation_time(transition_counts: list[list[int]], timestep: float) -> float:
    """Estimate the slowest rotational correlation time.
    
    Args:
        transition_counts: Matrix of transition counts.
        timestep: Time between frames in the sampling walk.
        
    Returns:
        The slowest correlation time.
    """
    times = rotamer_correlation_times(transition_counts, timestep)
    if len(times) > 0:
        return float(np.max(times))
    return 0.0


def reconstruct_rotamer_trajectory(
    lib: dict, 
    n_frames: int, 
    start_index: int | None = None,
    seed: int | None = None
) -> list[int]:
    """Generate a sequence of rotamer indices using transition probabilities.
    
    Args:
        lib: Rotamer library dict with 'transitions'
        n_frames: Length of desired trajectory
        start_index: Initial rotamer index (0-based)
        
    Returns:
        List of rotamer indices (0-based)
    """
    if "transitions" not in lib or lib["transitions"] is None:
        raise ValueError("Library does not contain transition data")
        
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
        
    p_matrix = rotamer_transition_matrix(lib["transitions"])
    n_states = p_matrix.shape[0]
    
    weights = np.array(lib["weight"], dtype=float)
    weights /= weights.sum() # Ensure exact sum to 1.0
    
    if start_index is None:
        # Sample starting state from Boltzmann weights
        start_index = np.random.choice(n_states, p=weights)
        
    traj = [start_index]
    current = start_index
    
    for _ in range(n_frames - 1):
        # Sample next state based on current row of P
        probs = p_matrix[current]
        probs /= probs.sum() # Ensure exact sum to 1.0
        nxt = np.random.choice(n_states, p=probs)
        traj.append(nxt)
        current = nxt
        
    return traj


# --------------------------------------------------------------------------
# langevin
# --------------------------------------------------------------------------
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


def make_langevin_simulator(
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
        from IMP.bff.scoring import build_dye_restraints
        from IMP.bff.cgdye.topology import dye_forcefield_system

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
        sites = self.system.sites
        if len(sites) != len(self.dye_particles):
            raise ValueError(f"MOL2 has {len(sites)} atoms, the dye hierarchy {len(self.dye_particles)}")
        self.site_particles = {s.id: p for s, p in zip(sites, self.dye_particles)}
        self.atom_names = tuple(s.atom_name for s in sites)
        anchor_ids = set(self.system.groups["dye_anchor"])
        self.mobile = [p for s, p in zip(sites, self.dye_particles) if s.id not in anchor_ids]
        self.fixed = [p for s, p in zip(sites, self.dye_particles) if s.id in anchor_ids]

        radii = [_VDW_RADIUS.get(_element(s.atom_name), 1.7) for s in sites]
        masses = [float(s.mass) for s in sites]
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
        self.simulator = make_langevin_simulator(
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


# --------------------------------------------------------------------------
# library_gen
# --------------------------------------------------------------------------
"""Core logic for generating rotamer libraries from stochastic sampling."""

class LinkerSampler:
    """Stochastic sampler for linker internal DOFs (angles and dihedrals)."""

    def __init__(self, dye_mol2, dye_pdb=None):
        self.dye_mol2 = dye_mol2
        self.dye_pdb = dye_pdb
        
        # Setup minimal model
        self.model = IMP.Model()
        self.hier = IMP.atom.read_mol2(dye_mol2, self.model)
        self.atoms = list(IMP.atom.get_by_type(self.hier, IMP.atom.ATOM_TYPE))
        
        # Setup coordinates
        self.idx_to_particle = {}
        self.serial_to_pos0 = {}
        for a in self.atoms:
            serial = IMP.atom.Atom(a).get_input_index()
            self.idx_to_particle[serial] = a
            c = IMP.core.XYZ(a).get_coordinates()
            self.serial_to_pos0[serial] = np.array([c[0], c[1], c[2]])
            
        # Parse topology
        self.mol2_atoms, self.mol2_bonds = parse_dye_mol2(dye_mol2, "dye")
        self.idx_to_name = {i: a["atom_name"] for i, a in self.mol2_atoms.items()}
        self.graph = build_graph(self.mol2_bonds)
        
        from .topology import ring_atoms_from_graph, _directed_bond_with_anchor, _directed_angle_with_anchor
        
        self.ring_atoms = ring_atoms_from_graph(self.graph, self.mol2_atoms)
        # CA is the labeling site anchor
        self.anchor_idx = next((i for i, n in self.idx_to_name.items() if n == "CA"), None)
        if self.anchor_idx is None:
            # Fallback to first serial if CA not found
            self.anchor_idx = sorted(self.idx_to_particle.keys())[0]
        
        # Find rotatable bonds
        self.rot_bonds = []
        anchor_names = {"N", "CA", "C"}
        for a, b in sorted(self.mol2_bonds):
            na, nb = self.idx_to_name[a], self.idx_to_name[b]
            if na in anchor_names or nb in anchor_names: continue
            if a in self.ring_atoms and b in self.ring_atoms: continue
            if na.startswith("H") or nb.startswith("H"): continue
            fixed_idx, moving_idx, moving_set = _directed_bond_with_anchor(a, b, self.graph, self.anchor_idx)
            if moving_set:
                moving_ids = [i for i in moving_set if i in self.idx_to_particle]
                if moving_ids:
                    self.rot_bonds.append((fixed_idx, moving_idx, moving_ids))

        # Find rotatable angles
        self.rot_angles = []
        for b in self.graph:
            neigh = sorted(self.graph[b])
            for i in range(len(neigh)):
                for j in range(i + 1, len(neigh)):
                    a, c = neigh[i], neigh[j]
                    if b in self.ring_atoms: continue
                    if self.idx_to_name[a].startswith("H") or self.idx_to_name[c].startswith("H"): continue
                    center_idx, moving_idx, moving_set = _directed_angle_with_anchor(a, b, c, self.graph, self.anchor_idx)
                    if moving_set:
                        moving_ids = [idx for i in moving_set if (idx := i) in self.idx_to_particle]
                        if moving_ids:
                            fixed_neighbor = a if moving_idx == c else c
                            self.rot_angles.append((center_idx, moving_idx, fixed_neighbor, moving_ids))

    def apply_config(self, cfg):
        for i, p in self.idx_to_particle.items():
            IMP.core.XYZ(p).set_coordinates(self.serial_to_pos0[i])
        
        n_dih = len(self.rot_bonds)
        dih_cfg = cfg[:n_dih]
        ang_cfg = cfg[n_dih:]

        # Apply torsion rotations sequentially
        for angle, (fixed_idx, moving_idx, moving_ids) in zip(dih_cfg, self.rot_bonds):
            pf, pm = self.idx_to_particle[fixed_idx], self.idx_to_particle[moving_idx]
            cf, cm = IMP.core.XYZ(pf).get_coordinates(), IMP.core.XYZ(pm).get_coordinates()
            axis = cm - cf
            if axis.get_magnitude() < 1e-8: continue
            rot = IMP.algebra.get_rotation_about_axis(axis, angle)
            tf = IMP.algebra.get_rotation_about_point(cf, rot)
            for mid in moving_ids:
                p = self.idx_to_particle[mid]
                c = IMP.core.XYZ(p).get_coordinates()
                IMP.core.XYZ(p).set_coordinates(tf.get_transformed(c))

        # Apply angle rotations sequentially
        for angle, (b_idx, c_idx, a_idx, moving_ids) in zip(ang_cfg, self.rot_angles):
            pb, pc, pa = self.idx_to_particle[b_idx], self.idx_to_particle[c_idx], self.idx_to_particle[a_idx]
            cb, cc, ca = [IMP.core.XYZ(p).get_coordinates() for p in [pb, pc, pa]]
            v_ba, v_bc = ca - cb, cc - cb
            axis = IMP.algebra.get_vector_product(v_ba, v_bc)
            if axis.get_magnitude() < 1e-8: continue
            rot = IMP.algebra.get_rotation_about_axis(axis, angle)
            tf = IMP.algebra.get_rotation_about_point(cb, rot)
            for mid in moving_ids:
                p = self.idx_to_particle[mid]
                c = IMP.core.XYZ(p).get_coordinates()
                IMP.core.XYZ(p).set_coordinates(tf.get_transformed(c))

    def sample(self, n_steps=1000, write_every=10, step_size_dih=0.1, step_size_ang=0.02, temperature=298.15, seed=None):
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)
            
        n_dih = len(self.rot_bonds)
        n_ang = len(self.rot_angles)
        current_cfg = [0.0] * (n_dih + n_ang)
        
        # Setup evaluator for Metropolis
        atoms_dict, bonds = parse_dye_mol2(self.dye_mol2, "dye")
        # bonded 1-2/1-3/1-4 pairs excluded from the LJ score
        evaluator = DyeInternalEnergyEvaluator(dye_internal_system(atoms_dict, bonds))
        
        self.apply_config(current_cfg)
        current_energy = evaluator.evaluate(self.get_coords())
        
        frames = []
        sorted_serials = sorted(self.idx_to_particle.keys())
        
        kb = 0.0019872041 # kcal/(mol*K)
        beta = 1.0 / (kb * temperature)
        
        for step in range(n_steps):
            proposal = []
            for i, v in enumerate(current_cfg):
                if i < n_dih: proposal.append(v + random.gauss(0, step_size_dih))
                else: proposal.append(v + random.gauss(0, step_size_ang))
            
            # Save current coords to restore if rejected
            old_coords = self.get_coords()
            
            self.apply_config(proposal)
            proposal_energy = evaluator.evaluate(self.get_coords())
            
            # Metropolis acceptance
            accept = False
            if proposal_energy <= current_energy:
                accept = True
            else:
                prob = math.exp(-beta * (proposal_energy - current_energy))
                if random.random() < prob:
                    accept = True
            
            if accept:
                current_cfg = proposal
                current_energy = proposal_energy
            else:
                # Restore old coords if needed (though apply_config uses pos0, so it's fine)
                pass

            if (step + 1) % write_every == 0:
                self.apply_config(current_cfg)
                frames.append(self.get_coords())
        return np.asarray(frames)

    def get_coords(self):
        sorted_serials = sorted(self.idx_to_particle.keys())
        coords = np.zeros((len(self.atoms), 3))
        for i, serial in enumerate(sorted_serials):
            p = self.idx_to_particle[serial]
            c = IMP.core.XYZ(p).get_coordinates()
            coords[i] = [c[0], c[1], c[2]]
        return coords


def generate_linker_rotamers(
    dye_mol2: str,
    n_steps: int = 10000,
    write_every: int = 10,
    step_size_dih: float = 0.1,
    step_size_ang: float = 0.02,
    cluster_threshold: float = 0.5,
    temperature: float = 298.15,
    seed: int = 42,
    protein_coords: np.ndarray | None = None,
    protein_elements: list[str] | None = None,
    mean_field_K: float = 1.0,
    mean_field_n_iter: int = 10,
) -> dict:
    """Generate a clustered rotamer library for a dye, including transitions.

    Args:
        dye_mol2:          Path to the dye MOL2 file.
        n_steps:           Metropolis steps.
        write_every:       Save every Nth step.
        step_size_dih:     Dihedral perturbation magnitude (rad).
        step_size_ang:     Bond-angle perturbation magnitude (rad).
        cluster_threshold: RMSD threshold for leader clustering (Å).
        temperature:       Sampling temperature in Kelvin.
        seed:              RNG seed.
        protein_coords:    Static protein atom coordinates, shape (N, 3).
                           If provided, Trick 2 (mean-field reweighting) is
                           applied after Boltzmann clustering.
        protein_elements:  Element symbols for each protein atom (len = N).
                           Required when protein_coords is provided.
        mean_field_K:      Energy scale factor K for Eq. 42 (default 1.0).
        mean_field_n_iter: Mean-field iterations (default 10).

    Returns:
        Dict with keys: 'weight', 'atom_names', 'coords', 'transitions'.
    """
    sampler = LinkerSampler(dye_mol2)
    all_coords = sampler.sample(
        n_steps=n_steps,
        write_every=write_every,
        step_size_dih=step_size_dih,
        step_size_ang=step_size_ang,
        temperature=temperature,
        seed=seed
    )

    # Setup evaluator for Boltzmann weights
    atoms_dict, bonds = parse_dye_mol2(dye_mol2, "dye")
    # bonded 1-2/1-3/1-4 pairs excluded from the LJ score (same as the sampler)
    evaluator = DyeInternalEnergyEvaluator(dye_internal_system(atoms_dict, bonds))

    # Trick 3: vectorized batch scoring (evaluate_batch is now NumPy-based)
    energies = evaluator.evaluate_batch(all_coords)

    # Clustering
    centers = cluster_frames_leader(all_coords, cluster_threshold)
    assignments = assign_frames_to_clusters(all_coords, centers)

    # Initial Boltzmann weights from internal (self) energy
    frame_weights = boltzmann_weights(energies, temperature)
    c_weights = rotamer_cluster_weights(assignments, frame_weights, len(centers))
    c_weights_arr = np.asarray(c_weights, dtype=np.float64)

    # Trick 2: mean-field reweighting by protein–dye interaction energies
    if protein_coords is not None and len(centers) > 0:
        import re
        # Derive dye element list from atom names
        dye_elements = []
        for a in sorted(atoms_dict.values(), key=lambda x: x["serial"]):
            name = a.get("atom_name", "C")
            m = re.match(r"([A-Za-z])", name)
            dye_elements.append(m.group(1).upper() if m else "C")

        prot_elems = protein_elements if protein_elements is not None else ["C"] * len(protein_coords)

        cluster_center_coords = np.array(
            [all_coords[c_idx] for c_idx in centers], dtype=np.float64
        )

        c_weights_arr = rotamer_mean_field_weights(
            rotamer_coords=cluster_center_coords,
            initial_weights=c_weights_arr,
            protein_coords=np.asarray(protein_coords, dtype=np.float64),
            dye_elements=dye_elements,
            protein_elements=list(prot_elems),
            K=mean_field_K,
            n_iter=mean_field_n_iter,
        )

    # Build Transition Matrix
    n_clusters = len(centers)
    transitions = np.zeros((n_clusters, n_clusters), dtype=int)
    for i in range(len(assignments) - 1):
        src = assignments[i]
        dst = assignments[i + 1]
        transitions[src, dst] += 1

    return {
        "weight": c_weights_arr.tolist(),
        "atom_names": [IMP.atom.Atom(a).get_name() for a in sampler.atoms],
        "coords": {i + 1: all_coords[c_idx] for i, c_idx in enumerate(centers)},
        "transitions": transitions.tolist()
    }


# --------------------------------------------------------------------------
# rrt
# --------------------------------------------------------------------------
"""Fluorescent Protein (FP) library accessor (``data/cgdye/fp_library.json``).

Copy of fpsim's FP registry (imp.bff has 4 FPs, fpsim has drifted); tracked in
chisurf PRD-96. Do not extend it here.
"""

# --------------------------------------------------------------------------
# fp_lib
# --------------------------------------------------------------------------
"""Fluorescent Protein (FP) library accessor (``data/cgdye/fp_library.json``).

Copy of fpsim's FP registry (imp.bff has 4 FPs, fpsim has drifted); tracked in
chisurf PRD-96. Do not extend it here.
"""

def load_fp_library() -> Dict[str, Dict]:
    """Load the FP library from the bundled JSON file."""
    from IMP.bff.tools import _data_root
    path = _data_root() / "fp_library.json"
    with open(path, "r") as f:
        return json.load(f)

def get_fp_names() -> List[str]:
    """Return a list of all supported FP names."""
    return list(load_fp_library().keys())

def get_fp_data(name: str) -> Dict:
    """Return the data for a specific FP."""
    lib = load_fp_library()
    if name not in lib:
        raise ValueError(f"FP '{name}' not found in library. Available: {', '.join(lib.keys())}")
    return lib[name]


# --------------------------------------------------------------------------
# rrt
# --------------------------------------------------------------------------
"""RRT samplers: rigid transforms (Transformation3D) and periodic torsion vectors."""

@dataclass
class IMPNode:
    node_id: int
    transform: IMP.algebra.Transformation3D
    parent_id: int | None


class IMPRRTTree:
    def __init__(self, root_transform):
        self.nodes = {0: IMPNode(0, root_transform, None)}
        self.next_id = 1

    def add_node(self, transform, parent_id):
        nid = self.next_id
        self.nodes[nid] = IMPNode(nid, transform, parent_id)
        self.next_id += 1
        return nid

    def nearest(self, transform, rot_weight=0.25):
        best_id = None
        best_d = float("inf")
        for nid, n in self.nodes.items():
            d = transform_distance(n.transform, transform, rot_weight=rot_weight)
            if d < best_d:
                best_d = d
                best_id = nid
        return best_id

    def path_to_root(self, node_id):
        path = []
        cur = node_id
        while cur is not None:
            path.append(cur)
            cur = self.nodes[cur].parent_id
        return list(reversed(path))


def make_transform(tx, ty, tz, rx, ry, rz):
    rot = IMP.algebra.get_rotation_from_fixed_xyz(rx, ry, rz)
    tr = IMP.algebra.Vector3D(tx, ty, tz)
    return IMP.algebra.Transformation3D(rot, tr)


def transform_to_tuple(t):
    tr = t.get_translation()
    e = IMP.algebra.get_fixed_xyz_from_rotation(t.get_rotation())
    return (tr[0], tr[1], tr[2], e.get_x(), e.get_y(), e.get_z())


def transform_distance(a, b, rot_weight=0.25):
    at = a.get_translation()
    bt = b.get_translation()
    dt = math.sqrt((at[0] - bt[0]) ** 2 + (at[1] - bt[1]) ** 2 + (at[2] - bt[2]) ** 2)

    ae = IMP.algebra.get_fixed_xyz_from_rotation(a.get_rotation())
    be = IMP.algebra.get_fixed_xyz_from_rotation(b.get_rotation())
    dr = math.sqrt(
        (ae.get_x() - be.get_x()) ** 2
        + (ae.get_y() - be.get_y()) ** 2
        + (ae.get_z() - be.get_z()) ** 2
    )
    return dt + rot_weight * dr


def sample_transform(bounds, rng):
    vals = [rng.uniform(lo, hi) for lo, hi in bounds]
    return make_transform(*vals)


def steer_transform(from_t, to_t, step_size, rot_weight=0.25):
    d = transform_distance(from_t, to_t, rot_weight=rot_weight)
    if d == 0:
        return from_t
    if d <= step_size:
        return to_t
    alpha = step_size / d
    fa = transform_to_tuple(from_t)
    ta = transform_to_tuple(to_t)
    vals = [fa[i] + alpha * (ta[i] - fa[i]) for i in range(6)]
    return make_transform(*vals)


def is_collision_sphere(transform, obstacles, probe_radius):
    tr = transform.get_translation()
    x, y, z = tr[0], tr[1], tr[2]
    for ox, oy, oz, orad in obstacles:
        d2 = (x - ox) ** 2 + (y - oy) ** 2 + (z - oz) ** 2
        md = probe_radius + orad
        if d2 < md * md:
            return True
    return False


def rrt_grow_step(tree, target_t, step_size, collision_fn, rot_weight=0.25):
    near_id = tree.nearest(target_t, rot_weight=rot_weight)
    near_t = tree.nodes[near_id].transform
    new_t = steer_transform(near_t, target_t, step_size, rot_weight=rot_weight)
    if collision_fn(new_t):
        return None
    return tree.add_node(new_t, near_id)


def run_rigid_body_rrt(
    start_transform,
    bounds,
    n_iter,
    step_size,
    collision_fn,
    goal_transform=None,
    goal_bias=0.1,
    goal_tolerance=1.0,
    rot_weight=0.25,
    seed=0,
):
    """Grow an RRT over rigid-body transforms (``IMP.algebra.Transformation3D``).

    ``bounds`` are six (lo, hi) pairs for (tx, ty, tz, rx, ry, rz);
    ``collision_fn(transform) -> bool`` gates growth. Returns ``(tree, goal_node_id)``.
    """
    rng = random.Random(seed)
    tree = IMPRRTTree(start_transform)
    goal_node_id = None

    for _ in range(n_iter):
        if goal_transform is not None and rng.random() < goal_bias:
            target = goal_transform
        else:
            target = sample_transform(bounds, rng)

        new_id = rrt_grow_step(
            tree, target, step_size, collision_fn, rot_weight=rot_weight
        )
        if new_id is None:
            continue
        if goal_transform is not None:
            d = transform_distance(
                tree.nodes[new_id].transform, goal_transform, rot_weight=rot_weight
            )
            if d <= goal_tolerance:
                goal_node_id = new_id
                break

    return tree, goal_node_id


# ---------------------------------------------------------------------------
# RRT over a vector of periodic torsion angles (internal linker DOFs)
# ---------------------------------------------------------------------------

def _wrap_angle(a):
    """Wrap an angle into [-pi, pi)."""
    return (a + math.pi) % (2.0 * math.pi) - math.pi


def torsion_distance(a, b):
    """Euclidean distance between two torsion vectors on the torus."""
    return math.sqrt(sum(_wrap_angle(x - y) ** 2 for x, y in zip(a, b)))


def steer_torsions(from_cfg, to_cfg, step_size):
    """Move from ``from_cfg`` towards ``to_cfg`` by at most ``step_size`` (rad)."""
    d = torsion_distance(from_cfg, to_cfg)
    if d == 0.0:
        return list(from_cfg)
    alpha = min(1.0, step_size / d)
    return [_wrap_angle(f + alpha * _wrap_angle(t - f)) for f, t in zip(from_cfg, to_cfg)]


class TorsionRRTTree:
    """RRT tree whose nodes are torsion vectors (list of angles in rad)."""

    def __init__(self, root_cfg):
        self.configs = [list(root_cfg)]
        self.parents = [None]

    def add_node(self, cfg, parent_id):
        self.configs.append(list(cfg))
        self.parents.append(parent_id)
        return len(self.configs) - 1

    def nearest(self, cfg):
        best_id, best_d = 0, float("inf")
        for i, c in enumerate(self.configs):
            d = torsion_distance(c, cfg)
            if d < best_d:
                best_id, best_d = i, d
        return best_id

    def __len__(self):
        return len(self.configs)


def run_torsion_rrt(
    n_dofs,
    n_iter,
    step_size,
    collision_fn,
    *,
    start_cfg=None,
    goal_cfg=None,
    goal_bias=0.1,
    goal_tolerance=0.1,
    seed=0,
):
    """Grow an RRT over ``n_dofs`` periodic torsions.

    ``collision_fn(cfg) -> bool`` returns True when the configuration clashes
    (it is expected to apply ``cfg`` to the molecule as a side effect if the
    caller needs the coordinates). Returns ``(tree, goal_node_id)``.
    """
    rng = random.Random(seed)
    root = list(start_cfg) if start_cfg is not None else [0.0] * n_dofs
    tree = TorsionRRTTree(root)
    goal_node_id = None
    for _ in range(n_iter):
        if goal_cfg is not None and rng.random() < goal_bias:
            target = list(goal_cfg)
        else:
            target = [rng.uniform(-math.pi, math.pi) for _ in range(n_dofs)]
        near_id = tree.nearest(target)
        new_cfg = steer_torsions(tree.configs[near_id], target, step_size)
        if collision_fn(new_cfg):
            continue
        new_id = tree.add_node(new_cfg, near_id)
        if goal_cfg is not None and torsion_distance(new_cfg, goal_cfg) <= goal_tolerance:
            goal_node_id = new_id
            break
    return tree, goal_node_id


# --------------------------------------------------------------------------
# segments
# --------------------------------------------------------------------------
"""Segmentation and FP detection ported from fpsim (Biopython-free copy).

This is a copy of fpsim's segmenter kept for ``dye label-fusion``; the
duplication and its drift are tracked in chisurf PRD-96 (fpsim <-> imp.bff).
Do not extend it here -- fpsim is the home of the algorithm.
"""

logger = logging.getLogger(__name__)

def _smith_waterman(query: str, templ: str,
                    match: float = 2.0, mismatch: float = -1.0,
                    gap_open: float = -5.0, gap_extend: float = -1.0):
    """Local alignment with affine gaps -- :func:`IMP.bff.smith_waterman`.

    Replaces Biopython's ``PairwiseAligner`` in local mode with the same
    scoring, so imp.bff needs nothing beyond what IMP brings. Gotoh's
    three-matrix formulation, which is what makes an affine penalty exact
    rather than a per-position approximation.

    Was 73 lines of Python over three ``(n+1, m+1)`` numpy planes, indexed
    element by element: 1,287 ms on a 900x900 pair against 8 ms here. Gated on
    300 randomised pairs, identical blocks on all of them.

    Returns
    -------
    tuple
        ``(query_blocks, templ_blocks)``, each a list of ``(start, end)``
        half-open index pairs, in the spelling Biopython's ``aligned`` uses.
    """
    blocks = IMP.bff.smith_waterman(
        query, templ, float(match), float(mismatch),
        float(gap_open), float(gap_extend))
    return ([(b.query_start, b.query_end) for b in blocks],
            [(b.template_start, b.template_end) for b in blocks])


def _align_best_window(query: str, templ: str) -> tuple[int, int, float]:
    """Local alignment; return (q_start_1based, q_end_1based, identity)."""
    q_blocks, t_blocks = _smith_waterman(query, templ)
    if not q_blocks:
        return (0, 0, 0.0)

    q_start = int(q_blocks[0][0])
    q_end = int(q_blocks[-1][1])

    matches = 0
    aligned_len = 0
    for (qs, qe), (ts, te) in zip(q_blocks, t_blocks):
        qseg = query[qs:qe]
        tseg = templ[ts:te]
        aligned_len += len(qseg)
        for qc, tc in zip(qseg, tseg):
            if qc == tc:
                matches += 1
    identity = (matches / aligned_len) if aligned_len else 0.0
    return (int(q_start + 1), int(q_end), float(identity))

def find_fp_domains(seq: str, min_identity: float = 0.35) -> List[Tuple[str, int, int, float]]:
    """Detect FP domains in a fusion protein sequence."""
    lib = load_fp_library()
    hits: List[Tuple[str, int, int, float]] = []

    for name, data in lib.items():
        templ = data["sequence"]
        motifs = data.get("motifs", [])
        
        best_seed = (-1, -1.0)
        for m in motifs:
            i = seq.find(m)
            if i >= 0 and len(m) > best_seed[1]:
                best_seed = (i, len(m))

        if best_seed[0] >= 0:
            i = best_seed[0]
            start = max(0, i - 50)
            end = min(len(seq), i + 300)
            qs, qe, idy = _align_best_window(seq[start:end], templ)
            if idy >= min_identity and qs > 0:
                hits.append((name, int(start + qs), int(start + qe), float(idy)))
            continue

        qs, qe, idy = _align_best_window(seq, templ)
        if idy >= min_identity and qs > 0:
            hits.append((name, int(qs), int(qe), float(idy)))

    hits.sort(key=lambda x: (-x[3], x[2] - x[1], x[1]))
    covered = set()
    non_overlap: List[Tuple[str, int, int, float]] = []
    for nm, s, e, idy in hits:
        if any(i in covered for i in range(s, e + 1)):
            continue
        for i in range(s, e + 1):
            covered.add(i)
        non_overlap.append((nm, s, e, idy))
    return non_overlap

def parse_plddt_from_pdb(pdb_path: str | Path, chain_id: str = "A") -> Dict[int, float]:
    """Parse pLDDT from AlphaFold PDB B-factor column."""
    plddt_sum: Dict[int, float] = {}
    plddt_cnt: Dict[int, int] = {}
    
    with open(pdb_path, "r") as f:
        for line in f:
            if not line.startswith("ATOM"): continue
            ch = line[21].strip()
            if chain_id and ch and ch != chain_id: continue
            try:
                resseq = int(line[22:26])
                b = float(line[60:66])
                plddt_sum[resseq] = plddt_sum.get(resseq, 0.0) + b
                plddt_cnt[resseq] = plddt_cnt.get(resseq, 0) + 1
            except ValueError: continue

    if not plddt_sum:
        raise ValueError(f"No pLDDT data found for chain '{chain_id}' in {pdb_path}")
    return {r: plddt_sum[r] / plddt_cnt[r] for r in sorted(plddt_sum.keys())}

def segments_from_plddt(
    seq_len: int,
    plddt: Dict[int, float],
    fp_domains: List[Tuple[str, int, int, float]],
    rigid_threshold: float = 70.0,
    min_rb_len: int = 12
) -> List[Dict]:
    """Create rigid/linker segments from pLDDT and FP detections."""
    labels = ['R' if plddt.get(i+1, 0.0) >= rigid_threshold else 'L' for i in range(seq_len)]
    
    # FP domains are always rigid
    for _, s, e, _ in fp_domains:
        for i in range(s-1, e):
            labels[i] = 'R'
            
    # Small rigid islands become linkers
    i = 0
    while i < seq_len:
        if labels[i] == 'R':
            j = i
            while j < seq_len and labels[j] == 'R': j += 1
            if (j - i) < min_rb_len:
                # Only if not part of an FP? Let's keep it simple for now
                for k in range(i, j): labels[k] = 'L'
            i = j
        else: i += 1
        
    segs = []
    i = 0
    while i < seq_len:
        curr = labels[i]
        j = i
        while j < seq_len and labels[j] == curr: j += 1
        
        kind = "core" if curr == 'R' else "linker"
        name = kind
        
        # Check if it's an FP
        for nm, s, e, _ in fp_domains:
            if max(i+1, s) <= min(j, e):
                if (min(j, e) - max(i+1, s) + 1) > 0.5 * (j - i):
                    kind = "fp"
                    name = nm
                    break
                    
        segs.append({"kind": kind, "name": name, "start": i+1, "end": j})
        i = j
    return segs
