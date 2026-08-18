"""Core logic for generating rotamer libraries from stochastic sampling."""

from __future__ import annotations

import math
import random
import numpy as np

import IMP
import IMP.algebra
import IMP.atom

from ..topology.builder import parse_dye_mol2, build_graph
from .clustering import cluster_frames_leader, assign_frames_to_clusters
from IMP.bff.scoring.boltzmann import boltzmann_weights, rotamer_cluster_weights
from IMP.bff.scoring.dye_lj import DyeInternalEnergyEvaluator, compute_lj_pair_sites, dye_internal_system
from IMP.bff.scoring.mean_field import rotamer_mean_field_weights


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
        
        from ..topology.builder import ring_atoms_from_graph, _directed_bond_with_anchor, _directed_angle_with_anchor
        
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
