#!/usr/bin/env python
"""Collision-gated Metropolis walk over a dye linker's internal DOFs (hGBP1 site 481).

A proposal (random torsion/angle perturbation) is accepted whenever the dye
does not clash with the protein: no forces, friction or temperature are
involved, so this is *not* Langevin dynamics. Real Langevin/Brownian dynamics
lives in ``IMP.bff.cgdye.sampling.langevin`` (``dye sample-langevin``).
"""

import math
import os
import random
import sys
from collections import deque
from pathlib import Path

import click
import IMP
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.core
import IMP.rmf
import RMF



from IMP.bff.cgdye.labeling.attachment import attach_dyes
from IMP.bff.cgdye.topology.builder import parse_mol2

def _structure(name):
    """A bundled input structure, wherever IMP.bff is installed.

    These paths used to be relative to the working directory, so a script only
    ran from one place -- and stopped running at all once the data became IMP
    module data under data/cgdye.
    """
    from IMP.bff.cgdye.utils import get_structure_dir
    return str(get_structure_dir(name))



def _protein_residue_groups(protein_hier):
    groups = {}
    for a in IMP.atom.get_by_type(protein_hier, IMP.atom.ATOM_TYPE):
        res_p = a.get_parent()
        if not IMP.atom.Residue.get_is_setup(res_p):
            continue
        res = IMP.atom.Residue(res_p)
        chain_p = res_p.get_parent()
        if not IMP.atom.Chain.get_is_setup(chain_p):
            continue
        chain = IMP.atom.Chain(chain_p)
        key = (chain.get_id(), int(res.get_index()))
        groups.setdefault(key, []).append(a)
    return groups


def _center_and_radius(coords):
    if not coords:
        return None
    cx = sum(c[0] for c in coords) / len(coords)
    cy = sum(c[1] for c in coords) / len(coords)
    cz = sum(c[2] for c in coords) / len(coords)
    r = 0.0
    for c in coords:
        dx = c[0] - cx
        dy = c[1] - cy
        dz = c[2] - cz
        d = math.sqrt(dx * dx + dy * dy + dz * dz)
        if d > r:
            r = d
    return (cx, cy, cz, r)


def _build_protein_obstacles(
    protein_hier,
    site_ca,
    resolution,
    interaction_sphere,
):
    groups = _protein_residue_groups(protein_hier)
    residues = sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1]))
    obstacles = []

    if resolution <= 0:
        for _key, atoms in residues:
            for a in atoms:
                c = IMP.core.XYZ(a).get_coordinates()
                dx = c[0] - site_ca[0]
                dy = c[1] - site_ca[1]
                dz = c[2] - site_ca[2]
                if (
                    dx * dx + dy * dy + dz * dz
                    > interaction_sphere * interaction_sphere
                ):
                    continue
                obstacles.append((c[0], c[1], c[2], 0.0))
        return obstacles

    by_chain = {}
    for (cid, ridx), atoms in residues:
        by_chain.setdefault(cid, []).append((ridx, atoms))

    for cid, chain_res in by_chain.items():
        chain_res.sort(key=lambda x: x[0])
        for i in range(0, len(chain_res), resolution):
            block = chain_res[i : i + resolution]
            coords = []
            for _ridx, atoms in block:
                for a in atoms:
                    c = IMP.core.XYZ(a).get_coordinates()
                    coords.append((c[0], c[1], c[2]))
            cr = _center_and_radius(coords)
            if cr is None:
                continue
            cx, cy, cz, rad = cr
            dx = cx - site_ca[0]
            dy = cy - site_ca[1]
            dz = cz - site_ca[2]
            if dx * dx + dy * dy + dz * dz > interaction_sphere * interaction_sphere:
                continue
            obstacles.append((cx, cy, cz, rad))
    return obstacles


def _ring_atoms(graph):
    ring = set()
    for start in graph:
        stack = [(start, [start])]
        while stack:
            cur, path = stack.pop()
            for nb in graph[cur]:
                if nb == start and len(path) >= 3:
                    ring.update(path)
                elif nb not in path and len(path) < 8:
                    stack.append((nb, path + [nb]))
    return ring


def _component_without_edge(graph, start, block_u, block_v):
    q = deque([start])
    seen = {start}
    while q:
        u = q.popleft()
        for nb in graph[u]:
            if (u == block_u and nb == block_v) or (u == block_v and nb == block_u):
                continue
            if nb not in seen:
                seen.add(nb)
                q.append(nb)
    return seen


def _directed_bond_with_anchor(a, b, graph, anchor_idx):
    comp_b = _component_without_edge(graph, b, a, b)
    if anchor_idx not in comp_b:
        return a, b, comp_b
    comp_a = _component_without_edge(graph, a, a, b)
    if anchor_idx not in comp_a:
        return b, a, comp_a
    return None, None, None


def _directed_angle_with_anchor(a, b, c, graph, anchor_idx):
    # For angle a-b-c, we rotate the part starting at c around axis at b.
    comp_c = _component_without_edge(graph, c, b, c)
    if anchor_idx not in comp_c:
        return b, c, comp_c
    comp_a = _component_without_edge(graph, a, a, b)
    if anchor_idx not in comp_a:
        return b, a, comp_a
    return None, None, None


def _pdb_linker_indices(dye_pdb):
    # Linker atoms are in residue C1R.
    anchors = {"N", "CA", "C", "O"}
    out = set()
    with open(dye_pdb) as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            serial = int(line[6:11])
            aname = line[12:16].strip()
            resname = line[17:20].strip()
            if resname != "C1R":
                continue
            if aname.startswith("H") or aname in anchors:
                continue
            out.add(serial)
    return out


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--protein-pdb", default=_structure("1DG3.pdb"), show_default=True)
@click.option(
    "--dye-pdb", default=_structure("alexa488_r48.pdb"), show_default=True
)
@click.option(
    "--dye-mol2", default=_structure("alexa488_r48.mol2"), show_default=True
)
@click.option("--chain", default="A", show_default=True)
@click.option("--residue", default=481, show_default=True, type=int)
@click.option("--n-steps", default=1000, show_default=True, type=int)
@click.option("--write-every", default=10, show_default=True, type=int)
@click.option(
    "--step-size",
    default=0.1,
    show_default=True,
    type=float,
    help="StdDev of Gaussian perturbation for dihedrals in radians.",
)
@click.option(
    "--step-size-angle",
    default=0.02,
    show_default=True,
    type=float,
    help="StdDev of Gaussian perturbation for bond angles in radians.",
)
@click.option("--seed", default=481, show_default=True, type=int)
@click.option("--collision-cutoff", default=1.2, show_default=True, type=float)
@click.option("--vdw-scale", default=0.8, show_default=True, type=float)
@click.option("--protein-resolution", default=1, show_default=True, type=int)
@click.option("--interaction-sphere", default=35.0, show_default=True, type=float)
@click.option(
    "--output-rmf",
    default="output/test_runs/hgbp1_1dg3_site481/dof_walk/trajectory.rmf3",
    show_default=True,
)
def main(
    protein_pdb,
    dye_pdb,
    dye_mol2,
    chain,
    residue,
    n_steps,
    write_every,
    step_size,
    step_size_angle,
    seed,
    collision_cutoff,
    vdw_scale,
    protein_resolution,
    interaction_sphere,
    output_rmf,
):
    """Collision-gated Metropolis walk over the linker's internal DOFs."""
    random.seed(seed)
    model = IMP.Model()
    protein = IMP.atom.read_pdb(protein_pdb, model, IMP.atom.NonWaterPDBSelector())
    dye = IMP.atom.read_pdb(dye_pdb, model, IMP.atom.AllPDBSelector())

    # Attach once; this anchors the dye to labeling site using backbone frame logic.
    attached = attach_dyes(protein, [(dye, chain, residue)], strip_site_sidechain=True)
    site = attached[0]["site"]
    site_ca = IMP.core.XYZ(site["CA"]).get_coordinates()

    # Ensure bonds are present for visualization
    IMP.atom.add_bonds(protein)

    dye_atoms = list(IMP.atom.get_by_type(dye, IMP.atom.ATOM_TYPE))
    prot_obstacles = _build_protein_obstacles(
        protein,
        site_ca,
        protein_resolution,
        interaction_sphere,
    )
    eff_cutoff = collision_cutoff * vdw_scale

    # Map MOL2 serial index -> IMP atom.
    idx_to_particle = {}
    idx_to_xyz0 = {}
    for i, a in enumerate(dye_atoms, start=1):
        idx_to_particle[i] = a
        c = IMP.core.XYZ(a).get_coordinates()
        idx_to_xyz0[i] = IMP.algebra.Vector3D(c[0], c[1], c[2])

    # Mol2 graph for internal connectivity.
    mol2_atoms, mol2_bonds = parse_mol2(dye_mol2, "dye")
    idx_to_name = {i: a["atom_name"] for i, a in mol2_atoms.items()}

    # Add dye bonds to model for visualization
    added_bonds = set()
    for a_idx, b_idx in mol2_bonds:
        if a_idx in idx_to_particle and b_idx in idx_to_particle:
            pair = tuple(sorted((a_idx, b_idx)))
            if pair not in added_bonds:
                p1, p2 = idx_to_particle[a_idx], idx_to_particle[b_idx]
                if not IMP.atom.Bonded.get_is_setup(p1):
                    b1 = IMP.atom.Bonded.setup_particle(p1)
                else:
                    b1 = IMP.atom.Bonded(p1)
                if not IMP.atom.Bonded.get_is_setup(p2):
                    b2 = IMP.atom.Bonded.setup_particle(p2)
                else:
                    b2 = IMP.atom.Bonded(p2)
                IMP.atom.create_bond(b1, b2, 1)
                added_bonds.add(pair)

    # Filter protein obstacles that overlap with the initial anchored configuration.
    prot_obstacles_filt = []
    baseline_dye = [idx_to_xyz0[i] for i in sorted(idx_to_xyz0)]
    for pcx, pcy, pcz, prad in prot_obstacles:
        keep = True
        for dc in baseline_dye:
            dx = pcx - dc[0]
            dy = pcy - dc[1]
            dz = pcz - dc[2]
            thr = 1.8 + prad
            if dx * dx + dy * dy + dz * dz < thr * thr:
                keep = False
                break
        if keep:
            prot_obstacles_filt.append((pcx, pcy, pcz, prad))

    linker_indices = _pdb_linker_indices(dye_pdb)
    graph = {i: set() for i in mol2_atoms}
    for a, b in mol2_bonds:
        graph[a].add(b)
        graph[b].add(a)

    ring_atoms = _ring_atoms(graph)
    anchor_idx = next((i for i, n in idx_to_name.items() if n == "CA"), None)
    if anchor_idx is None:
        raise ValueError("Anchor CA not found in dye MOL2")

    # Find rotatable bonds in the linker.
    anchor_names = {"N", "CA", "C"}
    rot_bonds = []
    for a, b in sorted(mol2_bonds):
        na, nb = idx_to_name[a], idx_to_name[b]
        if na in anchor_names or nb in anchor_names:
            continue
        if not (a in linker_indices and b in linker_indices):
            continue
        if a in ring_atoms and b in ring_atoms:
            continue
        if idx_to_name[a].startswith("H") or idx_to_name[b].startswith("H"):
            continue
        fixed_idx, moving_idx, moving_set = _directed_bond_with_anchor(
            a, b, graph, anchor_idx
        )
        if moving_set:
            moving_ids = [i for i in moving_set if i in idx_to_particle]
            if moving_ids:
                rot_bonds.append((fixed_idx, moving_idx, moving_ids))

    # Find rotatable bond angles in the linker.
    rot_angles = []
    for b in graph:
        neigh = sorted(graph[b])
        for i in range(len(neigh)):
            for j in range(i + 1, len(neigh)):
                a, c = neigh[i], neigh[j]
                if (
                    a not in linker_indices
                    or b not in linker_indices
                    or c not in linker_indices
                ):
                    continue
                if b in ring_atoms:
                    continue
                if idx_to_name[a].startswith("H") or idx_to_name[c].startswith("H"):
                    continue
                center_idx, moving_idx, moving_set = _directed_angle_with_anchor(
                    a, b, c, graph, anchor_idx
                )
                if moving_set:
                    moving_ids = [idx for i in moving_set if (idx := i) in idx_to_particle]
                    if moving_ids:
                        fixed_neighbor = a if moving_idx == c else c
                        rot_angles.append(
                            (center_idx, moving_idx, fixed_neighbor, moving_ids)
                        )

    def apply_config(cfg):
        # Reset to baseline
        for i, p in idx_to_particle.items():
            IMP.core.XYZ(p).set_coordinates(idx_to_xyz0[i])

        n_dih = len(rot_bonds)
        dih_cfg = cfg[:n_dih]
        ang_cfg = cfg[n_dih:]

        # Apply torsion rotations sequentially
        for angle, (fixed_idx, moving_idx, moving_ids) in zip(dih_cfg, rot_bonds):
            # Get CURRENT coordinates for axis
            pf, pm = idx_to_particle[fixed_idx], idx_to_particle[moving_idx]
            cf, cm = (
                IMP.core.XYZ(pf).get_coordinates(),
                IMP.core.XYZ(pm).get_coordinates(),
            )
            axis = cm - cf
            if axis.get_magnitude() < 1e-8:
                continue
            rot = IMP.algebra.get_rotation_about_axis(axis, angle)
            tf = IMP.algebra.get_rotation_about_point(cf, rot)
            for mid in moving_ids:
                p = idx_to_particle[mid]
                c = IMP.core.XYZ(p).get_coordinates()
                IMP.core.XYZ(p).set_coordinates(tf.get_transformed(c))

        # Apply angle rotations sequentially
        for angle, (b_idx, c_idx, a_idx, moving_ids) in zip(ang_cfg, rot_angles):
            # Get CURRENT coordinates for axis
            pb, pc, pa = (
                idx_to_particle[b_idx],
                idx_to_particle[c_idx],
                idx_to_particle[a_idx],
            )
            cb, cc, ca = [IMP.core.XYZ(p).get_coordinates() for p in [pb, pc, pa]]
            v_ba = ca - cb
            v_bc = cc - cb
            axis = IMP.algebra.get_vector_product(v_ba, v_bc)
            if axis.get_magnitude() < 1e-8:
                continue
            rot = IMP.algebra.get_rotation_about_axis(axis, angle)
            tf = IMP.algebra.get_rotation_about_point(cb, rot)
            for mid in moving_ids:
                p = idx_to_particle[mid]
                c = IMP.core.XYZ(p).get_coordinates()
                IMP.core.XYZ(p).set_coordinates(tf.get_transformed(c))

    def has_collision():
        for a in dye_atoms:
            c = IMP.core.XYZ(a).get_coordinates()
            x, y, z = c[0], c[1], c[2]
            for pcx, pcy, pcz, prad in prot_obstacles_filt:
                dx, dy, dz = x - pcx, y - pcy, z - pcz
                thr = eff_cutoff + prad
                if dx * dx + dy * dy + dz * dz < thr * thr:
                    return True
        return False

    current_cfg = [0.0] * (len(rot_bonds) + len(rot_angles))
    apply_config(current_cfg)

    out_dir = os.path.dirname(os.path.abspath(output_rmf))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model, "root"))
    root.add_child(protein)
    root.add_child(dye)
    fh = RMF.create_rmf_file(output_rmf)
    IMP.rmf.add_hierarchies(fh, [root])

    IMP.rmf.save_frame(fh, "init")
    n_frames = n_steps // write_every
    accepted = 0

    for step in range(n_steps):
        n_dih = len(rot_bonds)
        proposal = []
        for i, v in enumerate(current_cfg):
            if i < n_dih:
                proposal.append(v + random.gauss(0, step_size))
            else:
                proposal.append(v + random.gauss(0, step_size_angle))

        apply_config(proposal)
        if not has_collision():
            current_cfg = proposal
            accepted += 1
        else:
            apply_config(current_cfg)

        if (step + 1) % write_every == 0:
            IMP.rmf.save_frame(fh, str((step + 1) // write_every))

    click.echo(
        f"Internal-DOF Metropolis walk finished:\n"
        f"  Steps: {n_steps}, Frames: {n_frames + 1}\n"
        f"  Accepted: {accepted} ({100.0 * accepted / n_steps:.1f}%)\n"
        f"  Internal DOFs (linker torsions): {len(rot_bonds)}\n"
        f"  Internal DOFs (linker angles): {len(rot_angles)}\n"
        f"  Wrote RMF: {output_rmf}"
    )


if __name__ == "__main__":
    main()
