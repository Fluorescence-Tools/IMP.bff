#!/usr/bin/env python
"""Generate a cgdye template CIF from a MOL2 file.

This is a bootstrap tool that auto-detects chemistry features from the
MOL2 bond graph and writes a starter template CIF. The generated template
should be hand-edited to add region definitions, colors, and flex groups.

Usage:
    python scripts/gen_template.py --mol2 inputs/structures/cx4.mol2 --name CX4 \\
        --out cgdye/templates/cx4.template.cif
"""

import re
from collections import defaultdict
from pathlib import Path

import click
import ihm.format


def _parse_mol2(path):
    """Parse MOL2 file, return (atoms dict, bonds set, graph)."""
    atoms = {}
    bonds = []
    in_atom = False
    in_bond = False

    with open(path) as fh:
        for line in fh:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom = True
                in_bond = False
                continue
            if line.startswith("@<TRIPOS>BOND"):
                in_bond = True
                in_atom = False
                continue
            if line.startswith("@<TRIPOS>"):
                in_atom = False
                in_bond = False
                continue

            if in_atom:
                parts = line.split()
                if len(parts) < 5:
                    continue
                serial = int(parts[0])
                atom_name = parts[1]
                x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                atom_type = parts[5]
                atoms[serial] = {
                    "serial": serial,
                    "name": atom_name,
                    "type": atom_type,
                    "x": x,
                    "y": y,
                    "z": z,
                }

            if in_bond:
                parts = line.split()
                if len(parts) < 3:
                    continue
                a = int(parts[1])
                b = int(parts[2])
                bonds.append((a, b))

    graph = defaultdict(set)
    for a, b in bonds:
        graph[a].add(b)
        graph[b].add(a)

    return atoms, bonds, graph


def _element_from_atom_name(atom_name):
    m = re.match(r"([A-Za-z]+)", atom_name)
    if not m:
        return "C"
    return m.group(1)[0].upper()


def _find_cycles(graph, max_len=7):
    """Find all cycles up to max_len using DFS."""
    nodes = sorted(graph.keys())
    cycles = set()

    def norm(cyc):
        cyc = list(cyc)
        n = len(cyc)
        rots = [tuple(cyc[i:] + cyc[:i]) for i in range(n)]
        rc = list(reversed(cyc))
        rots += [tuple(rc[i:] + rc[:i]) for i in range(n)]
        return min(rots)

    def dfs(start, cur, visited, path):
        for nbr in graph.get(cur, set()):
            if nbr == start and len(path) >= 3:
                cycles.add(norm(path))
                continue
            if nbr in visited or nbr < start or len(path) >= max_len:
                continue
            dfs(start, nbr, visited | {nbr}, path + [nbr])

    for s in nodes:
        dfs(s, s, {s}, [s])
    return cycles


def _ring_atoms_from_graph(graph, atoms):
    """Find atoms in rings of size 5, 6, or 7."""
    cycles = _find_cycles(graph, max_len=7)
    rings = [c for c in cycles if len(c) in (5, 6, 7)]
    out = set()
    for cyc in rings:
        for s in cyc:
            atom = atoms.get(s, {})
            if _element_from_atom_name(atom.get("name", "")) != "H":
                out.add(s)
    return out


def _detect_ring_centers(atoms, graph):
    """Return set of atom serials that are in rings."""
    return _ring_atoms_from_graph(graph, atoms)


def _detect_pi_centers(atoms, graph):
    """Return set of atom serials that are sp2 (pi) centers."""
    out = set()
    for s, atom in atoms.items():
        el = _element_from_atom_name(atom.get("name", ""))
        if el not in {"C", "N"}:
            continue
        nbrs = list(graph.get(s, set()))
        heavy_nbrs = [
            n
            for n in nbrs
            if _element_from_atom_name(atoms.get(n, {}).get("name", "")) != "H"
        ]
        if len(heavy_nbrs) == 3:
            out.add(s)
    return out


def _detect_flat_centers(atoms, graph):
    """Return set of atom serials that are bonded to >=3 O atoms (e.g., sulfonates)."""
    out = set()
    for s, atom in atoms.items():
        if not atom.get("name", "").startswith("S"):
            continue
        nbrs = list(graph.get(s, set()))
        o_nbrs = [
            n
            for n in nbrs
            if _element_from_atom_name(atoms.get(n, {}).get("name", "")) == "O"
        ]
        if len(o_nbrs) >= 3:
            out.add(s)
    return out


def _detect_orient_centers(atoms, graph):
    """Return set of atom serials that are bonded to 1 C + >=2 O (sulfonate orientation)."""
    out = set()
    for s, atom in atoms.items():
        if not atom.get("name", "").startswith("S"):
            continue
        nbrs = list(graph.get(s, set()))
        c_nbrs = [
            n
            for n in nbrs
            if _element_from_atom_name(atoms.get(n, {}).get("name", "")) == "C"
        ]
        o_nbrs = [
            n
            for n in nbrs
            if _element_from_atom_name(atoms.get(n, {}).get("name", "")) == "O"
        ]
        if len(c_nbrs) == 1 and len(o_nbrs) >= 2:
            out.add(s)
    return out


def _serial_to_name(atoms):
    return {s: a["name"] for s, a in atoms.items()}


@click.command()
@click.option("--mol2", required=True, help="Input MOL2 file")
@click.option("--name", required=True, help="Component name (e.g. CX4, atto655)")
@click.option("--out", required=True, help="Output template CIF file")
def main(mol2, name, out):
    atoms, bonds, graph = _parse_mol2(mol2)
    serial_to_name = _serial_to_name(atoms)

    ring_centers = _detect_ring_centers(atoms, graph)
    pi_centers = _detect_pi_centers(atoms, graph)
    flat_centers = _detect_flat_centers(atoms, graph)
    orient_centers = _detect_orient_centers(atoms, graph)

    print(f"Parsed {len(atoms)} atoms, {len(bonds)} bonds from {mol2}")
    print(f"Detected ring centers: {len(ring_centers)}")
    print(f"Detected pi centers: {len(pi_centers)}")
    print(f"Detected flat centers: {len(flat_centers)}")
    print(f"Detected orient centers: {len(orient_centers)}")

    impropers = []
    for s in sorted(ring_centers):
        impropers.append({"center_atom": serial_to_name[s], "type": "ring"})
    for s in sorted(pi_centers):
        impropers.append({"center_atom": serial_to_name[s], "type": "pi"})
    for s in sorted(flat_centers):
        impropers.append({"center_atom": serial_to_name[s], "type": "flat"})
    for s in sorted(orient_centers):
        impropers.append({"center_atom": serial_to_name[s], "type": "orient"})

    template = {
        "name": name,
        "features": {
            "ring_core": {
                "rb": False,
                "md_fixed": False,
                "feature_type": "dof",
                "region_color": None,
                "atoms": [
                    {"name": serial_to_name[s], "occurrence": 1}
                    for s in sorted(ring_centers)
                ],
            }
        },
        "impropers": impropers,
    }

    Path(out).parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        w = ihm.format.CifWriter(fh)
        w.start_block(name)

        with w.category("_cgdye_template") as c:
            c.write(name=name)

        with w.loop(
            "_cgdye_feature",
            ["feature_id", "feature_type", "rb", "md_fixed", "region_color"],
        ) as l:
            for fid, spec in template["features"].items():
                l.write(
                    feature_id=fid,
                    feature_type=spec.get("feature_type", "dof"),
                    rb="YES" if spec.get("rb", False) else "NO",
                    md_fixed="YES" if spec.get("md_fixed", False) else "NO",
                    region_color=spec.get("region_color") or ".",
                )

        with w.loop(
            "_cgdye_feature_atom",
            ["feature_id", "atom_name", "occurrence"],
        ) as l:
            for fid, spec in template["features"].items():
                for a in spec.get("atoms", []):
                    l.write(
                        feature_id=fid,
                        atom_name=a.get("name"),
                        occurrence=int(a.get("occurrence", 1)),
                    )

        if template.get("impropers"):
            with w.loop(
                "_cgdye_improper",
                ["center_atom", "type"],
            ) as l:
                for imp in template["impropers"]:
                    l.write(
                        center_atom=imp.get("center_atom"),
                        type=imp.get("type"),
                    )

    print(f"Wrote template to {out}")
    print("NOTE: Edit the template to add regions, flex groups, and region colors.")


if __name__ == "__main__":
    main()
