#!/usr/bin/env python
"""Build compact FF mmCIF system from N-component MOL2 files.

Uses IMP.atom.read_mol2() for bond topology, then reads template CIFs for
feature definitions and improper centers.

Usage:
    python -m cgdye.topology.builder \\
        --component name=hgbp1,mol2=inputs/structures/1DG3.mol2,role=fixed \\
        --component name=alexa488,mol2=inputs/structures/alexa488_r48.mol2,role=mobile \\
        --output-cif output/systems/hgbp1_alexa488.system.cif
"""

import math
import os
import re
from collections import defaultdict
from pathlib import Path

import IMP
import IMP.algebra
import IMP.atom
import IMP.core

from ..io.cif import write_ff_system
from ..io.template_cif import read_cgdye_template
from IMP.bff.cgdye.utils import import_click

click = import_click()  # optional: only the CLI entry point needs it


def _atom_name_from_type(atom_type_string):
    return (
        atom_type_string.replace("HET: ", "")
        .replace("HET:", "")
        .replace("ATOM: ", "")
        .replace("ATOM:", "")
        .strip()
    )


def _element_from_atom_name(atom_name):
    m = re.match(r"([A-Za-z]+)", atom_name)
    if not m:
        return "C"
    return m.group(1)[0].upper()


def _read_mol2_atom_names(path) -> dict:
    names = {}
    in_atom = False
    with open(path) as fh:
        for line in fh:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom = True
                continue
            if line.startswith("@<TRIPOS>"):
                in_atom = False
                continue
            if not in_atom:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            names[int(parts[0])] = parts[1]
    return names


def parse_mol2(path, component):
    model = IMP.Model()
    old_level = IMP.get_log_level()
    try:
        IMP.set_log_level(IMP.SILENT)
        hier = IMP.atom.read_mol2(str(path), model)
    finally:
        IMP.set_log_level(old_level)

    mol2_atom_names = _read_mol2_atom_names(path)

    atoms = {}
    imp_particles = {}
    for p in IMP.atom.get_by_type(hier, IMP.atom.ATOM_TYPE):
        atom = IMP.atom.Atom(p)
        serial = atom.get_input_index()
        atom_name = mol2_atom_names.get(serial) or _atom_name_from_type(
            atom.get_atom_type().get_string()
        )
        residue = IMP.atom.get_residue(atom)
        coord = IMP.core.XYZ(p).get_coordinates()
        resname = residue.get_residue_type().get_string() if residue else "UNK"
        atoms[serial] = {
            "serial": serial,
            "component": component,
            "atom_name": atom_name,
            "resname": resname or "UNK",
            "element": _element_from_atom_name(atom_name),
            "x": float(coord[0]),
            "y": float(coord[1]),
            "z": float(coord[2]),
        }
        imp_particles[serial] = p

    bonds = set()
    for bond_p in IMP.atom.get_internal_bonds(hier):
        b = IMP.atom.Bond(bond_p)
        s1 = IMP.atom.Atom(b.get_bonded(0).get_particle()).get_input_index()
        s2 = IMP.atom.Atom(b.get_bonded(1).get_particle()).get_input_index()
        bonds.add(tuple(sorted((s1, s2))))

    return atoms, bonds


def distance(a, b):
    va = IMP.algebra.Vector3D(float(a["x"]), float(a["y"]), float(a["z"]))
    vb = IMP.algebra.Vector3D(float(b["x"]), float(b["y"]), float(b["z"]))
    return float(IMP.algebra.get_distance(va, vb))


def build_graph(bonds):
    g = defaultdict(set)
    for a, b in bonds:
        g[a].add(b)
        g[b].add(a)
    return g


def build_angles(graph):
    angles = set()
    for b, neigh in graph.items():
        nn = sorted(neigh)
        for i in range(len(nn)):
            for j in range(i + 1, len(nn)):
                angles.add((nn[i], b, nn[j]))
    return sorted(angles)


def build_dihedrals(graph):
    ds = set()
    for b, neigh in graph.items():
        for c in neigh:
            if b > c:
                continue
            left = [a for a in graph[b] if a != c]
            right = [d for d in graph[c] if d != b]
            for a in left:
                for d in right:
                    if len({a, b, c, d}) < 4:
                        continue
                    t1 = (a, b, c, d)
                    t2 = (d, c, b, a)
                    ds.add(min(t1, t2))
    return sorted(ds)


def find_cycles(graph, max_len=7):
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


def _component_without_edge(graph, start, block_u, block_v):
    from collections import deque
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


def ring_atoms_from_graph(graph, serial_to_atom, ring_sizes=(5, 6, 7)):
    cycles = [
        c for c in find_cycles(graph, max_len=max(ring_sizes)) if len(c) in ring_sizes
    ]
    out = set()
    for cyc in cycles:
        for s in cyc:
            a = serial_to_atom.get(s)
            if a and a.get("element", "C") != "H":
                out.add(s)
    return out


def angle_value(a, b, c):
    bax = float(a["x"] - b["x"])
    bay = float(a["y"] - b["y"])
    baz = float(a["z"] - b["z"])
    bcx = float(c["x"] - b["x"])
    bcy = float(c["y"] - b["y"])
    bcz = float(c["z"] - b["z"])
    ba2 = bax * bax + bay * bay + baz * baz
    bc2 = bcx * bcx + bcy * bcy + bcz * bcz
    if ba2 == 0.0 or bc2 == 0.0:
        return 1.910633
    dot = bax * bcx + bay * bcy + baz * bcz
    cosang = max(-1.0, min(1.0, dot / math.sqrt(ba2 * bc2)))
    return float(math.acos(cosang))


def sid(comp, atom_name):
    return f"{comp}/{atom_name}"


def _alpha_suffix(i):
    chars = []
    n = int(i)
    while True:
        chars.append(chr(ord("A") + (n % 26)))
        n = n // 26 - 1
        if n < 0:
            break
    return "".join(reversed(chars))


def _serial_to_site_atom_names(atoms):
    by_name = defaultdict(list)
    for serial, atom in sorted(atoms.items()):
        by_name[atom.get("atom_name", "")].append(serial)

    out = {}
    for atom_name, serials in by_name.items():
        if len(serials) == 1:
            out[serials[0]] = atom_name
            continue
        for i, serial in enumerate(serials):
            out[serial] = f"{atom_name}{_alpha_suffix(i)}"
    return out


def _resolve_template_feature_ids(
    template, feature_id, component, atoms, serial_to_site_name
):
    spec = template.get("features", {}).get(feature_id, {})
    entries = spec.get("atoms", [])
    by_name = defaultdict(list)
    for serial, atom in sorted(atoms.items()):
        by_name[atom.get("atom_name", "")].append(serial)

    out = []
    for e in entries:
        name = e.get("name")
        occurrence = int(e.get("occurrence", 1))
        idx = occurrence - 1
        serials = by_name.get(name, [])
        if 0 <= idx < len(serials):
            out.append(sid(component, serial_to_site_name[serials[idx]]))
    return sorted(set(out))


def _build_ring_impropers_from_template(graph, serial_to_atom, center_serials):
    ring_atoms = ring_atoms_from_graph(graph, serial_to_atom)
    impropers = []
    seen = set()
    all_serials = set(serial_to_atom.keys())
    for center in sorted(center_serials):
        if center not in ring_atoms:
            continue
        nbrs = sorted(n for n in graph.get(center, set()) if n in all_serials)
        if len(nbrs) < 3:
            continue
        ring_nbrs = [n for n in nbrs if n in ring_atoms]
        if len(ring_nbrs) < 2:
            continue
        n1, n2 = ring_nbrs[:2]
        n3_choices = [n for n in nbrs if n not in (n1, n2)]
        if not n3_choices:
            continue
        non_ring = [n for n in n3_choices if n not in ring_atoms]
        n3 = non_ring[0] if non_ring else n3_choices[0]
        key = (center, frozenset((n1, n2, n3)))
        if key in seen:
            continue
        seen.add(key)
        impropers.append((n1, center, n2, n3))
    return impropers


def _build_pi_impropers_from_template(graph, serial_to_atom, center_serials):
    impropers = []
    seen = set()
    all_serials = set(serial_to_atom.keys())
    for center in sorted(center_serials):
        atom = serial_to_atom.get(center)
        if not atom or atom.get("element") not in {"C", "N"}:
            continue
        nbrs = sorted(n for n in graph.get(center, set()) if n in all_serials)
        if len(nbrs) != 3:
            continue
        n1, n2, n3 = nbrs
        key = (center, frozenset((n1, n2, n3)))
        if key in seen:
            continue
        seen.add(key)
        impropers.append((n1, center, n2, n3))
    return impropers


def _build_flat_impropers_from_template(graph, serial_to_atom, center_serials):
    impropers = []
    for s in sorted(center_serials):
        atom = serial_to_atom.get(s, {})
        if not atom.get("atom_name", "").strip().upper().startswith("S"):
            continue
        o_nbrs = [
            n
            for n in sorted(graph.get(s, set()))
            if serial_to_atom.get(n, {}).get("element") == "O"
        ]
        if len(o_nbrs) < 3:
            continue
        impropers.append((o_nbrs[0], s, o_nbrs[1], o_nbrs[2]))
    return impropers


def _build_orient_impropers_from_template(graph, serial_to_atom, center_serials):
    impropers = []
    for s in sorted(center_serials):
        atom = serial_to_atom.get(s, {})
        if not atom.get("atom_name", "").strip().upper().startswith("S"):
            continue
        nbrs = list(sorted(graph.get(s, set())))
        c_nbrs = [n for n in nbrs if serial_to_atom.get(n, {}).get("element") == "C"]
        o_nbrs = [n for n in nbrs if serial_to_atom.get(n, {}).get("element") == "O"]
        if len(c_nbrs) != 1 or len(o_nbrs) < 2:
            continue
        impropers.append((c_nbrs[0], s, o_nbrs[0], o_nbrs[1]))
    return impropers


def _center_atom_serials_from_template(
    template, improper_type, component_atoms, site_names
):
    center_names = set()
    for imp in template.get("impropers", []):
        if imp.get("type") == improper_type:
            center_names.add(imp.get("center_atom"))
    center_serials = set()
    for serial, atom in component_atoms.items():
        if atom.get("atom_name") in center_names:
            center_serials.add(serial)
    return center_serials


def _parse_component_spec(spec_str):
    parts = {}
    for part in spec_str.split(","):
        if "=" in part:
            k, v = part.split("=", 1)
            parts[k.strip()] = v.strip()
    return parts


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--component",
    "components",
    multiple=True,
    required=True,
    help="Component spec: name=X,mol2=Y,template=Z,role=fixed|mobile (repeatable)",
)
@click.option("--output-cif", type=click.Path(path_type=str), required=True)
@click.option("--bond-k", type=float, default=2000.0, show_default=True)
@click.option("--angle-k", type=float, default=400.0, show_default=True)
@click.option("--pi-dihedral-k", type=float, default=12.0, show_default=True)
@click.option("--linker-dihedral-k", type=float, default=1.5, show_default=True)
@click.option("--ring-improper-k", type=float, default=40.0, show_default=True)
@click.option("--pi-improper-k", type=float, default=180.0, show_default=True)
@click.option("--flat-improper-k", type=float, default=120.0, show_default=True)
@click.option("--orient-improper-k", type=float, default=220.0, show_default=True)
@click.option("--n-steps", type=int, default=500000, show_default=True)
@click.option("--write-every", type=int, default=1000, show_default=True)
@click.option(
    "--default-radius",
    type=float,
    default=1.7,
    show_default=True,
    help="Default site radius in Angstrom.",
)
@click.option(
    "--default-mass",
    type=float,
    default=12.0,
    show_default=True,
    help="Default site mass in Da.",
)
@click.option("--nonbonded-k", type=float, default=5.0, show_default=True)
@click.option(
    "--nonbonded-cutoff",
    type=float,
    default=6.0,
    show_default=True,
    help="Nonbonded cutoff in Angstrom.",
)
@click.option("--minimize-steps", type=int, default=200, show_default=True)
def main(
    components,
    output_cif,
    bond_k,
    angle_k,
    pi_dihedral_k,
    linker_dihedral_k,
    ring_improper_k,
    pi_improper_k,
    flat_improper_k,
    orient_improper_k,
    n_steps,
    write_every,
    default_radius,
    default_mass,
    nonbonded_k,
    nonbonded_cutoff,
    minimize_steps,
):
    parsed_components = []
    for spec in components:
        parts = _parse_component_spec(spec)
        if "name" not in parts or "mol2" not in parts or "role" not in parts:
            raise click.ClickException(
                "Component spec must include name, mol2, template, role"
            )
        parsed_components.append(
            {
                "name": parts["name"],
                "mol2": parts["mol2"],
                "template": parts.get("template"),
                "role": parts["role"],
            }
        )

    fixed_comps = [c for c in parsed_components if c["role"] == "fixed"]
    mobile_comps = [c for c in parsed_components if c["role"] == "mobile"]

    if len(fixed_comps) != 1:
        raise click.ClickException("Exactly one fixed component required")

    fixed_name = fixed_comps[0]["name"]
    fixed_mol2 = fixed_comps[0]["mol2"]
    fixed_template_path = fixed_comps[0].get("template")

    fixed_atoms, fixed_bonds = parse_mol2(fixed_mol2, fixed_name)
    fixed_site_names = _serial_to_site_atom_names(fixed_atoms)
    print(f"Fixed ({fixed_name}): {len(fixed_atoms)} atoms, {len(fixed_bonds)} bonds")

    templates = {}
    if fixed_template_path:
        templates[fixed_name] = read_cgdye_template(fixed_template_path)

    all_atoms = {(fixed_name, k): v for k, v in fixed_atoms.items()}
    all_site_names = {fixed_name: fixed_site_names}
    all_graphs = {fixed_name: build_graph(fixed_bonds)}

    for comp in mobile_comps:
        name = comp["name"]
        mol2 = comp["mol2"]
        template_path = comp.get("template")

        atoms, bonds = parse_mol2(mol2, name)
        site_names = _serial_to_site_atom_names(atoms)

        print(f"Mobile ({name}): {len(atoms)} atoms, {len(bonds)} bonds")

        all_atoms.update({(name, k): v for k, v in atoms.items()})
        all_site_names[name] = site_names
        all_graphs[name] = build_graph(bonds)

        if template_path:
            templates[name] = read_cgdye_template(template_path)

    system_name_parts = [c["name"] for c in parsed_components]
    system_name = "_".join(system_name_parts)

    sites = []
    for comp_name, comp_atoms in [(fixed_name, fixed_atoms)] + [
        (c["name"], parse_mol2(c["mol2"], c["name"])[0]) for c in mobile_comps
    ]:
        for serial in sorted(comp_atoms.keys()):
            atom_name = all_site_names[comp_name][serial]
            sites.append(
                {
                    "id": sid(comp_name, atom_name),
                    "component": comp_name,
                    "atom_name": atom_name,
                    "site_serial": int(serial),
                    "radius": default_radius,
                    "mass": default_mass,
                }
            )

    bond_type = "B1"
    angle_type = "A1"
    torsion_type_pi = "T_PI"
    torsion_type_link = "T_LINK"
    improper_type_ring = "I_RING"
    improper_type_pi = "I_PI"
    improper_type_flat = "I_FLAT"
    improper_type_orient = "I_ORIENT"

    def sid_from_serial(component, serial):
        return sid(component, all_site_names[component][serial])

    bonds = []
    for comp in parsed_components:
        comp_name = comp["name"]
        comp_bonds = parse_mol2(comp["mol2"], comp_name)[1]
        for a, b in sorted(comp_bonds):
            bonds.append(
                [
                    sid_from_serial(comp_name, a),
                    sid_from_serial(comp_name, b),
                    distance(all_atoms[(comp_name, a)], all_atoms[(comp_name, b)]),
                    bond_type,
                ]
            )

    angles = []
    for comp in parsed_components:
        comp_name = comp["name"]
        comp_atoms = parse_mol2(comp["mol2"], comp_name)[0]
        g = build_graph(parse_mol2(comp["mol2"], comp_name)[1])
        for a, b, c in build_angles(g):
            angles.append(
                [
                    sid_from_serial(comp_name, a),
                    sid_from_serial(comp_name, b),
                    sid_from_serial(comp_name, c),
                    angle_value(
                        all_atoms[(comp_name, a)],
                        all_atoms[(comp_name, b)],
                        all_atoms[(comp_name, c)],
                    ),
                    angle_type,
                ]
            )

    dihedrals = []
    for comp in mobile_comps:
        comp_name = comp["name"]
        comp_atoms = parse_mol2(comp["mol2"], comp_name)[0]
        g = all_graphs[comp_name]

        pi_dihedrals = []
        linker_dihedrals = []
        for a, b, c, d in build_dihedrals(g):
            if all(s in comp_atoms for s in (a, b, c, d)):
                pi_dihedrals.append((a, b, c, d))

        for a, b, c, d in pi_dihedrals:
            atom_b = comp_atoms.get(b, {}).get("element", "")
            atom_c = comp_atoms.get(c, {}).get("element", "")
            if atom_b in {"C", "N"} and atom_c in {"C", "N"}:
                dihedrals.append(
                    [
                        sid_from_serial(comp_name, a),
                        sid_from_serial(comp_name, b),
                        sid_from_serial(comp_name, c),
                        sid_from_serial(comp_name, d),
                        torsion_type_pi,
                    ]
                )
            else:
                dihedrals.append(
                    [
                        sid_from_serial(comp_name, a),
                        sid_from_serial(comp_name, b),
                        sid_from_serial(comp_name, c),
                        sid_from_serial(comp_name, d),
                        torsion_type_link,
                    ]
                )

    impropers = []
    for comp in parsed_components:
        comp_name = comp["name"]
        comp_atoms = parse_mol2(comp["mol2"], comp_name)[0]
        g = all_graphs[comp_name]
        template = templates.get(comp_name)

        if template:
            center_ring = _center_atom_serials_from_template(
                template, "ring", comp_atoms, all_site_names[comp_name]
            )
            center_pi = _center_atom_serials_from_template(
                template, "pi", comp_atoms, all_site_names[comp_name]
            )
            center_flat = _center_atom_serials_from_template(
                template, "flat", comp_atoms, all_site_names[comp_name]
            )
            center_orient = _center_atom_serials_from_template(
                template, "orient", comp_atoms, all_site_names[comp_name]
            )

            for a, b, c, d in _build_ring_impropers_from_template(
                g, comp_atoms, center_ring
            ):
                impropers.append(
                    [
                        sid_from_serial(comp_name, a),
                        sid_from_serial(comp_name, b),
                        sid_from_serial(comp_name, c),
                        sid_from_serial(comp_name, d),
                        improper_type_ring,
                    ]
                )

            for a, b, c, d in _build_pi_impropers_from_template(
                g, comp_atoms, center_pi
            ):
                impropers.append(
                    [
                        sid_from_serial(comp_name, a),
                        sid_from_serial(comp_name, b),
                        sid_from_serial(comp_name, c),
                        sid_from_serial(comp_name, d),
                        improper_type_pi,
                    ]
                )

            for a, b, c, d in _build_flat_impropers_from_template(
                g, comp_atoms, center_flat
            ):
                impropers.append(
                    [
                        sid_from_serial(comp_name, a),
                        sid_from_serial(comp_name, b),
                        sid_from_serial(comp_name, c),
                        sid_from_serial(comp_name, d),
                        improper_type_flat,
                    ]
                )

            for a, b, c, d in _build_orient_impropers_from_template(
                g, comp_atoms, center_orient
            ):
                impropers.append(
                    [
                        sid_from_serial(comp_name, a),
                        sid_from_serial(comp_name, b),
                        sid_from_serial(comp_name, c),
                        sid_from_serial(comp_name, d),
                        improper_type_orient,
                    ]
                )

    groups = {}
    rb_groups = {}
    md_fixed_groups = {}

    for comp in parsed_components:
        comp_name = comp["name"]
        comp_atoms = parse_mol2(comp["mol2"], comp_name)[0]
        template = templates.get(comp_name)

        groups[f"{comp_name}_all"] = [
            sid_from_serial(comp_name, k) for k in sorted(comp_atoms.keys())
        ]

        if template:
            for fid, spec in template.get("features", {}).items():
                if spec.get("feature_type") == "dof":
                    ids = _resolve_template_feature_ids(
                        template, fid, comp_name, comp_atoms, all_site_names[comp_name]
                    )
                    if ids:
                        if spec.get("rb"):
                            rb_groups[f"{comp_name}_{fid}_rb"] = ids
                        if spec.get("md_fixed"):
                            md_fixed_groups[f"{comp_name}_{fid}_md_fixed"] = ids
                        if not spec.get("rb") and not spec.get("md_fixed"):
                            groups[f"{comp_name}_{fid}"] = ids

    fixed_groups = [f"{c['name']}_all" for c in fixed_comps]

    probes = []
    labeling_sites = []
    
    for i, c in enumerate(parsed_components, start=1):
        if c["role"] == "mobile":
            probes.append({
                "id": len(probes) + 1,
                "name": c["name"],
                "origin": "extrinsic",
                "link_type": "covalent"
            })
        elif c["role"] == "fixed":
            # For fixed components (proteins), we might want to record specific residues
            # as potential labeling sites. For now, let's assume if it's 'fixed' and
            # we know the residue (from previous steps/context), we add it.
            # Since this builder is general, we might need to pass labeling info via CLI.
            pass

    from .dye import CHARMM36_LJ
    lj_types = {}
    all_elements = set()
    for s in sites:
        ename = _element_from_atom_name(s["atom_name"])
        all_elements.add(ename)
    
    for elem in sorted(all_elements):
        params = CHARMM36_LJ.get(elem, CHARMM36_LJ["C"])
        lj_types[f"LJ_{elem}"] = {
            "element": elem,
            "rmin_half": params["rmin_half"],
            "epsilon": params["epsilon"],
        }

    out_dir = os.path.dirname(os.path.abspath(output_cif))
    system = {
        "name": system_name,
        "components": {
            c["name"]: {
                "mol2": os.path.relpath(os.path.abspath(c["mol2"]), out_dir),
                "role": c["role"],
            }
            for c in parsed_components
        },
        "probes": probes,
        "labeling_sites": labeling_sites,
        "sites": sites,
        "groups": groups,
        "rb_groups": rb_groups,
        "md_fixed_groups": md_fixed_groups,
        "fixed_groups": fixed_groups,
        "bond_types": {bond_type: {"k": bond_k}},
        "angle_types": {angle_type: {"k": angle_k}},
        "torsion_types": {
            torsion_type_pi: {
                "periodicity": 2,
                "phase_rad": math.pi,
                "k": pi_dihedral_k,
            },
            torsion_type_link: {
                "periodicity": 3,
                "phase_rad": 0.0,
                "k": linker_dihedral_k,
            },
        },
        "improper_types": {
            improper_type_ring: {
                "periodicity": 2,
                "phase_rad": 0.0,
                "k": ring_improper_k,
            },
            improper_type_pi: {"periodicity": 2, "phase_rad": 0.0, "k": pi_improper_k},
            improper_type_flat: {
                "periodicity": 2,
                "phase_rad": 0.0,
                "k": flat_improper_k,
            },
            improper_type_orient: {
                "periodicity": 2,
                "phase_rad": 0.0,
                "k": orient_improper_k,
            },
        },
        "bonds": bonds,
        "angles": angles,
        "dihedrals": dihedrals,
        "impropers": impropers,
        "nonbonded": {"enabled": True, "k": nonbonded_k, "cutoff_A": nonbonded_cutoff},
        "lj_types": lj_types,
        "sampling": {
            "n_steps": n_steps,
            "write_every": write_every,
            "minimize_steps": minimize_steps,
        },
    }

    os.makedirs(os.path.dirname(os.path.abspath(output_cif)), exist_ok=True)
    write_ff_system(output_cif, system)
    print(f"Wrote {output_cif}")
    print(
        f"sites={len(sites)} bonds={len(bonds)} angles={len(angles)} "
        f"dihedrals={len(dihedrals)} impropers={len(impropers)}"
    )


if __name__ == "__main__":
    main()
