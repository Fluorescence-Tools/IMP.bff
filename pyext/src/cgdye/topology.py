"""Topology builders for cgdye."""

from collections import defaultdict
import math
import os
import re


from IMP.bff.io.cif import read_component_template_cif, write_dye_forcefield_cif
from IMP.bff.scoring import CHARMM36_LJ, build_lj_type_table
import IMP
import IMP.algebra
import IMP.atom
import IMP.core

# --------------------------------------------------------------------------
# dye
# --------------------------------------------------------------------------
"""Dye topology builder utilities.

The Lennard-Jones table and kernel moved to :mod:`IMP.bff.scoring`
in the 2026-08-18 cleanup and are re-exported here, so the "one LJ source"
property this module established still holds: every scorer
(``sampling.scoring``, ``sampling.mean_field``,
``scoring.rotamer``) reads the same table, and the two that used
to sit side by side (``rmin_half``/``epsilon`` here, ``p_Rmin2``/``eps`` in the
rotamer scorer) still cannot drift apart.

They moved because ``rotamer/`` became part of ``representation`` and a core
domain must not import a legacy package to compute a steric term.
"""

#: CHARMM36 per-element LJ parameters: ``rmin_half`` (Rmin/2, Angstrom) and
#: ``epsilon`` (kcal/mol, negative as in the CHARMM parameter files).


















def build_dye_topology(atoms, bonds, template):
    """Bonds, angles, dihedrals and template-driven impropers of a dye from its MOL2 graph.

    ``atoms``/``bonds`` as returned by ``parse_dye_mol2``; ``template`` a dye
    template dict (``impropers`` entries with ``center_atom`` and ``type``).
    """
    graph = build_graph(bonds)
    result = {"bonds": [], "angles": [], "dihedrals": [], "impropers": []}

    for a, b in sorted(bonds):
        result["bonds"].append((a, b, _distance(atoms[a], atoms[b])))

    for a, b, c in build_angles(graph):
        result["angles"].append((a, b, c, _angle_value(atoms[a], atoms[b], atoms[c])))

    for a, b, c, d in build_dihedrals(graph):
        result["dihedrals"].append((a, b, c, d))

    result["impropers"] = _build_impropers(atoms, graph, template)
    return result


def _build_impropers(atoms, graph, template):
    out = []
    name_to_serial = defaultdict(list)
    for serial, atom in atoms.items():
        name_to_serial[atom["atom_name"]].append(serial)

    for imp in template.get("impropers", []):
        center_name = imp["center_atom"]
        imp_type = imp["type"]
        for center in name_to_serial.get(center_name, []):
            if imp_type == "ring":
                out.extend(_ring_impropers(center, graph))
            elif imp_type == "pi":
                out.extend(_pi_impropers(center, graph))
            elif imp_type == "flat":
                out.extend(_flat_impropers(center, graph, atoms))
            elif imp_type == "orient":
                out.extend(_orient_impropers(center, graph, atoms))
    return out


def _ring_impropers(center, graph):
    cycles = find_cycles(graph)
    ring_atoms = set()
    for cyc in cycles:
        ring_atoms.update(cyc)
    if center not in ring_atoms:
        return []
    nbrs = sorted(graph.get(center, set()))
    ring_nbrs = [n for n in nbrs if n in ring_atoms]
    if len(ring_nbrs) < 2:
        return []
    non_ring = [n for n in nbrs if n not in ring_nbrs]
    wing3 = non_ring[0] if non_ring else (nbrs[2] if len(nbrs) >= 3 else None)
    if wing3 is None:
        return []
    return [(ring_nbrs[0], center, ring_nbrs[1], wing3)]


def _pi_impropers(center, graph):
    nbrs = sorted(graph.get(center, set()))
    if len(nbrs) != 3:
        return []
    return [(nbrs[0], center, nbrs[1], nbrs[2])]


def _flat_impropers(center, graph, atoms):
    if atoms[center]["element"] != "S":
        return []
    o_nbrs = sorted(n for n in graph.get(center, set()) if atoms[n]["element"] == "O")
    if len(o_nbrs) < 3:
        return []
    return [(o_nbrs[0], center, o_nbrs[1], o_nbrs[2])]


def _orient_impropers(center, graph, atoms):
    if atoms[center]["element"] != "S":
        return []
    c_nbrs = sorted(n for n in graph.get(center, set()) if atoms[n]["element"] == "C")
    o_nbrs = sorted(n for n in graph.get(center, set()) if atoms[n]["element"] == "O")
    if not c_nbrs or len(o_nbrs) < 2:
        return []
    return [(c_nbrs[0], center, o_nbrs[0], o_nbrs[1])]


def _distance(a, b):
    return (
        (a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2
    ) ** 0.5


def _angle_value(a, b, c):
    ba = [a["x"] - b["x"], a["y"] - b["y"], a["z"] - b["z"]]
    bc = [c["x"] - b["x"], c["y"] - b["y"], c["z"] - b["z"]]
    dot = sum(x * y for x, y in zip(ba, bc))
    mag_ba = sum(x**2 for x in ba) ** 0.5
    mag_bc = sum(x**2 for x in bc) ** 0.5
    if mag_ba == 0 or mag_bc == 0:
        return 0.0
    cos_theta = max(-1.0, min(1.0, dot / (mag_ba * mag_bc)))
    return math.acos(cos_theta)


# --------------------------------------------------------------------------
# builder
# --------------------------------------------------------------------------
"""Build compact FF mmCIF system from N-component MOL2 files.

Uses IMP.atom.read_mol2() for bond topology, then reads template CIFs for
feature definitions and improper centers.

Usage:
    python -m cgdye.topology.builder         --component name=hgbp1,mol2=inputs/structures/1DG3.mol2,role=fixed         --component name=alexa488,mol2=inputs/structures/alexa488_r48.mol2,role=mobile         --output-cif output/systems/hgbp1_alexa488.system.cif
"""

# --------------------------------------------------------------------------
# builder
# --------------------------------------------------------------------------
"""Build compact FF mmCIF system from N-component MOL2 files.

Uses IMP.atom.read_mol2() for bond topology, then reads template CIFs for
feature definitions and improper centers.

Usage:
    python -m cgdye.topology.builder \
        --component name=hgbp1,mol2=inputs/structures/1DG3.mol2,role=fixed \
        --component name=alexa488,mol2=inputs/structures/alexa488_r48.mol2,role=mobile \
        --output-cif output/systems/hgbp1_alexa488.system.cif
"""

#!/usr/bin/env python






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


def parse_dye_mol2(path, component):
    """Read a MOL2 file into ``(atoms, bonds)``.

    ``atoms`` maps the MOL2 serial to a dict (serial, component, atom_name,
    resname, element, x, y, z); ``bonds`` is a set of sorted serial pairs.
    """
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


# The molecular-graph primitives -- ``build_graph``, ``build_angles``,
# ``build_dihedrals``, ``find_cycles`` -- existed **twice**, once in
# ``topology/dye.py`` and once in ``topology/builder.py``, with identical logic
# and different variable names. Merging the two modules surfaced that; only the
# second copy was ever reachable, because it shadowed the first.
#
# Checked before deleting one: on 400 random graphs the two agree exactly on
# the graph and the angles, and agree on the dihedrals **up to direction** --
# ``(1, 2, 4, 0)`` against ``(0, 4, 2, 1)`` is one torsion written backwards,
# and a torsion angle is the same either way. The copy kept is the one that was
# already winning.

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


# `_resolve_template_feature_ids` was here, identical to `_resolve_feature_ids`
# below but for a `defaultdict` where that one uses `setdefault`, and parameter
# names. The two builders in this file had one each.


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


# The click decorators for this moved to `bin/imp_bff` as
# `build-system`; what is left is the function they called, which is
# library code and now importable without click.
def build_system_from_specs(
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
            raise ValueError(
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
        raise ValueError("Exactly one fixed component required")

    fixed_name = fixed_comps[0]["name"]
    fixed_mol2 = fixed_comps[0]["mol2"]
    fixed_template_path = fixed_comps[0].get("template")

    fixed_atoms, fixed_bonds = parse_dye_mol2(fixed_mol2, fixed_name)
    fixed_site_names = _serial_to_site_atom_names(fixed_atoms)
    print(f"Fixed ({fixed_name}): {len(fixed_atoms)} atoms, {len(fixed_bonds)} bonds")

    templates = {}
    if fixed_template_path:
        templates[fixed_name] = read_component_template_cif(fixed_template_path)

    all_atoms = {(fixed_name, k): v for k, v in fixed_atoms.items()}
    all_site_names = {fixed_name: fixed_site_names}
    all_graphs = {fixed_name: build_graph(fixed_bonds)}

    for comp in mobile_comps:
        name = comp["name"]
        mol2 = comp["mol2"]
        template_path = comp.get("template")

        atoms, bonds = parse_dye_mol2(mol2, name)
        site_names = _serial_to_site_atom_names(atoms)

        print(f"Mobile ({name}): {len(atoms)} atoms, {len(bonds)} bonds")

        all_atoms.update({(name, k): v for k, v in atoms.items()})
        all_site_names[name] = site_names
        all_graphs[name] = build_graph(bonds)

        if template_path:
            templates[name] = read_component_template_cif(template_path)

    system_name_parts = [c["name"] for c in parsed_components]
    system_name = "_".join(system_name_parts)

    sites = []
    for comp_name, comp_atoms in [(fixed_name, fixed_atoms)] + [
        (c["name"], parse_dye_mol2(c["mol2"], c["name"])[0]) for c in mobile_comps
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
        comp_bonds = parse_dye_mol2(comp["mol2"], comp_name)[1]
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
        comp_atoms = parse_dye_mol2(comp["mol2"], comp_name)[0]
        g = build_graph(parse_dye_mol2(comp["mol2"], comp_name)[1])
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
        comp_atoms = parse_dye_mol2(comp["mol2"], comp_name)[0]
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
        comp_atoms = parse_dye_mol2(comp["mol2"], comp_name)[0]
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
        comp_atoms = parse_dye_mol2(comp["mol2"], comp_name)[0]
        template = templates.get(comp_name)

        groups[f"{comp_name}_all"] = [
            sid_from_serial(comp_name, k) for k in sorted(comp_atoms.keys())
        ]

        if template:
            for fid, spec in template.get("features", {}).items():
                if spec.get("feature_type") == "dof":
                    ids = _resolve_feature_ids(
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

    # (was: from .topology import ...) -- now in this module
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
    write_dye_forcefield_cif(output_cif, system)
    print(f"Wrote {output_cif}")
    print(
        f"sites={len(sites)} bonds={len(bonds)} angles={len(angles)} "
        f"dihedrals={len(dihedrals)} impropers={len(impropers)}"
    )




# --------------------------------------------------------------------------
# combined
# --------------------------------------------------------------------------
"""Build combined protein+dye systems."""

def _resolve_feature_ids(template, feature_id, comp_name, atoms, serial_to_site_name):
    entries = template.get("features", {}).get(feature_id, {}).get("atoms", [])
    by_name = {}
    for serial, atom in sorted(atoms.items()):
        by_name.setdefault(atom.get("atom_name", ""), []).append(serial)
    out = []
    for e in entries:
        name = e.get("name")
        occurrence = int(e.get("occurrence", 1))
        idx = occurrence - 1
        serials = by_name.get(name, [])
        if 0 <= idx < len(serials):
            out.append(sid(comp_name, serial_to_site_name[serials[idx]]))
    return sorted(set(out))


def build_dye_protein_system(
    protein_mol2,
    dye_mol2,
    protein_name,
    dye_name,
    protein_template=None,
    dye_template=None,
    default_radius=1.7,
    default_mass=12.0,
):
    """Build one force-field system for a protein and a dye from their MOL2 files.

    Bonded terms come from the MOL2 connectivity (and the optional templates),
    non-bonded LJ types from the CHARMM36 table; the result is the dict
    ``write_dye_forcefield_cif`` writes and ``sim.runner`` simulates.
    """
    p_atoms, p_bonds = parse_dye_mol2(protein_mol2, protein_name)
    d_atoms, d_bonds = parse_dye_mol2(dye_mol2, dye_name)
    p_names = _serial_to_site_atom_names(p_atoms)
    d_names = _serial_to_site_atom_names(d_atoms)
    p_graph = build_graph(p_bonds)
    d_graph = build_graph(d_bonds)

    p_template = read_component_template_cif(protein_template) if protein_template else {}
    d_template = read_component_template_cif(dye_template) if dye_template else {}

    sites = []
    for serial in sorted(p_atoms):
        aname = p_names[serial]
        sites.append(
            {
                "id": sid(protein_name, aname),
                "component": protein_name,
                "atom_name": aname,
                "site_serial": serial,
                "radius": default_radius,
                "mass": default_mass,
            }
        )
    for serial in sorted(d_atoms):
        aname = d_names[serial]
        sites.append(
            {
                "id": sid(dye_name, aname),
                "component": dye_name,
                "atom_name": aname,
                "site_serial": serial,
                "radius": default_radius,
                "mass": default_mass,
            }
        )

    bonds = []
    for a, b in sorted(p_bonds):
        bonds.append(
            [
                sid(protein_name, p_names[a]),
                sid(protein_name, p_names[b]),
                distance(p_atoms[a], p_atoms[b]),
                "B1",
            ]
        )
    for a, b in sorted(d_bonds):
        bonds.append(
            [
                sid(dye_name, d_names[a]),
                sid(dye_name, d_names[b]),
                distance(d_atoms[a], d_atoms[b]),
                "B1",
            ]
        )

    angles = []
    for comp_name, comp_atoms, comp_names, comp_graph in [
        (protein_name, p_atoms, p_names, p_graph),
        (dye_name, d_atoms, d_names, d_graph),
    ]:
        for a, b, c in build_angles(comp_graph):
            angles.append(
                [
                    sid(comp_name, comp_names[a]),
                    sid(comp_name, comp_names[b]),
                    sid(comp_name, comp_names[c]),
                    angle_value(comp_atoms[a], comp_atoms[b], comp_atoms[c]),
                    "A1",
                ]
            )

    dihedrals = []
    for a, b, c, d in build_dihedrals(d_graph):
        eb = d_atoms.get(b, {}).get("element", "")
        ec = d_atoms.get(c, {}).get("element", "")
        tt = "T_PI" if eb in {"C", "N"} and ec in {"C", "N"} else "T_LINK"
        dihedrals.append(
            [
                sid(dye_name, d_names[a]),
                sid(dye_name, d_names[b]),
                sid(dye_name, d_names[c]),
                sid(dye_name, d_names[d]),
                tt,
            ]
        )

    impropers = []

    groups = {
        f"{protein_name}_all": [sid(protein_name, p_names[k]) for k in sorted(p_atoms)],
        f"{dye_name}_all": [sid(dye_name, d_names[k]) for k in sorted(d_atoms)],
    }
    rb_groups = {}
    md_fixed_groups = {}

    for comp_name, templ, atoms, names in [
        (protein_name, p_template, p_atoms, p_names),
        (dye_name, d_template, d_atoms, d_names),
    ]:
        for fid, spec in templ.get("features", {}).items():
            ids = _resolve_feature_ids(templ, fid, comp_name, atoms, names)
            if not ids:
                continue
            if spec.get("rb"):
                rb_groups[f"{comp_name}_{fid}_rb"] = ids
            elif spec.get("md_fixed"):
                md_fixed_groups[f"{comp_name}_{fid}_md_fixed"] = ids
            else:
                groups[f"{comp_name}_{fid}"] = ids

    elements = set()
    for s in sites:
        m = re.match(r"([A-Za-z])", s["atom_name"])
        if m:
            elements.add(m.group(1).upper())

    from IMP.bff.io.cif import forcefield_system_from_dict
    return forcefield_system_from_dict(
    {
            "name": f"{protein_name}_{dye_name}",
            "components": {
                protein_name: {"mol2": protein_mol2, "role": "fixed"},
                dye_name: {"mol2": dye_mol2, "role": "mobile"},
            },
            "sites": sites,
            "groups": groups,
            "rb_groups": rb_groups,
            "md_fixed_groups": md_fixed_groups,
            "fixed_groups": [f"{protein_name}_all"],
            "bond_types": {"B1": {"k": 2000.0}},
            "angle_types": {"A1": {"k": 400.0}},
            "torsion_types": {
                "T_PI": {"periodicity": 2, "phase_rad": math.pi, "k": 12.0},
                "T_LINK": {"periodicity": 3, "phase_rad": 0.0, "k": 1.5},
            },
            "improper_types": {
                "I_RING": {"periodicity": 2, "phase_rad": 0.0, "k": 40.0},
                "I_PI": {"periodicity": 2, "phase_rad": 0.0, "k": 180.0},
                "I_FLAT": {"periodicity": 2, "phase_rad": 0.0, "k": 120.0},
                "I_ORIENT": {"periodicity": 2, "phase_rad": 0.0, "k": 220.0},
            },
            "lj_types": build_lj_type_table(elements),
            "bonds": bonds,
            "angles": angles,
            "dihedrals": dihedrals,
            "impropers": impropers,
            "nonbonded": {"enabled": True, "k": 5.0, "cutoff_A": 6.0},
            "sampling": {"n_steps": 500000, "write_every": 1000, "minimize_steps": 200},
        }
    )


def dye_forcefield_system(dye_mol2, dye_name="dye", dye_template=None, default_radius=1.7, default_mass=12.0):
    """A force-field system dict for one dye alone (no protein component).

    The same bonded/non-bonded terms and types as :func:`build_dye_protein_system`
    gives the dye; the dye's backbone-anchor atoms (``N``, ``CA``, ``C``, ``O``)
    are collected in ``fixed_groups`` so a sampler can hold them on the
    labelled residue. Site ids are ``<dye_name>:<atom_name>`` (unique via
    ``_serial_to_site_atom_names``).
    """
    d_atoms, d_bonds = parse_dye_mol2(dye_mol2, dye_name)
    d_names = _serial_to_site_atom_names(d_atoms)
    d_graph = build_graph(d_bonds)
    d_template = read_component_template_cif(dye_template) if dye_template else {}
    sites = []
    for serial in sorted(d_atoms):
        aname = d_names[serial]
        elem = d_atoms[serial].get("element", "C") or "C"
        sites.append({"id": sid(dye_name, aname), "component": dye_name, "atom_name": aname,
                      "site_serial": serial, "radius": default_radius, "mass": _ELEMENT_MASS.get(elem, default_mass),
                      "element": elem})
    bonds = [[sid(dye_name, d_names[a]), sid(dye_name, d_names[b]), distance(d_atoms[a], d_atoms[b]), "B1"]
             for a, b in sorted(d_bonds)]
    angles = [[sid(dye_name, d_names[a]), sid(dye_name, d_names[b]), sid(dye_name, d_names[c]),
               angle_value(d_atoms[a], d_atoms[b], d_atoms[c]), "A1"] for a, b, c in build_angles(d_graph)]
    dihedrals = []
    for a, b, c, d in build_dihedrals(d_graph):
        eb = d_atoms.get(b, {}).get("element", "")
        ec = d_atoms.get(c, {}).get("element", "")
        tt = "T_PI" if eb in {"C", "N"} and ec in {"C", "N"} else "T_LINK"
        dihedrals.append([sid(dye_name, d_names[a]), sid(dye_name, d_names[b]), sid(dye_name, d_names[c]), sid(dye_name, d_names[d]), tt])
    anchor = [sid(dye_name, d_names[k]) for k in sorted(d_atoms) if d_names[k].upper() in {"N", "CA", "C", "O"}]
    elements = {s["element"] for s in sites}
    from IMP.bff.io.cif import forcefield_system_from_dict
    return forcefield_system_from_dict(
    {
            "name": dye_name,
            "components": {dye_name: {"mol2": dye_mol2, "role": "mobile"}},
            "sites": sites,
            "groups": {f"{dye_name}_all": [s["id"] for s in sites], f"{dye_name}_anchor": anchor},
            "rb_groups": {}, "md_fixed_groups": {},
            "fixed_groups": [f"{dye_name}_anchor"],
            "bond_types": {"B1": {"k": 2000.0}},
            "angle_types": {"A1": {"k": 400.0}},
            "torsion_types": {
                "T_PI": {"periodicity": 2, "phase_rad": math.pi, "k": 12.0},
                "T_LINK": {"periodicity": 3, "phase_rad": 0.0, "k": 1.5},
            },
            "improper_types": {},
            "lj_types": build_lj_type_table(elements),
            "bonds": bonds, "angles": angles, "dihedrals": dihedrals, "impropers": [],
        }
    )


#: atomic masses (Da) of the elements a dye MOL2 carries
_ELEMENT_MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06, "P": 30.974, "F": 18.998, "CL": 35.45, "BR": 79.904}
