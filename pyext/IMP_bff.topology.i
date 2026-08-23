/*
 * cgdye/topology.py to %pythoncode: the dye topology, the graph helpers,
 * and the thin wrappers over the C++ system builder. `build_forcefield_
 * system` itself is C++ (TopologyBuild.h) -- it reads the MOL2 files, derives
 * the terms and groups, and builds the typed system through the one JSON
 * conversion path. What stays Python: the dict-shaped `parse_dye_mol2`, the
 * graph dict views, and the CLI-shaped wrappers that call the builder.
 */

%include "IMP/bff/TopologyBuild.h"

%pythoncode %{
# C++ free functions from earlier .i files
read_mol2_component = _IMP_bff.read_mol2_component
element_from_atom_name = _IMP_bff.element_from_atom_name

from collections import defaultdict
import json
import math
import os
import re


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
        result["angles"].append((a, b, c, angle_value(atoms[a], atoms[b], atoms[c])))

    for a, b, c, d in build_dihedrals(graph):
        result["dihedrals"].append((a, b, c, d))

    result["impropers"] = _build_impropers(atoms, graph, template)
    return result


def _build_impropers(atoms, graph, template):
    """The template's improper centres expanded against the bond graph.

    One expander family per kind, shared with `build_forcefield_system`. There
    were two -- a per-centre set used here and a per-centre-list set used
    there, differing in bounds checks and a dedup -- and they agreed on every
    kind for atto655, cx4 and alexa488_r48. The list form is the one kept: it
    is the one whose output ships.
    """
    by_name = defaultdict(list)
    for serial, atom in atoms.items():
        by_name[atom["atom_name"]].append(serial)

    expand = {
        "ring": _build_ring_impropers_from_template,
        "pi": _build_pi_impropers_from_template,
        "flat": _build_flat_impropers_from_template,
        "orient": _build_orient_impropers_from_template,
    }
    centers_by_kind = defaultdict(list)
    for imp in template.impropers:
        centers_by_kind[imp.type].extend(by_name.get(imp.center_atom, []))

    out = []
    for kind, centers in centers_by_kind.items():
        if kind in expand:
            out.extend(expand[kind](graph, atoms, sorted(set(centers))))
    return out


def _distance(a, b):
    return (
        (a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2
    ) ** 0.5


# `_angle_value` was here: the same formula as `angle_value` below, differing
# only in what it returns when an arm has zero length -- 0.0 against that one's
# 1.910633 rad, the tetrahedral angle. Neither fires: 6,364 angles across the
# shipped structures have no zero-length arm. The tetrahedral fallback is the
# one kept, because a harmonic minimum at a collapsed angle is not a default
# anyone wants if it ever does.


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






def parse_dye_mol2(path, component):
    """Read a MOL2 file into ``(atoms, bonds)``.

    :func:`IMP.bff.read_mol2_component`. ``atoms`` maps the MOL2 serial to a
    dict (serial, component, atom_name, resname, element, x, y, z); ``bonds``
    is a set of sorted serial pairs.

    The atom *name* comes from column 2 of ``@<TRIPOS>ATOM`` rather than the
    TRIPOS type, because IMP's ``read_mol2`` maps ``C.3`` to ``C3`` and loses
    ``C12``/``N1`` -- the names a template's features and impropers refer to.

    Gated against the Python it replaced on atto655, cx4, alexa488_r48 and the
    4,698-atom 1DG3: identical atoms, identical bonds, every field equal.
    """
    component_data = IMP.bff.read_mol2_component(str(path), component)
    atoms = {
        a.serial: {
            "serial": a.serial,
            "component": a.component,
            "atom_name": a.atom_name,
            "resname": a.resname,
            "element": a.element,
            "x": a.x,
            "y": a.y,
            "z": a.z,
        }
        for a in component_data.atoms
    }
    return atoms, {tuple(b) for b in component_data.bonds}


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


def _as_int_graph(pairs):
    """A :class:`IMP.bff.MolecularGraph` over *pairs*, plus the way back.

    The graph is C++ and its nodes are ints; callers here have MOL2 serials
    (already ints) but `scoring` has site-id strings. Numbering the nodes in
    first-appearance order keeps one C++ implementation for both, and keeps the
    numbering stable so derived terms come out in the order the Python's did.

    The four derivations were Python -- adjacency, angles, torsions and rings,
    here and again in `cgdye.sim` and `scoring`. One implementation now, gated
    against the Python on atto655, cx4 and alexa488_r48: identical angles,
    torsions and rings on all three.
    """
    index, nodes = {}, []
    for a, b in pairs:
        for n in (a, b):
            if n not in index:
                index[n] = len(nodes)
                nodes.append(n)
    edges = [(index[a], index[b]) for a, b in pairs]
    return IMP.bff.MolecularGraph(edges), nodes


def _graph_of(graph):
    """*graph* as ``(MolecularGraph, nodes)``, numbered in sorted node order.

    Sorted, not first-appearance: the C++ emits an angle's ends ascending by
    node number and a torsion in whichever direction is smaller, so numbering
    by sorted node is what makes those come back ascending by *node* -- which
    is what the Python did. Numbering by first appearance gave the same 138
    angles on atto655 with 8 of them reversed.
    """
    pairs = set()
    for node, nbrs in graph.items():
        for nb in nbrs:
            pairs.add((node, nb) if node <= nb else (nb, node))
    nodes = sorted({n for pair in pairs for n in pair})
    index = {n: i for i, n in enumerate(nodes)}
    return IMP.bff.MolecularGraph([(index[a], index[b]) for a, b in sorted(pairs)]), nodes


def build_graph(bonds):
    """{node: set of bonded nodes} -- kept as the dict shape callers index.

    The connectivity itself is :class:`IMP.bff.MolecularGraph`; this is the
    Python view of it.
    """
    # Key order is first-appearance in *bonds*, not ascending: this was a
    # `defaultdict(set)` filled by iterating the bonds, and `LinkerSampler`
    # walks the result, so the order reaches a seeded sampler and its pinned
    # weights. Returning C++'s ascending node order moved them.
    g, nodes = _as_int_graph([tuple(b) for b in bonds])
    return {nodes[i]: {nodes[j] for j in g.get_neighbors(i)}
            for i in range(len(nodes))}


def build_angles(graph):
    """Every ``(a, b, c)`` with *b* bonded to both -- see
    :meth:`IMP.bff.MolecularGraph.get_angles`."""
    g, nodes = _graph_of(graph)
    return sorted(tuple(nodes[i] for i in a) for a in g.get_angles())


def build_dihedrals(graph):
    """Every proper torsion, each once -- see
    :meth:`IMP.bff.MolecularGraph.get_dihedrals`."""
    g, nodes = _graph_of(graph)
    return sorted(tuple(nodes[i] for i in d) for d in g.get_dihedrals())


def find_cycles(graph, max_len=7):
    """Simple cycles of at most *max_len* nodes -- see
    :meth:`IMP.bff.MolecularGraph.get_rings`."""
    g, nodes = _graph_of(graph)
    return {tuple(nodes[i] for i in c) for c in g.get_rings(max_len)}


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


def _expand_impropers(kind, graph, serial_to_atom, center_serials):
    """Improper quadruples of *kind* about *center_serials*.

    :meth:`IMP.bff.MolecularGraph.expand_impropers`. There were four Python
    functions here -- ring, pi, flat, orient -- reading the same graph and the
    same per-atom element and name. Gated against them on atto655, cx4 and
    alexa488_r48: identical quadruples for every kind on all three.

    `flat` and `orient` test the *atom name* starting with S rather than the
    element, which is what the Python did and is deliberately preserved: a
    MOL2's types are less reliable than its names.
    """
    nodes = sorted(serial_to_atom)
    g, mapped = _as_int_graph(sorted(
        (a, b) for a, nbrs in graph.items() for b in nbrs if a < b))
    index = {n: i for i, n in enumerate(mapped)}
    keep = [n for n in nodes if n in index]
    quads = g.expand_impropers(
        kind,
        [index[n] for n in sorted(set(center_serials)) if n in index],
        [index[n] for n in keep],
        [serial_to_atom[n].get("element", "") for n in keep],
        [serial_to_atom[n].get("atom_name", "") for n in keep],
        8,
    )
    return [tuple(mapped[i] for i in q) for q in quads]


def _build_ring_impropers_from_template(graph, serial_to_atom, center_serials):
    return _expand_impropers("ring", graph, serial_to_atom, center_serials)


def _build_pi_impropers_from_template(graph, serial_to_atom, center_serials):
    return _expand_impropers("pi", graph, serial_to_atom, center_serials)


def _build_flat_impropers_from_template(graph, serial_to_atom, center_serials):
    return _expand_impropers("flat", graph, serial_to_atom, center_serials)


def _build_orient_impropers_from_template(graph, serial_to_atom, center_serials):
    return _expand_impropers("orient", graph, serial_to_atom, center_serials)


def _center_atom_serials_from_template(
    template, improper_type, component_atoms, site_names
):
    center_names = set()
    for imp in template.impropers:
        if imp.type == improper_type:
            center_names.add(imp.center_atom)
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
def build_forcefield_system(
    components,
    bond_k=2000.0,
    angle_k=400.0,
    pi_dihedral_k=12.0,
    linker_dihedral_k=1.5,
    ring_improper_k=40.0,
    pi_improper_k=180.0,
    flat_improper_k=120.0,
    orient_improper_k=220.0,
    n_steps=20000,
    write_every=100,
    default_radius=1.7,
    default_mass=12.0,
    nonbonded_k=5.0,
    nonbonded_cutoff=6.0,
    minimize_steps=200,
    relative_to=None,
):
    """Build a :class:`IMP.bff.DyeForceFieldSystem` from component specs.

    *components* is a list of ``name=X,mol2=Y,template=Z,role=fixed|mobile``
    strings, or of dicts with those keys. The builder itself is C++
    (`TopologyBuild.h`): it reads the MOL2s, derives the terms and groups and
    builds the typed system through the one JSON conversion path. *relative_to*
    is a directory the recorded mol2 paths are made relative to -- the writer
    passes the output file's directory.
    """
    specs = []
    for spec in components:
        parts = spec if isinstance(spec, dict) else _parse_component_spec(spec)
        specs.append({
            "name": parts["name"],
            "mol2": parts["mol2"],
            "template": parts.get("template"),
            "role": parts["role"],
        })
    return _IMP_bff.build_forcefield_system(
        json.dumps(specs), float(bond_k), float(angle_k),
        float(pi_dihedral_k), float(linker_dihedral_k),
        float(ring_improper_k), float(pi_improper_k),
        float(flat_improper_k), float(orient_improper_k),
        int(n_steps), int(write_every), float(default_radius),
        float(default_mass), float(nonbonded_k), float(nonbonded_cutoff),
        int(minimize_steps), relative_to or "")


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
    """`build_forcefield_system` written to *output_cif* -- the `build-system`
    command's body, and nothing more than the write."""
    system = build_forcefield_system(
        components, bond_k=bond_k, angle_k=angle_k, pi_dihedral_k=pi_dihedral_k,
        linker_dihedral_k=linker_dihedral_k, ring_improper_k=ring_improper_k,
        pi_improper_k=pi_improper_k, flat_improper_k=flat_improper_k,
        orient_improper_k=orient_improper_k, n_steps=n_steps,
        write_every=write_every, default_radius=default_radius,
        default_mass=default_mass, nonbonded_k=nonbonded_k,
        nonbonded_cutoff=nonbonded_cutoff, minimize_steps=minimize_steps,
        relative_to=os.path.dirname(os.path.abspath(output_cif)))
    os.makedirs(os.path.dirname(os.path.abspath(output_cif)), exist_ok=True)
    write_dye_forcefield_cif(output_cif, system)
    print(f"Wrote {output_cif}")
    print(
        f"sites={len(system.sites)} bonds={len(system.bonds)} "
        f"angles={len(system.angles)} dihedrals={len(system.dihedrals)} "
        f"impropers={len(system.impropers)}"
    )




# --------------------------------------------------------------------------
# combined
# --------------------------------------------------------------------------
"""Build combined protein+dye systems."""

def _resolve_feature_ids(template, feature_id, comp_name, atoms, serial_to_site_name):
    entries = (template.features[feature_id].atoms
               if feature_id in template.features else [])
    by_name = {}
    for serial, atom in sorted(atoms.items()):
        by_name.setdefault(atom.get("atom_name", ""), []).append(serial)
    out = []
    for e in entries:
        name = e.name
        occurrence = int(e.occurrence)
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
    """One force-field system for a protein and a dye, both from MOL2.

    A two-component :func:`build_forcefield_system`. It used to be a second
    implementation of that function -- 167 lines that agreed with it on sites,
    bonds, angles, dihedrals, groups and every type table, and differed only in
    hard-coding ``impropers = []``. Systems built here now carry the impropers
    the templates declare: 81 of them on the two shipped components, where this
    produced none. See okf/validation/impropers_are_dropped.md.
    """
    return build_forcefield_system(
        [
            {"name": protein_name, "mol2": protein_mol2,
             "template": protein_template, "role": "fixed"},
            {"name": dye_name, "mol2": dye_mol2,
             "template": dye_template, "role": "mobile"},
        ],
        default_radius=default_radius,
        default_mass=default_mass,
    )


def dye_forcefield_system(
    dye_mol2,
    dye_name="dye",
    dye_template=None,
    default_radius=1.7,
    default_mass=12.0,
):
    """A force-field system for one dye alone, with no protein component.

    A one-component :func:`build_forcefield_system`, and formerly a third
    implementation of it.

    The one thing that is this function's own: the dye's backbone-anchor atoms
    (``N``, ``CA``, ``C``, ``O``) are collected into a ``<dye_name>_anchor``
    group and named in ``fixed_groups``, so a sampler can hold them on the
    labelled residue. With no protein component there is nothing else fixed.
    """
    system = build_forcefield_system(
        [{"name": dye_name, "mol2": dye_mol2,
          "template": dye_template, "role": "mobile"}],
        default_radius=default_radius,
        default_mass=default_mass,
    )

    anchor = [site.id for site in system.sites
              if site.atom_name.upper() in {"N", "CA", "C", "O"}]
    groups = {k: list(v) for k, v in system.groups.items()}
    groups[f"{dye_name}_anchor"] = sorted(anchor)
    system.groups = groups
    system.fixed_groups = [f"{dye_name}_anchor"]
    return system


#: atomic masses (Da) of the elements a dye MOL2 carries
_ELEMENT_MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06, "P": 30.974, "F": 18.998, "CL": 35.45, "BR": 79.904}

%}
