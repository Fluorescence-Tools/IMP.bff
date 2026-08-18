"""Dye topology builder utilities.

The Lennard-Jones table and kernel moved to :mod:`IMP.bff.tools.lennard_jones`
in the 2026-08-18 cleanup and are re-exported here, so the "one LJ source"
property this module established still holds: every scorer
(``sampling.scoring``, ``sampling.mean_field``,
``representation.rotamer.scoring``) reads the same table, and the two that used
to sit side by side (``rmin_half``/``epsilon`` here, ``p_Rmin2``/``eps`` in the
rotamer scorer) still cannot drift apart.

They moved because ``rotamer/`` became part of ``representation`` and a core
domain must not import a legacy package to compute a steric term.
"""

# Sterics live in a tool, not here -- see above.
from IMP.bff.tools.lennard_jones import (  # noqa: F401
    CHARMM36_LJ,
    lj_cross,
    lj_params,
    lj_energy,
    lj_parameter_arrays,
)

import math
from collections import defaultdict

import numpy as np


#: CHARMM36 per-element LJ parameters: ``rmin_half`` (Rmin/2, Angstrom) and
#: ``epsilon`` (kcal/mol, negative as in the CHARMM parameter files).










def build_graph(bonds):
    g = defaultdict(set)
    for a, b in bonds:
        g[a].add(b)
        g[b].add(a)
    return g


def build_angles(graph):
    angles = set()
    for b, neighbors in graph.items():
        nlist = sorted(neighbors)
        for i in range(len(nlist)):
            for j in range(i + 1, len(nlist)):
                angles.add((nlist[i], b, nlist[j]))
    return sorted(angles)


def build_dihedrals(graph):
    dihedrals = set()
    for b, nb in graph.items():
        for c in nb:
            if b >= c:
                continue
            for a in nb:
                if a == c:
                    continue
                for d in graph[c]:
                    if d == b or d == a:
                        continue
                    t = (a, b, c, d)
                    rev = (d, c, b, a)
                    if t not in dihedrals and rev not in dihedrals:
                        dihedrals.add(t)
    return sorted(dihedrals)


def find_cycles(graph, max_len=7):
    visited = set()
    cycles = set()

    def dfs(start, current, path):
        if len(path) > max_len:
            return
        for nb in graph[current]:
            if nb == start and len(path) >= 3:
                cycles.add(tuple(sorted(path)))
                continue
            if nb in visited:
                continue
            visited.add(nb)
            dfs(start, nb, path + [nb])
            visited.discard(nb)

    for node in graph:
        visited.add(node)
        dfs(node, node, [node])
        visited.discard(node)
    return cycles


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


def torsion_cosine(torsion_type):
    """The ``IMP.core.Cosine`` for a torsion type stored in the CHARMM convention.

    cgdye's ``torsion_types`` (``k``, ``periodicity`` n, ``phase_rad`` δ) mean
    the CHARMM/AMBER form ``V = k (1 + cos(n φ − δ))``: ``T_PI`` (n = 2,
    δ = π) is minimal at the planar 0°/180°, ``T_LINK`` (n = 3, δ = 0) at
    the staggered ±60°/180°. ``IMP.core.Cosine(k, n, δ')`` scores
    ``k (1 − cos(n φ − δ'))`` — the *opposite* sign — so δ' = δ + π. Passing
    δ straight through (as the code did) put every conjugated torsion's
    minimum at 90° and every linker torsion's at the eclipsed 0°/±120°.
    """
    import IMP.core
    return IMP.core.Cosine(float(torsion_type["k"]), int(torsion_type["periodicity"]),
                           float(torsion_type["phase_rad"]) + math.pi)

