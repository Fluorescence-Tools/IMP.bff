"""Dye topology builder utilities and the one CHARMM36 Lennard-Jones source.

``CHARMM36_LJ`` / ``lj_cross`` / ``lj_energy`` are the single LJ table and
kernel of cgdye. Every scorer (``sampling.scoring``, ``sampling.mean_field``,
``rotamer.scoring``) derives its parameters and energies from here so the two
tables that used to live side by side (``rmin_half``/``epsilon`` here,
``p_Rmin2``/``eps`` in ``rotamer.scoring``) cannot drift apart again.
"""

import math
from collections import defaultdict

import numpy as np


#: CHARMM36 per-element LJ parameters: ``rmin_half`` (Rmin/2, Angstrom) and
#: ``epsilon`` (kcal/mol, negative as in the CHARMM parameter files).
CHARMM36_LJ = {
    "C": {"rmin_half": 2.02446316, "epsilon": -0.06394724},
    "N": {"rmin_half": 1.89285714, "epsilon": -0.15428571},
    "O": {"rmin_half": 1.693, "epsilon": -0.12642017},
    "S": {"rmin_half": 2.1, "epsilon": -0.47},
    "H": {"rmin_half": 0.98357778, "epsilon": -0.03466645},
}


def lj_params(element):
    return CHARMM36_LJ.get(element, CHARMM36_LJ["C"])


def lj_cross(elem_i, elem_j):
    """Lorentz-Berthelot cross parameters ``(rmin, eps)`` for two elements.

    ``rmin = rmin_half_i + rmin_half_j`` (Angstrom), ``eps = sqrt(eps_i eps_j)``
    (kcal/mol, positive well depth).
    """
    pi = lj_params(elem_i)
    pj = lj_params(elem_j)
    rmin = pi["rmin_half"] + pj["rmin_half"]
    eps = (pi["epsilon"] * pj["epsilon"]) ** 0.5
    return rmin, eps


def lj_parameter_arrays(elements):
    """``(rmin_half, epsilon)`` arrays for a sequence of element symbols.

    Unknown elements fall back to carbon, as ``lj_params`` does.
    """
    params = [lj_params(e) for e in elements]
    rmin_half = np.array([p["rmin_half"] for p in params], dtype=np.float64)
    epsilon = np.array([p["epsilon"] for p in params], dtype=np.float64)
    return rmin_half, epsilon


def lj_energy(r, rmin, eps, *, repulsive_only=False, cutoff=None, r_floor=0.01):
    """12-6 Lennard-Jones energy ``eps * ((rmin/r)^12 - 2 (rmin/r)^6)``.

    Vectorised over numpy arrays (``r``, ``rmin``, ``eps`` broadcast against
    each other); scalars work too. ``r`` is clamped to ``r_floor`` to keep the
    energy finite. ``repulsive_only=True`` zeroes the attractive tail
    (``r >= rmin``), ``cutoff`` zeroes pairs beyond that distance. With
    positive ``eps`` (the ``lj_cross`` convention) the well depth at ``rmin``
    is ``-eps``.
    """
    r_arr = np.asarray(r, dtype=np.float64)
    r_safe = np.maximum(r_arr, r_floor)
    ratio6 = np.power(np.asarray(rmin, dtype=np.float64) / r_safe, 6)
    energy = np.asarray(eps, dtype=np.float64) * (ratio6 * ratio6 - 2.0 * ratio6)
    if repulsive_only:
        energy = np.where(r_arr < rmin, energy, 0.0)
    if cutoff is not None:
        energy = np.where(r_arr < cutoff, energy, 0.0)
    if np.ndim(energy) == 0:
        return float(energy)
    return energy


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
