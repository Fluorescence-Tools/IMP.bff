#!/usr/bin/env python
"""IMP MD runner from compact FF mmCIF with md23-style hybrid sampling."""

import json
import math
import os
import random
from collections import defaultdict
from pathlib import Path


import IMP
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.core
import IMP.pmi.dof
import IMP.pmi.samplers
import IMP.pmi.tools
import IMP.rmf
import RMF

from IMP.bff.io.cif import read_dye_forcefield_cif, write_dye_forcefield_cif
from ..system import fixed_components, mobile_components
from IMP.bff.scoring import torsion_cosine
from IMP.bff.tools import import_click

click = import_click()  # optional: only the CLI entry point needs it


def _read_mol2_quiet(path, model):
    old = IMP.get_log_level()
    try:
        IMP.set_log_level(IMP.SILENT)
        return IMP.atom.read_mol2(str(path), model)
    finally:
        IMP.set_log_level(old)


def _read_mol2_atom_names(path: str) -> dict:
    """Return {serial: atom_name} from @<TRIPOS>ATOM section.

    IMP's read_mol2 maps TRIPOS types (C.3 → C3) losing the actual column-2
    atom name (e.g. C12, N1).  Read them directly from the file by serial.
    """
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


def _atom_name(p):
    return (
        IMP.atom.Atom(p)
        .get_atom_type()
        .get_string()
        .replace("HET: ", "")
        .replace("HET:", "")
        .strip()
    )


def _validate_system(system):
    if not system.get("components"):
        raise ValueError("system requires components")
    if not system.get("sites"):
        raise ValueError("system requires sites")

    sids = [s["id"] for s in system["sites"]]
    if len(sids) != len(set(sids)):
        raise ValueError("duplicate site ids")
    sidset = set(sids)

    for s in system["sites"]:
        if s["component"] not in system["components"]:
            raise ValueError(f"site {s['id']} unknown component {s['component']}")
        if s.get("site_serial") is None:
            raise ValueError(f"site {s['id']} missing site_serial")

    for a, b, *_ in system.get("bonds", []):
        if a not in sidset or b not in sidset:
            raise ValueError("bond references unknown site")
    for a, b, c, *_ in system.get("angles", []):
        if a not in sidset or b not in sidset or c not in sidset:
            raise ValueError("angle references unknown site")
    for block in ("dihedrals", "impropers"):
        for a, b, c, d, *_ in system.get(block, []):
            if a not in sidset or b not in sidset or c not in sidset or d not in sidset:
                raise ValueError(f"{block} references unknown site")


def _group_ids(system, group_name):
    g = system.get("groups", {})
    if group_name not in g:
        raise ValueError(f"unknown group {group_name}")
    return g[group_name]


def _component_for_group(system, group_name):
    ids = set(_group_ids(system, group_name))
    comps = {s["component"] for s in system.get("sites", []) if s.get("id") in ids}
    if not comps:
        raise ValueError(f"group {group_name} has no sites")
    if len(comps) != 1:
        raise ValueError(
            f"group {group_name} spans multiple components: {sorted(comps)}"
        )
    return next(iter(comps))


def _infer_fixed_component(system):
    fixed = fixed_components(system)
    if fixed:
        return fixed[0]
    raise ValueError("Could not infer fixed component; define role: fixed in system")


def _infer_mobile_group(system):
    mobile = mobile_components(system)
    if not mobile:
        raise ValueError("No component with role 'mobile' found in system")
    comp = mobile[0]
    # Check for '{comp}_all' group
    g = system.get("groups", {})
    target = f"{comp}_all"
    if target in g:
        return target
    # Fallback: check if the component name itself is a group
    if comp in g:
        return comp
    raise ValueError(f"Could not infer mobile group for {comp}; define '{comp}_all' in system")


def _load_hierarchies(model, system, system_cif, only_components=None):
    base = os.path.dirname(os.path.abspath(system_cif))
    root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model, "root"))

    component_hiers = {}
    component_atom_serial = {}
    component_atom_name_serial = {}
    component_atom_name_particles = {}

    for comp, spec in system["components"].items():
        if only_components is not None and comp not in only_components:
            continue
        # Support both new mol2 key and legacy pdb key
        struct_path = spec.get("mol2") or spec.get("pdb")
        if struct_path is None:
            raise ValueError(f"Component {comp} has no mol2 or pdb path in system CIF")
        if not os.path.isabs(struct_path):
            struct_path = os.path.normpath(os.path.join(base, struct_path))
        if struct_path.endswith(".mol2"):
            hier = _read_mol2_quiet(struct_path, model)
        else:
            # Legacy PDB fallback
            old = IMP.get_log_level()
            try:
                IMP.set_log_level(IMP.SILENT)
                hier = IMP.atom.read_pdb(struct_path, model, IMP.atom.AllPDBSelector())
            finally:
                IMP.set_log_level(old)
        hier.set_name(comp)
        root.add_child(hier)
        component_hiers[comp] = hier

        # For MOL2 files, read true atom names from the file (IMP maps TRIPOS
        # types C.3->C3 etc., losing the actual name like C12, N1).
        mol2_names = (
            _read_mol2_atom_names(struct_path) if struct_path.endswith(".mol2") else {}
        )

        atom_map = {}
        atom_name_map = {}
        atom_name_particles = defaultdict(list)
        for a in IMP.atom.get_by_type(hier, IMP.atom.ATOM_TYPE):
            idx = IMP.atom.Atom(a).get_input_index()
            atom_map[idx] = a
            # Use the real MOL2 column-2 name; fall back to IMP atom type
            anm = mol2_names.get(idx) or _atom_name(a)
            atom_name_map[idx] = anm
            atom_name_particles[anm].append(a)
        component_atom_serial[comp] = atom_map
        component_atom_name_serial[comp] = atom_name_map
        component_atom_name_particles[comp] = atom_name_particles

    site_particles = {}
    site_atom_names = {}
    for s in system["sites"]:
        sid = s["id"]
        comp = s["component"]
        if only_components is not None and comp not in only_components:
            continue
        serial = s.get("site_serial")
        atom_name = s.get("atom_name")
        amap = component_atom_serial.get(comp, {})
        p = None
        if serial is not None:
            serial = int(serial)
            if serial not in amap:
                raise ValueError(
                    f"{sid}: serial {serial} missing from component {comp} pdb"
                )
            p = amap[serial]
        else:
            by_name = component_atom_name_particles.get(comp, {})
            candidates = by_name.get(atom_name, []) if atom_name else []
            if len(candidates) == 1:
                p = candidates[0]
            elif len(candidates) > 1:
                raise ValueError(
                    f"{sid}: atom name {atom_name} is ambiguous in component {comp}; provide site_serial"
                )
            else:
                raise ValueError(
                    f"{sid}: atom name {atom_name} missing from component {comp} pdb"
                )
        if hasattr(p, "get_particle"):
            p = p.get_particle()

        xyzr = IMP.core.XYZR(p)
        xyzr.set_radius(float(s.get("radius", xyzr.get_radius())))
        if not IMP.atom.Mass.get_is_setup(p):
            IMP.atom.Mass.setup_particle(p, float(s.get("mass", 12.0)))
        else:
            IMP.atom.Mass(p).set_mass(float(s.get("mass", IMP.atom.Mass(p).get_mass())))

        if not IMP.atom.Bonded.get_is_setup(p):
            IMP.atom.Bonded.setup_particle(p)

        site_particles[sid] = p
        site_atom_names[sid] = component_atom_name_serial.get(comp, {}).get(serial, "")

    return root, site_particles, component_hiers, site_atom_names


def _dist(p1, p2):
    c1 = IMP.core.XYZ(p1).get_coordinates()
    c2 = IMP.core.XYZ(p2).get_coordinates()
    dx = c1[0] - c2[0]
    dy = c1[1] - c2[1]
    dz = c1[2] - c2[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _angle(p1, p2, p3):
    a = IMP.core.XYZ(p1).get_coordinates()
    b = IMP.core.XYZ(p2).get_coordinates()
    c = IMP.core.XYZ(p3).get_coordinates()
    bax = a[0] - b[0]
    bay = a[1] - b[1]
    baz = a[2] - b[2]
    bcx = c[0] - b[0]
    bcy = c[1] - b[1]
    bcz = c[2] - b[2]
    nba = math.sqrt(bax * bax + bay * bay + baz * baz)
    nbc = math.sqrt(bcx * bcx + bcy * bcy + bcz * bcz)
    if nba == 0 or nbc == 0:
        return 1.910633
    dot = bax * bcx + bay * bcy + baz * bcz
    cosang = max(-1.0, min(1.0, dot / (nba * nbc)))
    return math.acos(cosang)


def _build_restraints(model, system, site_particles):
    restraints = []
    softsphere_restraints = []

    bt = system.get("bond_types", {})
    at = system.get("angle_types", {})
    tt = system.get("torsion_types", {})
    it = system.get("improper_types", {})

    sp = site_particles  # shorthand; only contains loaded components

    for a, b, length, tid in system.get("bonds", []):
        if a not in sp or b not in sp:
            continue
        p1, p2 = sp[a], sp[b]
        k = float(bt[tid]["k"])
        d0 = _dist(p1, p2) if length is None else float(length)
        restraints.append(
            IMP.core.DistanceRestraint(model, IMP.core.Harmonic(d0, k), p1, p2)
        )

    for a, b, c, theta, tid in system.get("angles", []):
        if a not in sp or b not in sp or c not in sp:
            continue
        p1, p2, p3 = sp[a], sp[b], sp[c]
        k = float(at[tid]["k"])
        t0 = _angle(p1, p2, p3) if theta is None else float(theta)
        restraints.append(
            IMP.core.AngleRestraint(model, IMP.core.Harmonic(t0, k), p1, p2, p3)
        )

    for a, b, c, d, tid in system.get("dihedrals", []):
        if a not in sp or b not in sp or c not in sp or d not in sp:
            continue
        p1, p2, p3, p4 = sp[a], sp[b], sp[c], sp[d]
        # CHARMM-convention type -> IMP.core.Cosine (sign flip, see topology.dye.torsion_cosine)
        restraints.append(IMP.core.DihedralRestraint(model, torsion_cosine(tt[tid]), p1, p2, p3, p4))

    for a, b, c, d, tid in system.get("impropers", []):
        if a not in sp or b not in sp or c not in sp or d not in sp:
            continue
        p1, p2, p3, p4 = sp[a], sp[b], sp[c], sp[d]
        t = it[tid]
        theta0 = IMP.core.get_dihedral(
            IMP.core.XYZ(p1), IMP.core.XYZ(p2), IMP.core.XYZ(p3), IMP.core.XYZ(p4)
        )
        fun = IMP.core.Harmonic(theta0, float(t["k"]))
        restraints.append(IMP.core.DihedralRestraint(model, fun, p1, p2, p3, p4))

    nb = system.get("nonbonded", {})
    if nb.get("enabled", True):
        sids = sorted(site_particles.keys())
        idx_pairs = []
        excl = _derive_exclusions(system)
        for i in range(len(sids)):
            for j in range(i + 1, len(sids)):
                a, b = sids[i], sids[j]
                if (a, b) in excl:
                    continue
                idx_pairs.append(
                    (site_particles[a].get_index(), site_particles[b].get_index())
                )
        lpc = IMP.container.ListPairContainer(model, idx_pairs)
        pair_score = IMP.core.SoftSpherePairScore(float(nb.get("k", 5.0)))
        rnb = IMP.container.PairsRestraint(pair_score, lpc)
        restraints.append(rnb)
        softsphere_restraints.append(rnb)

    return restraints, softsphere_restraints


def _derive_exclusions(system):
    excl = set()

    def add_pair(a, b):
        excl.add((a, b) if a <= b else (b, a))

    explicit = system.get("exclusions", [])
    for a, b in explicit:
        add_pair(a, b)

    if explicit:
        return excl

    for a, b, *_ in system.get("bonds", []):
        add_pair(a, b)
    for a, _, c, *_ in system.get("angles", []):
        add_pair(a, c)
    for a, _, _, d, *_ in system.get("dihedrals", []):
        add_pair(a, d)

    return excl


def _build_site_graph(system):
    g = defaultdict(set)
    for a, b, *_ in system.get("bonds", []):
        g[a].add(b)
        g[b].add(a)
    return g


def _find_cycles(graph, max_len=8):
    nodes = sorted(graph.keys())
    cycles = set()

    def canonical(cyc):
        cyc = list(cyc)
        n = len(cyc)
        rots = [tuple(cyc[i:] + cyc[:i]) for i in range(n)]
        rc = list(reversed(cyc))
        rots.extend(tuple(rc[i:] + rc[:i]) for i in range(n))
        return min(rots)

    def dfs(start, cur, visited, path):
        for nbr in graph.get(cur, set()):
            if nbr == start and len(path) >= 3:
                cycles.add(canonical(path))
                continue
            if nbr in visited:
                continue
            if nbr < start:
                continue
            if len(path) >= max_len:
                continue
            dfs(start, nbr, visited | {nbr}, path + [nbr])

    for start in nodes:
        dfs(start, start, {start}, [start])
    return cycles


def _mobile_ring_site_ids(
    system, site_atom_names, mobile_component, include_attached_h=True
):
    full_graph = _build_site_graph(system)
    mobile_ids = {
        s["id"]
        for s in system.get("sites", [])
        if s.get("component") == mobile_component
    }
    mobile_graph = {
        sid: {nb for nb in full_graph.get(sid, set()) if nb in mobile_ids}
        for sid in mobile_ids
    }

    cycles = _find_cycles(mobile_graph, max_len=8)
    ring_ids = set()
    for cyc in cycles:
        if len(cyc) in (5, 6, 7):
            ring_ids.update(cyc)

    if include_attached_h:
        extra = set()
        for sid in ring_ids:
            for nb in mobile_graph.get(sid, set()):
                if site_atom_names.get(nb, "").startswith("H"):
                    extra.add(nb)
        ring_ids.update(extra)

    return ring_ids


def _topo_distance_leq(graph, src, dst, max_depth):
    if src == dst:
        return True
    seen = {src}
    frontier = {src}
    depth = 0
    while frontier and depth < max_depth:
        depth += 1
        nxt = set()
        for node in frontier:
            for nb in graph.get(node, set()):
                if nb == dst:
                    return True
                if nb in seen:
                    continue
                seen.add(nb)
                nxt.add(nb)
        frontier = nxt
    return False


def _fixed_flex_ids(system, fixed_flex_mode, fixed_name):
    """Return the set of fixed component site IDs that are allowed to move.

    ``fixed_flex_mode`` may be:
    - ``"static"``   – fixed component is fully rigid (default)
    - ``"flex"``     – atoms listed in the ``{name}_flex`` group of the
                        system CIF are released; all other fixed atoms stay fixed
    """
    if fixed_flex_mode != "flex":
        return set()
    groups = system.get("groups", {})
    flex_group = f"{fixed_name}_flex"
    return set(groups.get(flex_group, []))


def _build_go_restraints(
    model,
    system,
    site_particles,
    site_atom_names,
    mobile_component,
    fixed_component,
    fixed_flex_mode,
    go_mobile_k,
    go_fixed_k,
    cutoff,
):
    graph = _build_site_graph(system)
    s2c = {s["id"]: s["component"] for s in system.get("sites", [])}

    def is_heavy(sid):
        nm = site_atom_names.get(sid, "")
        return not nm.startswith("H")

    def build_contacts(ids):
        ids = sorted(ids)
        contacts = []
        for i, s1 in enumerate(ids):
            p1 = site_particles[s1]
            for s2 in ids[i + 1 :]:
                if _topo_distance_leq(graph, s1, s2, 2):
                    continue
                p2 = site_particles[s2]
                d0 = _dist(p1, p2)
                if d0 <= cutoff:
                    contacts.append((s1, s2, d0))
        return contacts

    mobile_heavy = [
        sid for sid, comp in s2c.items() if comp == mobile_component and is_heavy(sid)
    ]
    fixed_heavy = [
        sid for sid, comp in s2c.items() if comp == fixed_component and is_heavy(sid)
    ]
    fixed_movable = _fixed_flex_ids(system, fixed_flex_mode, fixed_component)

    restraints = []
    for s1, s2, d0 in build_contacts(mobile_heavy):
        r = IMP.core.DistanceRestraint(
            model,
            IMP.core.Harmonic(max(d0, 1.0), go_mobile_k),
            site_particles[s1],
            site_particles[s2],
        )
        r.set_name(f"go_mobile_{s1}_{s2}")
        restraints.append(r)

    for s1, s2, d0 in build_contacts(fixed_heavy):
        if s1 not in fixed_movable and s2 not in fixed_movable:
            continue
        r = IMP.core.DistanceRestraint(
            model,
            IMP.core.Harmonic(max(d0, 1.0), go_fixed_k),
            site_particles[s1],
            site_particles[s2],
        )
        r.set_name(f"go_fixed_{s1}_{s2}")
        restraints.append(r)

    return restraints


def _build_component_atom_maps(component_hiers, system=None, system_cif=None):
    """Build {comp: {by_name, by_serial}} maps.

    For MOL2-backed components the true atom names are read directly from the
    file (IMP's read_mol2 maps TRIPOS types, losing the real column-2 name).
    """
    base = os.path.dirname(os.path.abspath(system_cif)) if system_cif else None
    out = {}
    for comp, hier in component_hiers.items():
        # Resolve real atom names from MOL2 if available
        mol2_names = {}
        if system and base:
            spec = system.get("components", {}).get(comp, {})
            struct_path = spec.get("mol2") or spec.get("pdb") or ""
            if struct_path and not os.path.isabs(struct_path):
                struct_path = os.path.normpath(os.path.join(base, struct_path))
            if struct_path.endswith(".mol2") and os.path.exists(struct_path):
                mol2_names = _read_mol2_atom_names(struct_path)

        by_name = defaultdict(list)
        by_serial = {}
        atoms = IMP.atom.get_by_type(hier, IMP.atom.ATOM_TYPE)
        for a in atoms:
            idx = IMP.atom.Atom(a).get_input_index()
            nm = mol2_names.get(idx) or _atom_name(a)
            by_name[nm].append(a)
            by_serial[idx] = a
        out[comp] = {"by_name": by_name, "by_serial": by_serial}
    return out


def _resolve_position_particle(
    pos_name,
    pos_def,
    guest_component,
    atom_maps,
    host_component=None,
):
    comp = pos_def.get("component_id")
    if comp not in atom_maps:
        comp = None

    if comp is None:
        serial = pos_def.get("atom_serial")
        if serial is not None:
            candidates = [
                c
                for c, cmap in atom_maps.items()
                if int(serial) in cmap.get("by_serial", {})
            ]
            if len(candidates) == 1:
                comp = candidates[0]

    if comp is None:
        anm = pos_def.get("atom_name")
        if anm:
            candidates = [
                c for c, cmap in atom_maps.items() if cmap.get("by_name", {}).get(anm)
            ]
            if len(candidates) == 1:
                comp = candidates[0]

    if comp is None:
        comp = guest_component

    cmap = atom_maps.get(comp, {})
    by_serial = cmap.get("by_serial", {})
    by_name = cmap.get("by_name", {})

    serial = pos_def.get("atom_serial")
    if serial is not None and int(serial) in by_serial:
        return by_serial[int(serial)]

    anm = pos_def.get("atom_name")
    if anm and anm in by_name and by_name[anm]:
        return by_name[anm][0]

    return None




def _center(ps):
    c = IMP.algebra.Vector3D(0.0, 0.0, 0.0)
    for p in ps:
        c += IMP.core.XYZ(p).get_coordinates()
    return c / max(1, len(ps))


def _rand_unit(rng):
    z = rng.uniform(-1.0, 1.0)
    t = rng.uniform(0.0, 2.0 * math.pi)
    rxy = math.sqrt(max(0.0, 1.0 - z * z))
    return IMP.algebra.Vector3D(rxy * math.cos(t), rxy * math.sin(t), z)


def _score_based_initial_placement(
    placement_sf,
    host_ps,
    guest_ps,
    distance_a,
    n_trials,
    seed,
):
    rng = random.Random(seed)

    host_center = _center(host_ps)
    guest_template = [IMP.core.XYZ(p).get_coordinates() for p in guest_ps]
    guest_center = _center(guest_ps)

    best_coords = list(guest_template)
    best_score = None

    for _ in range(max(1, n_trials)):
        axis = _rand_unit(rng)
        angle = rng.uniform(0.0, 2.0 * math.pi)
        rot = IMP.algebra.get_rotation_about_axis(axis, angle)
        direction = _rand_unit(rng)
        target_center = host_center + direction * distance_a

        trial_coords = []
        for c in guest_template:
            rel = c - guest_center
            rr = rot.get_rotated(rel)
            trial_coords.append(target_center + rr)

        for p, c in zip(guest_ps, trial_coords):
            IMP.core.XYZ(p).set_coordinates(c)

        sc = placement_sf.evaluate(False)
        if best_score is None or sc < best_score:
            best_score = sc
            best_coords = list(trial_coords)

    for p, c in zip(guest_ps, best_coords):
        IMP.core.XYZ(p).set_coordinates(c)

    return best_score


def _configure_movable(
    system,
    site_particles,
    site_atom_names,
    fixed_flex_mode,
    md_fixed_site_ids,
    explicit_movable_ids=None,
):
    fixed = set()
    for gid in system.get("fixed_groups", []):
        fixed.update(_group_ids(system, gid))

    # In flex mode, release atoms listed in any {component}_flex group that
    # belongs to a fixed component.  Those atoms were included in the broad
    # fixed_groups entry but should be free to move during MD.
    if fixed_flex_mode == "flex":
        for gid, ids in system.get("groups", {}).items():
            if gid.endswith("_flex"):
                fixed.difference_update(ids)

    if explicit_movable_ids is not None:
        explicit = set(explicit_movable_ids)
        md_fixed = set(md_fixed_site_ids)
        movable = []
        for sid, p in site_particles.items():
            xyz = IMP.core.XYZ(p)
            opt = sid in explicit and sid not in md_fixed
            xyz.set_coordinates_are_optimized(opt)
            if opt:
                movable.append(p)
        return movable

    fixed.update(md_fixed_site_ids)

    movable = []
    for sid, p in site_particles.items():
        xyz = IMP.core.XYZ(p)
        opt = sid not in fixed
        xyz.set_coordinates_are_optimized(opt)
        if opt:
            movable.append(p)
    return movable


def _absolutize_component_pdbs(system, system_cif):
    sys_out = json.loads(json.dumps(system))
    base_dir = Path(system_cif).resolve().parent
    for comp, spec in sys_out.get("components", {}).items():
        for key in ("mol2", "pdb"):
            p = spec.get(key)
            if p and not os.path.isabs(p):
                sys_out["components"][comp][key] = str((base_dir / p).resolve())
    return sys_out


def _write_system_copy(system, system_cif, out_dir, output_root):
    sys_out = _absolutize_component_pdbs(system, system_cif)

    out_path = Path(out_dir)
    copy_path = out_path / "system.cif"
    write_dye_forcefield_cif(str(copy_path), sys_out)

    out_root = Path(output_root)
    trajs_anchor = None
    for p in [out_root] + list(out_root.parents):
        if p.name == "trajs":
            trajs_anchor = p
            break
    systems_dir = (
        trajs_anchor.parent / "systems"
        if trajs_anchor is not None
        else out_root / "systems"
    )
    systems_dir.mkdir(parents=True, exist_ok=True)
    root_copy_path = systems_dir / f"{system.get('name', 'ff_system')}.system.cif"
    write_dye_forcefield_cif(str(root_copy_path), sys_out)


def _run_simple_md(
    model,
    root,
    sf,
    out_dir,
    movable,
    n_steps,
    write_every,
    temperature_k,
    friction_ps,
    timestep_fs,
    log_every_frames,
):
    md = IMP.atom.MolecularDynamics(model)
    md.set_scoring_function(sf)
    md.set_maximum_time_step(timestep_fs)
    md.assign_velocities(temperature_k)

    if movable:
        thermo = IMP.atom.LangevinThermostatOptimizerState(
            model, movable, temperature_k, friction_ps
        )
        thermo.set_period(1)
        md.add_optimizer_state(thermo)

    init_fh = RMF.create_rmf_file(os.path.join(out_dir, "initial.0.rmf3"))
    IMP.rmf.add_hierarchies(init_fh, [root])
    IMP.rmf.save_frame(init_fh, "init")
    del init_fh

    rmf_fh = RMF.create_rmf_file(os.path.join(out_dir, "rmfs", "0.rmf3"))
    IMP.rmf.add_hierarchies(rmf_fh, [root])

    n_frames = max(1, n_steps // max(1, write_every))
    

    with open(os.path.join(out_dir, "stat.0.out"), "w") as stat:


        stat.write("frame\tscore\tkinetic_energy\n")
        IMP.rmf.save_frame(rmf_fh, "init")
        init_score = sf.evaluate(False)
        stat.write(f"init\t{init_score:.6f}\t{md.get_kinetic_energy():.6f}\n")

        for frame in range(n_frames):
            md.optimize(write_every)
            score = sf.evaluate(False)
            ke = md.get_kinetic_energy()
            IMP.rmf.save_frame(rmf_fh, str(frame))
            stat.write(f"{frame}\t{score:.6f}\t{ke:.6f}\n")

            if frame % max(1, int(log_every_frames)) == 0:
                print(f"  frame {frame:5d}/{n_frames} score={score:.2f} KE={ke:.2f}")


def _rb_mc_step(
    rb_particles_per_group, mc_sf, mc_temp, mc_steps, max_trans, max_rot_rad
):
    """Manual rigid-body MC: translate+rotate all atoms in each group as a body.

    Returns (n_accepted, n_proposed).
    IMP's MolecularDynamics cannot integrate NonRigidMember particles, so we
    avoid IMP RigidBody membership entirely and manage the rigid move ourselves.
    """
    import random as _rnd

    kT = mc_temp  # kcal/mol, same units as score
    t_scale = max_trans / 3.0
    r_scale = max_rot_rad / 3.0
    n_acc = 0
    n_prop = 0
    for _ in range(mc_steps):
        for rb_ps in rb_particles_per_group:
            if not rb_ps:
                continue
            # Save current coords
            old_coords = [IMP.core.XYZ(p).get_coordinates() for p in rb_ps]
            # Compute COM
            n = len(old_coords)
            com = IMP.algebra.Vector3D(
                sum(c[0] for c in old_coords) / n,
                sum(c[1] for c in old_coords) / n,
                sum(c[2] for c in old_coords) / n,
            )
            # Random translation ~ N(0, t_scale) per axis
            trans = IMP.algebra.Vector3D(
                _rnd.gauss(0, t_scale),
                _rnd.gauss(0, t_scale),
                _rnd.gauss(0, t_scale),
            )
            # Random rotation around random axis by angle ~ N(0, r_scale)
            angle = _rnd.gauss(0, r_scale)
            axis = IMP.algebra.get_random_vector_on_unit_sphere()
            half = angle / 2.0
            s = math.sin(half)
            rot = IMP.algebra.Rotation3D(
                IMP.algebra.Vector4D(
                    math.cos(half),
                    axis[0] * s,
                    axis[1] * s,
                    axis[2] * s,
                )
            )
            score_before = mc_sf.evaluate(False)
            # Apply: rotate around COM then translate
            for p, old_c in zip(rb_ps, old_coords):
                local = old_c - com
                new_c = rot.get_rotated(local) + com + trans
                IMP.core.XYZ(p).set_coordinates(new_c)
            score_after = mc_sf.evaluate(False)
            delta = score_after - score_before
            n_prop += 1
            if delta <= 0 or _rnd.random() < math.exp(-delta / kT):
                n_acc += 1
            else:
                # Reject: restore
                for p, old_c in zip(rb_ps, old_coords):
                    IMP.core.XYZ(p).set_coordinates(old_c)
    return n_acc, n_prop


def _run_alternating_rb_mc_md(
    model,
    mc_sf,
    full_sf,
    root,
    out_dir,
    site_particles,
    movable,
    rb_site_groups,
    md_steps,
    md_steps_per_block,
    mc_steps_per_block,
    temperature_k,
    friction_ps,
    timestep_fs,
    rb_max_trans,
    rb_max_rot_deg,
    mc_temp,
    log_every_frames,
    fixed_flex_ids,
    mobile_ring_particles=None,
    mobile_linker_particles=None,
):
    # Manual RB-MC: no IMP RigidBody membership so MD can freely integrate all
    # movable particles.  Each MC step applies a random rigid translation+rotation
    # to all atoms in a site-group and accepts/rejects via Metropolis.
    rb_site_groups = [list(g) for g in rb_site_groups if g]
    if not rb_site_groups:
        raise RuntimeError("No rb_groups in system mmCIF")

    rb_particles_per_group = [
        [site_particles[sid] for sid in g] for g in rb_site_groups
    ]

    md = IMP.atom.MolecularDynamics(model)
    md.set_scoring_function(full_sf)
    md.set_maximum_time_step(timestep_fs)
    md.assign_velocities(temperature_k)

    if movable:
        thermo = IMP.atom.LangevinThermostatOptimizerState(
            model, movable, temperature_k, friction_ps
        )
        thermo.set_period(1)
        md.add_optimizer_state(thermo)


    init_fh = RMF.create_rmf_file(os.path.join(out_dir, "initial.0.rmf3"))
    IMP.rmf.add_hierarchies(init_fh, [root])
    IMP.rmf.save_frame(init_fh, "init")
    del init_fh

    rmf_fh = RMF.create_rmf_file(os.path.join(out_dir, "rmfs", "0.rmf3"))
    IMP.rmf.add_hierarchies(rmf_fh, [root])

    n_frames = max(1, md_steps // max(1, md_steps_per_block))

    rb_particles = list({site_particles[sid] for g in rb_site_groups for sid in g})
    rb_com0 = (
        _center(rb_particles) if rb_particles else IMP.algebra.Vector3D(0.0, 0.0, 0.0)
    )

    fixed_flex_particles = [
        site_particles[sid] for sid in fixed_flex_ids if sid in site_particles
    ]
    fixed_flex_coords0 = [
        IMP.core.XYZ(p).get_coordinates() for p in fixed_flex_particles
    ]

    # Internal linker displacement: position of each linker atom relative to ring COM
    _ring_ps = list(mobile_ring_particles) if mobile_ring_particles else []
    _linker_ps = list(mobile_linker_particles) if mobile_linker_particles else []

    def _ring_com():
        if not _ring_ps:
            return None
        return _center(_ring_ps)

    _ring_com0 = _ring_com()
    _linker_rel0 = []
    if _ring_com0 is not None and _linker_ps:
        for p in _linker_ps:
            _linker_rel0.append(IMP.core.XYZ(p).get_coordinates() - _ring_com0)

    acc_total = 0
    prop_total = 0

    print(
        f"    Alternating RB-MC/MD: frames={n_frames}, mc_steps={mc_steps_per_block}, md_steps={md_steps_per_block}, rb_groups={len(rb_site_groups)}",
        flush=True,
    )

    with open(os.path.join(out_dir, "stat.0.out"), "w") as stat:


        stat.write("frame\tscore\tkinetic_energy\n")
        IMP.rmf.save_frame(rmf_fh, "init")
        stat.write(
            f"init\t{full_sf.evaluate(False):.6f}\t{md.get_kinetic_energy():.6f}\n"
        )

        for frame in range(n_frames):
            n_acc, n_prop = _rb_mc_step(
                rb_particles_per_group,
                mc_sf,
                mc_temp,
                mc_steps_per_block,
                rb_max_trans,
                math.radians(rb_max_rot_deg),
            )
            acc_total += n_acc
            prop_total += n_prop
            md.optimize(md_steps_per_block)
            score = full_sf.evaluate(False)
            ke = md.get_kinetic_energy()
            rb_com = _center(rb_particles) if rb_particles else rb_com0
            rb_com_disp = IMP.algebra.get_distance(rb_com, rb_com0)

            fixed_flex_disp = 0.0
            if fixed_flex_particles:
                vals = []
                for p, c0 in zip(fixed_flex_particles, fixed_flex_coords0):
                    c = IMP.core.XYZ(p).get_coordinates()
                    vals.append(IMP.algebra.get_distance(c, c0))
                fixed_flex_disp = sum(vals) / len(vals)

            # Mean internal linker displacement relative to ring COM
            linker_int_disp = 0.0
            if _ring_com0 is not None and _linker_ps and _linker_rel0:
                rc = _ring_com()
                vals = []
                for p, rel0 in zip(_linker_ps, _linker_rel0):
                    rel_now = IMP.core.XYZ(p).get_coordinates() - rc
                    vals.append(IMP.algebra.get_distance(rel_now, rel0))
                linker_int_disp = sum(vals) / len(vals)

            IMP.rmf.save_frame(rmf_fh, str(frame))
            stat.write(f"{frame}\t{score:.6f}\t{ke:.6f}\n")
            if frame % max(1, int(log_every_frames)) == 0:
                acc_txt = f"{n_acc}/{n_prop}"
                print(
                    f"  frame {frame:5d}/{n_frames} score={score:.2f} KE={ke:.2f} "
                    f"rbCOM={rb_com_disp:.3f}A linkerInt={linker_int_disp:.3f}A "
                    f"so3Disp={fixed_flex_disp:.3f}A mcAcc={acc_txt}",
                    flush=True,
                )


def _run(
    system,
    system_cif,
    output_root,
    md_steps,
    write_every_override,
    sampling_mode,
    mc_steps,
    mc_temperature,
    rb_max_trans,
    rb_max_rot_deg,
    mobile_group,
    fixed_flex_mode,
    freeze_mobile_rings,
    init_placement,
    init_distance_a,
    init_trials,
    init_seed,
    com_pull_k,
    go_mobile_k,
    go_fixed_k,
    go_cutoff,
    log_every_frames,
    n_restarts,
):
    model = IMP.Model()
    if mobile_group is None:
        mobile_group = _infer_mobile_group(system)
        print(f"Inferred mobile group: {mobile_group}")
    mobile_component = _component_for_group(system, mobile_group)
    fixed_component = _infer_fixed_component(system)
    root, site_particles, component_hiers, site_atom_names = _load_hierarchies(
        model,
        system,
        system_cif,
        only_components={fixed_component, mobile_component},
    )
    restraints, softsphere_restraints = _build_restraints(model, system, site_particles)

    # Output dir is named after this specific run: {fixed}_{mobile}_imp
    # This makes each mobile's trajectory independently addressable even when
    # the system CIF contains all components.
    run_name = f"{fixed_component}_{mobile_component}"
    out_dir = os.path.join(output_root, f"{run_name}_imp")
    os.makedirs(os.path.join(out_dir, "rmfs"), exist_ok=True)
    go_restraints = _build_go_restraints(
        model,
        system,
        site_particles,
        site_atom_names,
        mobile_component,
        fixed_component,
        fixed_flex_mode,
        go_mobile_k,
        go_fixed_k,
        go_cutoff,
    )

    com_pull_restraint = None
    if com_pull_k > 0:
        fixed_ids = (
            _group_ids(system, f"{fixed_component}_all")
            if f"{fixed_component}_all" in system.get("groups", {})
            else []
        )
        mobile_ids = _group_ids(system, mobile_group)
        fixed_ps = [site_particles[sid] for sid in fixed_ids if sid in site_particles]
        mobile_ps = [site_particles[sid] for sid in mobile_ids if sid in site_particles]
        if fixed_ps and mobile_ps:
            c1 = IMP.atom.CenterOfMass.setup_particle(IMP.Particle(model), fixed_ps)
            c2 = IMP.atom.CenterOfMass.setup_particle(IMP.Particle(model), mobile_ps)
            com_pull_restraint = IMP.core.DistanceRestraint(
                model, IMP.core.Harmonic(0.0, com_pull_k), c1, c2
            )
            com_pull_restraint.set_name("com_pull")

    all_restraints = restraints + go_restraints
    if com_pull_restraint:
        all_restraints.append(com_pull_restraint)

    for r in all_restraints:
        IMP.pmi.tools.add_restraint_to_model(model, r)

    sf = IMP.core.RestraintsScoringFunction(all_restraints)
    mc_restraints = list(softsphere_restraints)
    if not mc_restraints:
        mc_restraints = all_restraints
    mc_sf = IMP.core.RestraintsScoringFunction(mc_restraints)

    if init_placement:
        fixed_ids = (
            _group_ids(system, f"{fixed_component}_all")
            if f"{fixed_component}_all" in system.get("groups", {})
            else []
        )
        mobile_ids_for_init = _group_ids(system, mobile_group)
        fixed_ps = (
            [site_particles[sid] for sid in fixed_ids]
            if fixed_ids
            else [
                p
                for sid, p in site_particles.items()
                if sid not in set(mobile_ids_for_init)
            ]
        )
        mobile_ps = [site_particles[sid] for sid in mobile_ids_for_init]
        seed = init_seed + sum((i + 1) * ord(ch) for i, ch in enumerate(run_name))
        best = _score_based_initial_placement(
            mc_sf,
            fixed_ps,
            mobile_ps,
            init_distance_a,
            init_trials,
            seed,
        )
        print(
            f"Initial score-based placement: distance={init_distance_a:.1f}A trials={init_trials} best_score={best:.2f}",
            flush=True,
        )

    sampling = dict(system.get("sampling", {}))
    if md_steps is not None:
        sampling["n_steps"] = md_steps
    if write_every_override is not None:
        sampling["write_every"] = write_every_override

    n_steps = int(sampling.get("n_steps", 100000))
    write_every = int(sampling.get("write_every", 1000))
    temperature_k = float(sampling.get("temperature_K", 300.0))
    friction_ps = float(sampling.get("friction_ps", 10.0))
    timestep_fs = float(sampling.get("timestep_fs", 0.25))
    minimize_steps = int(sampling.get("minimize_steps", 200))

    mobile_all_ids = _group_ids(system, mobile_group)
    ring_ids = _mobile_ring_site_ids(
        system, site_atom_names, mobile_component, include_attached_h=True
    )
    # The RB always covers the entire mobile component so that every atom is
    # displaced as a rigid body during MC moves.  md_fixed_groups control which
    # atoms are frozen during the MD sub-steps (ring core); flex/linker atoms
    # remain MD-movable but still ride the RB during MC.
    default_rb_site_groups = [list(mobile_all_ids)]

    md_fixed_groups_cfg = system.get("md_fixed_groups", {})
    mobile_md_fixed = sorted(
        sid
        for sids in md_fixed_groups_cfg.values()
        for sid in sids
        if sid in set(mobile_all_ids)
    )
    default_md_fixed_site_ids = mobile_md_fixed
    if not default_md_fixed_site_ids and freeze_mobile_rings:
        default_md_fixed_site_ids = sorted(set(mobile_all_ids) & set(ring_ids))
    known_site_ids = set(site_particles.keys())

    rb_site_groups = [list(g) for g in default_rb_site_groups]
    md_fixed_site_ids = list(default_md_fixed_site_ids)
    movable = _configure_movable(
        system,
        site_particles,
        site_atom_names,
        fixed_flex_mode,
        md_fixed_site_ids,
    )

    rb_site_groups = [g for g in rb_site_groups if g]
    if not rb_site_groups:
        rb_site_groups = [list(mobile_all_ids)]
    rb_site_ids = sorted({sid for g in rb_site_groups for sid in g})

    unknown_md_fixed = sorted(set(md_fixed_site_ids) - known_site_ids)
    if unknown_md_fixed:
        raise ValueError(
            f"mmCIF md_fixed ids not in system: {', '.join(unknown_md_fixed)}"
        )

    movable = [
        p
        for sid, p in site_particles.items()
        if IMP.core.XYZ(p).get_coordinates_are_optimized()
    ]

    fixed_flex_ids = _fixed_flex_ids(system, fixed_flex_mode, fixed_component)
    movable_ids = {
        sid
        for sid, p in site_particles.items()
        if IMP.core.XYZ(p).get_coordinates_are_optimized()
    }
    mobile_linker_movable = [
        sid for sid in movable_ids if sid in mobile_all_ids and sid not in set(ring_ids)
    ]
    print(
        f"DOF: movable={len(movable)} / {len(site_particles)} (freeze_mobile_rings={freeze_mobile_rings}, fixed_flex_mode={fixed_flex_mode})",
        flush=True,
    )
    print(
        f"DOF detail: rb_groups={len(rb_site_groups)} rb_sites={len(rb_site_ids)} md_fixed_sites={len(md_fixed_site_ids)} mobile_linker_movable={len(mobile_linker_movable)} fixed_flex_movable={len(fixed_flex_ids)}",
        flush=True,
    )

    if minimize_steps > 0 and movable:
        cg = IMP.core.ConjugateGradients(model)
        cg.set_scoring_function(sf)
        cg.optimize(minimize_steps)

    _write_system_copy(system, system_cif, out_dir, output_root)

    print(
        f"\n{run_name}: {n_steps} steps, write_every={write_every}, mode={sampling_mode}"
    )
    print(
        f"Restraints: total={len(all_restraints)} softsphere={len(softsphere_restraints)} go={len(go_restraints)}",
        flush=True,
    )

    # Particles for internal displacement tracking
    mobile_ring_ps = [
        site_particles[sid] for sid in md_fixed_site_ids if sid in site_particles
    ]
    mobile_linker_ps = [
        site_particles[sid]
        for sid in mobile_all_ids
        if sid in site_particles
        and sid not in set(md_fixed_site_ids)
        and not site_atom_names.get(sid, "").startswith("H")
    ]

    if sampling_mode == "hybrid_md_mc":
        _run_alternating_rb_mc_md(
            model,
            mc_sf,
            sf,
            root,
            out_dir,
            site_particles,
            movable,
            rb_site_groups,
            n_steps,
            write_every,
            mc_steps,
            mc_temperature,
            rb_max_trans,
            rb_max_rot_deg,
            temperature_k,
            friction_ps,
            timestep_fs,
            log_every_frames,
            fixed_flex_ids,
            mobile_ring_ps,
            mobile_linker_ps,
        )
    elif sampling_mode == "multi_restart":
        _run_multi_restart(
            model,
            mc_sf,
            sf,
            root,
            out_dir,
            site_particles,
            movable,
            rb_site_groups,
            n_restarts,
            n_steps,
            write_every,
            mc_steps,
            mc_temperature,
            rb_max_trans,
            rb_max_rot_deg,
            temperature_k,
            friction_ps,
            timestep_fs,
            log_every_frames,
            fixed_flex_ids,
            mobile_ring_ps,
            mobile_linker_ps,
            system,
            mobile_group,
            init_distance_a,
            init_trials,
            init_seed,
            fixed_component,
        )
    else:
        _run_simple_md(
            model,
            root,
            sf,
            out_dir,
            movable,
            n_steps,
            write_every,
            temperature_k,
            friction_ps,
            timestep_fs,
            log_every_frames,
        )


def _run_multi_restart(
    model,
    mc_sf,
    full_sf,
    root,
    out_dir,
    site_particles,
    movable,
    rb_site_groups,
    n_restarts,
    md_steps,
    md_steps_per_block,
    mc_steps_per_block,
    mc_temp,
    rb_max_trans,
    rb_max_rot_deg,
    temperature_k,
    friction_ps,
    timestep_fs,
    log_every_frames,
    fixed_flex_ids,
    mobile_ring_particles,
    mobile_linker_particles,
    system,
    guest_group,
    init_distance_a,
    init_trials,
    init_seed,
    host_component,
):
    """Run multiple independent restarts of the hybrid MC/MD simulation and collect final states."""
    print(f"Running {n_restarts} independent restarts for multi_restart mode.")
    os.makedirs(out_dir, exist_ok=True)

    # Collective RMF for final frames
    final_rmf_path = os.path.join(out_dir, "rmfs", "0.rmf3")
    os.makedirs(os.path.dirname(final_rmf_path), exist_ok=True)
    final_rh = RMF.create_rmf_file(final_rmf_path)
    IMP.rmf.add_hierarchies(final_rh, [root])

    # Collective stat file
    stat_path = os.path.join(out_dir, "stat.0.out")
    
    # Save initial state to restore for each restart
    all_ps = [model.get_particle(i) for i in model.get_particle_indexes()]
    initial_coordinates = [IMP.core.XYZ(p).get_coordinates() for p in all_ps]
    vk = [IMP.FloatKey(f"velocity_{i}") for i in range(3)]
    initial_velocities = [
        IMP.atom.get_velocity(p) if p.has_attribute(vk[0]) else IMP.algebra.Vector3D(0, 0, 0)
        for p in all_ps
    ]

    converged = False

    with open(stat_path, "w") as fstats:
        fstats.write("frame\tscore\tkinetic_energy\n")

        # Reuse MD and Thermostat objects
        md = IMP.atom.MolecularDynamics(model)
        md.set_scoring_function(full_sf)
        md.set_maximum_time_step(timestep_fs)

        if movable:
            thermo = IMP.atom.LangevinThermostatOptimizerState(
                model, movable, temperature_k, friction_ps
            )
            thermo.set_period(1)
            md.add_optimizer_state(thermo)


        rb_particles_per_group = [
            [site_particles[sid] for sid in g] for g in rb_site_groups
        ]

        for i in range(n_restarts):
            # Restore initial state
            for p, coords, vel in zip(all_ps, initial_coordinates, initial_velocities):
                IMP.core.XYZ(p).set_coordinates(coords)
                if p.has_attribute(vk[0]):
                    IMP.atom.set_velocity(p, vel)

            # Re-run initial placement for each restart
            mobile_ids_for_init = _group_ids(system, guest_group)
            fixed_ids = (
                _group_ids(system, f"{host_component}_all")
                if f"{host_component}_all" in system.get("groups", {})
                else []
            )
            fixed_ps = (
                [site_particles[sid] for sid in fixed_ids]
                if fixed_ids
                else [
                    p
                    for sid, p in site_particles.items()
                    if sid not in set(mobile_ids_for_init)
                ]
            )
            mobile_ps = [site_particles[sid] for sid in mobile_ids_for_init]
            
            seed = init_seed + i
            _score_based_initial_placement(
                mc_sf,
                fixed_ps,
                mobile_ps,
                init_distance_a,
                init_trials,
                seed,
            )
            
            # Re-assign velocities after placement
            md.assign_velocities(temperature_k)

            # Short hybrid run
            n_blocks = max(1, md_steps // max(1, md_steps_per_block))
            for b in range(n_blocks):
                _rb_mc_step(
                    rb_particles_per_group,
                    mc_sf,
                    mc_temp,
                    mc_steps_per_block,
                    rb_max_trans,
                    math.radians(rb_max_rot_deg),
                )
                md.optimize(md_steps_per_block)
                if (b + 1) % 10 == 0:
                    IMP.rmf.save_frame(final_rh, f"restart_{i}_block_{b}")
                    if (b + 1) % 50 == 0:
                        score = full_sf.evaluate(False)
                        print(f"    [Restart {i+1} Block {b+1}/{n_blocks}] score={score:.2f}", flush=True)

            # Save final frame to collective RMF
            score = full_sf.evaluate(False)
            ke = md.get_kinetic_energy()
            
            IMP.rmf.save_frame(final_rh, str(i))
            fstats.write(f"{i}\t{score:.6f}\t{ke:.6f}\n")
            fstats.flush()
            
            if (i+1) % 100 == 0 or i == 0:
                print(f"  [Restart {i+1}/{n_restarts}] score={score:.2f} KE={ke:.2f}", flush=True)


    if converged:
        print(f"Simulation stopped early due to convergence.")
    print(f"Meta-sampling complete. Final frames saved to {final_rmf_path}")


def _parse_paths(system_cif, systems_dir):
    if system_cif:
        return [system_cif]
    if not systems_dir:
        raise click.ClickException("Provide --system-cif or --systems-dir")
    paths = []
    for fn in sorted(os.listdir(systems_dir)):
        if fn.endswith(".cif"):
            paths.append(os.path.join(systems_dir, fn))
    if not paths:
        raise click.ClickException(f"No .cif files found in {systems_dir}")
    return paths


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--system-cif", type=click.Path(path_type=str), default=None)
@click.option("--systems-dir", type=click.Path(path_type=str), default=None)
@click.option(
    "--output-root",
    type=click.Path(path_type=str),
    default=None,
    show_default=True,
    help="Root directory for trajectory output. Defaults to 'output/trajs' relative to the current working directory.",
)
@click.option(
    "--md-steps",
    type=int,
    default=20000,
    show_default=True,
    help="Total MD steps after MC pre-stage.",
)
@click.option("--write-every", type=int, default=500, show_default=True)
@click.option(
    "--sampling-mode",
    type=click.Choice(["simple_md", "hybrid_md_mc", "multi_restart"]),
    default="hybrid_md_mc",
    show_default=True,
)
@click.option("--mc-steps", type=int, default=8, show_default=True)
@click.option("--mc-temperature", type=float, default=2.0, show_default=True)
@click.option("--rb-max-translation", type=float, default=0.35, show_default=True)
@click.option("--rb-max-rotation-deg", type=float, default=4.0, show_default=True)
@click.option(
    "--friction-ps",
    type=float,
    default=None,
    show_default=True,
    help="Langevin friction coefficient in ps^-1. Overrides system CIF value.",
)
@click.option(
    "--temperature-k",
    type=float,
    default=None,
    show_default=True,
    help="MD temperature in K. Overrides system CIF value.",
)
@click.option("--mobile-group", type=str, default=None, show_default=True)
@click.option(
    "--fixed-flex-mode",
    type=click.Choice(["static", "flex"]),
    default="static",
    show_default=True,
    help="'static': fixed component fully rigid. 'flex': atoms in the {name}_flex group are released.",
)
@click.option(
    "--freeze-mobile-rings/--no-freeze-mobile-rings",
    default=True,
    show_default=True,
    help="Fallback only: use ring detection as MD-fixed set when mmCIF has no _ff_dof_md_fixed_member.",
)
@click.option("--init-placement/--no-init-placement", default=True, show_default=True)
@click.option("--init-distance-a", type=float, default=62.5, show_default=True)
@click.option("--init-trials", type=int, default=200, show_default=True)
@click.option("--init-seed", type=int, default=42, show_default=True)
@click.option(
    "--com-pull-k", type=float, default=0.0, help="Weak harmonic force pulling guest COM to host COM."
)
@click.option("--go-mobile-k", type=float, default=3.0, show_default=True)
@click.option("--go-fixed-k", type=float, default=2.0, show_default=True)
@click.option("--go-cutoff", type=float, default=6.0, show_default=True)
@click.option(
    "--log-every-frames",
    type=int,
    default=10,
    show_default=True,
    help="Progress print frequency in written frames.",
)
@click.option("--n-restarts", type=int, default=1, show_default=True, help="Number of independent restarts (for multi_restart mode).")
def main(
    system_cif,
    systems_dir,
    output_root,
    md_steps,
    write_every,
    sampling_mode,
    mc_steps,
    mc_temperature,
    rb_max_translation,
    rb_max_rotation_deg,
    friction_ps,
    temperature_k,
    mobile_group,
    fixed_flex_mode,
    freeze_mobile_rings,
    init_placement,
    init_distance_a,
    init_trials,
    init_seed,
    com_pull_k,
    go_mobile_k,
    go_fixed_k,
    go_cutoff,
    log_every_frames,
    n_restarts,
):
    IMP.setup_from_argv([os.path.basename(__file__)], "IMP FF mmCIF MD (hybrid)")
    IMP.set_check_level(IMP.NONE)

    if output_root is None:
        output_root = str(Path.cwd() / "output" / "trajs")


    for path in _parse_paths(system_cif, systems_dir):
        system = read_dye_forcefield_cif(path)
        _validate_system(system)
        # CLI overrides for sampling params not settable per system CIF
        if friction_ps is not None:
            system.setdefault("sampling", {})["friction_ps"] = friction_ps
        if temperature_k is not None:
            system.setdefault("sampling", {})["temperature_K"] = temperature_k
        _run(
            system,
            path,
            output_root,
            md_steps,
            write_every,
            sampling_mode,
            mc_steps,
            mc_temperature,
            rb_max_translation,
            rb_max_rotation_deg,
            mobile_group,
            fixed_flex_mode,
            freeze_mobile_rings,
            init_placement,
            init_distance_a,
            init_trials,
            init_seed,
            com_pull_k,
            go_mobile_k,
            go_fixed_k,
            go_cutoff,
            log_every_frames,
            n_restarts=n_restarts,
        )


if __name__ == "__main__":
    main()
