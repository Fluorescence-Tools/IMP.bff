"""Build combined protein+dye systems."""

import math
import re

from IMP.bff.cgdye.io.template_cif import read_cgdye_template
from IMP.bff.cgdye.sampling.scoring import build_lj_type_table
from IMP.bff.cgdye.topology.builder import (
    _serial_to_site_atom_names,
    angle_value,
    build_angles,
    build_dihedrals,
    build_graph,
    distance,
    parse_mol2,
    sid,
)


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


def build_combined_system(
    protein_mol2,
    dye_mol2,
    protein_name,
    dye_name,
    protein_template=None,
    dye_template=None,
    default_radius=1.7,
    default_mass=12.0,
):
    p_atoms, p_bonds = parse_mol2(protein_mol2, protein_name)
    d_atoms, d_bonds = parse_mol2(dye_mol2, dye_name)
    p_names = _serial_to_site_atom_names(p_atoms)
    d_names = _serial_to_site_atom_names(d_atoms)
    p_graph = build_graph(p_bonds)
    d_graph = build_graph(d_bonds)

    p_template = read_cgdye_template(protein_template) if protein_template else {}
    d_template = read_cgdye_template(dye_template) if dye_template else {}

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

    return {
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
