"""Read/write compact FF topology mmCIF tables."""

import IMP.bff
import ihm.format
import os
import re

import numpy as np

# --------------------------------------------------------------------------
# forcefield_cif
# --------------------------------------------------------------------------
"""Read/write compact FF topology mmCIF tables."""

#!/usr/bin/env python




class _BaseHandler:
    not_in_file = object()
    omitted = object()
    unknown = object()


def _as_str(v, h):
    if v in (h.not_in_file, h.omitted, h.unknown):
        return None
    return str(v)


def _as_int(v, h):
    s = _as_str(v, h)
    return None if s is None else int(s)


def _as_float(v, h):
    s = _as_str(v, h)
    return None if s is None else float(s)


_SITE_ID_RE = re.compile(r"^(?P<prefix>.+?)(?P<serial>\d+)$")


def _split_site_id(site_id):
    m = _SITE_ID_RE.match(site_id)
    if not m:
        return None, None
    return m.group("prefix"), int(m.group("serial"))


def _expand_site_range(start_id, end_id):
    p1, i1 = _split_site_id(start_id)
    p2, i2 = _split_site_id(end_id)
    if i1 is None or i2 is None:
        if start_id == end_id:
            return [start_id]
        raise ValueError(f"invalid group range: {start_id}..{end_id}")
    if p1 is None or p2 is None or p1 != p2:
        if start_id == end_id:
            return [start_id]
        raise ValueError(f"invalid group range: {start_id}..{end_id}")
    start = int(i1)
    end = int(i2)
    if end < start:
        start, end = end, start
    return [f"{p1}{i}" for i in range(start, end + 1)]


def _compress_int_ranges(nos):
    """Compress a list of integers into sorted contiguous (start, end) pairs."""
    nos = sorted(set(nos))
    if not nos:
        return []
    ranges = []
    start = prev = nos[0]
    for cur in nos[1:]:
        if cur == prev + 1:
            prev = cur
        else:
            ranges.append((start, prev))
            start = prev = cur
    ranges.append((start, prev))
    return ranges


def _compress_site_ranges(site_ids):
    parseable = {}
    literal = set()
    for sid in set(site_ids):
        prefix, serial = _split_site_id(sid)
        if prefix is None:
            literal.add(sid)
            continue
        parseable.setdefault(prefix, []).append(serial)

    ranges = []
    for prefix in sorted(parseable.keys()):
        serials = sorted(set(parseable[prefix]))
        if not serials:
            continue
        start = serials[0]
        prev = serials[0]
        for cur in serials[1:]:
            if cur == prev + 1:
                prev = cur
                continue
            ranges.append((f"{prefix}{start}", f"{prefix}{prev}"))
            start = cur
            prev = cur
        ranges.append((f"{prefix}{start}", f"{prefix}{prev}"))

    for sid in sorted(literal):
        ranges.append((sid, sid))
    return ranges


def _parse_mol2_atom_site_rows(mol2_path, asym_id, entity_id):
    """Parse @<TRIPOS>ATOM section of a MOL2 file into _atom_site rows."""
    rows = []
    in_atom = False
    with open(mol2_path) as fh:
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
            if len(parts) < 6:
                continue
            serial = int(parts[0])
            atom_name = parts[1]
            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
            atype = parts[5]  # TRIPOS type, e.g. "C.3", "N.am"
            comp_id = parts[7] if len(parts) > 7 else "UNK"
            # derive element from TRIPOS type (before first '.')
            element = atype.split(".")[0].capitalize()
            if not element:
                m = re.match(r"[A-Za-z]+", atom_name)
                element = m.group(0)[0].upper() if m else "C"
            rows.append(
                {
                    "group_PDB": "HETATM",
                    "id": serial,
                    "type_symbol": element,
                    "label_atom_id": atom_name,
                    "label_comp_id": comp_id,
                    "label_asym_id": asym_id,
                    "label_entity_id": int(entity_id),
                    "label_seq_id": 1,
                    "Cartn_x": x,
                    "Cartn_y": y,
                    "Cartn_z": z,
                    "occupancy": 1.0,
                    "B_iso_or_equiv": 0.0,
                    "auth_asym_id": asym_id,
                    "auth_comp_id": comp_id,
                    "auth_seq_id": 1,
                    "auth_atom_id": atom_name,
                    "pdbx_PDB_model_num": 1,
                }
            )
    return rows


def _parse_struct_atom_site_rows(struct_path, asym_id, entity_id):
    """Dispatch to MOL2 or PDB parser based on file extension."""
    if struct_path.endswith(".mol2"):
        return _parse_mol2_atom_site_rows(struct_path, asym_id, entity_id)
    return _parse_pdb_atom_site_rows(struct_path, asym_id, entity_id)


def _parse_pdb_atom_site_rows(pdb_path, asym_id, entity_id):
    rows = []
    with open(pdb_path) as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            serial_txt = line[6:11].strip()
            atom_name = line[12:16].strip()
            comp_id = line[17:20].strip() or "UNK"
            auth_asym = line[21:22].strip() or asym_id
            seq_txt = line[22:26].strip()
            x_txt = line[30:38].strip()
            y_txt = line[38:46].strip()
            z_txt = line[46:54].strip()
            occ_txt = line[54:60].strip()
            b_txt = line[60:66].strip()
            element = line[76:78].strip()
            if not element:
                # PDB fallback: infer from atom name
                m = re.match(r"([A-Za-z]+)", atom_name)
                element = (m.group(1)[0] if m else atom_name[:1]).upper()

            if not serial_txt or not x_txt or not y_txt or not z_txt:
                continue
            serial = int(serial_txt)
            seq_id = int(seq_txt) if seq_txt else 1
            rows.append(
                {
                    "group_PDB": "ATOM" if line.startswith("ATOM") else "HETATM",
                    "id": serial,
                    "type_symbol": element,
                    "label_atom_id": atom_name,
                    "label_comp_id": comp_id,
                    "label_asym_id": asym_id,
                    "label_entity_id": int(entity_id),
                    "label_seq_id": seq_id,
                    "Cartn_x": float(x_txt),
                    "Cartn_y": float(y_txt),
                    "Cartn_z": float(z_txt),
                    "occupancy": float(occ_txt) if occ_txt else 1.0,
                    "B_iso_or_equiv": float(b_txt) if b_txt else 0.0,
                    "auth_asym_id": auth_asym,
                    "auth_comp_id": comp_id,
                    "auth_seq_id": seq_id,
                    "auth_atom_id": atom_name,
                }
            )
    return rows


def read_dye_forcefield_cif(path):
    """Read a compact force-field system (sites, bonds, angles, dihedrals, LJ types, ...) from mmCIF."""
    system = {
        "name": os.path.splitext(os.path.basename(path))[0],
        "components": {},
        "sites": [],
        "groups": {},
        "rb_groups": {},
        "md_fixed_groups": {},
        "fixed_groups": [],
        "bond_types": {},
        "angle_types": {},
        "torsion_types": {},
        "improper_types": {},
        "bonds": [],
        "angles": [],
        "dihedrals": [],
        "impropers": [],
        "nonbonded": {"enabled": True, "k": 5.0, "cutoff_A": 6.0},
        "lj_types": {},
        "sampling": {
            "temperature_K": 300.0,
            "friction_ps": 10.0,
            "timestep_fs": 0.25,
            "n_steps": 500000,
            "write_every": 1000,
            "minimize_steps": 200,
        },
    }

    def _as_site_token(v, h):
        s = _as_str(v, h)
        if s is None:
            return None
        if s.isdigit():
            return int(s)
        return s

    def _pick_value(preferred, fallback, h):
        if preferred in (None, h.not_in_file, h.omitted, h.unknown):
            return fallback
        return preferred

    class SystemH(_BaseHandler):
        def __call__(self, name):
            s = _as_str(name, self)
            if s:
                system["name"] = s

    class ComponentH(_BaseHandler):
        def __call__(self, id, mol2_path=None, pdb_path=None, role=None):
            cid = _as_str(id, self)
            mol2 = _as_str(mol2_path, self) if mol2_path is not None else None
            pdb = _as_str(pdb_path, self) if pdb_path is not None else None
            rv = _as_str(role, self) if role is not None else None
            # prefer mol2; fall back to pdb for legacy CIFs
            spec = {"mol2": mol2} if mol2 else {"pdb": pdb}
            if rv:
                spec["role"] = rv
            system["components"][cid] = spec

    class SiteH(_BaseHandler):
        def __call__(
            self,
            site_id,
            component_id,
            site_no=None,
            site_serial=None,
            atom_name=None,
            radius_A=None,
            mass_Da=None,
        ):
            system["sites"].append(
                {
                    "site_no": _as_int(site_no, self),
                    "id": _as_str(site_id, self),
                    "component": _as_str(component_id, self),
                    "atom_name": _as_str(atom_name, self),
                    "site_serial": _as_int(site_serial, self),
                    "radius": _as_float(radius_A, self) or 1.7,
                    "mass": _as_float(mass_Da, self) or 12.0,
                }
            )

    class GroupMemberH(_BaseHandler):
        def __call__(
            self,
            group_id,
            site_id=None,
            site_id_start=None,
            site_id_end=None,
            n=None,
            n_start=None,
            n_end=None,
        ):
            gid = _as_str(group_id, self)
            ni = _as_int(n, self)
            ns = _as_int(n_start, self)
            ne = _as_int(n_end, self)
            if ni is not None:
                system["groups"].setdefault(gid, []).append(ni)
                return
            if ns is not None or ne is not None:
                ns = int(ns if ns is not None else ne)  # type: ignore[arg-type]
                ne = int(ne if ne is not None else ns)
                system["groups"].setdefault(gid, []).extend(range(ns, ne + 1))
                return
            sid = _as_str(site_id, self)
            sid_start = _as_str(site_id_start, self)
            sid_end = _as_str(site_id_end, self)
            if sid:
                system["groups"].setdefault(gid, []).append(sid)
                return
            if sid_start is None and sid_end is None:
                return
            if sid_start is None:
                sid_start = sid_end
            if sid_end is None:
                sid_end = sid_start
            system["groups"].setdefault(gid, []).extend(
                _expand_site_range(sid_start, sid_end)
            )

    class FixedGroupH(_BaseHandler):
        def __call__(self, group_id):
            gid = _as_str(group_id, self)
            if gid:
                system["fixed_groups"].append(gid)

    class RbMemberH(_BaseHandler):
        def __call__(
            self,
            rb_id,
            site_id=None,
            site_id_start=None,
            site_id_end=None,
            n=None,
            n_start=None,
            n_end=None,
        ):
            rid = _as_str(rb_id, self)
            if not rid:
                return
            ni = _as_int(n, self)
            ns = _as_int(n_start, self)
            ne = _as_int(n_end, self)
            if ni is not None:
                system["rb_groups"].setdefault(rid, []).append(ni)
                return
            if ns is not None or ne is not None:
                ns = ns if ns is not None else ne
                ne = ne if ne is not None else ns
                system["rb_groups"].setdefault(rid, []).extend(range(ns, ne + 1))
                return
            sid = _as_str(site_id, self)
            sid_start = _as_str(site_id_start, self)
            sid_end = _as_str(site_id_end, self)
            if sid:
                system["rb_groups"].setdefault(rid, []).append(sid)
                return
            if sid_start is None and sid_end is None:
                return
            if sid_start is None:
                sid_start = sid_end
            if sid_end is None:
                sid_end = sid_start
            system["rb_groups"].setdefault(rid, []).extend(
                _expand_site_range(sid_start, sid_end)
            )

    class MdFixedMemberH(_BaseHandler):
        def __call__(
            self,
            md_fixed_id,
            site_id=None,
            site_id_start=None,
            site_id_end=None,
            n=None,
            n_start=None,
            n_end=None,
        ):
            fid = _as_str(md_fixed_id, self)
            if not fid:
                return
            ni = _as_int(n, self)
            ns = _as_int(n_start, self)
            ne = _as_int(n_end, self)
            if ni is not None:
                system["md_fixed_groups"].setdefault(fid, []).append(ni)
                return
            if ns is not None or ne is not None:
                ns = ns if ns is not None else ne
                ne = ne if ne is not None else ns
                system["md_fixed_groups"].setdefault(fid, []).extend(range(ns, ne + 1))
                return
            sid = _as_str(site_id, self)
            sid_start = _as_str(site_id_start, self)
            sid_end = _as_str(site_id_end, self)
            if sid:
                system["md_fixed_groups"].setdefault(fid, []).append(sid)
                return
            if sid_start is None and sid_end is None:
                return
            if sid_start is None:
                sid_start = sid_end
            if sid_end is None:
                sid_end = sid_start
            system["md_fixed_groups"].setdefault(fid, []).extend(
                _expand_site_range(sid_start, sid_end)
            )

    class SamplingH(_BaseHandler):
        def __call__(
            self,
            temperature_K,
            friction_ps,
            timestep_fs,
            n_steps,
            write_every,
            minimize_steps,
        ):
            for k, v, conv in [
                ("temperature_K", temperature_K, _as_float),
                ("friction_ps", friction_ps, _as_float),
                ("timestep_fs", timestep_fs, _as_float),
                ("n_steps", n_steps, _as_int),
                ("write_every", write_every, _as_int),
                ("minimize_steps", minimize_steps, _as_int),
            ]:
                vv = conv(v, self)
                if vv is not None:
                    system["sampling"][k] = vv

    class NonbondedH(_BaseHandler):
        def __call__(self, enabled, k, cutoff_A):
            ev = _as_str(enabled, self)
            if ev is not None:
                system["nonbonded"]["enabled"] = ev.upper() in {"1", "TRUE", "YES"}
            kv = _as_float(k, self)
            cv = _as_float(cutoff_A, self)
            if kv is not None:
                system["nonbonded"]["k"] = kv
            if cv is not None:
                system["nonbonded"]["cutoff_A"] = cv

    class BondTypeH(_BaseHandler):
        def __call__(self, type_id, k_kcal_mol_A2):
            tid = _as_str(type_id, self)
            system["bond_types"][tid] = {"k": _as_float(k_kcal_mol_A2, self)}

    class AngleTypeH(_BaseHandler):
        def __call__(self, type_id, k_kcal_mol_rad2):
            tid = _as_str(type_id, self)
            system["angle_types"][tid] = {"k": _as_float(k_kcal_mol_rad2, self)}

    class TorsionTypeH(_BaseHandler):
        def __call__(self, type_id, periodicity, phase_rad, k_kcal_mol):
            tid = _as_str(type_id, self)
            system["torsion_types"][tid] = {
                "periodicity": _as_int(periodicity, self),
                "phase_rad": _as_float(phase_rad, self),
                "k": _as_float(k_kcal_mol, self),
            }

    class ImproperTypeH(_BaseHandler):
        def __call__(self, type_id, periodicity, phase_rad, k_kcal_mol):
            tid = _as_str(type_id, self)
            system["improper_types"][tid] = {
                "periodicity": _as_int(periodicity, self),
                "phase_rad": _as_float(phase_rad, self),
                "k": _as_float(k_kcal_mol, self),
            }

    class LjTypeH(_BaseHandler):
        def __call__(self, type_id, element, rmin_half_A, epsilon_kcal_mol):
            tid = _as_str(type_id, self)
            system["lj_types"][tid] = {
                "element": _as_str(element, self),
                "rmin_half": _as_float(rmin_half_A, self),
                "epsilon": _as_float(epsilon_kcal_mol, self),
            }

    class BondH(_BaseHandler):
        def __call__(
            self,
            site_id_1=None,
            site_id_2=None,
            length_A=None,
            type_id=None,
            n1=None,
            n2=None,
        ):
            a = _pick_value(n1, site_id_1, self)
            b = _pick_value(n2, site_id_2, self)
            system["bonds"].append(
                [
                    _as_site_token(a, self),
                    _as_site_token(b, self),
                    _as_float(length_A, self),
                    _as_str(type_id, self),
                ]
            )

    class AngleH(_BaseHandler):
        def __call__(
            self,
            site_id_1=None,
            site_id_2=None,
            site_id_3=None,
            theta_rad=None,
            type_id=None,
            n1=None,
            n2=None,
            n3=None,
        ):
            a = _pick_value(n1, site_id_1, self)
            b = _pick_value(n2, site_id_2, self)
            c = _pick_value(n3, site_id_3, self)
            system["angles"].append(
                [
                    _as_site_token(a, self),
                    _as_site_token(b, self),
                    _as_site_token(c, self),
                    _as_float(theta_rad, self),
                    _as_str(type_id, self),
                ]
            )

    class TorsionH(_BaseHandler):
        def __call__(
            self,
            site_id_1=None,
            site_id_2=None,
            site_id_3=None,
            site_id_4=None,
            type_id=None,
            n1=None,
            n2=None,
            n3=None,
            n4=None,
        ):
            a = _pick_value(n1, site_id_1, self)
            b = _pick_value(n2, site_id_2, self)
            c = _pick_value(n3, site_id_3, self)
            d = _pick_value(n4, site_id_4, self)
            system["dihedrals"].append(
                [
                    _as_site_token(a, self),
                    _as_site_token(b, self),
                    _as_site_token(c, self),
                    _as_site_token(d, self),
                    _as_str(type_id, self),
                ]
            )

    class ImproperH(_BaseHandler):
        def __call__(
            self,
            site_id_1=None,
            site_id_2=None,
            site_id_3=None,
            site_id_4=None,
            type_id=None,
            n1=None,
            n2=None,
            n3=None,
            n4=None,
        ):
            a = _pick_value(n1, site_id_1, self)
            b = _pick_value(n2, site_id_2, self)
            c = _pick_value(n3, site_id_3, self)
            d = _pick_value(n4, site_id_4, self)
            system["impropers"].append(
                [
                    _as_site_token(a, self),
                    _as_site_token(b, self),
                    _as_site_token(c, self),
                    _as_site_token(d, self),
                    _as_str(type_id, self),
                ]
            )

    handlers = {
        "_ff_system": SystemH(),
        "_ff_component": ComponentH(),
        "_ff_site": SiteH(),
        "_ff_group_member": GroupMemberH(),
        "_ff_dof_fixed_group": FixedGroupH(),
        "_ff_dof_rb_member": RbMemberH(),
        "_ff_dof_md_fixed_member": MdFixedMemberH(),
        "_ff_sampling": SamplingH(),
        "_ff_nonbonded": NonbondedH(),
        "_ff_bond_type": BondTypeH(),
        "_ff_angle_type": AngleTypeH(),
        "_ff_torsion_type": TorsionTypeH(),
        "_ff_improper_type": ImproperTypeH(),
        "_ff_lj_type": LjTypeH(),
        "_ff_bond": BondH(),
        "_ff_angle": AngleH(),
        "_ff_torsion": TorsionH(),
        "_ff_improper": ImproperH(),
    }

    with open(path) as fh:
        r = ihm.format.CifReader(fh, handlers)
        r.read_file()

    site_no_to_id = {}
    next_no = 1
    for s in system["sites"]:
        sid = s["id"]
        num = s["site_no"]
        if num is None:
            while next_no in site_no_to_id:
                next_no += 1
            num = next_no
            next_no += 1
            s["site_no"] = num
        site_no_to_id[int(num)] = sid

    def _resolve_site_ref(v):
        if isinstance(v, int):
            if v not in site_no_to_id:
                raise ValueError(f"unknown site number {v} in {path}")
            return site_no_to_id[v]
        return v

    system["bonds"] = [
        [_resolve_site_ref(a), _resolve_site_ref(b), length, tid]
        for a, b, length, tid in system["bonds"]
    ]
    system["angles"] = [
        [_resolve_site_ref(a), _resolve_site_ref(b), _resolve_site_ref(c), theta, tid]
        for a, b, c, theta, tid in system["angles"]
    ]
    system["dihedrals"] = [
        [
            _resolve_site_ref(a),
            _resolve_site_ref(b),
            _resolve_site_ref(c),
            _resolve_site_ref(d),
            tid,
        ]
        for a, b, c, d, tid in system["dihedrals"]
    ]
    system["impropers"] = [
        [
            _resolve_site_ref(a),
            _resolve_site_ref(b),
            _resolve_site_ref(c),
            _resolve_site_ref(d),
            tid,
        ]
        for a, b, c, d, tid in system["impropers"]
    ]
    system["groups"] = {
        gid: [_resolve_site_ref(v) for v in sids]
        for gid, sids in system["groups"].items()
    }
    system["rb_groups"] = {
        rid: [_resolve_site_ref(v) for v in sids]
        for rid, sids in system["rb_groups"].items()
    }
    system["md_fixed_groups"] = {
        fid: [_resolve_site_ref(v) for v in sids]
        for fid, sids in system["md_fixed_groups"].items()
    }

    return forcefield_system_from_dict(system)


def write_dye_forcefield_cif(path, system):
    """Write a force-field system (see :func:`read_dye_forcefield_cif`) to mmCIF.

    Accepts a :class:`IMP.bff.DyeForceFieldSystem` or the plain dictionary that
    preceded it, normalising once here so everything below sees one type. The
    dictionary form is what callers assembling a system by hand still write.
    """
    system = as_forcefield_system(system)
    with open(path, "w") as out:
        w = ihm.format.CifWriter(out)
        w.start_block((system.name or "ff_system"))
        out_dir = os.path.dirname(os.path.abspath(path))

        with w.category("_ff_system") as c:
            c.write(name=(system.name or "ff_system"))

        with w.loop("_ff_component", ["id", "mol2_path", "pdb_path", "role"]) as l:
            for cid, spec in system.components.items():
                l.write(
                    id=cid,
                    # `None`, not `""`: the type uses an empty string for an
                    # absent path, and mmCIF spells absent as an omitted value.
                    # Writing "" round-trips into a path that is falsy but not
                    # None, which every `is None` check downstream misses.
                    mol2_path=spec.mol2_path or None,
                    pdb_path=spec.pdb_path or None,
                    role=spec.role or None,
                )

        if system.probes:
            with w.loop("_flr_probe_list", ["probe_id", "chromophore_name", "probe_origin", "probe_link_type"]) as l:
                for p in system.probes:
                    l.write(
                        probe_id=p.id,
                        chromophore_name=p.name,
                        probe_origin=p.origin,
                        probe_link_type=p.link_type
                    )

        # `labeling_sites` was **deleted**, not typed. It was initialised to an
        # empty list in the topology builder, never appended to anywhere, and
        # written here into a `_flr_poly_probe_position` loop that therefore
        # could never emit a row. It survived because a dictionary key that is
        # always empty looks exactly like one that is sometimes populated.

        atom_site_rows = []
        asym_ord = 0
        for ent_id, (cid, spec) in enumerate(
            system.components.items(), start=1
        ):
            struct_path = spec.mol2_path or spec.pdb_path
            if not struct_path:
                continue
            if not os.path.isabs(struct_path):
                struct_path = os.path.abspath(os.path.join(out_dir, struct_path))
            if not os.path.exists(struct_path):
                continue
            asym_ord += 1
            asym_id = chr(ord("A") + ((asym_ord - 1) % 26))
            atom_site_rows.extend(
                _parse_struct_atom_site_rows(struct_path, asym_id, ent_id)
            )

        if atom_site_rows:
            atom_site_cols = [
                "group_PDB",
                "id",
                "type_symbol",
                "label_atom_id",
                "label_comp_id",
                "label_asym_id",
                "label_entity_id",
                "label_seq_id",
                "Cartn_x",
                "Cartn_y",
                "Cartn_z",
                "occupancy",
                "B_iso_or_equiv",
                "pdbx_PDB_model_num",
                "auth_asym_id",
                "auth_comp_id",
                "auth_seq_id",
                "auth_atom_id",
            ]
            with w.loop("_atom_site", atom_site_cols) as l:
                atom_id = 1
                for r in atom_site_rows:
                    l.write(
                        group_PDB=r["group_PDB"],
                        id=atom_id,
                        type_symbol=r["type_symbol"],
                        label_atom_id=r["label_atom_id"],
                        label_comp_id=r["label_comp_id"],
                        label_asym_id=r["label_asym_id"],
                        label_entity_id=r["label_entity_id"],
                        label_seq_id=r["label_seq_id"],
                        Cartn_x=r["Cartn_x"],
                        Cartn_y=r["Cartn_y"],
                        Cartn_z=r["Cartn_z"],
                        occupancy=r["occupancy"],
                        B_iso_or_equiv=r["B_iso_or_equiv"],
                        pdbx_PDB_model_num=1,
                        auth_asym_id=r["auth_asym_id"],
                        auth_comp_id=r["auth_comp_id"],
                        auth_seq_id=r["auth_seq_id"],
                        auth_atom_id=r["auth_atom_id"],
                    )
                    atom_id += 1

        sites = system.sites
        compact_site_loop = all(
            abs(float(s.radius) - 1.7) < 1e-12
            and abs(float(s.mass) - 12.0) < 1e-12
            for s in sites
        )

        def _site_atom_name(site):
            if site.atom_name:
                return site.atom_name
            sid = str(site.id)
            if "/" in sid:
                return sid.split("/", 1)[1]
            if ":" in sid:
                return sid.split(":", 1)[1]
            return sid

        site_rows = []
        site_id_to_no = {}
        used_site_nos = set()
        next_no = 1
        for s in sites:
            sid = s.id
            num = s.site_no
            if num is not None:
                num = int(num)
                if num in used_site_nos:
                    num = None
            if num is None:
                while next_no in used_site_nos:
                    next_no += 1
                num = next_no
                next_no += 1
            used_site_nos.add(num)
            site_rows.append((s, num))
            site_id_to_no[sid] = num

        if compact_site_loop:
            with w.loop(
                "_ff_site",
                ["site_no", "site_id", "component_id", "atom_name", "site_serial"],
            ) as l:
                for s, num in site_rows:
                    l.write(
                        site_no=num,
                        site_id=s.id,
                        component_id=s.component,
                        atom_name=_site_atom_name(s),
                        site_serial=s.site_serial,
                    )
        else:
            with w.loop(
                "_ff_site",
                [
                    "site_no",
                    "site_id",
                    "component_id",
                    "atom_name",
                    "site_serial",
                    "radius_A",
                    "mass_Da",
                ],
            ) as l:
                for s, num in site_rows:
                    l.write(
                        site_no=num,
                        site_id=s.id,
                        component_id=s.component,
                        atom_name=_site_atom_name(s),
                        site_serial=s.site_serial,
                        radius_A=s.radius,
                        mass_Da=s.mass,
                    )

        def _write_no_ranges(loop_name, id_col, member_col, mapping):
            """Write group membership as integer site-number ranges."""
            rows = []
            for gid, sids in mapping.items():
                nos = sorted(site_id_to_no[s] for s in sids)
                for start, end in _compress_int_ranges(nos):
                    rows.append((gid, start, end))
            has_ranges = any(s != e for _, s, e in rows)
            if has_ranges:
                with w.loop(loop_name, [id_col, "n_start", "n_end"]) as l:
                    for gid, start, end in rows:
                        l.write(**{id_col: gid, "n_start": start, "n_end": end})
            else:
                with w.loop(loop_name, [id_col, "n"]) as l:
                    for gid, start, _ in rows:
                        l.write(**{id_col: gid, "n": start})

        _write_no_ranges("_ff_group_member", "group_id", "n", system.groups)

        with w.loop("_ff_dof_fixed_group", ["group_id"]) as l:
            for gid in system.fixed_groups:
                l.write(group_id=gid)

        _write_no_ranges("_ff_dof_rb_member", "rb_id", "n", system.rb_groups)
        _write_no_ranges(
            "_ff_dof_md_fixed_member",
            "md_fixed_id",
            "n",
            system.md_fixed_groups,
        )

        samp = system.sampling
        with w.category("_ff_sampling") as c:
            c.write(
                temperature_K=samp.temperature_K,
                friction_ps=samp.friction_ps,
                timestep_fs=samp.timestep_fs,
                n_steps=samp.n_steps,
                write_every=samp.write_every,
                minimize_steps=samp.minimize_steps,
            )

        nb = system.nonbonded
        with w.category("_ff_nonbonded") as c:
            c.write(
                enabled="YES" if nb.enabled else "NO",
                k=nb.k,
                cutoff_A=nb.cutoff,
            )

        with w.loop("_ff_bond_type", ["type_id", "k_kcal_mol_A2"]) as l:
            for tid, k in system.bond_types.items():
                l.write(type_id=tid, k_kcal_mol_A2=k)

        with w.loop("_ff_angle_type", ["type_id", "k_kcal_mol_rad2"]) as l:
            for tid, k in system.angle_types.items():
                l.write(type_id=tid, k_kcal_mol_rad2=k)

        with w.loop(
            "_ff_torsion_type", ["type_id", "periodicity", "phase_rad", "k_kcal_mol"]
        ) as l:
            for tid, t in system.torsion_types.items():
                l.write(
                    type_id=tid,
                    periodicity=t.periodicity,
                    phase_rad=t.phase,
                    k_kcal_mol=t.k,
                )

        with w.loop(
            "_ff_improper_type", ["type_id", "periodicity", "phase_rad", "k_kcal_mol"]
        ) as l:
            for tid, t in system.improper_types.items():
                l.write(
                    type_id=tid,
                    periodicity=t.periodicity,
                    phase_rad=t.phase,
                    k_kcal_mol=t.k,
                )

        if system.lj_types:
            with w.loop(
                "_ff_lj_type",
                ["type_id", "element", "rmin_half_A", "epsilon_kcal_mol"],
            ) as l:
                for tid, t in system.lj_types.items():
                    l.write(
                        type_id=tid,
                        element=t.element,
                        rmin_half_A=t.rmin_half,
                        epsilon_kcal_mol=t.epsilon,
                    )

        def _no(sid):
            return site_id_to_no[sid]

        with w.loop("_ff_bond", ["n1", "n2", "length_A", "type_id"]) as l:
            for bd in system.bonds:
                l.write(n1=_no(bd.site_a), n2=_no(bd.site_b),
                        length_A=bd.length, type_id=bd.type_id)

        with w.loop("_ff_angle", ["n1", "n2", "n3", "theta_rad", "type_id"]) as l:
            for an in system.angles:
                l.write(n1=_no(an.site_a), n2=_no(an.site_b), n3=_no(an.site_c),
                        theta_rad=an.theta, type_id=an.type_id)

        with w.loop("_ff_torsion", ["n1", "n2", "n3", "n4", "type_id"]) as l:
            for to in system.dihedrals:
                l.write(n1=_no(to.site_a), n2=_no(to.site_b), n3=_no(to.site_c),
                        n4=_no(to.site_d), type_id=to.type_id)

        with w.loop("_ff_improper", ["n1", "n2", "n3", "n4", "type_id"]) as l:
            for to in system.impropers:
                l.write(n1=_no(to.site_a), n2=_no(to.site_b), n3=_no(to.site_c),
                        n4=_no(to.site_d), type_id=to.type_id)


# --------------------------------------------------------------------------
# nmr_cif -- retired 2026-08-19
#
# The NMR-restraint dialect and its reader/writer are in
# ``junk/nmr_restraints/``. Their only consumer was the cgdye MD runner's
# optional distance-restraint path, which was retired with them.
# --------------------------------------------------------------------------
"""Read/write NMR restraints in official mmCIF categories."""

#!/usr/bin/env python




def _distance_type_to_ihm(distance_type):
    if distance_type == "AtomLowerBound":
        return "lower bound"
    if distance_type == "AtomUpperBound":
        return "upper bound"
    return "lower and upper bound"


def _ihm_to_distance_type(restraint_type):
    if restraint_type == "lower bound":
        return "AtomLowerBound"
    if restraint_type == "upper bound":
        return "AtomUpperBound"
    return "AtomDistance"


def _set_order(k):
    order = {
        "subunit_prefix": 0,
        "global_prefix": 1,
        "common_prefix": 2,
        "single_s1_prefix": 3,
    }
    return (order.get(k, 100), k)


def _group_prefix_pairs(restraint_sets):
    pairs = []
    for key, prefix in restraint_sets.items():
        if not key.endswith("_prefix"):
            continue
        p = "" if prefix is None else str(prefix)
        if not p:
            continue
        pairs.append((key, p))
    pairs.sort(key=lambda kv: (len(kv[1]) * -1, _set_order(kv[0])))
    return pairs


def _infer_group_for_name(name, key_prefix_pairs):
    for key, prefix in key_prefix_pairs:
        if name.startswith(prefix):
            return key
    return None


def _constraint_subtype(set_key):
    t = set_key.lower()
    if "common" in t or "csp" in t:
        return "CSP"
    return "NOE"


def _prefix_for_set_key(set_key, subtype, used_prefixes):
    defaults = {
        "subunit_prefix": "SUB_",
        "global_prefix": "GLOB_",
        "common_prefix": "CSP_",
        "single_s1_prefix": "S1ONLY_",
        "single_s1_common_prefix": "",
    }
    if set_key in defaults:
        return defaults[set_key]
    if subtype == "CSP":
        base = "CSP_"
    else:
        # Use a sanitized version of the set_key or a unique ID
        base = f"{set_key}_" if set_key else f"G{len(used_prefixes) + 1}_"
    if base in used_prefixes:
        i = 2
        while f"{base}{i}_" in used_prefixes:
            i += 1
        return f"{base}{i}_"
    return base








# --------------------------------------------------------------------------
# rotamer_cif
# --------------------------------------------------------------------------
"""Read/write rotamer libraries.

Coordinates stored as numpy .npy files (shape: n_rotamers x n_atoms x 3).
Weights stored as plain text (one per line). Atom names stored as text.
"""

#!/usr/bin/env python




def _base_path(rmf_or_npy_path):
    if rmf_or_npy_path.endswith(".rmf"):
        return rmf_or_npy_path[:-4]
    if rmf_or_npy_path.endswith(".npy"):
        return rmf_or_npy_path[:-4]
    return rmf_or_npy_path


def read_rotamer_library(path):
    """Read a rotamer library from numpy/text files.

    ``path`` can be the .npy, .rmf, or base path (without extension).
    Files expected: ``{base}_coords.npy``, ``{base}_weights.txt``,
    ``{base}_atoms.txt``.

    Returns dict with:
        - id: list of rotamer IDs (1-indexed ints)
        - weight: list of weights (float)
        - atom_names: list of atom name strings
        - coords: dict {rotamer_id: ndarray shape (n_atoms, 3)}
    """
    base = _base_path(path)
    coords = np.load(base + "_coords.npy")
    n_rotamers = coords.shape[0]

    weights_file = base + "_weights.txt"
    if os.path.exists(weights_file):
        with open(weights_file) as fh:
            weights = [float(line.strip()) for line in fh if line.strip()]
    else:
        weights = [1.0 / n_rotamers] * n_rotamers

    atoms_file = base + "_atoms.txt"
    if os.path.exists(atoms_file):
        with open(atoms_file) as fh:
            atom_names = [line.strip() for line in fh if line.strip()]
    else:
        atom_names = [f"AT{i}" for i in range(coords.shape[1])]

    data = {
        "id": list(range(1, n_rotamers + 1)),
        "weight": weights,
        "atom_names": atom_names,
        "coords": {i + 1: coords[i] for i in range(n_rotamers)},
    }
    return data


def write_rotamer_library(path, library):
    """Write a rotamer library to numpy/text files.

    ``library`` must have:
        - weight: list of floats
        - atom_names: list of atom name strings
        - coords: dict {rotamer_id: list/array of [x, y, z] per atom}
    """
    base = _base_path(path)
    atom_names = library["atom_names"]
    n_atoms = len(atom_names)
    ids = sorted(library["coords"].keys())

    coord_array = np.zeros((len(ids), n_atoms, 3))
    for i, rid in enumerate(ids):
        c = library["coords"][rid]
        arr = np.asarray(c)
        if arr.ndim == 1 and arr.shape[0] == 3:
            arr = arr.reshape(1, 3)
        coord_array[i, : arr.shape[0]] = arr[:n_atoms]

    np.save(base + "_coords.npy", coord_array)

    with open(base + "_weights.txt", "w") as fh:
        for w in library["weight"]:
            fh.write(f"{w}\n")

    with open(base + "_atoms.txt", "w") as fh:
        for name in atom_names:
            fh.write(f"{name}\n")


def normalize_weights(library):
    """Normalize rotamer weights to sum to 1.0 in-place."""
    total = sum(library["weight"])
    if total > 0:
        library["weight"] = [w / total for w in library["weight"]]


# --------------------------------------------------------------------------
# template_cif
# --------------------------------------------------------------------------
"""Read/write cgdye feature templates in mmCIF format."""

#!/usr/bin/env python




def _as_bool(v, h):
    s = _as_str(v, h)
    if s is None:
        return None
    return s.upper() in {"1", "YES", "TRUE"}


def read_component_template_cif(path, with_dye_metadata=False):
    """Read a cgdye component template from mmCIF.

    Named regions and per-atom features, plus the impropers the template
    declares. With *with_dye_metadata* -- which :func:`read_dye_template_cif`
    asks for by name -- the ``_cgdye_metadata`` category is read too: the
    centre atom, the two dipole atoms, and the charged-atom lists. It is off
    by default so this function returns exactly the three keys it always did.

    This was two functions. The dye one was the component one plus that
    handler and its five keys, copied; 27 lines of 122 differed, and of those
    27 all but three were the copy itself.
    """

    template = {
        "name": os.path.splitext(os.path.basename(path))[0],
        "features": {},
        "impropers": [],
        "center_atom": None,
        "dipole_atom_1": None,
        "dipole_atom_2": None,
        "positive_atoms": [],
        "negative_atoms": [],
    }

    class TemplateH(_BaseHandler):
        def __call__(self, name):
            n = _as_str(name, self)
            if n:
                template["name"] = n

    _metadata_map = {
        "center_atom": lambda v: v if v else None,
        "dipole_atom_1": lambda v: v if v else None,
        "dipole_atom_2": lambda v: v if v else None,
        "positive_atoms": lambda v: [a.strip() for a in v.split(",")] if v else [],
        "negative_atoms": lambda v: [a.strip() for a in v.split(",")] if v else [],
    }

    class MetadataH(_BaseHandler):
        def __call__(self, key, value):
            k = _as_str(key, self)
            v = _as_str(value, self)
            if not k or not v or k not in _metadata_map:
                return
            result = _metadata_map[k](v)
            if k in ("positive_atoms", "negative_atoms"):
                template[k].extend(result)
            else:
                template[k] = result

    class FeatureH(_BaseHandler):
        def __call__(self, feature_id, rb, md_fixed, feature_type, region_color):
            fid = _as_str(feature_id, self)
            if not fid:
                return
            template["features"].setdefault(
                fid,
                {
                    "rb": False,
                    "md_fixed": False,
                    "feature_type": "dof",
                    "region_color": None,
                    "atoms": [],
                },
            )
            rbv = _as_bool(rb, self)
            mdv = _as_bool(md_fixed, self)
            ftv = _as_str(feature_type, self)
            rcv = _as_str(region_color, self)
            if rbv is not None:
                template["features"][fid]["rb"] = rbv
            if mdv is not None:
                template["features"][fid]["md_fixed"] = mdv
            if ftv is not None:
                template["features"][fid]["feature_type"] = ftv
            if rcv is not None and rcv != ".":
                template["features"][fid]["region_color"] = rcv

    class FeatureAtomH(_BaseHandler):
        def __call__(self, feature_id, atom_name, occurrence):
            fid = _as_str(feature_id, self)
            anm = _as_str(atom_name, self)
            occ = _as_int(occurrence, self)
            if not fid or not anm:
                return
            template["features"].setdefault(
                fid,
                {
                    "rb": False,
                    "md_fixed": False,
                    "feature_type": "dof",
                    "region_color": None,
                    "atoms": [],
                },
            )
            template["features"][fid]["atoms"].append(
                {
                    "name": anm,
                    "occurrence": 1 if occ is None else occ,
                }
            )

    class ImproperH(_BaseHandler):
        def __call__(self, center_atom, type):
            ca = _as_str(center_atom, self)
            it = _as_str(type, self)
            if not ca or not it:
                return
            template["impropers"].append(
                {
                    "center_atom": ca,
                    "type": it,
                }
            )

    handlers = {
        "_cgdye_template": TemplateH(),
        "_cgdye_feature": FeatureH(),
        "_cgdye_feature_atom": FeatureAtomH(),
        "_cgdye_improper": ImproperH(),
    }

    if with_dye_metadata:
        handlers["_cgdye_metadata"] = MetadataH()
    else:
        for key in ("center_atom", "dipole_atom_1", "dipole_atom_2",
                    "positive_atoms", "negative_atoms"):
            template.pop(key)

    with open(path) as fh:
        r = ihm.format.CifReader(fh, handlers)
        r.read_file()

    return template


def _write_cif_safe(path, template_name, write_fn):
    import io

    buf = io.StringIO()
    w = ihm.format.CifWriter(buf)
    w.start_block(template_name)
    write_fn(w)
    raw = buf.getvalue()
    lines = raw.split("\n")
    fixed = []
    for line in lines:
        stripped = line.lstrip()
        if stripped and not stripped.startswith(("_", "loop_", "data_", "#")):
            parts = stripped.split()
            new_parts = []
            for p in parts:
                if p.startswith("#"):
                    p = f'"{p}"'
                new_parts.append(p)
            line = " " * (len(line) - len(stripped)) + " ".join(new_parts)
        fixed.append(line)
    with open(path, "w") as out:
        out.write("\n".join(fixed))


def write_component_template_cif(path, template):
    """Write a component feature template (see :func:`read_component_template_cif`) to mmCIF."""
    def _write(w):
        _write_header(w, template)
        _write_features(w, template)
        _write_feature_atoms(w, template)
        _write_impropers(w, template)

    _write_cif_safe(path, template.get("name", "cgdye_template"), _write)


def _write_header(w, template):
    with w.category("_cgdye_template") as c:
        c.write(name=template.get("name", "cgdye_template"))


def _write_features(w, template):
    with w.loop(
        "_cgdye_feature",
        ["feature_id", "feature_type", "rb", "md_fixed", "region_color"],
    ) as l:
        for fid, spec in template.get("features", {}).items():
            l.write(
                feature_id=fid,
                feature_type=spec.get("feature_type", "dof"),
                rb="YES" if spec.get("rb", False) else "NO",
                md_fixed="YES" if spec.get("md_fixed", False) else "NO",
                region_color=spec.get("region_color") or ".",
            )


def _write_feature_atoms(w, template):
    with w.loop(
        "_cgdye_feature_atom",
        ["feature_id", "atom_name", "occurrence"],
    ) as l:
        for fid, spec in template.get("features", {}).items():
            for a in spec.get("atoms", []):
                l.write(
                    feature_id=fid,
                    atom_name=a.get("name"),
                    occurrence=int(a.get("occurrence", 1)),
                )


def _write_impropers(w, template):
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


def region_features(template):
    """Return dict of {feature_id: region_color} for features with non-null region_color."""
    out = {}
    for fid, spec in template.get("features", {}).items():
        rc = spec.get("region_color")
        if rc:
            out[fid] = rc
    return out


# ============================================================================
# Dye-specific template handlers (for IMP rotamer library system)
# ============================================================================


def read_dye_template_cif(path):
    """A component template read with its dye metadata -- see
    :func:`read_component_template_cif`, which this names."""
    return read_component_template_cif(path, with_dye_metadata=True)


def write_dye_template_cif(path, template):
    """Write a dye template to mmCIF format."""

    def _write(w):
        _write_header(w, template)

        metadata_entries = []
        for key in ("center_atom", "dipole_atom_1", "dipole_atom_2"):
            val = template.get(key)
            if val:
                metadata_entries.append((key, val))
        for key in ("positive_atoms", "negative_atoms"):
            val = template.get(key, [])
            if val:
                metadata_entries.append((key, ", ".join(val)))
        if metadata_entries:
            with w.loop("_cgdye_metadata", ["key", "value"]) as l:
                for k, v in metadata_entries:
                    l.write(key=k, value=v)

        _write_features(w, template)
        _write_feature_atoms(w, template)
        _write_impropers(w, template)

    _write_cif_safe(path, template.get("name", "dye_template"), _write)


def as_forcefield_system(system):
    """A system, whichever way it was given.

    The typed object passes through; a dictionary is converted. Public entry
    points that take a system call this once, so everything below them sees
    :class:`IMP.bff.DyeForceFieldSystem` and nothing has to ask.

    The dictionary spelling is kept because assembling a system by hand -- a
    handful of sites and one bond, in a test or a script -- is genuinely easier
    as a literal than as a dozen constructor calls. What is not kept is the
    dictionary reaching any of the code that *uses* a system.
    """
    if isinstance(system, dict):
        return forcefield_system_from_dict(system)
    return system


def forcefield_system_from_dict(d):
    """A parsed or built dictionary as a typed :class:`IMP.bff.DyeForceFieldSystem`.

    **One conversion, used by both producers.** The mmCIF reader and
    ``cgdye.topology.build_dye_protein_system`` both built "the system", and
    they did not agree: the reader's sites carried ``site_no`` and no element,
    the builder's carried an element and no ``site_no``, and the writer coped
    with ``s.site_no``. Funnelling both through here is what makes the
    shape single.

    **Group members are canonicalised to site id tokens.** The mmCIF schema
    spells a member two ways -- by index (``n``, or ``n_start``/``n_end`` for a
    run) or by id (``site_id``, ...) -- and the dictionary carried whichever
    the file used, so a consumer could receive ``4`` or ``"CX4/S1"`` for the
    same site with no way to tell which. An index is resolved here against the
    site list; anything already a string is kept.
    """
    sites_in = list(d.get("sites") or [])
    by_no = {}
    for i, row in enumerate(sites_in):
        by_no[i] = str(row.get("id") or "")
        n = row.get("site_no")
        if n is not None:
            by_no[int(n)] = str(row.get("id") or "")

    def _token(x):
        if isinstance(x, str):
            return x
        return by_no.get(int(x), str(x))

    sys_ = IMP.bff.DyeForceFieldSystem(str(d.get("name") or ""))

    comps = IMP.bff.FFComponentMap()
    for cid, spec in (d.get("components") or {}).items():
        c = IMP.bff.FFComponent()
        # Two spellings, because two producers: the mmCIF reader stores
        # `mol2_path` (the column name) and the topology builder stores `mol2`.
        # Accepting both here is what lets the one type serve both.
        c.mol2_path = str(spec.get("mol2_path") or spec.get("mol2") or "")
        c.pdb_path = str(spec.get("pdb_path") or spec.get("pdb") or "")
        c.role = str(spec.get("role") or "")
        comps[str(cid)] = c
    sys_.components = comps

    sites = IMP.bff.FFSiteVector()
    for i, row in enumerate(sites_in):
        st = IMP.bff.FFSite()
        st.id = str(row.get("id") or "")
        st.component = str(row.get("component") or "")
        st.atom_name = str(row.get("atom_name") or "")
        st.element = str(row.get("element") or "")
        st.site_no = int(row.get("site_no") if row.get("site_no") is not None else i)
        st.site_serial = int(row.get("site_serial") or 0)
        st.radius = float(row.get("radius") if row.get("radius") is not None else 1.7)
        st.mass = float(row.get("mass") if row.get("mass") is not None else 12.0)
        sites.append(st)
    sys_.sites = sites

    def _groups(key):
        m = IMP.bff.MapStringVectorString()
        for gid, members in (d.get(key) or {}).items():
            v = IMP.bff.VectorString()
            for x in members:
                v.append(_token(x))
            m[str(gid)] = v
        return m
    sys_.groups = _groups("groups")
    sys_.rb_groups = _groups("rb_groups")
    sys_.md_fixed_groups = _groups("md_fixed_groups")

    fg = IMP.bff.VectorString()
    for g in d.get("fixed_groups") or []:
        fg.append(str(g))
    sys_.fixed_groups = fg

    def _scalar_types(key):
        m = IMP.bff.MapStringDouble()
        for tid, params in (d.get(key) or {}).items():
            m[str(tid)] = float((params or {}).get("k") or 0.0)
        return m
    sys_.bond_types = _scalar_types("bond_types")
    sys_.angle_types = _scalar_types("angle_types")

    def _torsion_types(key):
        m = IMP.bff.FFTorsionTypeMap()
        for tid, params in (d.get(key) or {}).items():
            params = params or {}
            t = IMP.bff.FFTorsionType()
            t.periodicity = int(params.get("periodicity") or 1)
            t.phase = float(params.get("phase_rad") or 0.0)
            t.k = float(params.get("k") or 0.0)
            m[str(tid)] = t
        return m
    sys_.torsion_types = _torsion_types("torsion_types")
    sys_.improper_types = _torsion_types("improper_types")

    lj = IMP.bff.FFLJTypeMap()
    for tid, params in (d.get("lj_types") or {}).items():
        params = params or {}
        t = IMP.bff.FFLJType()
        t.element = str(params.get("element") or "")
        t.rmin_half = float(params.get("rmin_half") or 0.0)
        t.epsilon = float(params.get("epsilon") or 0.0)
        lj[str(tid)] = t
    sys_.lj_types = lj

    bonds = IMP.bff.FFBondVector()
    for row in d.get("bonds") or []:
        a, b, length, tid = row
        r = IMP.bff.FFBond()
        r.site_a, r.site_b = _token(a), _token(b)
        r.length = float(length or 0.0)
        r.type_id = str(tid or "")
        bonds.append(r)
    sys_.bonds = bonds

    angles = IMP.bff.FFAngleVector()
    for row in d.get("angles") or []:
        a, b, c, theta, tid = row
        r = IMP.bff.FFAngle()
        r.site_a, r.site_b, r.site_c = _token(a), _token(b), _token(c)
        r.theta = float(theta or 0.0)
        r.type_id = str(tid or "")
        angles.append(r)
    sys_.angles = angles

    def _torsions(key):
        v = IMP.bff.FFTorsionVector()
        for row in d.get(key) or []:
            a, b, c, dd, tid = row
            r = IMP.bff.FFTorsion()
            r.site_a, r.site_b = _token(a), _token(b)
            r.site_c, r.site_d = _token(c), _token(dd)
            r.type_id = str(tid or "")
            v.append(r)
        return v
    sys_.dihedrals = _torsions("dihedrals")
    sys_.impropers = _torsions("impropers")

    probes = IMP.bff.FFProbeVector()
    for row in d.get("probes") or []:
        pr = IMP.bff.FFProbe()
        pr.id = int(row.get("id") or 0)
        pr.name = str(row.get("name") or "")
        pr.origin = str(row.get("origin") or "extrinsic")
        pr.link_type = str(row.get("link_type") or "covalent")
        probes.append(pr)
    sys_.probes = probes

    nb = IMP.bff.FFNonbonded()
    src = d.get("nonbonded") or {}
    nb.enabled = bool(src.get("enabled", True))
    nb.k = float(src.get("k", 5.0))
    nb.cutoff = float(src.get("cutoff_A", 6.0))
    sys_.nonbonded = nb

    sp = IMP.bff.FFSampling()
    src = d.get("sampling") or {}
    sp.temperature_K = float(src.get("temperature_K", 300.0))
    sp.friction_ps = float(src.get("friction_ps", 10.0))
    sp.timestep_fs = float(src.get("timestep_fs", 0.25))
    sp.n_steps = int(src.get("n_steps", 500000))
    sp.write_every = int(src.get("write_every", 1000))
    sp.minimize_steps = int(src.get("minimize_steps", 200))
    sys_.sampling = sp
    return sys_
