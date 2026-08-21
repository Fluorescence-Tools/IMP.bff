/*
 * io/cif.py to C++: the FF writer, template CIF, rotamer library IO, and
 * string utilities. The FF reader is already in ForceFieldCIF.h.
 *
 * The writer writes mmCIF text directly (the ihm C library is read-only).
 * The template reader uses the ihm C reader (same as ForceFieldCIF.cpp).
 * The rotamer library IO reads/writes numpy .npy + text files directly.
 *
 * forcefield_system_from_dict and as_forcefield_system stay as %pythoncode
 * because they bridge Python dicts to the typed C++ DyeForceFieldSystem --
 * the dict is a Python literal convenience, not a kernel.
 */

%include "IMP/bff/CifIO.h"

// Only the FF writer and rotamer library IO are C++ here. The template CIF
// reader/writer stays as %pythoncode because it returns dicts (consumers use
// template["impropers"], template["features"]["f1"]["rb"]) and uses
// ihm.format.CifReader/CifWriter (Python-only, deferred import).

IMP_SWIG_VALUE(IMP::bff, RotamerLibraryData, RotamerLibraryDatas);

%pythoncode %{
import os as _os
import re as _re
import numpy as _np


def read_dye_forcefield_cif(path):
    """Read a force-field system from mmCIF."""
    return IMP.bff.read_forcefield_cif(str(path))


def as_forcefield_system(system):
    """A system, whichever way it was given."""
    if isinstance(system, dict):
        return forcefield_system_from_dict(system)
    return system


def forcefield_system_from_dict(d):
    """A parsed or built dictionary as a typed DyeForceFieldSystem."""
    sys_ = IMP.bff.DyeForceFieldSystem(str(d.get("name") or ""))

    comps = IMP.bff.FFComponentMap()
    for cid, spec in (d.get("components") or {}).items():
        c = IMP.bff.FFComponent()
        c.mol2_path = str(spec.get("mol2_path") or spec.get("mol2") or "")
        c.pdb_path = str(spec.get("pdb_path") or spec.get("pdb") or "")
        c.role = str(spec.get("role") or "")
        comps[str(cid)] = c
    sys_.components = comps

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


def write_dye_forcefield_cif(path, system):
    """Write a force-field system to mmCIF. Accepts DyeForceFieldSystem or dict."""
    system = as_forcefield_system(system)
    IMP.bff._write_dye_forcefield_cif(str(path), system)


# -- Template CIF reader/writer: stays Python (returns dicts, uses ihm.format) --

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


def _as_bool(v, h):
    s = _as_str(v, h)
    if s is None:
        return None
    return s.upper() in {"1", "YES", "TRUE"}


def read_component_template_cif(path, with_dye_metadata=False):
    """Read a cgdye component template from mmCIF."""
    import ihm.format

    template = {
        "name": _os.path.splitext(_os.path.basename(path))[0],
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
                {"rb": False, "md_fixed": False, "feature_type": "dof",
                 "region_color": None, "atoms": []})
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
                {"rb": False, "md_fixed": False, "feature_type": "dof",
                 "region_color": None, "atoms": []})
            template["features"][fid]["atoms"].append(
                {"name": anm, "occurrence": 1 if occ is None else occ})

    class ImproperH(_BaseHandler):
        def __call__(self, center_atom, type):
            ca = _as_str(center_atom, self)
            it = _as_str(type, self)
            if not ca or not it:
                return
            template["impropers"].append({"center_atom": ca, "type": it})

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
    import ihm.format

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


def _write_header(w, template):
    with w.category("_cgdye_template") as c:
        c.write(name=template.get("name", "cgdye_template"))


def _write_features(w, template):
    with w.loop("_cgdye_feature",
                ["feature_id", "feature_type", "rb", "md_fixed", "region_color"]) as l:
        for fid, spec in template.get("features", {}).items():
            l.write(feature_id=fid, feature_type=spec.get("feature_type", "dof"),
                    rb="YES" if spec.get("rb", False) else "NO",
                    md_fixed="YES" if spec.get("md_fixed", False) else "NO",
                    region_color=spec.get("region_color") or ".")


def _write_feature_atoms(w, template):
    with w.loop("_cgdye_feature_atom", ["feature_id", "atom_name", "occurrence"]) as l:
        for fid, spec in template.get("features", {}).items():
            for a in spec.get("atoms", []):
                l.write(feature_id=fid, atom_name=a.get("name"),
                        occurrence=int(a.get("occurrence", 1)))


def _write_impropers(w, template):
    if template.get("impropers"):
        with w.loop("_cgdye_improper", ["center_atom", "type"]) as l:
            for imp in template["impropers"]:
                l.write(center_atom=imp.get("center_atom"), type=imp.get("type"))


def write_component_template_cif(path, template):
    """Write a component feature template to mmCIF."""
    def _write(w):
        _write_header(w, template)
        _write_features(w, template)
        _write_feature_atoms(w, template)
        _write_impropers(w, template)
    _write_cif_safe(path, template.get("name", "cgdye_template"), _write)


def read_dye_template_cif(path):
    """A component template read with its dye metadata."""
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


def region_features(template):
    """Return dict of {feature_id: region_color} for features with non-null region_color."""
    out = {}
    for fid, spec in template.get("features", {}).items():
        rc = spec.get("region_color")
        if rc:
            out[fid] = rc
    return out


# -- Rotamer library IO: stays Python (tests pass dicts, not C++ structs) ---

def _base_path(rmf_or_npy_path):
    if rmf_or_npy_path.endswith(".rmf"):
        return rmf_or_npy_path[:-4]
    if rmf_or_npy_path.endswith(".npy"):
        return rmf_or_npy_path[:-4]
    return rmf_or_npy_path


def read_rotamer_library(path):
    """Read a rotamer library from numpy/text files."""
    base = _base_path(path)
    coords = _np.load(base + "_coords.npy")
    n_rotamers = coords.shape[0]
    weights_file = base + "_weights.txt"
    if _os.path.exists(weights_file):
        with open(weights_file) as fh:
            weights = [float(line.strip()) for line in fh if line.strip()]
    else:
        weights = [1.0 / n_rotamers] * n_rotamers
    atoms_file = base + "_atoms.txt"
    if _os.path.exists(atoms_file):
        with open(atoms_file) as fh:
            atom_names = [line.strip() for line in fh if line.strip()]
    else:
        atom_names = [f"AT{i}" for i in range(coords.shape[1])]
    return {
        "id": list(range(1, n_rotamers + 1)),
        "weight": weights,
        "atom_names": atom_names,
        "coords": {i + 1: coords[i] for i in range(n_rotamers)},
    }


def write_rotamer_library(path, library):
    """Write a rotamer library to numpy/text files."""
    base = _base_path(path)
    atom_names = library["atom_names"]
    n_atoms = len(atom_names)
    ids = sorted(library["coords"].keys())
    coord_array = _np.zeros((len(ids), n_atoms, 3))
    for i, rid in enumerate(ids):
        c = library["coords"][rid]
        arr = _np.asarray(c)
        if arr.ndim == 1 and arr.shape[0] == 3:
            arr = arr.reshape(1, 3)
        coord_array[i, : arr.shape[0]] = arr[:n_atoms]
    _np.save(base + "_coords.npy", coord_array)
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


# -- Atom-site row parsers (for _atom_site in FF writer) -------------------

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
            atype = parts[5]
            comp_id = parts[7] if len(parts) > 7 else "UNK"
            element = atype.split(".")[0].capitalize()
            if not element:
                m = _re.match(r"[A-Za-z]+", atom_name)
                element = m.group(0)[0].upper() if m else "C"
            rows.append({
                "group_PDB": "HETATM", "id": serial, "type_symbol": element,
                "label_atom_id": atom_name, "label_comp_id": comp_id,
                "label_asym_id": asym_id, "label_entity_id": int(entity_id),
                "label_seq_id": 1, "Cartn_x": x, "Cartn_y": y, "Cartn_z": z,
                "occupancy": 1.0, "B_iso_or_equiv": 0.0,
                "auth_asym_id": asym_id, "auth_comp_id": comp_id,
                "auth_seq_id": 1, "auth_atom_id": atom_name, "pdbx_PDB_model_num": 1,
            })
    return rows


def _parse_pdb_atom_site_rows(pdb_path, asym_id, entity_id):
    """Parse ATOM/HETATM records from a PDB file into _atom_site rows."""
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
                m = _re.match(r"([A-Za-z]+)", atom_name)
                element = (m.group(1)[0] if m else atom_name[:1]).upper()
            if not serial_txt or not x_txt or not y_txt or not z_txt:
                continue
            serial = int(serial_txt)
            seq_id = int(seq_txt) if seq_txt else 1
            rows.append({
                "group_PDB": "ATOM" if line.startswith("ATOM") else "HETATM",
                "id": serial, "type_symbol": element,
                "label_atom_id": atom_name, "label_comp_id": comp_id,
                "label_asym_id": asym_id, "label_entity_id": int(entity_id),
                "label_seq_id": seq_id,
                "Cartn_x": float(x_txt), "Cartn_y": float(y_txt), "Cartn_z": float(z_txt),
                "occupancy": float(occ_txt) if occ_txt else 1.0,
                "B_iso_or_equiv": float(b_txt) if b_txt else 0.0,
                "auth_asym_id": auth_asym, "auth_comp_id": comp_id,
                "auth_seq_id": seq_id, "auth_atom_id": atom_name,
            })
    return rows


def _parse_struct_atom_site_rows(struct_path, asym_id, entity_id):
    """Dispatch to MOL2 or PDB parser based on file extension."""
    if struct_path.endswith(".mol2"):
        return _parse_mol2_atom_site_rows(struct_path, asym_id, entity_id)
    return _parse_pdb_atom_site_rows(struct_path, asym_id, entity_id)
%}
