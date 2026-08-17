#!/usr/bin/env python
"""Read/write cgdye feature templates in mmCIF format."""

import os

import ihm.format


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


def read_component_template_cif(path):
    """Read a component feature template (regions/features per component) from mmCIF.

    The generic cgdye template: named regions and per-atom features used by the
    topology builder and the density analysis. Dye-specific templates with
    dipole/charge metadata are read by :func:`read_dye_template_cif`.
    """
    template = {
        "name": os.path.splitext(os.path.basename(path))[0],
        "features": {},
        "impropers": [],
    }

    class TemplateH(_BaseHandler):
        def __call__(self, name):
            n = _as_str(name, self)
            if n:
                template["name"] = n

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


def flex_features(template):
    """Return list of feature IDs that are DOF features with no rb/md_fixed (i.e., flex)."""
    out = []
    for fid, spec in template.get("features", {}).items():
        if (
            spec.get("feature_type") == "dof"
            and not spec.get("rb")
            and not spec.get("md_fixed")
        ):
            out.append(fid)
    return out


# ============================================================================
# Dye-specific template handlers (for IMP rotamer library system)
# ============================================================================


def read_dye_template_cif(path):
    """Read a dye template with additional metadata fields.

    Returns dict with: name, features, impropers, center_atom, dipole_atom_1,
                       dipole_atom_2, positive_atoms, negative_atoms
    """
    from collections import defaultdict

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
        "_cgdye_metadata": MetadataH(),
        "_cgdye_feature": FeatureH(),
        "_cgdye_feature_atom": FeatureAtomH(),
        "_cgdye_improper": ImproperH(),
    }

    with open(path) as fh:
        r = ihm.format.CifReader(fh, handlers)
        r.read_file()

    return template


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


def get_dye_metadata(template):
    """Extract dye-specific metadata from template.

    Returns dict with center_atom, dipole_atoms (tuple), positive_atoms, negative_atoms.
    """
    return {
        "center_atom": template.get("center_atom"),
        "dipole_atoms": (template.get("dipole_atom_1"), template.get("dipole_atom_2")),
        "positive_atoms": template.get("positive_atoms", []),
        "negative_atoms": template.get("negative_atoms", []),
    }
