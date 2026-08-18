#!/usr/bin/env python
"""Read/write NMR restraints in official mmCIF categories."""

import json
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


def _as_float(v, h):
    s = _as_str(v, h)
    return None if s is None else float(s)


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


def read_nmr_restraints(path, fixed_component_name=None):
    """Read NMR restraints from mmCIF (or legacy JSON).

    Parameters
    ----------
    fixed_component_name : str, optional
        The residue/component name that identifies the fixed component (e.g. the ``comp_id``
        written in the ``_ihm_non_poly_feature`` loop).  When provided, features
        whose ``comp_id`` matches this string receive ``component_id = "fixed"``
        and all others receive ``component_id = "mobile"``.  When *not* provided
        the ``component_id`` field in each feature is set to the raw ``comp_id``
        string so the caller can classify features itself.
    """
    if not path or not os.path.exists(path):
        return {
            "default_restraint_set": "",
            "Distances": {},
            "Positions": {},
            "restraint_sets": {},
        }

    if str(path).lower().endswith(".json"):
        with open(path) as fh:
            return json.load(fh)

    data = {
        "default_restraint_set": "",
        "Distances": {},
        "Positions": {},
        "restraint_sets": {},
    }
    feature_name = {}
    feature_pos = {}
    rows = []
    group_meta = {}

    class FeatureH(_BaseHandler):
        def __call__(self, feature_id, feature_type, entity_type, details):
            fid = _as_int(feature_id, self)
            if fid is None:
                return
            feature_name[fid] = _as_str(details, self) or f"feature_{fid}"

    class NonPolyH(_BaseHandler):
        def __call__(
            self, ordinal_id, feature_id, entity_id, asym_id, comp_id, atom_id
        ):
            fid = _as_int(feature_id, self)
            if fid is None:
                return
            comp = _as_str(comp_id, self)
            if fixed_component_name is not None:
                component_id = "fixed" if comp == fixed_component_name else "mobile"
            else:
                component_id = comp or "guest"
            feature_pos[fid] = {
                "simulation_type": "Atom",
                "atom_name": _as_str(atom_id, self),
                "residue_name": comp,
                "residue_seq_number": 0,
                "chain_identifier": _as_str(asym_id, self) or "",
                "component_id": component_id,
            }

    class DerivedH(_BaseHandler):
        def __call__(
            self,
            id,
            group_id,
            dataset_list_id,
            feature_id_1,
            feature_id_2,
            restraint_type,
            group_conditionality,
            probability,
            mic_value,
            distance_lower_limit,
            distance_upper_limit,
        ):
            rid = _as_int(id, self)
            if rid is None:
                return
            rows.append(
                {
                    "id": rid,
                    "group_id": _as_int(group_id, self),
                    "f1": _as_int(feature_id_1, self),
                    "f2": _as_int(feature_id_2, self),
                    "rtype": _as_str(restraint_type, self),
                    "lower": _as_float(distance_lower_limit, self),
                    "upper": _as_float(distance_upper_limit, self),
                }
            )

    class ConstraintFileH(_BaseHandler):
        def __call__(
            self,
            entry_id,
            constraint_filename,
            constraint_type,
            constraint_subtype,
            constraint_number,
        ):
            gid = _as_int(entry_id, self)
            if gid is None:
                return
            fname = _as_str(constraint_filename, self) or ""
            set_key = fname[:-4] if fname.endswith(".tbl") else fname
            group_meta[gid] = {
                "set_key": set_key,
                "subtype": _as_str(constraint_subtype, self) or "NOE",
            }

    handlers = {
        "_ihm_feature_list": FeatureH(),
        "_ihm_non_poly_feature": NonPolyH(),
        "_ihm_derived_distance_restraint": DerivedH(),
        "_pdbx_nmr_constraint_file": ConstraintFileH(),
    }

    with open(path) as fh:
        r = ihm.format.CifReader(fh, handlers)
        r.read_file()

    for fid, pos in sorted(feature_pos.items()):
        name = feature_name.get(fid, f"feature_{fid}")
        data["Positions"][name] = pos

    used_prefixes = set()
    for gid, meta in sorted(group_meta.items()):
        set_key = meta.get("set_key") or f"group_{gid}_prefix"
        subtype = meta.get("subtype")
        prefix = _prefix_for_set_key(set_key, subtype, used_prefixes)
        used_prefixes.add(prefix)
        data["restraint_sets"][set_key] = prefix

    if "common_prefix" not in data["restraint_sets"]:
        data["restraint_sets"]["common_prefix"] = ""
    if "subunit_prefix" not in data["restraint_sets"]:
        data["restraint_sets"]["subunit_prefix"] = "SUB_"
    if "global_prefix" not in data["restraint_sets"]:
        data["restraint_sets"]["global_prefix"] = "GLOB_"
    if "single_s1_prefix" not in data["restraint_sets"]:
        data["restraint_sets"]["single_s1_prefix"] = "S1ONLY_"
    if "single_s1_common_prefix" not in data["restraint_sets"]:
        data["restraint_sets"]["single_s1_common_prefix"] = ""
    data["default_restraint_set"] = "subunit"

    gid_to_prefix = {}
    for set_key, prefix in data["restraint_sets"].items():
        if not set_key.endswith("_prefix"):
            continue
        for gid, meta in group_meta.items():
            if meta.get("set_key") == set_key:
                gid_to_prefix[gid] = prefix

    for row in sorted(rows, key=lambda x: x["id"]):
        p1 = feature_name.get(row["f1"], f"feature_{row['f1']}")
        p2 = feature_name.get(row["f2"], f"feature_{row['f2']}")
        rtype = _ihm_to_distance_type(row.get("rtype"))
        lower = row.get("lower")
        upper = row.get("upper")
        if rtype == "AtomUpperBound":
            target = lower if lower is not None else upper if upper is not None else 0.0
            err_pos = 1.0 if upper is None else max(1e-6, float(upper) - float(target))
            err_neg = 1.0
        elif rtype == "AtomLowerBound":
            target = upper if upper is not None else lower if lower is not None else 0.0
            err_neg = 1.0 if lower is None else max(1e-6, float(target) - float(lower))
            err_pos = 1.0
        else:
            if lower is None and upper is None:
                target = 0.0
                err_pos = 1.0
                err_neg = 1.0
            else:
                lo = float(lower if lower is not None else upper)
                up = float(upper if upper is not None else lower)
                target = 0.5 * (lo + up)
                err_pos = max(1e-6, up - target)
                err_neg = max(1e-6, target - lo)

        prefix = gid_to_prefix.get(row.get("group_id"), "")
        dname = f"{prefix}RID{row['id']}" if prefix else f"RID{row['id']}"
        data["Distances"][dname] = {
            "distance": float(target),
            "error_pos": float(err_pos),
            "error_neg": float(err_neg),
            "position1_name": p1,
            "position2_name": p2,
            "distance_type": rtype,
        }

    return data


def write_nmr_restraints(path, data):
    """Write NMR restraints to an mmCIF file using official categories."""
    restraint_sets = dict(data.get("restraint_sets", {}))
    distances = dict(data.get("Distances", {}))
    positions = dict(data.get("Positions", {}))

    key_prefix_pairs = _group_prefix_pairs(restraint_sets)
    if not key_prefix_pairs:
        key_prefix_pairs = [
            ("subunit_prefix", "SUB_"),
            ("global_prefix", "GLOB_"),
            ("common_prefix", "CSP_"),
            ("single_s1_prefix", "S1ONLY_"),
        ]

    group_key_to_id = {}
    for i, (key, _prefix) in enumerate(
        sorted(key_prefix_pairs, key=lambda kv: _set_order(kv[0]))
    ):
        group_key_to_id[key] = i + 1

    distance_group_key = {}
    for dname in distances:
        gk = _infer_group_for_name(dname, key_prefix_pairs)
        if gk is None:
            gk = key_prefix_pairs[0][0]
        distance_group_key[dname] = gk

    with open(path, "w") as out:
        w = ihm.format.CifWriter(out)
        w.start_block("nmr_restraints")

        feature_id_by_name = {
            pname: i + 1 for i, (pname, _pos) in enumerate(sorted(positions.items()))
        }

        with w.loop(
            "_ihm_feature_list",
            ["feature_id", "feature_type", "entity_type", "details"],
        ) as l:
            for pname, _ in sorted(positions.items()):
                l.write(
                    feature_id=feature_id_by_name[pname],
                    feature_type="atom",
                    entity_type=None,
                    details=pname,
                )

        entity_id_by_component = {}
        next_entity_id = 1
        with w.loop(
            "_ihm_non_poly_feature",
            ["ordinal_id", "feature_id", "entity_id", "asym_id", "comp_id", "atom_id"],
        ) as l:
            ordinal = 1
            for pname, pos in sorted(positions.items()):
                chain = pos.get("chain_identifier", "")
                asym_id = chain if chain else "B"
                comp_id = pos.get("residue_name") or (
                    "DYE" if asym_id == "B" else "UNK"
                )
                component = pos.get("component_id")
                if component is None:
                    component = "guest"
                if component not in entity_id_by_component:
                    entity_id_by_component[component] = next_entity_id
                    next_entity_id += 1
                l.write(
                    ordinal_id=ordinal,
                    feature_id=feature_id_by_name[pname],
                    entity_id=entity_id_by_component[component],
                    asym_id=asym_id,
                    comp_id=comp_id,
                    atom_id=pos.get("atom_name"),
                )
                ordinal += 1

        with w.loop(
            "_ihm_derived_distance_restraint",
            [
                "id",
                "group_id",
                "dataset_list_id",
                "feature_id_1",
                "feature_id_2",
                "restraint_type",
                "group_conditionality",
                "probability",
                "mic_value",
                "distance_lower_limit",
                "distance_upper_limit",
            ],
        ) as l:
            rid = 1
            for dname, d in sorted(distances.items()):
                p1 = d.get("position1_name")
                p2 = d.get("position2_name")
                if p1 not in feature_id_by_name or p2 not in feature_id_by_name:
                    continue
                target = float(d.get("distance", 0.0))
                err_pos = (
                    float(d.get("error_pos", 1.0))
                    if d.get("error_pos") is not None
                    else 1.0
                )
                err_neg = (
                    float(d.get("error_neg", 1.0))
                    if d.get("error_neg") is not None
                    else 1.0
                )
                dtype = d.get("distance_type", "AtomUpperBound")
                rtype = _distance_type_to_ihm(dtype)

                if rtype == "upper bound":
                    lower = target
                    upper = target + err_pos
                elif rtype == "lower bound":
                    lower = target - err_neg
                    upper = target
                else:
                    lower = target - err_neg
                    upper = target + err_pos

                gkey = distance_group_key.get(dname, key_prefix_pairs[0][0])
                l.write(
                    id=rid,
                    group_id=group_key_to_id[gkey],
                    dataset_list_id=None,
                    feature_id_1=feature_id_by_name[p1],
                    feature_id_2=feature_id_by_name[p2],
                    restraint_type=rtype,
                    group_conditionality="ALL",
                    probability=1.0,
                    mic_value=None,
                    distance_lower_limit=lower,
                    distance_upper_limit=upper,
                )
                rid += 1

        counts = {k: 0 for k in group_key_to_id}
        for dname in distances:
            counts[distance_group_key[dname]] += 1

        with w.loop(
            "_pdbx_nmr_constraint_file",
            [
                "entry_id",
                "constraint_filename",
                "constraint_type",
                "constraint_subtype",
                "constraint_number",
            ],
        ) as l:
            for set_key, gid in sorted(group_key_to_id.items(), key=lambda kv: kv[1]):
                l.write(
                    entry_id=gid,
                    constraint_filename=f"{set_key}.tbl",
                    constraint_type="distance",
                    constraint_subtype=_constraint_subtype(set_key),
                    constraint_number=counts.get(set_key, 0),
                )


def convert_nmr_json_to_cif(json_path, cif_path):
    """Convert a legacy JSON restraints file to mmCIF."""
    with open(json_path) as fh:
        data = json.load(fh)

    for _pname, pos in data.get("Positions", {}).items():
        if "component_id" not in pos:
            pos["component_id"] = "guest"

    write_nmr_restraints(cif_path, data)
