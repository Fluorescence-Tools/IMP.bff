"""The legacy C# FPS ``.txt`` formats — **read only**.

A format nobody can still read is data that has been lost, and a decade of
measurements live in these files. So they are readable, and deliberately not
writable: nothing should produce another one.

Split out of ``IMP.bff.fret.io`` by PRD-113 stage 7.
"""

from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Tuple
import json
import os

import numpy as np

# (was: import fps_schema) -- now in this module

__all__ = [
    'AV_SIMULATION_TYPES',
    'DISTANCE_FIELDS',
    'DISTANCE_TYPES',
    'POSITION_FIELDS',
    'SCHEMA_VERSION',
    'SCORE_SET_FIELDS',
    'SIMULATION_TYPES',
    'fps_positions_for_docking',
    'read_evaluators_json',
    'read_fps_json',
    'read_old_distances_txt',
    'read_old_lps_txt',
    'to_json_schema',
    'validate',
    'validate_distance',
    'validate_position',
    'write_evaluators_json',
    'write_fps_json',
]

# --------------------------------------------------------------------------
# fps_legacy
# --------------------------------------------------------------------------
"""The legacy C# FPS ``.txt`` formats — **read only**.

A format nobody can still read is data that has been lost, and a decade of
measurements live in these files. So they are readable, and deliberately not
writable: nothing should produce another one.

Split out of ``IMP.bff.fret.io`` by PRD-113 stage 7.
"""

def read_old_lps_txt(
    path: str | os.PathLike,
    pdb_paths: Optional[List[str]] = None,
) -> Tuple[Dict, List[str]]:
    """Read a C# FPS format labeling positions (.txt) file.

    Parameters
    ----------
    path : str or PathLike
        Path to the labeling positions text file.
    pdb_paths : list of str, optional
        Paths to PDB files to resolve atom IDs to chain/residue info.

    Returns
    -------
    positions : dict
        Standard positions dict.
    molecules : list of str
        The unique list of molecule names in order of appearance.
    """
    if pdb_paths is None:
        pdb_paths = []

    if not pdb_paths:
        dir_name = os.path.dirname(str(path))
        try:
            if dir_name:
                pdb_paths = [os.path.join(dir_name, f) for f in os.listdir(dir_name) if f.endswith(".pdb")]
            else:
                pdb_paths = [f for f in os.listdir(".") if f.endswith(".pdb")]
        except Exception:
            pass

    # Map molecule name (basename without ext) to its full PDB path
    pdb_map = {}
    for p in pdb_paths:
        base = os.path.splitext(os.path.basename(p))[0].lower()
        pdb_map[base] = p

    positions = {}
    molecules = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue

            name = parts[0]
            mol = parts[1]
            dye = parts[2]
            av_type = parts[3]

            if mol not in molecules:
                molecules.append(mol)
            body_id = molecules.index(mol)

            # Resolve atom ID to chain, residue sequence number, and atom name
            chain = ""
            resseq = 1
            aname = "CA"
            atom_id = None

            # Look up PDB file if available
            pdb_file = pdb_map.get(mol.lower())

            # Parse parameters
            r2 = 0.0
            r3 = 0.0
            if av_type == "AV1" and len(parts) >= 8:
                ll = float(parts[4])
                lw = float(parts[5])
                r1 = float(parts[6])
                atom_id = int(parts[7])
            elif av_type == "AV3" and len(parts) >= 10:
                ll = float(parts[4])
                lw = float(parts[5])
                r1 = float(parts[6])
                r2 = float(parts[7])
                r3 = float(parts[8])
                atom_id = int(parts[9])
            elif av_type == "XYZ" and len(parts) >= 7:
                # Fixed coordinates
                x = float(parts[4])
                y = float(parts[5])
                z = float(parts[6])
                positions[name] = {
                    "simulation_type": "XYZ",
                    "x": x,
                    "y": y,
                    "z": z,
                    "body_id": body_id,
                }
                continue
            else:
                continue

            if atom_id is not None and pdb_file and os.path.exists(pdb_file):
                try:
                    with open(pdb_file, "r") as pf:
                        for pline in pf:
                            if pline.startswith(("ATOM  ", "HETATM")):
                                try:
                                    cur_id = int(pline[6:11].strip())
                                except ValueError:
                                    continue
                                if cur_id == atom_id:
                                    chain = pline[21].strip()
                                    try:
                                        resseq = int(pline[22:26].strip())
                                    except ValueError:
                                        resseq = 1
                                    aname = pline[12:16].strip()
                                    break
                except Exception:
                    pass
            elif atom_id is not None:
                # Fallback proxy if PDB not parsed
                resseq = atom_id

            positions[name] = {
                "chain_identifier": chain,
                "residue_seq_number": resseq,
                "atom_name": aname,
                "linker_length": ll,
                "linker_width": lw,
                "radius1": r1,
                "radius2": r2,
                "radius3": r3,
                "simulation_grid_resolution": 1.5,
                "simulation_type": av_type,
                "body_id": body_id,
            }

    return positions, molecules


def read_old_distances_txt(
    path: str | os.PathLike,
) -> Dict:
    """Read a C# FPS format experimental distances (.txt) file.

    Parameters
    ----------
    path : str or PathLike

    Returns
    -------
    distances : dict
        Standard distances dict.
    """
    distances = {}
    distance_type = "RDAMean"

    with open(path, "r") as f:
        # Check first line for distance type
        first_line = f.readline().strip()
        if first_line and not first_line.split()[0].isalnum():
            # e.g. "RDAMeanE"
            pass
        elif first_line:
            parts = first_line.split()
            if len(parts) == 1:
                distance_type = parts[0]
            else:
                # Re-wind or process as first distance
                f.seek(0)

        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue

            pos1 = parts[0]
            pos2 = parts[1]
            dist = float(parts[2])
            err_neg = float(parts[3])
            err_pos = float(parts[4])
            forster = float(parts[5]) if len(parts) >= 6 else 52.0

            dname = f"{pos1}_{pos2}"
            distances[dname] = {
                "position1_name": pos1,
                "position2_name": pos2,
                "distance": dist,
                "error_neg": err_neg,
                "error_pos": err_pos,
                "distance_type": distance_type,
                "Forster_radius": forster,
            }

    return distances


# --------------------------------------------------------------------------
# fps
# --------------------------------------------------------------------------
"""fps.json — the one reader and writer.

fps.json describes labelling positions and distance measurements on a
structure. It had three readers and no definition until PRD-97; the definition
is :mod:`IMP.bff.io.fps_schema` and this is the only code that reads or writes
the format. Files are checked against the schema rather than against another
parser.

Moved from ``IMP.bff.fret.io`` by PRD-113 stage 7, and split there: the legacy
C# text formats are in :mod:`IMP.bff.io.fps_legacy` and the structure readers in
:mod:`IMP.bff.io.structure`. One file held all three, which is how "the fps
reader" came to own PDB writing.
"""

#: The ``simulation_type`` values that describe an accessible volume, as opposed
#: to a fixed point. Lives here rather than in the schema because it is what the
#: *readers* accept, and the schema records what a file may say.
AV_SIMULATION_TYPES = ("AV1", "AV3", "XYZ")

# ---------------------------------------------------------------------------
# fps.json reader / writer
# ---------------------------------------------------------------------------


def read_fps_json(
    path: str | os.PathLike,
    pdb_paths: Optional[List[str]] = None,
    validate: bool = False,
) -> Tuple[Dict, Dict, Dict, Dict]:
    """Load an fps.json labeling file, or old C# .txt format files.

    If path is a .txt file, it reads positions from path and looks for
    a Distances.txt in the same directory.

    Parameters
    ----------
    path : str or os.PathLike
        Path to an fps.json (or legacy C# FPS ``.txt``) file.
    pdb_paths : list of str, optional
        PDB paths used to resolve atom ids in the legacy ``.txt`` format.
    validate : bool
        Check the payload against the fps.json schema
        (:func:`IMP.bff.io.fps_schema.validate`) and raise ``ValueError``
        listing the violations when it does not conform.

    Returns
    -------
    positions : dict
        {name: {param: value}}  — may be empty.
    distances : dict
        {name: {param: value}}  — may be empty.
    score_sets : dict
        {name: {"distances": [...], ...}} — from the ``"χ²"`` key.
    extra : dict
        Any top-level keys not in the above (e.g. ``Evaluators``).
    """
    path_str = str(path)
    if not path_str.endswith(".json"):
        # Old C# txt format
        positions, molecules = read_old_lps_txt(path_str, pdb_paths=pdb_paths)

        # Look for Distances.txt in the same directory
        dir_name = os.path.dirname(path_str)
        dist_path = os.path.join(dir_name, "Distances.txt")
        if os.path.exists(dist_path):
            distances = read_old_distances_txt(dist_path)
        else:
            distances = {}
        return positions, distances, {}, {}

    with open(path) as f:
        payload = json.load(f)

    if validate:
        errors, _warnings = validate_fps(payload)
        if errors:
            raise ValueError(
                f"{path_str} does not conform to the fps.json schema:\n  "
                + "\n  ".join(errors))

    positions = payload.pop("Positions", {})
    distances = payload.pop("Distances", {})
    score_sets = payload.pop("χ²", {})
    return positions, distances, score_sets, payload


def write_fps_json(
    path: str | os.PathLike,
    positions: Dict,
    distances: Dict,
    score_sets: Optional[Dict] = None,
    extra: Optional[Dict] = None,
    validate: bool = False,
    **kwargs,
) -> None:
    """Write an fps.json file.

    Parameters
    ----------
    path : str or os.PathLike
        Output JSON path.
    positions : dict
        Position entries.
    distances : dict
        Distance entries.
    score_sets : dict, optional
        Optional chi-square score sets (the ``"χ²"`` section).
    extra : dict, optional
        Optional extra top-level payload.
    validate : bool
        Check the assembled payload against the fps.json schema before
        writing and raise ``ValueError`` when it does not conform.
    **kwargs
        Additional ``json.dump`` options.
    """
    payload: Dict = {}
    if extra is not None:
        payload.update(extra)
    payload["Positions"] = positions
    payload["Distances"] = distances
    if score_sets:
        payload["χ²"] = score_sets
    if validate:
        errors, _warnings = validate_fps(payload)
        if errors:
            raise ValueError(
                "refusing to write a non-conforming fps.json:\n  "
                + "\n  ".join(errors))
    with open(path, "w") as f:
        json.dump(payload, f, indent=2, **kwargs)


def fps_positions_for_docking(
    positions: Dict,
    distances: Optional[Dict] = None,
) -> Tuple[Dict, Dict]:
    """Keep only the positions (and distances between them) the C++ AV scorer understands.

    Rotamer-ensemble positions (``simulation_type == "R1"``, PRD-108) are
    Python-only; ``IMP::bff::AVNetworkRestraint`` would score them as AV1
    with default parameters. Returns ``(positions, distances)`` restricted to
    AV1/AV3/XYZ positions and to distances whose two ends survive.
    """
    kept = {
        name: pos for name, pos in positions.items()
        if str((pos or {}).get("simulation_type", "AV1")) in AV_SIMULATION_TYPES
    }
    if distances is None:
        return kept, {}
    kept_d = {
        name: d for name, d in distances.items()
        if d.get("position1_name") in kept and d.get("position2_name") in kept
    }
    return kept, kept_d


def write_evaluators_json(path: str | os.PathLike, evaluators: List[Any]) -> None:
    """Append or overwrite the ``Evaluators`` key of an fps.json file.

    Parameters
    ----------
    path : str or PathLike
        Path to the json file.
    evaluators : list
        Evaluator objects with a ``to_dict()`` method, or plain dicts.
    """
    if os.path.exists(path):
        with open(path) as f:
            try:
                payload = json.load(f)
            except Exception:
                payload = {}
    else:
        payload = {}

    payload["Evaluators"] = [
        ev.to_dict() if hasattr(ev, "to_dict") else dict(ev)
        for ev in evaluators
    ]
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)


def read_evaluators_json(
    path: str | os.PathLike,
    factory: Optional[Callable[[Dict], Any]] = None,
) -> List[Any]:
    """Read the ``Evaluators`` list from an fps.json file.

    Parameters
    ----------
    path : str or PathLike
        Path to the json file.
    factory : callable, optional
        Called on each evaluator dict to instantiate an application object
        (e.g. ChiSurf passes its ``evaluators.from_dict``). Entries the
        factory raises on are skipped. Without a factory the raw dicts are
        returned — this module defines the *format*; what an evaluator *is*
        belongs to the application.

    Returns
    -------
    list
        Instantiated evaluators, or raw dicts when no factory is given.
    """
    if not os.path.exists(path):
        return []
    with open(path) as f:
        try:
            payload = json.load(f)
        except Exception:
            return []

    ev_list = payload.get("Evaluators", [])
    if factory is None:
        return [dict(d) for d in ev_list if isinstance(d, dict)]
    res = []
    for d in ev_list:
        try:
            res.append(factory(d))
        except Exception:
            pass
    return res


# --------------------------------------------------------------------------
# fps_schema
# --------------------------------------------------------------------------
"""The fps.json schema — the single authored definition (PRD-97 stage 0).

Moved from ``IMP.bff.fret.fps_schema`` by PRD-113 stage 7.

fps.json describes labelling positions and distance measurements on a
structure. Historically it had three readers (the C++
``IMP::bff::AVNetworkRestraint``/``FPSReaderWriter``, a small Python reader,
and ChiSurf's elaborate one) and **no definition**: two dialects existed whose
field sets differed, and a file written by one was not fully readable by the
other. This module is the definition. It follows mmfdb's discipline — the
schema is *one authored artifact* and everything else derives from it:

* :data:`POSITION_FIELDS` / :data:`DISTANCE_FIELDS` name every field either
  dialect uses, with type, default, and — where the flrCIF dictionary defines
  the concept — the canonical flrCIF item name. Where flrCIF has no item, the
  entry is an authored local definition, visible here rather than implied by a
  parser.
* :func:`to_json_schema` derives a JSON-Schema document from the field tables;
  the shipped ``data/fps_json_schema.json`` is generated by it and a test
  holds the two in agreement (drift guardrail, like mmfdb's).
* :func:`validate` checks a parsed fps.json payload against the definition —
  readers and writers are checked against the schema, not against each other.

Top-level layout of an fps.json (the "network" dialect)::

    {
      "Positions":  {"<name>": {<position fields>}, ...},
      "Distances":  {"<name>": {<distance fields>}, ...},
      "χ²":         {"<set name>": {"distances": ["<name>", ...]}, ...},
      "Evaluators": [ {...}, ... ],          # optional, application-defined
      ...                                     # extra keys are preserved
    }

The second dialect ("flat", the template files in ``data/``) is a single
position object with no wrapper — the same position fields, validated by
:func:`validate_position`.

References: the flrCIF extension dictionary (mmCIF ``ihm_flr`` extension),
categories ``_flr_FPS_AV_parameter``, ``_flr_FPS_global_parameter``,
``_flr_FPS_mean_probe_position`` and ``_flr_fret_distance_restraint``.
"""

SCHEMA_VERSION = "1.0"

#: ``simulation_type`` values. AV1: one dye radius; AV3: three; XYZ: a fixed
#: mean position (no volume simulation), carried by legacy C# FPS conversions;
#: R1: a rotamer ensemble (PRD-108) -- a FRETpredict-style 1:1 rotamer
#: library placed in the residue's backbone frame and Boltzmann-screened
#: against the structure, per rotamer centre + dipole + weight (R2 = sampled
#: library, R3 = mixture are reserved). R1 positions are Python-only:
#: ``IMP::bff::AVNetworkRestraint`` never reads ``simulation_type`` and would
#: score such a position as an AV1 with its AV parameters (the C++ side warns);
#: filter with :func:`IMP.bff.io.fps.fps_positions_for_docking` first.
SIMULATION_TYPES = ("AV1", "AV3", "XYZ", "R1")

#: ``distance_type`` values and their flrCIF ``distance_type`` spellings.
#: RDAMean = mean inter-dye distance <R_DA>; RDAMeanE = FRET-averaged distance
#: <R_DA>_E; Rmp = distance between mean dye positions R_mp.
DISTANCE_TYPES = {
    "RDAMean": "<R_DA>",
    "RDAMeanE": "<R_DA>_E",
    "Rmp": "R_mp",
}

# ---------------------------------------------------------------------------
# Field tables
#
# Each entry: name -> {
#   "type":     python type or tuple of accepted types,
#   "default":  value a reader may assume when the key is absent
#               (None = no default; the key is then required for the dialects
#                listed in "required_for"),
#   "flrcif":   canonical flrCIF item name, or None when flrCIF does not
#               define the concept (then "authored" documents the local
#               definition),
#   "dialects": which historical dialect(s) carried the key
#               ("network" = ChiSurf Positions/Distances files,
#                "flat" = imp.bff template files, "csfps" = values produced by
#                converting legacy C# FPS .txt input),
# }
# ---------------------------------------------------------------------------

_S, _F, _I, _B = str, float, int, bool

POSITION_FIELDS: Dict[str, Dict[str, Any]] = {
    # --- labelling site --------------------------------------------------
    "chain_identifier": dict(
        type=_S, default="", flrcif=None,
        authored="Chain (asym) id of the attachment atom; empty matches any.",
        dialects=("network", "flat")),
    "residue_seq_number": dict(
        type=_I, default=0, flrcif=None,
        authored="Author residue sequence number of the attachment atom.",
        dialects=("network", "flat")),
    "residue_name": dict(
        type=_S, default="", flrcif=None,
        authored="Residue name of the attachment atom (informational).",
        dialects=("flat",)),
    "atom_name": dict(
        type=_S, default="CA", flrcif=None,
        authored="Atom name the dye linker attaches to (e.g. CA, CB).",
        dialects=("network", "flat")),
    # --- AV geometry (flrCIF FPS AV parameters) --------------------------
    "linker_length": dict(
        type=_F, default=20.0, flrcif="_flr_FPS_AV_parameter.linker_length",
        dialects=("network", "flat")),
    "linker_width": dict(
        type=_F, default=0.5, flrcif="_flr_FPS_AV_parameter.linker_width",
        dialects=("network", "flat")),
    "radius1": dict(
        type=_F, default=3.5, flrcif="_flr_FPS_AV_parameter.probe_radius_1",
        dialects=("network", "flat")),
    "radius2": dict(
        type=_F, default=0.0, flrcif="_flr_FPS_AV_parameter.probe_radius_2",
        dialects=("network", "flat")),
    "radius3": dict(
        type=_F, default=0.0, flrcif="_flr_FPS_AV_parameter.probe_radius_3",
        dialects=("network", "flat")),
    "simulation_type": dict(
        type=_S, default="AV1", flrcif=None, enum=SIMULATION_TYPES,
        authored="Label model: AV1 (one radius), AV3 (three radii), XYZ "
                 "(fixed mean position, no simulation), or R1 (rotamer "
                 "ensemble: screened rotamer library, Python-only).",
        dialects=("network", "flat", "csfps")),
    # --- rotamer-ensemble positions (R1, imp.bff dialect, PRD-108) --------
    "rotamer_library": dict(
        type=_S, default="", flrcif=None,
        authored="Rotamer library name (IMP.bff registry / FRETpredict "
                 "spelling, e.g. 'AlexaFluor 488 C1R cutoff30'). Required "
                 "for simulation_type R1.",
        dialects=("flat",)),
    "dye_name": dict(
        type=_S, default="", flrcif=None,
        authored="Dye name for R0 from spectra (e.g. 'AlexaFluor 488').",
        dialects=("flat",)),
    "temperature": dict(
        type=_F, default=298.15, flrcif=None,
        authored="Screening temperature (K) of an R1 ensemble.",
        dialects=("flat",)),
    "electrostatic": dict(
        type=_B, default=False, flrcif=None,
        authored="Add Debye-Hueckel electrostatics to the R1 screening.",
        dialects=("flat",)),
    "potential": dict(
        type=_S, default="lj", flrcif=None, enum=("lj", "gauss"),
        authored="Screening potential of an R1 ensemble: 'lj' or 'gauss'.",
        dialects=("flat",)),
    "simulation_grid_resolution": dict(
        type=_F, default=1.5, flrcif=None,
        authored="AV grid spacing in Angstrom. Related to (but not the same "
                 "item as) _flr_FPS_global_parameter.AV_min_grid_A, which is "
                 "the *minimum* grid constant of FPS's relative-grid scheme.",
        dialects=("network", "flat")),
    # --- accessible-contact-volume extensions (imp.bff dialect) ----------
    "allowed_sphere_radius": dict(
        type=_F, default=1.5,
        flrcif="_flr_FPS_global_parameter.AV_allowed_sphere",
        dialects=("flat",)),
    "contact_volume_thickness": dict(
        type=_F, default=0.0, flrcif=None,
        authored="Thickness of the contact layer above the molecular surface "
                 "for accessible-contact-volume (ACV) weighting; 0 disables.",
        dialects=("flat",)),
    "contact_volume_trapped_fraction": dict(
        type=_F, default=-1.0, flrcif=None,
        authored="Fraction of dye density trapped in the contact volume; "
                 "negative disables ACV re-weighting.",
        dialects=("flat",)),
    "min_sphere_volume_fraction": dict(
        type=_F, default=0.0, flrcif=None,
        authored="Minimum accessible fraction of the free-dye sphere volume "
                 "below which a position is flagged as buried.",
        dialects=("flat",)),
    "anchor_atoms": dict(
        type=_S, default="", flrcif=None,
        authored="Explicit anchor atom names overriding the default linker "
                 "anchor search (comma separated).",
        dialects=("flat",)),
    "strip_mask": dict(
        type=_S, default="", flrcif=None,
        authored="Selection mask of atoms removed as obstacles before the AV "
                 "simulation, in the PyMOL dialect 'chain <id> and resid <n> "
                 "and [not] name A+B+...' ('+'-separated lists, as in PyMOL "
                 "itself). Empty means the default strip: the attachment "
                 "residue's side chain minus the attachment atom.",
        dialects=("flat",)),
    "chain_weighting": dict(
        type=_B, default=False, flrcif=None,
        authored="Weight AV grid points by linker-chain statistics instead "
                 "of uniformly.",
        dialects=("flat",)),
    # --- rigid-body / multi-body assignment (ChiSurf dialect) ------------
    "body_id": dict(
        type=_I, default=0, flrcif=None,
        authored="Index of the rigid body (input PDB) this position labels; "
                 "docking moves bodies relative to each other.",
        dialects=("network", "csfps")),
    # --- fixed positions (XYZ simulation_type, legacy C# conversions) ----
    "x": dict(
        type=_F, default=None, flrcif="_flr_FPS_mean_probe_position.mpp_xcoord",
        dialects=("csfps",)),
    "y": dict(
        type=_F, default=None, flrcif="_flr_FPS_mean_probe_position.mpp_ycoord",
        dialects=("csfps",)),
    "z": dict(
        type=_F, default=None, flrcif="_flr_FPS_mean_probe_position.mpp_zcoord",
        dialects=("csfps",)),
}

DISTANCE_FIELDS: Dict[str, Dict[str, Any]] = {
    "position1_name": dict(
        type=_S, default=None,
        flrcif="_flr_fret_distance_restraint.sample_probe_id_1",
        required=True, dialects=("network", "flat")),
    "position2_name": dict(
        type=_S, default=None,
        flrcif="_flr_fret_distance_restraint.sample_probe_id_2",
        required=True, dialects=("network", "flat")),
    "distance": dict(
        type=_F, default=None,
        flrcif="_flr_fret_distance_restraint.distance",
        required=True, dialects=("network", "flat")),
    "error_neg": dict(
        type=_F, default=None,
        flrcif="_flr_fret_distance_restraint.distance_error_minus",
        required=True, dialects=("network", "flat")),
    "error_pos": dict(
        type=_F, default=None,
        flrcif="_flr_fret_distance_restraint.distance_error_plus",
        required=True, dialects=("network", "flat")),
    "distance_type": dict(
        type=_S, default="RDAMean",
        flrcif="_flr_fret_distance_restraint.distance_type",
        enum=tuple(DISTANCE_TYPES), dialects=("network", "flat")),
    "Forster_radius": dict(
        type=_F, default=52.0,
        flrcif="_flr_fret_forster_radius.forster_radius",
        dialects=("network", "flat")),
}

#: A score set (the "χ²" top-level section) selects a named subset of the
#: distances to score together. flrCIF groups restraints via
#: ``_flr_fret_distance_restraint.group_id``; the named-set spelling here is
#: the authored fps.json form of the same concept.
SCORE_SET_FIELDS: Dict[str, Dict[str, Any]] = {
    "distances": dict(
        type=list, default=None,
        flrcif="_flr_fret_distance_restraint.group_id",
        required=True,
        authored="Names of the Distances entries belonging to this set."),
    "maximum_NaNs_allowed": dict(
        type=_I, default=None, flrcif=None,
        authored="Olga screening: maximum number of NaN model distances "
                 "tolerated before a structure's score set is rejected.",
        dialects=("network",)),
    "penalty_NaN": dict(
        type=_F, default=None, flrcif=None,
        authored="Olga screening: chi-square penalty added per NaN model "
                 "distance.",
        dialects=("network",)),
}


# ---------------------------------------------------------------------------
# Derivation: JSON Schema
# ---------------------------------------------------------------------------

_JSON_TYPES = {_S: "string", _F: "number", _I: "integer", _B: "boolean",
               list: "array"}


def _field_to_property(spec: Dict[str, Any]) -> Dict[str, Any]:
    """Derive one JSON-Schema property from a field-table entry."""
    prop: Dict[str, Any] = {"type": _JSON_TYPES[spec["type"]]}
    # JSON has one number type; integral floats are fine for "number" but an
    # "integer" field must stay integral.
    if spec.get("enum"):
        prop["enum"] = list(spec["enum"])
    if spec.get("default") is not None:
        prop["default"] = spec["default"]
    flr = spec.get("flrcif")
    if flr:
        prop["x-flrcif-item"] = flr
    if spec.get("authored"):
        prop["description"] = spec["authored"]
    prop["x-dialects"] = list(spec.get("dialects", ()))
    return prop


def _object_schema(fields: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    required = sorted(n for n, s in fields.items() if s.get("required"))
    schema: Dict[str, Any] = {
        "type": "object",
        "properties": {n: _field_to_property(s) for n, s in fields.items()},
        # Unknown keys are preserved by readers/writers but are not silently
        # part of the format: additionalProperties stays true so old files
        # validate, and validate() reports unknown keys as warnings.
        "additionalProperties": True,
    }
    if required:
        schema["required"] = required
    return schema


def to_json_schema() -> Dict[str, Any]:
    """Derive the JSON-Schema document for the network-dialect fps.json.

    The shipped ``data/fps_json_schema.json`` is this function's output; a
    test regenerates it and fails on drift, so the Python tables here stay the
    single authored definition.
    """
    position = _object_schema(POSITION_FIELDS)
    distance = _object_schema(DISTANCE_FIELDS)
    score_set = _object_schema(SCORE_SET_FIELDS)
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://integrativemodeling.org/schemas/bff/fps.json",
        "title": "fps.json — FRET labelling positions and distances",
        "x-schema-version": SCHEMA_VERSION,
        "type": "object",
        "properties": {
            "Positions": {"type": "object",
                          "additionalProperties": position},
            "Distances": {"type": "object",
                          "additionalProperties": distance},
            "χ²": {"type": "object",
                   "additionalProperties": score_set},
            "Evaluators": {"type": "array", "items": {"type": "object"}},
        },
        "additionalProperties": True,
        "$defs": {
            "position": position,
            "distance": distance,
            "score_set": score_set,
        },
    }


# ---------------------------------------------------------------------------
# Validation (self-contained; no jsonschema dependency)
# ---------------------------------------------------------------------------

def _check_fields(
    obj: Dict[str, Any],
    fields: Dict[str, Dict[str, Any]],
    where: str,
    errors: List[str],
    warnings: List[str],
) -> None:
    if not isinstance(obj, dict):
        errors.append(f"{where}: expected an object, got {type(obj).__name__}")
        return
    for name, spec in fields.items():
        if name not in obj:
            if spec.get("required"):
                errors.append(f"{where}: missing required field '{name}'")
            continue
        value = obj[name]
        expected = spec["type"]
        ok = isinstance(value, expected)
        # bool is an int subclass; a bool where a number is expected is a bug.
        if expected in (_F, _I) and isinstance(value, bool):
            ok = False
        # ints are acceptable where floats are expected.
        if expected is _F and isinstance(value, int) and not isinstance(value, bool):
            ok = True
        if not ok:
            errors.append(
                f"{where}.{name}: expected {expected.__name__}, "
                f"got {type(value).__name__} ({value!r})")
            continue
        enum = spec.get("enum")
        if enum and value not in enum:
            errors.append(
                f"{where}.{name}: {value!r} is not one of {sorted(enum)}")
    for key in obj:
        if key not in fields:
            warnings.append(f"{where}: unknown field '{key}'")


def validate_position(
    position: Dict[str, Any], name: str = "position"
) -> Tuple[List[str], List[str]]:
    """Validate one position object. Returns ``(errors, warnings)``."""
    errors: List[str] = []
    warnings: List[str] = []
    _check_fields(position, POSITION_FIELDS, name, errors, warnings)
    if isinstance(position, dict):
        stype = position.get("simulation_type", "AV1")
        if stype == "R1":
            if not str(position.get("rotamer_library", "") or "").strip():
                errors.append(
                    f"{name}: simulation_type R1 requires 'rotamer_library'")
            for k in ("linker_length", "linker_width", "radius1"):
                if k in position:
                    warnings.append(
                        f"{name}: '{k}' is an AV parameter; an R1 position "
                        "ignores it (Python side) or is scored as AV1 by "
                        "AVNetworkRestraint -- filter with fps_positions_for_docking()")
        if stype == "XYZ":
            for k in ("x", "y", "z"):
                if k not in position:
                    errors.append(
                        f"{name}: simulation_type XYZ requires '{k}'")
        elif stype == "AV3":
            for k in ("radius2", "radius3"):
                if float(position.get(k, 0.0) or 0.0) <= 0.0:
                    warnings.append(
                        f"{name}: simulation_type AV3 with non-positive "
                        f"'{k}' behaves like AV1")
    return errors, warnings


def validate_distance(
    distance: Dict[str, Any],
    name: str = "distance",
    position_names: Optional[set] = None,
) -> Tuple[List[str], List[str]]:
    """Validate one distance object. Returns ``(errors, warnings)``."""
    errors: List[str] = []
    warnings: List[str] = []
    _check_fields(distance, DISTANCE_FIELDS, name, errors, warnings)
    if position_names is not None and isinstance(distance, dict):
        for key in ("position1_name", "position2_name"):
            pname = distance.get(key)
            if pname and pname not in position_names:
                errors.append(
                    f"{name}.{key}: '{pname}' is not a defined position")
    return errors, warnings


def validate_fps(payload: Dict[str, Any]) -> Tuple[List[str], List[str]]:
    """Validate a parsed network-dialect fps.json payload.

    Parameters
    ----------
    payload : dict
        The full JSON document (``Positions``/``Distances``/``χ²``/...).

    Returns
    -------
    (errors, warnings) : (list of str, list of str)
        ``errors`` violate the schema; ``warnings`` flag unknown fields and
        suspicious-but-legal values. An empty ``errors`` list means the file
        conforms.
    """
    errors: List[str] = []
    warnings: List[str] = []
    if not isinstance(payload, dict):
        return [f"payload: expected an object, got {type(payload).__name__}"], []

    positions = payload.get("Positions", {})
    if not isinstance(positions, dict):
        errors.append("Positions: expected an object")
        positions = {}
    for pname, pos in positions.items():
        e, w = validate_position(pos, f"Positions[{pname!r}]")
        errors += e
        warnings += w

    distances = payload.get("Distances", {})
    if not isinstance(distances, dict):
        errors.append("Distances: expected an object")
        distances = {}
    for dname, dist in distances.items():
        e, w = validate_distance(
            dist, f"Distances[{dname!r}]", position_names=set(positions))
        errors += e
        warnings += w

    score_sets = payload.get("χ²", {})
    if not isinstance(score_sets, dict):
        errors.append("χ²: expected an object")
        score_sets = {}
    for sname, sset in score_sets.items():
        e, w = [], []
        _check_fields(sset, SCORE_SET_FIELDS, f"χ²[{sname!r}]", e, w)
        errors += e
        warnings += w
        if isinstance(sset, dict):
            for ref in sset.get("distances", []) or []:
                if ref not in distances:
                    errors.append(
                        f"χ²[{sname!r}]: '{ref}' is not a defined distance")

    evaluators = payload.get("Evaluators")
    if evaluators is not None and not isinstance(evaluators, list):
        errors.append("Evaluators: expected an array")

    return errors, warnings
