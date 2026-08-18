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

from __future__ import annotations

import json
import os
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

from IMP.bff.io import fps_schema
# ``read_fps_json`` dispatches on the file: a legacy ``LPs.txt`` next to a PDB is
# converted on the way in, so the one entry point reads both eras. That edge is
# one-way -- fps_legacy knows nothing about fps.json.
from IMP.bff.io.fps_legacy import read_old_distances_txt, read_old_lps_txt

__all__ = [
    "read_fps_json",
    "write_fps_json",
    "fps_positions_for_docking",
    "AV_SIMULATION_TYPES",
    "read_evaluators_json",
    "write_evaluators_json",
]

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
        errors, _warnings = fps_schema.validate(payload)
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
        errors, _warnings = fps_schema.validate(payload)
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
