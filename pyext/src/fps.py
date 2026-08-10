"""Read and write fps.json labeling files."""

from __future__ import annotations

import json
import os
from typing import Any, Optional


def read_fps_json(
    path: str | os.PathLike,
    pdb_paths: Optional[list[str]] = None,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Load an fps.json labeling file.

    Parameters
    ----------
    path : str or os.PathLike
        Path to a JSON labeling file.
    pdb_paths : list of str, optional
        Optional PDB path list retained for compatibility with legacy readers.

    Returns
    -------
    tuple[dict, dict, dict, dict]
        Position, distance, score-set, and extra payload dictionaries.
    """
    path_str = str(path)
    if not path_str.endswith(".json"):
        return {}, {}, {}, {}

    with open(path) as handle:
        payload = json.load(handle)

    positions = payload.pop("Positions", {})
    distances = payload.pop("Distances", {})
    score_sets = payload.pop("χ²", {})
    return positions, distances, score_sets, payload


def write_fps_json(
    path: str | os.PathLike,
    positions: dict[str, Any],
    distances: dict[str, Any],
    score_sets: Optional[dict[str, Any]] = None,
    extra: Optional[dict[str, Any]] = None,
    **kwargs: Any,
) -> None:
    """Write an fps.json labeling file.

    Parameters
    ----------
    path : str or os.PathLike
        Output JSON path.
    positions : dict
        Position entries.
    distances : dict
        Distance entries.
    score_sets : dict, optional
        Optional chi-square score sets.
    extra : dict, optional
        Optional extra top-level payload.
    **kwargs : Any
        Additional ``json.dump`` options.
    """
    payload: dict[str, Any] = {}
    if extra is not None:
        payload.update(extra)
    payload["Positions"] = positions
    payload["Distances"] = distances
    if score_sets:
        payload["χ²"] = score_sets
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2, **kwargs)
