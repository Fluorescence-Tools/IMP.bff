"""fps.json and legacy FPS file I/O — the one reader (PRD-97 stage 1).

This module is the single Python reader/writer for fps.json labelling files
and the legacy C# FPS ``.txt`` formats. It absorbed ChiSurf's
``fret/core/io.py`` (the only implementation that round-trips) and the old
``IMP.bff.fps`` module (a read-only subset, now deleted). Files are checked
against the authored schema in :mod:`.fps_schema` rather than against another
parser.

Legacy C# FPS support is **read only**: a format nobody can still read is data
that has been lost, and a decade of measurements live in these files — but
nothing should write them again.
"""

from __future__ import annotations

import json
import os
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

from . import fps_schema

__all__ = [
    "read_fps_json",
    "write_fps_json",
    "read_evaluators_json",
    "write_evaluators_json",
    "read_old_lps_txt",
    "read_old_distances_txt",
    "load_structure",
    "load_structure_with_particles",
    "write_pdb",
    "write_rmf",
    "compute_rmsd",
]


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
        (:func:`IMP.bff.fret.fps_schema.validate`) and raise ``ValueError``
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


# ---------------------------------------------------------------------------
# Evaluators section
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Legacy C# FPS formats (read only)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# PDB writing (pure Python, no IMP)
# ---------------------------------------------------------------------------

def write_pdb(
    atoms: np.ndarray,
    path: str | os.PathLike,
    chain: str = "A",
    res_name: str = "ALA",
    transform: Optional[np.ndarray] = None,
    model_index: int = 0,
) -> None:
    """Write (N,3) or (N,4) coordinates as a multi-model PDB.

    Parameters
    ----------
    atoms : (N, 3) or (N, 4) ndarray
        If (N,4), column 4 is ignored (element hint).
    transform : (4,4), (3,3), or (3,) optional
        Homogeneous transform, rotation matrix, or translation.
    """
    if atoms.ndim != 2 or atoms.shape[1] < 3:
        raise ValueError(f"Expected (N,3) or (N,4) array, got {atoms.shape}")
    coords = atoms[:, :3].copy()
    if transform is not None:
        coords = _apply_transform(coords, transform)

    lines = [f"MODEL     {model_index:>4d}"]
    for i in range(coords.shape[0]):
        x, y, z = coords[i]
        name = "CA"
        name_fmt = f"{name:>4s}" if len(name) >= 4 else f" {name:<3s}"
        lines.append(
            f"ATOM  {i + 1:>5d} {name_fmt}{res_name:<3s} {chain}"
            f"{i + 1:>4d}    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
            f"{1.0:>6.2f}{0.0:>6.2f}"
        )
    lines.append("ENDMDL")
    lines.append("END")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _apply_transform(coords: np.ndarray, transform: np.ndarray) -> np.ndarray:
    """Apply a (4,4) homogeneous, (3,3) rotation, or (3,) translation."""
    t = np.asarray(transform, dtype=np.float64)
    if t.shape == (4, 4):
        return coords @ t[:3, :3].T + t[:3, 3]
    if t.shape == (3, 3):
        return coords @ t.T
    if t.shape == (3,):
        return coords + t
    raise ValueError(f"Unexpected transform shape {t.shape}")


# ---------------------------------------------------------------------------
# RMF writing (IMP.rmf directly — no application writer involved)
# ---------------------------------------------------------------------------

def write_rmf(
    atoms: np.ndarray,
    path: str | os.PathLike,
    model_name: str = "structure",
    transform: Optional[np.ndarray] = None,
    metadata: Optional[Dict] = None,
    radius: float = 1.5,
) -> None:
    """Write ``(N, 3)`` coordinates to an RMF file via ``IMP.rmf``.

    Each coordinate becomes an ``XYZR`` ball under a single hierarchy named
    ``model_name``; ``metadata`` (when given) is stored as the RMF file
    description. This is a self-contained replacement for the application
    RMF writer the pre-move implementation borrowed.

    Parameters
    ----------
    atoms : (N, 3) or (N, 4) ndarray
        Coordinates; a 4th column is ignored.
    path : str or PathLike
        Output ``.rmf``/``.rmf3`` path.
    model_name : str
        Name of the root hierarchy node.
    transform : (4,4), (3,3), or (3,) optional
        Applied to the coordinates before writing.
    metadata : dict, optional
        JSON-serialised into the RMF file description.
    radius : float
        Ball radius (Angstrom) given to every particle.
    """
    import IMP
    import IMP.atom
    import IMP.core
    import IMP.rmf
    import RMF

    coords = np.asarray(atoms[:, :3], dtype=np.float64).copy()
    if transform is not None:
        coords = _apply_transform(coords, transform)

    model = IMP.Model()
    root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model, str(model_name)))
    for i in range(coords.shape[0]):
        p = IMP.Particle(model, f"p{i}")
        IMP.core.XYZR.setup_particle(
            p, IMP.algebra.Sphere3D(IMP.algebra.Vector3D(*coords[i]), radius))
        IMP.atom.Mass.setup_particle(p, 1.0)
        root.add_child(IMP.atom.Hierarchy.setup_particle(p))

    fh = RMF.create_rmf_file(str(path))
    try:
        if metadata:
            fh.set_description(json.dumps(metadata, default=str))
        IMP.rmf.add_hierarchy(fh, root)
        IMP.rmf.save_frame(fh, "frame_0")
    finally:
        del fh


# ---------------------------------------------------------------------------
# Structure loading (IMP.atom)
# ---------------------------------------------------------------------------

def load_structure(pdb_path: str | os.PathLike) -> np.ndarray:
    """Load a PDB and return (N, 3) coordinates.

    Requires ``IMP.atom``.
    """
    coords, _leaves, _model, _hier = load_structure_with_particles(pdb_path)
    return coords


def load_structure_with_particles(pdb_path: str | os.PathLike):
    """Load a PDB and return (coords, leaves, model, hierarchy).

    Requires ``IMP.atom``.
    """
    import IMP
    import IMP.atom
    import IMP.core

    model = IMP.Model()
    hier = IMP.atom.read_pdb(str(pdb_path), model, IMP.atom.NonWaterPDBSelector())
    leaves = IMP.atom.get_leaves(hier)
    n = len(leaves)
    coords = np.empty((n, 3), dtype=np.float64)
    for i, p in enumerate(leaves):
        xyz = IMP.core.XYZ(p)
        coords[i, 0] = xyz.get_x()
        coords[i, 1] = xyz.get_y()
        coords[i, 2] = xyz.get_z()
    return coords, leaves, model, hier


def compute_rmsd(
    coords_a: np.ndarray,
    coords_b: np.ndarray,
    selection_mask: Optional[np.ndarray] = None,
    superpose: bool = False,
) -> float:
    """Compute RMSD between two coordinate arrays, optionally with superposition.

    Parameters
    ----------
    coords_a : (N, 3) ndarray
        First coordinate array.
    coords_b : (N, 3) ndarray
        Second coordinate array.
    selection_mask : (N,) bool ndarray, optional
        If provided, only the selected atoms are used to calculate the RMSD
        (and the superposition transform, if enabled).
    superpose : bool
        If True, apply Kabsch rotation to align coords_a onto coords_b
        before computing RMSD.

    Returns
    -------
    rmsd : float
        The root-mean-square deviation.
    """
    if coords_a.shape != coords_b.shape:
        raise ValueError(f"Shape mismatch: {coords_a.shape} vs {coords_b.shape}")

    if selection_mask is not None:
        a_sel = coords_a[selection_mask]
        b_sel = coords_b[selection_mask]
    else:
        a_sel = coords_a
        b_sel = coords_b

    if len(a_sel) == 0:
        return 0.0

    if superpose:
        # Compute centroids
        centroid_a = a_sel.mean(axis=0)
        centroid_b = b_sel.mean(axis=0)

        # Center coordinates
        a_centered = a_sel - centroid_a
        b_centered = b_sel - centroid_b

        # SVD of covariance matrix
        c = a_centered.T @ b_centered
        v, s, wt = np.linalg.svd(c)

        # Det for reflection detection
        d = np.linalg.det(v @ wt)
        s_diag = np.array([1.0, 1.0, np.sign(d) if d != 0 else 1.0])
        r = v @ np.diag(s_diag) @ wt

        # Rotate entire coords_a
        coords_a_aligned = (coords_a - centroid_a) @ r + centroid_b
        if selection_mask is not None:
            diff = coords_a_aligned[selection_mask] - coords_b[selection_mask]
        else:
            diff = coords_a_aligned - coords_b
    else:
        diff = a_sel - b_sel

    return float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))
