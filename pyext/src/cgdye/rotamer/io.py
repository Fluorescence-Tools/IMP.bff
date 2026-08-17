"""Input/output helpers for rotamer libraries and protein frames."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import IMP
import IMP.atom
import IMP.core
import IMP.rmf
import numpy as np
import RMF
import json
from IMP.bff.cgdye.io.rotamer_rmf import read_rotamer_library_rmf
from IMP.bff.cgdye.rotamer.scoring import _selector_resnames
from IMP.bff.cgdye.utils import get_template_dir

_LIBRARY_REGISTRY: dict[str, dict[str, Any]] | None = None


def _registry_path() -> Path:
    """Return the bundled FRETpredict rotamer-library registry path.

    Returns
    -------
    pathlib.Path
        Path to ``libraries.json``.
    """
    import IMP.bff
    return Path(IMP.bff.get_data_path("rotamer_library")) / "libraries.json"


def get_library_registry() -> dict[str, dict[str, Any]]:
    """Load the bundled FRETpredict rotamer-library registry.

    Returns
    -------
    dict
        Registry keyed by FRETpredict library name.
    """
    global _LIBRARY_REGISTRY
    if _LIBRARY_REGISTRY is None:
        with _registry_path().open() as handle:
            _LIBRARY_REGISTRY = json.load(handle)
    return _LIBRARY_REGISTRY


def normalize_library_name(library_name: str) -> str:
    """Return the registry key for a rotamer-library name.

    Parameters
    ----------
    library_name : str
        User-facing library name, for example
        ``AlexaFluor 488 C1R cutoff30``.

    Returns
    -------
    str
        Registry key without the cutoff suffix.
    """
    return re.sub(r"\s+cutoff\d+$", "", str(library_name).strip())


def get_library_metadata(library_name: str) -> dict[str, Any]:
    """Return metadata for a FRETpredict-style rotamer library.

    Parameters
    ----------
    library_name : str
        User-facing library name.

    Returns
    -------
    dict
        Library metadata from ``libraries.json``.
    """
    key = normalize_library_name(library_name)
    registry = get_library_registry()
    if key not in registry:
        raise ValueError(f"Unknown rotamer library {library_name!r}")
    metadata = dict(registry[key])
    metadata["name"] = key
    metadata["library_name"] = library_name
    metadata["cutoff"] = _cutoff_from_library_name(library_name)
    return metadata


def _cutoff_from_library_name(library_name: str) -> int | None:
    """Extract the cutoff value from a library name.

    Parameters
    ----------
    library_name : str
        User-facing library name.

    Returns
    -------
    int or None
        Cutoff value if present.
    """
    match = re.search(r"cutoff(\d+)", str(library_name))
    return int(match.group(1)) if match else None


def _library_filename(metadata: dict[str, Any], cutoff: int | None = None) -> str:
    """Return the FRETpredict library file stem.

    Parameters
    ----------
    metadata : dict
        Library metadata.
    cutoff : int, optional
        Cutoff value.

    Returns
    -------
    str
        File stem without extension.
    """
    base = str(metadata["filename"])
    if cutoff is None:
        return base
    # The registry's filename carries FRETpredict's default cutoff
    # (``A48_C1R_cutoff30``); a name that asks for another cutoff replaces it.
    stem = re.sub(r"_?cutoff\d+$", "", base)
    return f"{stem}_cutoff{int(cutoff)}"


def resolve_rotamer_library_path(library_name: str, lib_dir: str | Path | None = None) -> Path:
    """Resolve a FRETpredict-style library name to an IMP-native RMF file.

    Parameters
    ----------
    library_name : str
        Library name or explicit path.
    lib_dir : pathlib.Path or str, optional
        Directory containing RMF rotamer libraries.

    Returns
    -------
    pathlib.Path
        Existing RMF library path.
    """
    candidate = Path(str(library_name))
    if candidate.exists():
        return candidate

    metadata = get_library_metadata(library_name)
    cutoff = metadata.get("cutoff")
    filename = _library_filename(metadata, cutoff)
    stem = filename.split("_cutoff")[0]

    # The FRETpredict library files (module data, data/rotamer_library) are the
    # canonical libraries: <stem>.pdb + <stem>_cutoff<N>.dcd (+ weights) for
    # each cutoff. They are tried first so that the *requested cutoff* is the
    # one loaded. The RMF templates under templates/rotamer hold only the
    # cutoff-30 clustering, so resolving every name to <stem>.rmf3 silently
    # returned the wrong library for cutoff10/cutoff20 names.
    if lib_dir is None:
        dcd = _registry_path().parent / f"{filename}.dcd"
        if dcd.exists() and dcd.with_name(f"{stem}.pdb").exists():
            return dcd

    template_dir = Path(lib_dir) if lib_dir is not None else get_template_dir("rotamer")
    candidates = [
        template_dir / f"{filename}.dcd",
        template_dir / f"{filename}.rmf3",
        template_dir / f"{stem}.rmf3",
        template_dir / f"{stem}.pdb",
        template_dir / f"{filename}.pdb",
    ]
    for path in candidates:
        if path.exists():
            if path.suffix.lower() == ".rmf3" and cutoff not in (None, 30) and path.stem == stem:
                raise FileNotFoundError(
                    f"{library_name!r}: only the cutoff-30 RMF template {path.name} is available; "
                    f"the cutoff-{cutoff} library needs {filename}.dcd next to {stem}.pdb")
            return path
    raise FileNotFoundError(f"No rotamer library found for {library_name!r}")


def _coords_to_array(library: dict[str, Any]) -> np.ndarray:
    """Convert RMF rotamer coordinates to an array.

    Parameters
    ----------
    library : dict
        Rotamer library returned by ``read_rotamer_library_rmf``.

    Returns
    -------
    numpy.ndarray
        Array with shape ``(n_rotamers, n_atoms, 3)``.
    """
    ids = sorted(library["coords"].keys())
    return np.stack([library["coords"][rid] for rid in ids], axis=0).astype(np.float64)


def _metadata_from_path(path: Path) -> dict[str, Any]:
    """Infer FRETpredict metadata from an RMF template file stem.

    Parameters
    ----------
    path : pathlib.Path
        RMF template path.

    Returns
    -------
    dict
        Metadata dictionary.
    """
    stem = path.stem
    for key, metadata in get_library_registry().items():
        base = str(metadata["filename"]).split("_cutoff")[0]
        if stem == base:
            result = dict(metadata)
            result["name"] = key
            result["library_name"] = key
            result["cutoff"] = None
            return result
    return {}


def _resnames_from_pdb(path: Path) -> list[str] | None:
    """Return residue names from a PDB atom block.

    Parameters
    ----------
    path : pathlib.Path
        PDB path.

    Returns
    -------
    list of str or None
        Residue names if the file contains atom records.
    """
    resnames: list[str] = []
    with path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            resnames.append(line[17:20].strip())
            if len(resnames) >= 10000:
                break
    return resnames or None


def _pdb_path_for_metadata(metadata: dict[str, Any]) -> Path | None:
    """Return a bundled PDB template path for rotamer metadata.

    Parameters
    ----------
    metadata : dict
        Rotamer metadata.

    Returns
    -------
    pathlib.Path or None
        Matching PDB path if available.
    """
    filename = str(metadata.get("filename", ""))
    stem = filename.split("_cutoff")[0]
    if not stem:
        return None
    candidates = [
        _registry_path().parent / f"{stem}.pdb",
        _registry_path().parent / f"{filename}.pdb",
    ]
    for path in candidates:
        if path.exists():
            return path
    return None


def _infer_rotamer_resnames(atom_names: list[str], metadata: dict[str, Any]) -> list[str] | None:
    """Infer rotamer residue names when a sidecar PDB is unavailable.

    Parameters
    ----------
    atom_names : list of str
        Atom names.
    metadata : dict
        Rotamer metadata.

    Returns
    -------
    list of str or None
        Residue names, or None when no inference is possible.
    """
    selector_resnames = _selector_resnames_from_metadata(metadata)
    if not selector_resnames:
        return None
    dye_resname = next(iter(selector_resnames))
    match = re.match(r".+_(?P<linker>[A-Z]\d?[A-Z]R)$", str(metadata.get("name", "")))
    linker_resname = match.group("linker") if match else None
    if linker_resname is None:
        return [dye_resname] * len(atom_names)
    linker_atoms = {"CA", "HA", "C", "O", "C6", "H10", "H11", "S1", "C7", "C8", "H12", "C9", "O3", "N3", "C10", "O4", "C11", "H13", "H14", "C12", "H15", "H16", "C13", "H17", "H18", "C14", "H19", "H20", "C15", "H21", "H22", "N99", "H23", "N", "H", "HX2", "HX3"}
    return [linker_resname if name in linker_atoms else dye_resname for name in atom_names]


def _selector_resnames_from_metadata(metadata: dict[str, Any]) -> set[str]:
    """Return selector residue names from rotamer metadata.

    Parameters
    ----------
    metadata : dict
        Rotamer metadata.

    Returns
    -------
    set[str]
        Residue names referenced by selectors.
    """
    resnames: set[str] = set()
    for key in ("mu", "r", "positive", "negative"):
        resnames.update(_selector_resnames(metadata.get(key)))
    return resnames


def load_rotamer_library(
    library_name: str,
    lib_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Load an IMP-native rotamer library.

    Parameters
    ----------
    library_name : str
        FRETpredict-style library name or explicit RMF path.
    lib_dir : pathlib.Path or str, optional
        Directory containing RMF rotamer libraries.

    Returns
    -------
    dict
        Library dictionary with ``coords``, ``weights``, ``atom_names``, and
        ``metadata``.
    """
    path = resolve_rotamer_library_path(library_name, lib_dir=lib_dir)
    explicit_path = Path(str(library_name)).exists()
    metadata = get_library_metadata(library_name) if not explicit_path else _metadata_from_path(path)
    suffix = path.suffix.lower()
    if suffix == ".dcd":
        # FRETpredict library set: <stem>.pdb (names, residues) + DCD frames +
        # per-rotamer weights, read with the in-tree DCD reader.
        from IMP.bff.cgdye.sampling.rotamer import load_reference_rotamers
        stem = path.stem.split("_cutoff")[0]
        pdb_path = path.with_name(f"{stem}.pdb")
        weights_path = path.with_name(f"{path.stem}_weights.txt")
        ref = load_reference_rotamers(pdb_path, path, weights_path if weights_path.exists() else None)
        coords = np.asarray(ref["coords"], dtype=np.float64)
        weights = np.asarray(ref["weights"], dtype=np.float64)
        library = {
            "id": list(range(1, coords.shape[0] + 1)),
            "atom_names": [str(n) for n in ref["atom_names"]],
            "transitions": None,
        }
    elif suffix == ".rmf3":
        library = read_rotamer_library_rmf(str(path))
        weights = np.asarray(library["weight"], dtype=np.float64)
        if weights.size == 0 or np.sum(weights) <= 0:
            weights = np.ones(library["coords"][1].shape[0], dtype=np.float64)
        weights = weights / np.sum(weights)
        coords = _coords_to_array(library)
        pdb_path = path.with_suffix(".pdb")
    else:
        raise ValueError(f"Unsupported rotamer library file: {path}")
    library["path"] = str(path)
    library["coords"] = coords
    library["weight"] = weights
    library["weights"] = weights
    pdb_resnames = _resnames_from_pdb(pdb_path) if pdb_path.exists() else None
    if pdb_resnames is None:
        bundled_pdb = _pdb_path_for_metadata(metadata)
        pdb_resnames = _resnames_from_pdb(bundled_pdb) if bundled_pdb is not None else None
    library["resnames"] = pdb_resnames or _infer_rotamer_resnames(library["atom_names"], metadata)
    library["metadata"] = metadata
    return library


def _atom_particles(hierarchy) -> list[Any]:
    """Return atom-like leaf particles from a hierarchy.

    Parameters
    ----------
    hierarchy
        IMP hierarchy.

    Returns
    -------
    list
        Leaf particles with XYZ decorators.
    """
    leaves = IMP.atom.get_leaves(hierarchy)
    return [p for p in leaves if IMP.core.XYZ.get_is_setup(p)]


def _collect_frame(hierarchy) -> tuple[np.ndarray, list[str], list[str], list[str], list[str], list[int]]:
    """Collect coordinates and atom metadata from one IMP frame.

    Parameters
    ----------
    hierarchy
        IMP hierarchy.

    Returns
    -------
    tuple
        ``(coords, atom_names, atom_types, resnames, chain_ids, residue_indices)``.
    """
    particles = _atom_particles(hierarchy)
    coords = np.zeros((len(particles), 3), dtype=np.float64)
    atom_names: list[str] = []
    atom_types: list[str] = []
    resnames: list[str] = []
    chain_ids: list[str] = []
    residue_indices: list[int] = []
    for i, p in enumerate(particles):
        xyz = IMP.core.XYZ(p)
        coords[i] = [xyz.get_x(), xyz.get_y(), xyz.get_z()]
        if IMP.atom.Atom.get_is_setup(p):
            atom = IMP.atom.Atom(p)
            name = atom.get_name()
            parts = name.split()
            atom_names.append(parts[1] if len(parts) > 1 else parts[0])
            atom_types.append(atom.get_atom_type().get_string())
            res_p = p.get_parent()
            if IMP.atom.Residue.get_is_setup(res_p):
                resnames.append(str(IMP.atom.Residue(res_p).get_residue_type().get_string()))
                chain_p = res_p.get_parent()
                chain_ids.append(str(IMP.atom.Chain(chain_p).get_id()) if IMP.atom.Chain.get_is_setup(chain_p) else "")
                residue_indices.append(int(IMP.atom.Residue(res_p).get_index()))
            else:
                resnames.append("")
                chain_ids.append("")
                residue_indices.append(-1)
        else:
            name = p.get_name()
            atom_names.append(name)
            atom_types.append(name[0] if name else "C")
            resnames.append("")
            chain_ids.append("")
            residue_indices.append(-1)
    return coords, atom_names, atom_types, resnames, chain_ids, residue_indices


def load_protein_frames(
    protein: str | Path,
    max_frames: int | None = None,
) -> list[dict[str, Any]]:
    """Load protein frames from PDB or RMF.

    Parameters
    ----------
    protein : pathlib.Path or str
        PDB or RMF path.
    max_frames : int, optional
        Maximum number of frames to load.

    Returns
    -------
    list of dict
        Frame dictionaries containing ``coords``, ``atom_names``,
        ``atom_types``, and ``resnames``.
    """
    path = Path(protein)
    suffix = path.suffix.lower()
    if suffix in {".pdb", ".ent"}:
        model = IMP.Model()
        # A multi-MODEL PDB (a trajectory written as models) yields one frame
        # per model; a plain PDB yields one frame.
        hierarchies = IMP.atom.read_multimodel_pdb(str(path), model, IMP.atom.NonWaterPDBSelector())
        if max_frames is not None:
            hierarchies = hierarchies[:max_frames]
        frames = []
        for hierarchy in hierarchies:
            coords, atom_names, atom_types, resnames, chain_ids, residue_indices = _collect_frame(hierarchy)
            frames.append({"coords": coords, "atom_names": atom_names, "atom_types": atom_types, "resnames": resnames, "chain_ids": chain_ids, "residue_indices": residue_indices})
        return frames

    if suffix in {".rmf", ".rmf3"}:
        model = IMP.Model()
        fh = RMF.open_rmf_file_read_only(str(path))
        hierarchies = IMP.rmf.create_hierarchies(fh, model)
        frames: list[dict[str, Any]] = []
        n_frames = fh.get_number_of_frames()
        limit = n_frames if max_frames is None else min(n_frames, max_frames)
        for frame_id in range(limit):
            IMP.rmf.load_frame(fh, RMF.FrameID(frame_id))
            coords, atom_names, atom_types, resnames, chain_ids, residue_indices = _collect_frame(hierarchies[0])
            frames.append(
                {
                    "coords": coords,
                    "atom_names": atom_names,
                    "atom_types": atom_types,
                    "resnames": resnames,
                    "chain_ids": chain_ids,
                    "residue_indices": residue_indices,
                }
            )
        return frames

    raise ValueError(f"Unsupported protein file type: {path}")
