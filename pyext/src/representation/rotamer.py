"""Rotamer-based cgdye utilities."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Sequence
import RMF
import json
import logging
import math
import re

import numpy as np

from IMP.bff import forster_radius_from_spectra
from IMP.bff.io.fps import read_fps_json
from IMP.bff.io.structure import read_rotamer_library_rmf
from IMP.bff.representation.distance import fret_pair_efficiencies, fret_pair_geometry
from IMP.bff import States
from IMP.bff.scoring import compute_rotamer_score, selector_resnames
from IMP.bff import get_template_dir
import IMP
import IMP.atom
import IMP.core
import IMP.bff
import IMP.rmf

__all__ = [
    'RotamerDistance',
    'RotamerEnsemble',
    'RotamerFRET',
    'RotamerPosition',
    'SIMULATION_TYPE_R1',
    'backbone_rotation',
    'read_rotamer_fps',
    'resolve_backbone_site',
    'rotamer_ensembles_from_fps',
    'rotamer_fret_from_fps',
    'selector_atom_indices',
    'transform_library_to_site',
]

# --------------------------------------------------------------------------
# fps
# --------------------------------------------------------------------------
"""fps.json helpers for rotamer FRET."""

@dataclass
class RotamerPosition:
    """Rotamer labeling position from an fps.json file.

    Attributes
    ----------
    name : str
        Position name.
    chain : str or None
        Chain identifier.
    residue : int
        Residue number.
    atom_name : str
        Anchor atom name.
    dye : str or None
        Optional dye name.
    library : str or None
        Optional rotamer library name.
    role : str or None
        Optional label role, such as donor or acceptor.
    """

    name: str
    chain: str | None
    residue: int
    atom_name: str = "CA"
    dye: str | None = None
    library: str | None = None
    role: str | None = None

    @classmethod
    def from_payload(cls, name: str, payload: dict[str, Any]) -> "RotamerPosition":
        """Create a rotamer position from an fps.json position entry.

        Parameters
        ----------
        name : str
            Position name.
        payload : dict
            Position entry.

        Returns
        -------
        RotamerPosition
            Parsed position.
        """
        return cls(
            name=name,
            chain=_first(payload, "chain_identifier", "chain", "segid"),
            residue=int(_first(payload, "residue_seq_number", "residue", "resid")),
            atom_name=str(_first(payload, "atom_name", default="CA")),
            dye=_first(payload, "dye_name", "dye", "chromophore"),
            library=_first(payload, "rotamer_library", "library_name", "libname", "library"),
            role=_first(payload, "role", "label_type"),
        )


@dataclass
class RotamerDistance:
    """FRET distance between two rotamer positions.

    Attributes
    ----------
    name : str
        Distance name.
    donor_position : str
        Donor position name.
    acceptor_position : str
        Acceptor position name.
    donor : str or None
        Optional donor dye name.
    acceptor : str or None
        Optional acceptor dye name.
    libname_1 : str or None
        Optional donor rotamer library name.
    libname_2 : str or None
        Optional acceptor rotamer library name.
    """

    name: str
    donor_position: str
    acceptor_position: str
    donor: str | None = None
    acceptor: str | None = None
    libname_1: str | None = None
    libname_2: str | None = None

    @classmethod
    def from_payload(
        cls,
        name: str,
        payload: dict[str, Any],
        positions: dict[str, Any],
    ) -> "RotamerDistance":
        """Create a rotamer distance from an fps.json distance entry.

        Parameters
        ----------
        name : str
            Distance name.
        payload : dict
            Distance entry.
        positions : dict
            Position entries.

        Returns
        -------
        RotamerDistance
            Parsed distance.
        """
        donor_position = _first(payload, "position1_name", "donor_position", "donor_position_name", "dye1_position")
        acceptor_position = _first(payload, "position2_name", "acceptor_position", "acceptor_position_name", "dye2_position")
        donor = _first(payload, "donor", "donor_dye", "dye1", "dye_1")
        acceptor = _first(payload, "acceptor", "acceptor_dye", "dye2", "dye_2")
        if donor is None:
            donor = positions.get(donor_position, {}).get("dye")
        if acceptor is None:
            acceptor = positions.get(acceptor_position, {}).get("dye")
        libname_1 = _first(payload, "libname_1", "donor_library", "donor_rotamer_library", "library_1")
        libname_2 = _first(payload, "libname_2", "acceptor_library", "acceptor_rotamer_library", "library_2")
        if libname_1 is None:
            libname_1 = positions.get(donor_position, {}).get("rotamer_library") or positions.get(donor_position, {}).get("library")
        if libname_2 is None:
            libname_2 = positions.get(acceptor_position, {}).get("rotamer_library") or positions.get(acceptor_position, {}).get("library")
        return cls(
            name=name,
            donor_position=donor_position,
            acceptor_position=acceptor_position,
            donor=donor,
            acceptor=acceptor,
            libname_1=libname_1,
            libname_2=libname_2,
        )


def _first(payload: dict[str, Any], *keys: str, default: Any = None) -> Any:
    """Return the first non-empty value from a payload.

    Parameters
    ----------
    payload : dict
        Mapping to search.
    *keys : str
        Candidate keys.
    default : Any, optional
        Default value.

    Returns
    -------
    Any
        First non-empty value.
    """
    for key in keys:
        value = payload.get(key)
        if value not in (None, ""):
            return value
    return default


def read_rotamer_fps(
    path: str | Path,
    distance_name: str | None = None,
) -> tuple[RotamerPosition, RotamerPosition, RotamerDistance, dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Read the first rotamer FRET distance from an fps.json file.

    Parameters
    ----------
    path : str or pathlib.Path
        fps.json path.
    distance_name : str, optional
        Distance entry to select.

    Returns
    -------
    tuple[RotamerPosition, RotamerPosition, RotamerDistance, dict, dict, dict]
        Donor position, acceptor position, distance, positions, distances,
        and extra payload.
    """
    positions, distances, score_sets, extra = read_fps_json(path)
    if not distances:
        raise ValueError(f"No distances found in {path}")
    if distance_name is None:
        distance_name = next(iter(distances))
    if distance_name not in distances:
        raise ValueError(f"Distance {distance_name!r} not found in {path}")
    distance = RotamerDistance.from_payload(distance_name, distances[distance_name], positions)
    if distance.donor_position not in positions:
        raise ValueError(f"Donor position {distance.donor_position!r} not found in {path}")
    if distance.acceptor_position not in positions:
        raise ValueError(f"Acceptor position {distance.acceptor_position!r} not found in {path}")
    donor = RotamerPosition.from_payload(distance.donor_position, positions[distance.donor_position])
    acceptor = RotamerPosition.from_payload(distance.acceptor_position, positions[distance.acceptor_position])
    return donor, acceptor, distance, positions, distances, score_sets | extra


def rotamer_fret_from_fps(
    fps_path: str | Path,
    protein_path: str | Path,
    distance_name: str | None = None,
    **kwargs: Any,
) -> Any:
    """Build a ``RotamerFRET`` object from an fps.json file.

    Parameters
    ----------
    fps_path : str or pathlib.Path
        fps.json path.
    protein_path : str or pathlib.Path
        Protein PDB or RMF path.
    distance_name : str, optional
        Distance entry to select.
    **kwargs : Any
        Additional ``RotamerFRET`` options.

    Returns
    -------
    RotamerFRET
        Configured rotamer FRET object.
    """
    # (was: from .fret import ...) -- now in this module

    donor, acceptor, distance, _positions, _distances, _extra = read_rotamer_fps(fps_path, distance_name=distance_name)
    chains = kwargs.get("chains")
    if chains is None:
        chains = [donor.chain, acceptor.chain]
    if kwargs.get("donor") is None:
        kwargs["donor"] = distance.donor or donor.dye
    if kwargs.get("acceptor") is None:
        kwargs["acceptor"] = distance.acceptor or acceptor.dye
    if kwargs.get("libname_1") is None:
        kwargs["libname_1"] = distance.libname_1 or donor.library
    if kwargs.get("libname_2") is None:
        kwargs["libname_2"] = distance.libname_2 or acceptor.library
    return RotamerFRET(
        protein_path,
        [donor.residue, acceptor.residue],
        chains=chains,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# Writing: rotamer ensembles -> fps.json R1 positions and predicted distances
# ---------------------------------------------------------------------------

def rotamer_position_payload(
    chain: str | None,
    residue: int,
    library: str,
    *,
    atom_name: str = "CA",
    dye: str | None = None,
    temperature: float | None = None,
    electrostatic: bool | None = None,
    potential: str | None = None,
) -> dict[str, Any]:
    """The fps.json entry of a rotamer-ensemble position (``simulation_type`` ``R1``)."""
    payload: dict[str, Any] = {
        "chain_identifier": chain or "",
        "residue_seq_number": int(residue),
        "atom_name": atom_name,
        "simulation_type": "R1",
        "rotamer_library": str(library),
    }
    if dye:
        payload["dye_name"] = str(dye)
    if temperature is not None:
        payload["temperature"] = float(temperature)
    if electrostatic is not None:
        payload["electrostatic"] = bool(electrostatic)
    if potential is not None:
        payload["potential"] = str(potential)
    return payload


def rotamer_ensemble_payload(ensemble: Any, atom_name: str = "CA") -> dict[str, Any]:
    """The fps.json ``R1`` entry describing an existing ``RotamerEnsemble``."""
    params = dict(getattr(ensemble, "params", {}) or {})
    return rotamer_position_payload(
        getattr(ensemble, "chain", params.get("chain")),
        getattr(ensemble, "residue", params.get("residue")),
        getattr(ensemble, "library", params.get("library")),
        atom_name=atom_name,
        temperature=params.get("temperature"),
        electrostatic=params.get("electrostatic"),
        potential=params.get("potential"),
    )


def distances_from_ensembles(
    ensembles: dict[str, Any],
    pairs: list[tuple[str, str]],
    forster_radius: float,
    *,
    distance_type: str = "RDAMeanE",
    error: float | None = None,
    error_fraction: float = 0.05,
    kappa2: str = "isotropic",
) -> dict[str, dict[str, Any]]:
    """Predicted fps.json distance entries between rotamer ensembles.

    ``distance_type`` selects what the ensembles predict: ``RDAMean`` (⟨R_DA⟩),
    ``RDAMeanE`` (FRET-averaged ⟨R_DA⟩_E) or ``Rmp`` (distance between mean
    positions), computed from the full pair matrix (no sampling).
    ``forster_radius`` in Å for κ² = 2/3. ⟨R_DA⟩_E follows the fps.json /
    AVNetworkRestraint convention -- κ² = 2/3 for every pair
    (``kappa2="isotropic"``), so the number is comparable with an AV's;
    ``kappa2="dipoles"`` uses the ensembles' per-pair κ² instead (the
    orientation-resolved efficiency ``RotamerEnsemble.fret_efficiencies``
    reports). Errors are ``error`` or ``error_fraction`` × distance.
    """
    if kappa2 not in ("isotropic", "dipoles"):
        raise ValueError("kappa2 must be 'isotropic' or 'dipoles'")
    from IMP.bff.representation.distance import fret_pair_geometry, fret_pair_efficiencies

    if distance_type not in ("RDAMean", "RDAMeanE", "Rmp"):
        raise ValueError(f"unknown distance_type {distance_type!r}")
    out: dict[str, dict[str, Any]] = {}
    for name1, name2 in pairs:
        e1, e2 = ensembles[name1], ensembles[name2]
        mu1 = getattr(e1, "mu", None) if kappa2 == "dipoles" else None
        mu2 = getattr(e2, "mu", None) if kappa2 == "dipoles" else None
        geometry = fret_pair_geometry(e1.points[:, :3], e1.points[:, 3], e2.points[:, :3], e2.points[:, 3], mu1, mu2)
        w = geometry["weight"]
        if distance_type == "RDAMean":
            value = float(np.sum(geometry["R"] * w))
        elif distance_type == "Rmp":
            value = float(np.linalg.norm(e1.mean_position - e2.mean_position))
        else:
            eff = fret_pair_efficiencies(geometry, forster_radius)
            mean_e = eff["static"]
            if mean_e <= 0:
                value = float(np.sum(geometry["R"] * w))
            elif mean_e >= 1:
                value = 0.0
            else:
                value = float(forster_radius * (1.0 / mean_e - 1.0) ** (1.0 / 6.0))
        err = float(error) if error is not None else float(error_fraction * value)
        out[f"{name1}_{name2}"] = {
            "position1_name": name1,
            "position2_name": name2,
            "distance_type": distance_type,
            "distance": value,
            "error_neg": err,
            "error_pos": err,
            "Forster_radius": float(forster_radius),
        }
    return out


def write_rotamer_fps(
    path: str | Path,
    positions: dict[str, dict[str, Any]],
    distances: dict[str, dict[str, Any]] | None = None,
    *,
    merge_into: str | Path | None = None,
    validate: bool = True,
) -> None:
    """Write (or merge into) an fps.json with rotamer (``R1``) positions.

    ``positions``/``distances`` are fps.json entries (see
    :func:`rotamer_position_payload`, :func:`distances_from_ensembles`). With
    ``merge_into`` the existing file's positions and distances are kept and
    the new ones added (same names overwrite). The payload is validated
    against the fps.json schema before writing.
    """
    from IMP.bff.io.fps import read_fps_json, write_fps_json

    all_positions: dict[str, Any] = {}
    all_distances: dict[str, Any] = {}
    score_sets: dict[str, Any] = {}
    extra: dict[str, Any] = {}
    if merge_into is not None:
        p0, d0, score_sets, extra = read_fps_json(merge_into)
        all_positions.update(p0)
        all_distances.update(d0)
    all_positions.update(positions)
    all_distances.update(distances or {})
    write_fps_json(path, all_positions, all_distances, score_sets or None, extra or None, validate=validate)


# --------------------------------------------------------------------------
# io
# --------------------------------------------------------------------------
"""Input/output helpers for rotamer libraries and protein frames."""

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


def rotamer_library_registry() -> dict[str, dict[str, Any]]:
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


def rotamer_library_metadata(library_name: str) -> dict[str, Any]:
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
    registry = rotamer_library_registry()
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

    metadata = rotamer_library_metadata(library_name)
    cutoff = metadata.get("cutoff")
    filename = _library_filename(metadata, cutoff)
    stem = filename.split("_cutoff")[0]

    # The FRETpredict library files (module data, data/rotamer_library) are the
    # canonical libraries: <stem>.pdb + <stem>_cutoff<N>.bcif (+ weights) for
    # each cutoff. They are tried first so that the *requested cutoff* is the
    # one loaded. The RMF templates under templates/rotamer hold only the
    # cutoff-30 clustering, so resolving every name to <stem>.rmf3 silently
    # returned the wrong library for cutoff10/cutoff20 names.
    if lib_dir is None:
        dcd = _registry_path().parent / f"{filename}.bcif"
        if dcd.exists() and dcd.with_name(f"{stem}.pdb").exists():
            return dcd

    template_dir = Path(lib_dir) if lib_dir is not None else Path(get_template_dir("rotamer"))
    candidates = [
        template_dir / f"{filename}.bcif",
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
                    f"the cutoff-{cutoff} library needs {filename}.bcif next to {stem}.pdb")
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
    for key, metadata in rotamer_library_registry().items():
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
    selector_resnames = selector_resnames_from_metadata(metadata)
    if not selector_resnames:
        return None
    dye_resname = next(iter(selector_resnames))
    match = re.match(r".+_(?P<linker>[A-Z]\d?[A-Z]R)$", str(metadata.get("name", "")))
    linker_resname = match.group("linker") if match else None
    if linker_resname is None:
        return [dye_resname] * len(atom_names)
    linker_atoms = {"CA", "HA", "C", "O", "C6", "H10", "H11", "S1", "C7", "C8", "H12", "C9", "O3", "N3", "C10", "O4", "C11", "H13", "H14", "C12", "H15", "H16", "C13", "H17", "H18", "C14", "H19", "H20", "C15", "H21", "H22", "N99", "H23", "N", "H", "HX2", "HX3"}
    return [linker_resname if name in linker_atoms else dye_resname for name in atom_names]


def selector_resnames_from_metadata(metadata: dict[str, Any]) -> set[str]:
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
        resnames.update(selector_resnames(metadata.get(key)))
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
    metadata = rotamer_library_metadata(library_name) if not explicit_path else _metadata_from_path(path)
    suffix = path.suffix.lower()
    if suffix in (".bcif", ".dcd"):
        # FRETpredict library set: <stem>.pdb (names, residues) + frames +
        # per-rotamer weights. The shipped libraries are BinaryCIF as of
        # 2026-08-19; .dcd still reads, because a user's own library may be
        # one, but it is not what this package stores.
        from IMP.bff.sampling import load_rotamer_library_dcd
        stem = path.stem.split("_cutoff")[0]
        pdb_path = path.with_name(f"{stem}.pdb")
        weights_path = path.with_name(f"{path.stem}_weights.txt")
        ref = load_rotamer_library_dcd(pdb_path, path, weights_path if weights_path.exists() else None)
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
    # Two C++ calls instead of about a dozen SWIG crossings per atom. The loop
    # this replaces built a decorator, read three coordinates, looked up a type
    # and walked two parents for every leaf -- tens of thousands of crossings
    # per frame to move a few kilobytes, and it was the single largest cost in
    # loading a rotamer library.
    packed = np.asarray(IMP.bff.hierarchy_atom_coordinates(hierarchy),
                        dtype=np.float64).reshape(-1, 4)
    coords = np.ascontiguousarray(packed[:, :3])
    residue_indices = packed[:, 3].astype(int).tolist()
    meta = list(IMP.bff.hierarchy_atom_metadata(hierarchy))
    atom_names = meta[0::4]
    atom_types = meta[1::4]
    resnames = meta[2::4]
    chain_ids = meta[3::4]
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


# --------------------------------------------------------------------------
# ensemble
# --------------------------------------------------------------------------
"""``RotamerEnsemble``: a screened rotamer library at a site, usable as an AV.

The fps.json position type ``R1`` (PRD-108): a FRETpredict-style 1:1 rotamer
library transformed into the residue's backbone frame and Boltzmann-screened
against the protein. Per rotamer the chromophore centre, the transition
dipole and the weight are kept -- and every atom, so nothing downstream has
to re-derive them -- which is what a FRET *rate distribution* over a pair of
labels needs (R_ij, κ²_ij, w_i·w_j), not just a mean position.

It subclasses :class:`IMP.bff.AccessibleVolume` with ``points`` =
(N, 4) centre + weight, so every AV helper in ``fret`` (``av_pair_statistics``,
``histogram_rda``, ``mean_fret_distance``, ...) accepts it unchanged; the pair
kernels live in :mod:`IMP.bff.representation.distance` (``fret_pair_geometry``,
``fret_pair_efficiencies``).
"""

#: fps.json ``simulation_type`` of a screened 1:1 rotamer library.
SIMULATION_TYPE_R1 = "R1"


# ---------------------------------------------------------------------------
# frame / library helpers (module level; RotamerFRET builds on them too)
# ---------------------------------------------------------------------------

def resolve_backbone_site(frame: dict[str, Any], chain: Optional[str], residue: int):
    """CA, N, C coordinates of ``(chain, residue)`` in a protein frame dict."""
    coords = np.asarray(frame["coords"], dtype=np.float64)
    atom_names = [str(name).upper() for name in frame["atom_names"]]
    chain_id = (chain or "").upper()
    chain_ids = [str(v).upper() for v in frame.get("chain_ids", [""] * len(atom_names))]
    residue_indices = [int(v) for v in frame.get("residue_indices", [-1] * len(atom_names))]
    matches: dict[str, np.ndarray] = {}
    for coord, atom_name, frame_chain, frame_residue in zip(coords, atom_names, chain_ids, residue_indices):
        if atom_name not in {"CA", "N", "C"} or atom_name in matches:
            continue
        if chain_id and frame_chain and frame_chain != chain_id:
            continue
        if frame_residue != -1 and frame_residue != residue:
            continue
        matches[atom_name] = coord
    missing = [name for name in ("CA", "N", "C") if name not in matches]
    if missing:
        raise ValueError(f"Missing backbone atoms {missing} for chain {chain or 'A'} residue {residue}")
    return matches["CA"], matches["N"], matches["C"]


def backbone_rotation(ca, n, c) -> np.ndarray:
    """Rows are the site frame axes: x along CA→N, y in the N–CA–C plane, z = x × y."""
    ca_v = np.asarray(ca, dtype=np.float64)
    x = np.asarray(n, dtype=np.float64) - ca_v
    x /= np.linalg.norm(x)
    yt = np.asarray(c, dtype=np.float64) - ca_v
    yt /= np.linalg.norm(yt)
    z = np.cross(x, yt)
    z /= np.linalg.norm(z)
    y = np.cross(z, x)
    return np.vstack([x, y, z])


def transform_library_to_site(coords: np.ndarray, ca, n, c) -> np.ndarray:
    """Library coordinates (n_rot, n_atoms, 3) into the backbone frame at CA."""
    rotation = backbone_rotation(ca, n, c)
    return np.tensordot(np.asarray(coords, dtype=np.float64), rotation, axes=([2], [0])) + np.asarray(ca, dtype=np.float64)


def selector_atom_indices(
    atom_names: Sequence[str],
    selector,
    resnames: Optional[Sequence[str]] = None,
) -> list[int]:
    """Indices of the atoms named by a FRETpredict selector (``'C7 and resname A48'``).

    The ``and resname X`` clause is **honoured** when ``resnames`` is supplied.
    It has to be: atom names repeat between the dye residue and its linker --
    ``C13`` is in both ``A48`` and ``C1R``, ``C9`` in both ``A35``/``T48`` and
    theirs -- and a selector that ignores the residue takes whichever comes
    first in the atom ordering. That is the dye today only by luck of the
    ordering; a library written linker-first would silently resolve the
    transition dipole to two linker atoms and raise nothing.

    Without ``resnames`` the clause cannot be checked and the first name match is
    returned, as before. Callers that have the residue names should pass them.

    A selector that matches nothing now raises naming **which** item failed. The
    previous behaviour raised only when *every* item failed, so a two-atom
    selector with one bad name returned a one-element list, and
    :meth:`RotamerEnsemble.from_site` then fell back to atoms 0 and 1 -- a
    silently wrong dipole.
    """
    items = [selector] if isinstance(selector, str) else list(selector)
    if not items:
        return []
    upper = [str(n).upper() for n in atom_names]
    upper_res = [str(r).upper() for r in resnames] if resnames is not None else None

    indices: list[int] = []
    for item in items:
        text = str(item)
        head, _, tail = text.partition(" and ")
        wanted = head.strip().upper()
        want_res = None
        match = re.search(r"\bresname\s+(\S+)", tail, flags=re.IGNORECASE)
        if match:
            want_res = match.group(1).strip().upper()

        found = None
        for i, name in enumerate(upper):
            if name != wanted:
                continue
            if want_res is not None and upper_res is not None:
                if i >= len(upper_res) or upper_res[i] != want_res:
                    continue
            found = i
            break
        if found is None:
            raise ValueError(
                f"Atom selector {text!r} matched no atom in the rotamer library"
                + (f" (residue {want_res} not found with that atom name)"
                   if want_res and upper_res is not None else "")
            )
        indices.append(found)
    return indices


def _library(library) -> dict:
    return library if isinstance(library, dict) else load_rotamer_library(str(library))


def _frame(structure, frame_index: int) -> dict:
    if isinstance(structure, dict):
        return structure
    frames = load_protein_frames(structure, max_frames=frame_index + 1)
    if frame_index >= len(frames):
        raise ValueError(f"{structure} has {len(frames)} frame(s), frame {frame_index} requested")
    return frames[frame_index]


# ---------------------------------------------------------------------------
# the ensemble
# ---------------------------------------------------------------------------

class RotamerEnsemble(States):
    """A rotamer library placed and screened at one labelling site (fps ``R1``).

    A **sibling** of :class:`IMP.bff.AccessibleVolume`, not a subclass of it. It
    used to inherit from the concrete AV, which is why distance code worked for
    rotamers -- by inheritance rather than by design -- and it then had to carry
    a grid it does not have, filling ``density``, ``grid_step`` and
    ``grid_shape`` with empty placeholders. No consumer ever read them: they all
    used ``points``, ``mean_position``, ``n_points`` or ``has_volume``, which is
    the :class:`IMP.bff.States` surface both representations share.

    From ``States``, which is C++: ``points`` (N, 4) chromophore centre +
    weight, ``attachment_point`` (CA), ``orientations`` (the per-rotamer
    transition dipoles), ``position_name``, ``params`` (carries
    ``simulation_type='R1'``, library, temperature, ...). Adds ``atoms``
    (N, n_atoms, 3) in the protein frame, ``atom_names``, ``resnames``,
    ``energies``, ``partition`` (Z) and the ``library`` name.

    ``mu`` is the name the rotamer code uses for the dipoles, and it *is*
    ``States.orientations`` -- a property over the one array, so a caller cannot
    set one and read a stale other. It is why this is a plain class rather than
    a dataclass: the states are a C++ value, and its constructor is the one that
    has to run.
    """

    def __init__(
        self,
        *,
        points,
        attachment_point,
        position_name: str = "",
        params: Optional[dict] = None,
        mu=None,
        atoms=None,
        atom_names: tuple = (),
        resnames: tuple = (),
        energies=None,
        partition: float = 0.0,
        library: str = "",
        chain: str = "",
        residue: int = 0,
    ):
        States.__init__(
            self, points=points, attachment_point=attachment_point,
            orientations=mu, position_name=position_name, params=params or {})
        self.atoms = (np.zeros((0, 0, 3)) if atoms is None
                      else np.asarray(atoms, dtype=np.float64))
        self.atom_names = tuple(atom_names)
        self.resnames = tuple(resnames)
        self.energies = (np.zeros(0) if energies is None
                         else np.asarray(energies, dtype=np.float64))
        self.partition = float(partition)
        self.library = str(library)
        self.chain = str(chain)
        self.residue = int(residue)

    @property
    def mu(self) -> np.ndarray:
        """(N, 3) per-rotamer transition dipoles -- ``States.orientations``."""
        return self.orientations

    @mu.setter
    def mu(self, value):
        self.orientations = value

    def __repr__(self) -> str:
        return (f"RotamerEnsemble({self.n_rotamers} rotamers, "
                f"{self.library!r}, Z={self.partition:.3g})")

    # -- construction ------------------------------------------------------
    @classmethod
    def from_site(
        cls,
        structure,
        chain: Optional[str],
        residue: int,
        library,
        *,
        temperature: float = 298.15,
        electrostatic: bool = False,
        potential: str = "lj",
        ignore_h: bool = True,
        sigma_scaling: float = 0.5,
        epsilon_scaling: float = 1.0,
        frame_index: int = 0,
        position_name: str = "",
    ) -> "RotamerEnsemble":
        """Place ``library`` on ``(chain, residue)`` of ``structure`` and screen it.

        ``structure`` is a PDB / multi-MODEL PDB / RMF path or a frame dict
        from ``load_protein_frames``; ``library`` a registry name
        (``'AlexaFluor 488 C1R cutoff30'``), a path, or a loaded library dict.
        Scoring is FRETpredict's (LJ or Gauss, optional Debye–Hückel; the
        labelled residue and hydrogens are not obstacles).
        """
        frame = _frame(structure, frame_index)
        lib = _library(library)
        ca, n, c = resolve_backbone_site(frame, chain, residue)
        rotamers = transform_library_to_site(lib["coords"], ca, n, c)
        score = compute_rotamer_score(
            rotamers,
            frame["coords"],
            frame["atom_names"],
            frame["resnames"],
            lib["atom_names"],
            lib.get("metadata", {}),
            lib.get("resnames"),
            protein_residue_indices=frame.get("residue_indices"),
            protein_chain_ids=frame.get("chain_ids"),
            site_residue=residue,
            site_chain=chain,
            rotamer_weights=lib.get("weights"),
            temperature=temperature,
            ignore_h=ignore_h,
            electrostatic=electrostatic,
            potential=potential,
            sigma_scaling=sigma_scaling,
            epsilon_scaling=epsilon_scaling,
        )
        metadata = dict(lib.get("metadata", {}) or {})
        names = list(lib["atom_names"])
        lib_res = lib.get("resnames")
        centre_idx = selector_atom_indices(names, metadata.get("r", []), lib_res)[0]
        mu_idx = selector_atom_indices(names, metadata.get("mu", []), lib_res)
        if len(mu_idx) >= 2:
            mu = rotamers[:, mu_idx[1], :] - rotamers[:, mu_idx[0], :]
        else:
            mu = rotamers[:, 1, :] - rotamers[:, 0, :]
        mu = mu / np.linalg.norm(mu, axis=1, keepdims=True)
        centres = rotamers[:, centre_idx, :]
        points = np.hstack([centres, score.weights[:, None]]).astype(np.float64)
        ca_v = np.asarray(ca, dtype=np.float64)
        return cls(
            points=points,
            attachment_point=ca_v.copy(),
            position_name=position_name or f"{chain or ''}{residue}",
            params={
                "simulation_type": SIMULATION_TYPE_R1,
                "library": metadata.get("library_name", metadata.get("name", str(library))),
                "chain": chain or "",
                "residue": int(residue),
                "temperature": float(temperature),
                "electrostatic": bool(electrostatic),
                "potential": potential,
                "ignore_h": bool(ignore_h),
                "sigma_scaling": float(sigma_scaling),
                "epsilon_scaling": float(epsilon_scaling),
                "partition": float(score.partition),
            },
            mu=mu,
            atoms=rotamers,
            atom_names=tuple(names),
            resnames=tuple(lib.get("resnames") or ()),
            energies=np.asarray(score.energies, dtype=np.float64),
            partition=float(score.partition),
            library=str(metadata.get("library_name", metadata.get("name", str(library)))),
            chain=chain or "",
            residue=int(residue),
        )

    # -- accessors ---------------------------------------------------------
    @property
    def centres(self) -> np.ndarray:
        """(N, 3) chromophore centres in the protein frame (Å)."""
        return self.points[:, :3]

    @property
    def weights(self) -> np.ndarray:
        """(N,) normalised Boltzmann × library weights."""
        return self.points[:, 3]

    @property
    def n_rotamers(self) -> int:
        return int(self.points.shape[0])

    # -- pair physics ------------------------------------------------------
    def pair_geometry(self, other: "RotamerEnsemble") -> dict:
        """R_ij, κ²_ij and w_i·w_j against another ensemble (or any AV: κ² = 2/3)."""
        mu_other = getattr(other, "mu", None)
        if mu_other is not None and np.asarray(mu_other).shape[0] != other.points.shape[0]:
            mu_other = None
        return fret_pair_geometry(self.centres, self.weights, other.points[:, :3], other.points[:, 3], self.mu, mu_other)

    def pair_distribution(self, other: "RotamerEnsemble", forster_radius: float, tau0: Optional[float] = None) -> dict:
        """The pair's FRET rate distribution and averages (see ``fret.distance.fret_pair_efficiencies``).

        ``forster_radius`` in Å for κ² = 2/3.
        """
        geometry = self.pair_geometry(other)
        out = fret_pair_efficiencies(geometry, forster_radius, tau0)
        out["kappa2"] = geometry["kappa2"]
        return out

    def fret_efficiencies(
        self,
        other: "RotamerEnsemble",
        forster_radius: Optional[float] = None,
        *,
        donor: Optional[str] = None,
        acceptor: Optional[str] = None,
    ) -> dict:
        """E_static, E_dynamic1, E_dynamic2 and ⟨κ²⟩ for this (donor) and ``other`` (acceptor).

        Give ``forster_radius`` (Å, κ² = 2/3) or the dye names -- then R0 is
        computed from the spectra at the pair's ⟨κ²⟩, as FRETpredict does.
        """
        geometry = self.pair_geometry(other)
        if forster_radius is None:
            if donor is None or acceptor is None:
                raise ValueError("give forster_radius (A) or donor and acceptor names")
            r0_nm = forster_radius_from_spectra(donor, acceptor, geometry["kappa2_avg"])
            # FRETpredict applies its k2-dependent R0 with the isotropic formula
            # 1/(1 + (2/3/k2)(r/R0)^6); fret_pair_efficiencies expects R0 at k2 = 2/3
            forster_radius = r0_nm * 10.0
            out = fret_pair_efficiencies(geometry, forster_radius)
            out["forster_radius_nm"] = r0_nm
            return out
        return fret_pair_efficiencies(geometry, forster_radius)


def rotamer_ensembles_from_fps(
    fps_json,
    structure,
    library_map: Optional[Dict[str, str]] = None,
    *,
    frame_index: int = 0,
    **kwargs,
) -> Dict[str, RotamerEnsemble]:
    """One :class:`RotamerEnsemble` per fps.json position that names a library.

    A position's ``rotamer_library`` / ``library`` field selects the library;
    ``library_map`` (position name → library name) overrides or supplies it
    for AV-only files. Positions without a library are skipped. ``kwargs`` go
    to :meth:`RotamerEnsemble.from_site`.
    """
    from IMP.bff.io.fps import read_fps_json
    # (was: from .fps import ...) -- now in this module

    positions, _distances, _score_sets, _extra = read_fps_json(fps_json)
    frame = _frame(structure, frame_index)
    library_map = dict(library_map or {})
    out: Dict[str, RotamerEnsemble] = {}
    for name, payload in positions.items():
        pos = RotamerPosition.from_payload(name, payload)
        lib_name = library_map.get(name) or pos.library
        if not lib_name:
            continue
        out[name] = RotamerEnsemble.from_site(
            frame, pos.chain, pos.residue, lib_name, position_name=name, **kwargs)
    return out


# --------------------------------------------------------------------------
# fret
# --------------------------------------------------------------------------
"""Rotamer-based FRET prediction compatible with FRETpredict workflows."""

_log = logging.getLogger(__name__)


@dataclass
class FRETFrameResult:
    """Per-frame FRET calculation result.

    Attributes
    ----------
    z : tuple[float, float]
        Donor and acceptor partition functions.
    k2 : float
        Weighted average orientation factor.
    estatic : float
        Static-regime efficiency.
    edynamic1 : float
        Dynamic1-regime efficiency.
    edynamic2 : float
        Dynamic2-regime efficiency.
    """

    z: tuple[float, float]
    k2: float
    estatic: float
    edynamic1: float
    edynamic2: float


def _weighted_average_sd_se(values: np.ndarray, weights: np.ndarray) -> tuple[float, float, float]:
    """Weighted mean, standard deviation and standard error.

    :func:`IMP.bff.weighted_average_sd_se`. Non-finite values and their weights
    are dropped and the rest renormalised, so a frame where the dye could not
    be placed contributes nothing rather than poisoning the mean.
    """
    out = IMP.bff.weighted_average_sd_se(
        np.ascontiguousarray(values, dtype=np.float64),
        np.ascontiguousarray(weights, dtype=np.float64))
    return float(out[0]), float(out[1]), float(out[2])


def _calculate_ws(z_values: np.ndarray) -> np.ndarray:
    """Per-frame weights from a pair of partition functions.

    :func:`IMP.bff.rotamer_frame_weights`. Uniform when every product is zero:
    a frame in which neither dye has an accessible conformer says nothing about
    the others.
    """
    z_values = np.asarray(z_values, dtype=np.float64)
    if z_values.shape == (2,):
        return np.array([1.0], dtype=np.float64)
    if z_values.ndim != 2 or z_values.shape[1] != 2:
        raise ValueError(f"Expected Z array with shape (n_frames, 2), got {z_values.shape}")
    return np.asarray(IMP.bff.rotamer_frame_weights(z_values), dtype=np.float64)


def _effective_fraction(weights: np.ndarray) -> float:
    """The effective number of contributing frames --
    :func:`IMP.bff.effective_frame_fraction`."""
    return float(IMP.bff.effective_frame_fraction(
        np.ascontiguousarray(weights, dtype=np.float64)))


class RotamerFRET:
    """Predict FRET efficiencies from rotamer libraries.

    Parameters
    ----------
    protein : pathlib.Path or str
        Protein PDB or RMF path.
    residues : list[int]
        Placement residue numbers.
    donor : str
        Donor dye name for R0 calculation.
    acceptor : str
        Acceptor dye name for R0 calculation.
    libname_1 : str
        Donor rotamer library name.
    libname_2 : str
        Acceptor rotamer library name.
    chains : list[str], optional
        Placement chain IDs.
    temperature : float
        Temperature in K.
    electrostatic : bool
        Include Debye-Huckel electrostatics.
    output_prefix : str
        Prefix for output files.
    fixed_R0 : bool
        Use a fixed R0 value.
    r0 : float
        Fixed R0 value in nm when ``fixed_R0`` is true.
    r0lib : pathlib.Path or str, optional
        Optional R0 data directory.
    z_cutoff : float
        Partition-function cutoff.
    calc_distr : bool
        Save distance and k2 distributions.
    verbose : bool
        Enable debug logging.

    Examples
    --------
    >>> fret = RotamerFRET("openHsp90.pdb", [452, 637], donor="AlexaFluor 594", acceptor="AlexaFluor 568")
    >>> fret.run()
    """

    def __init__(self, protein: str | Path, residues: list[int], **kwargs: Any) -> None:
        self.protein_path = str(protein)
        self.residues = list(residues)
        if len(self.residues) != 2:
            raise ValueError("The residue_list must contain exactly 2 residue numbers")

        self.chains = kwargs.get("chains", [None, None])
        self.donor = kwargs.get("donor", "AlexaFluor 488")
        self.acceptor = kwargs.get("acceptor", "AlexaFluor 594")
        self.libname_1 = kwargs.get("libname_1", "AlexaFluor 488 C1R cutoff30")
        self.libname_2 = kwargs.get("libname_2", "AlexaFluor 594 C1R cutoff30")
        self.r0lib = kwargs.get("r0lib", None)
        self.z_cutoff = float(kwargs.get("z_cutoff", 0.05))
        self.fixed_R0 = bool(kwargs.get("fixed_R0", False))
        self.r0 = float(kwargs.get("r0", 5.4) or 5.4)
        self.temperature = float(kwargs.get("temperature", 300.0))
        self.electrostatic = bool(kwargs.get("electrostatic", False))
        self.potential = kwargs.get("potential", "lj")
        self.sigma_scaling = float(kwargs.get("sigma_scaling", 0.5))
        self.epsilon_scaling = float(kwargs.get("epsilon_scaling", 1.0))
        self.ignore_h = bool(kwargs.get("ign_H", True))
        self.output_prefix = kwargs.get("output_prefix", "res")
        self.weights = kwargs.get("weights", None)
        self.user_weights = kwargs.get("user_weights", None)
        self.filter_stdev = float(kwargs.get("filter_stdev", 0.02))
        self.verbose = bool(kwargs.get("verbose", False))
        self.calc_distr = bool(kwargs.get("calc_distr", False))
        self.max_frames = kwargs.get("max_frames", None)

        self.dr = 0.05
        self.rmin = -5.0
        self.rmax = float(kwargs.get("rmax", 20)) * 2.0 - self.rmin
        self.nr = int(round((self.rmax - self.rmin) / self.dr, 0) + 1)
        self.rax = np.linspace(self.rmin, self.rmax, self.nr)

        # Module logger; nothing is written unless the caller asks for a
        # log file (``log_file=``). The previous ``logging.basicConfig(filename="log")``
        # silently created a ``log`` file in the working directory on every
        # construction and reconfigured the root logger of the host process.
        log_file = kwargs.get("log_file", None)
        _log.setLevel(logging.DEBUG if self.verbose else logging.INFO)
        if log_file:
            handler = logging.FileHandler(str(log_file))
            handler.setFormatter(logging.Formatter("%(levelname)s:%(name)s:%(message)s"))
            _log.addHandler(handler)

        self.lib_1 = load_rotamer_library(self.libname_1)
        self.lib_2 = load_rotamer_library(self.libname_2)
        self.frames = load_protein_frames(self.protein_path, max_frames=self.max_frames)
        self.z_values: np.ndarray | None = None
        self.k2_values: np.ndarray | None = None
        self.estatic_values: np.ndarray | None = None
        self.edynamic1_values: np.ndarray | None = None
        self.edynamic2_values: np.ndarray | None = None
        self.distance_distributions: np.ndarray | None = None

    def _resolve_site(self, frame: dict[str, Any], chain: str | None, residue: int) -> dict[str, np.ndarray]:
        """CA/N/C of the site (see ``ensemble.resolve_backbone_site``)."""
        ca, n, c = resolve_backbone_site(frame, chain, residue)
        return {"CA": ca, "N": n, "C": c}

    def _transform_library(self, library: dict[str, Any], frame: dict[str, Any], chain: str | None, residue: int) -> np.ndarray:
        """Library coordinates in the site's backbone frame (``ensemble.transform_library_to_site``)."""
        ca, n, c = resolve_backbone_site(frame, chain, residue)
        return transform_library_to_site(library["coords"], ca, n, c)

    def _ensemble(self, library: dict[str, Any], frame: dict[str, Any], chain: str | None, residue: int) -> RotamerEnsemble:
        """The screened :class:`RotamerEnsemble` of one library at one site of one frame."""
        return RotamerEnsemble.from_site(
            frame, chain, residue, library,
            temperature=self.temperature, electrostatic=self.electrostatic,
            potential=self.potential, ignore_h=self.ignore_h,
            sigma_scaling=self.sigma_scaling, epsilon_scaling=self.epsilon_scaling)

    def _frame_fret(self, frame: dict[str, Any]) -> FRETFrameResult:
        """FRET quantities of one protein frame from the two screened ensembles.

        Static / dynamic1 / dynamic2 come from ``fret.distance.fret_pair_efficiencies``
        with R0 either fixed or computed from the spectra at this frame's ⟨κ²⟩
        (FRETpredict's convention).
        """
        donor = self._ensemble(self.lib_1, frame, self.chains[0], self.residues[0])
        acceptor = self._ensemble(self.lib_2, frame, self.chains[1], self.residues[1])
        geometry = donor.pair_geometry(acceptor)
        k2_avg = geometry["kappa2_avg"]
        if not self.fixed_R0:
            self.r0 = forster_radius_from_spectra(self.donor, self.acceptor, k2_avg,
                                                  library_cif="" if self.r0lib is None else str(self.r0lib))
            if self.r0 == 0:
                return FRETFrameResult((donor.partition, acceptor.partition), float("nan"), float("nan"), float("nan"), float("nan"))
        eff = fret_pair_efficiencies(geometry, float(self.r0) * 10.0)   # r0 in nm, geometry in A
        return FRETFrameResult(
            (donor.partition, acceptor.partition),
            k2_avg,
            eff["static"],
            eff["dynamic1"],
            eff["dynamic2"],
        )

    def trajectory_analysis(self) -> None:
        """Calculate FRET efficiencies for all protein frames.

        Returns
        -------
        None
        """
        n_frames = len(self.frames)
        z_values = np.empty((n_frames, 2), dtype=np.float64)
        k2_values = np.full(n_frames, np.nan, dtype=np.float64)
        estatic_values = np.full(n_frames, np.nan, dtype=np.float64)
        edynamic1_values = np.full(n_frames, np.nan, dtype=np.float64)
        edynamic2_values = np.full(n_frames, np.nan, dtype=np.float64)
        distance_distributions = np.zeros((n_frames, self.nr), dtype=np.float64) if self.calc_distr else None

        for frame_index, frame in enumerate(self.frames):
            result = self._frame_fret(frame)
            z_values[frame_index] = result.z
            if result.z[0] <= self.z_cutoff or result.z[1] <= self.z_cutoff:
                continue
            k2_values[frame_index] = result.k2
            estatic_values[frame_index] = result.estatic
            edynamic1_values[frame_index] = result.edynamic1
            edynamic2_values[frame_index] = result.edynamic2

        self.z_values = z_values
        self.k2_values = k2_values
        self.estatic_values = estatic_values
        self.edynamic1_values = edynamic1_values
        self.edynamic2_values = edynamic2_values
        self.distance_distributions = distance_distributions

    def save(self, reweight_output_prefix: str | None = None) -> None:
        """Save calculated FRET quantities to files.

        Parameters
        ----------
        reweight_output_prefix : str, optional
            Optional output prefix for summary files.

        Returns
        -------
        None
        """
        if self.z_values is None or self.k2_values is None:
            self.trajectory_analysis()
        assert self.z_values is not None
        assert self.k2_values is not None
        assert self.estatic_values is not None
        assert self.edynamic1_values is not None
        assert self.edynamic2_values is not None

        prefix = reweight_output_prefix or self.output_prefix
        r1, r2 = self.residues
        np.savetxt(f"{prefix}-Z-{r1}-{r2}.dat", self.z_values)
        np.savetxt(f"{prefix}-w_s-{r1}-{r2}.dat", _calculate_ws(self.z_values))
        np.savetxt(f"{prefix}-k2-{r1}-{r2}.dat", self.k2_values)
        np.savetxt(f"{prefix}-Es-{r1}-{r2}.dat", self.estatic_values)
        np.savetxt(f"{prefix}-Ed1-{r1}-{r2}.dat", self.edynamic1_values)
        np.savetxt(f"{prefix}-Ed2-{r1}-{r2}.dat", self.edynamic2_values)

        weights = np.ones_like(self.k2_values, dtype=np.float64)
        if self.user_weights is not None:
            weights = np.asarray(self.user_weights, dtype=np.float64)
        if weights.size != self.k2_values.size:
            raise ValueError(f"Weights array has size {weights.size} whereas the number of frames is {self.k2_values.size}")
        weights = weights / np.sum(weights)

        finite = np.isfinite(self.k2_values)
        labels = ["k2", "Estatic", "Edynamic1", "Edynamic2"]
        if self.k2_values.size == 1:
            rows = [
                (self.k2_values[0], np.nan, np.nan),
                (self.estatic_values[0], np.nan, np.nan),
                (self.edynamic1_values[0], np.nan, np.nan),
                (self.edynamic2_values[0], np.nan, np.nan),
            ]
        else:
            rows = [
                _weighted_average_sd_se(self.k2_values[finite], weights[finite]),
                _weighted_average_sd_se(self.estatic_values[finite], weights[finite]),
                _weighted_average_sd_se(self.edynamic1_values[finite], weights[finite]),
                _weighted_average_sd_se(self.edynamic2_values[finite], weights[finite]),
            ]
        # Written as a labelled text table rather than a pickled DataFrame: a
        # .pkl is unreadable without the library that wrote it, and IMP.bff
        # carries no dependency beyond what IMP itself brings.
        summary = np.asarray(rows, dtype=np.float64)
        np.savetxt(
            f"{prefix}-data-{r1}-{r2}.dat",
            summary,
            header="quantity Average SD SE\n" + " ".join(labels),
            comments="# ",
        )

    def reweight(self, **kwargs: Any) -> None:
        """Reweight saved FRET quantities.

        Parameters
        ----------
        boltzmann_weights : bool
            Use partition-function weights.
        user_weights : numpy.ndarray, optional
            User-provided per-frame weights.
        reweight_output_prefix : str
            Output prefix.
        **kwargs : dict
            Reweighting options.

        Returns
        -------
        None
        """
        prefix = kwargs.get("reweight_output_prefix", self.output_prefix)
        r1, r2 = self.residues
        if kwargs.get("boltzmann_weights", False):
            z_values = np.loadtxt(f"{self.output_prefix}-Z-{r1}-{r2}.dat")
            self.weights = _calculate_ws(z_values)
        elif kwargs.get("user_weights") is not None:
            self.user_weights = np.asarray(kwargs["user_weights"], dtype=np.float64)

        k2 = np.atleast_1d(np.loadtxt(f"{self.output_prefix}-k2-{r1}-{r2}.dat"))
        estatic = np.atleast_1d(np.loadtxt(f"{self.output_prefix}-Es-{r1}-{r2}.dat"))
        edynamic1 = np.atleast_1d(np.loadtxt(f"{self.output_prefix}-Ed1-{r1}-{r2}.dat"))
        edynamic2 = np.atleast_1d(np.loadtxt(f"{self.output_prefix}-Ed2-{r1}-{r2}.dat"))
        weights = np.ones_like(k2, dtype=np.float64)
        if self.weights is not None:
            weights = np.asarray(self.weights, dtype=np.float64)
        if self.user_weights is not None:
            user_weights = np.asarray(self.user_weights, dtype=np.float64)
            if user_weights.size != k2.size:
                raise ValueError(f"Weights array has size {user_weights.size} whereas the number of frames is {k2.size}")
            weights = weights * user_weights
        weights = weights / np.sum(weights)
        finite = np.isfinite(k2)
        rows = [
            _weighted_average_sd_se(k2[finite], weights[finite]),
            _weighted_average_sd_se(estatic[finite], weights[finite]),
            _weighted_average_sd_se(edynamic1[finite], weights[finite]),
            _weighted_average_sd_se(edynamic2[finite], weights[finite]),
        ]
        np.savetxt(
            f"{prefix}-data-{r1}-{r2}.dat",
            np.asarray(rows, dtype=np.float64),
            header="quantity Average SD SE\nk2 Estatic Edynamic1 Edynamic2",
            comments="# ",
        )

    def run(self) -> None:
        """Run trajectory analysis and save output files.

        Returns
        -------
        None
        """
        self.trajectory_analysis()
        self.save()
        _log.debug("Done")


