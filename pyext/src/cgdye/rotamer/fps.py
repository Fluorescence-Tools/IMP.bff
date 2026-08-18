"""fps.json helpers for rotamer FRET."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from IMP.bff.io.fps import read_fps_json


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
    from IMP.bff.cgdye.rotamer.fret import RotamerFRET

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

