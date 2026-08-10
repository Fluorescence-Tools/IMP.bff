"""fps.json helpers for rotamer FRET."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from IMP.bff.fps import read_fps_json


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
