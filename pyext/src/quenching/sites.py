"""Per-residue geometry for the PET model: where a residue sticks and quenches.

A residue enters the quenching model at two different points, and they are not
the same point in space:

*slow centre*
    The coarse side-chain position -- CB, else CA -- that drives the unspecific
    stickiness. Stickiness is a bulk property of the side chain, so its centroid
    is the right handle.
*quench centre*
    The centroid of the residue's PET-active atoms. Electron transfer happens at
    the redox-active moiety, and for tryptophan the indole ring sits ~3.3 A from
    CB -- which matters when contact distances are a few Angstrom.

Moved here from QuEst (``quest/core/dye_diffusion.py``) by PRD-109.
"""

from __future__ import annotations

from collections import OrderedDict
from typing import NamedTuple, Sequence

import numpy as np

from .pet import QUENCHER_ATOMS, normalize_amino_acid_quenching

__all__ = [
    "ResidueSites",
    "residue_sites",
    "slow_factors_for_residues",
    "quenching_rates_for_residues",
    "quench_radii_for_residues",
]

_FALLBACK = {
    "slow_factor": 1.0,
    "kQ": 0.0,
    "quench_radius": None,
    "quench_atoms": ("CB",),
}


class ResidueSites(NamedTuple):
    """One slow centre and one quench centre per residue, plus its type."""

    slow_centers: np.ndarray
    quench_centers: np.ndarray
    residue_names: list

    def __len__(self) -> int:
        return len(self.residue_names)


def _name(value) -> str:
    if isinstance(value, bytes):
        return value.decode("ascii", errors="ignore").strip().upper()
    return str(value).strip().upper()


def residue_sites(atoms, quenching_table=None) -> ResidueSites:
    """Group *atoms* by residue and locate its slow and quench centres.

    :param atoms: structured array with ``chain``, ``res_id``, ``res_name``,
        ``atom_name`` and ``coord`` fields.
    :param quenching_table: the per-residue interaction table, whose
        ``quench_atoms`` decide which atoms define each quench centre. The
        defaults from :data:`IMP.bff.QUENCHER_ATOMS` are used when omitted.

    Residues are keyed by ``(chain, res_id, res_name)``. **Keying on ``res_id``
    alone is wrong** and was a real defect in QuEst: residue numbers restart per
    chain, so in a homodimer every number occurs twice and two residues' atoms
    were folded into one centre.
    """
    table = normalize_amino_acid_quenching(quenching_table)

    def wanted_atoms(residue_name):
        params = table.get(residue_name)
        if params is None:
            return frozenset(QUENCHER_ATOMS.get(residue_name, ("CB",)))
        return frozenset(params.get("quench_atoms") or ("CB",))

    by_residue = OrderedDict()
    for index, atom in enumerate(atoms):
        key = (_name(atom["chain"]), int(atom["res_id"]), _name(atom["res_name"]))
        by_residue.setdefault(key, []).append(index)

    slow_centers = []
    quench_centers = []
    residue_names = []
    for (_chain, _res_id, res_name), indices in by_residue.items():
        block = atoms[np.asarray(indices, dtype=np.int64)]
        atom_names = [_name(n) for n in block["atom_name"]]
        coords = np.asarray(block["coord"], dtype=np.float64)
        residue_name = _name(res_name)

        if "CB" in atom_names:
            selected = atom_names.index("CB")
        elif "CA" in atom_names:
            selected = atom_names.index("CA")
        else:
            selected = 0
        slow_center = coords[selected]

        matched = [i for i, n in enumerate(atom_names) if n in wanted_atoms(residue_name)]
        quench_center = coords[matched].mean(axis=0) if matched else slow_center

        slow_centers.append(slow_center)
        quench_centers.append(quench_center)
        residue_names.append(residue_name)

    return ResidueSites(
        np.asarray(slow_centers, dtype=np.float64).reshape(-1, 3),
        np.asarray(quench_centers, dtype=np.float64).reshape(-1, 3),
        residue_names,
    )


def _lookup(table, residue_name):
    return table.get(_name(residue_name), _FALLBACK)


def slow_factors_for_residues(residue_names: Sequence[str], table) -> np.ndarray:
    """The stickiness factor of each residue, in ``residue_names`` order."""
    return np.asarray(
        [_lookup(table, r)["slow_factor"] for r in residue_names], dtype=np.float64
    )


def quenching_rates_for_residues(residue_names: Sequence[str], table) -> np.ndarray:
    """The quenching rate (1/ns) of each residue."""
    return np.asarray(
        [_lookup(table, r)["kQ"] for r in residue_names], dtype=np.float64
    )


def quench_radii_for_residues(
    residue_names: Sequence[str], table, critical_distance: float = 0.0
) -> np.ndarray:
    """The contact radius of each residue, inheriting *critical_distance*.

    A ``quench_radius`` of ``None`` in the table means "use the model-wide
    critical distance", which is how a project sets one radius for everything
    and overrides it per residue type where it matters.
    """
    default = float(critical_distance or 0.0)
    radii = []
    for residue in residue_names:
        value = _lookup(table, residue).get("quench_radius")
        radii.append(default if value is None else float(value))
    return np.asarray(radii, dtype=np.float64)
