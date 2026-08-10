"""Rotamer-protein interaction scoring."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

LJ_PARAMETERS: dict[str, dict[str, float]] = {
    "C": {"p_Rmin2": 2.02446316, "eps": -0.06394724},
    "N": {"p_Rmin2": 1.89285714, "eps": -0.15428571},
    "O": {"p_Rmin2": 1.693, "eps": -0.12642017},
    "S": {"p_Rmin2": 2.1, "eps": -0.47},
    "H": {"p_Rmin2": 0.98357778, "eps": -0.03466645},
}

_GAS_CONSTANT = 1.9858775e-3


@dataclass
class RotamerScoreResult:
    """Result of rotamer-protein scoring.

    Attributes
    ----------
    weights : numpy.ndarray
        Normalized Boltzmann weights.
    partition : float
        Partition function before normalization.
    energies : numpy.ndarray
        Raw potential plus electrostatic energies.
    """

    weights: np.ndarray
    partition: float
    energies: np.ndarray


def _atom_type(atom_name: str) -> str:
    """Return a coarse atom type from an atom name.

    Parameters
    ----------
    atom_name : str
        Atom name.

    Returns
    -------
    str
        One-letter atom type.
    """
    name = str(atom_name).strip().upper()
    for atom in ("S", "N", "O", "H", "C"):
        if name.startswith(atom):
            return atom
    return "C"


def _selector_matches(selector: str, atom_name: str, resname: str | None = None) -> bool:
    """Return whether one selector matches an atom.

    Parameters
    ----------
    selector : str
        FRETpredict-style selector.
    atom_name : str
        Atom name.
    resname : str, optional
        Atom residue name.

    Returns
    -------
    bool
        True if the selector matches the atom.
    """
    parts = [part.strip() for part in str(selector).split(" and ")]
    atom = parts[0].strip().upper() if parts else str(selector).strip().upper()
    if atom.startswith("name "):
        atom = atom.split(None, 1)[1].upper()
    if atom != str(atom_name).upper():
        return False
    selector_resname = None
    for part in parts[1:]:
        tokens = part.split()
        if len(tokens) >= 2 and tokens[0].lower() == "resname":
            selector_resname = tokens[1].upper()
    return selector_resname is None or (resname is not None and str(resname).upper() == selector_resname)


def _selector_resnames(selector: str | list[str] | None) -> set[str]:
    """Extract residue names from FRETpredict-style selectors.

    Parameters
    ----------
    selector : str or list of str or None
        Selector such as ``C7 and resname T48``.

    Returns
    -------
    set[str]
        Residue names appearing in the selector.
    """
    if selector is None:
        return set()
    if isinstance(selector, str):
        selectors = [selector]
    else:
        selectors = list(selector)
    resnames: set[str] = set()
    for item in selectors:
        for part in str(item).split(" and "):
            tokens = part.strip().split()
            if len(tokens) >= 2 and tokens[0].lower() == "resname":
                resnames.add(tokens[1].upper())
    return resnames


def _selector_atom_names(selector: str | list[str] | None, atom_names: list[str] | None = None, resnames: list[str] | None = None) -> set[str]:
    """Extract atom names from FRETpredict-style selectors.

    Parameters
    ----------
    selector : str or list of str or None
        Selector such as ``C7 and resname T48``.
    atom_names : list of str, optional
        Atom names aligned with ``resnames``.
    resnames : list of str, optional
        Residue names aligned with atom names. When provided, selectors are
        matched against their ``resname`` terms.

    Returns
    -------
    set[str]
        Atom names.
    """
    if selector is None:
        return set()
    if isinstance(selector, str):
        selectors = [selector]
    else:
        selectors = list(selector)
    selector_resnames = _selector_resnames(selectors)
    if selector_resnames and len(selector_resnames) > 1:
        return set()
    names: set[str] = set()
    if resnames is None:
        for item in selectors:
            parts = [part.strip() for part in str(item).split(" and ")]
            atom = parts[0].strip().upper() if parts else str(item).strip().upper()
            if atom.startswith("name "):
                atom = atom.split(None, 1)[1].upper()
            names.add(atom)
    else:
        if atom_names is None:
            raise ValueError("atom_names is required when resnames is provided")
        for atom, resname in zip(atom_names, resnames):
            if any(_selector_matches(item, atom, resname) for item in selectors):
                names.add(str(atom).upper())
    return names


def _site_mask(
    atom_names: list[str],
    residue_indices: list[int] | None = None,
    site_residue: list[int] | int | None = None,
    site_chain: str | None = None,
    chain_ids: list[str] | None = None,
    mask_backbone: bool = True,
) -> np.ndarray:
    """Return a mask for placement-residue atoms.

    Parameters
    ----------
    atom_names : list of str
        Atom names.
    residue_indices : list of int, optional
        Residue indices aligned with atom names.
    site_residue : list of int, int, or None
        Placement residue indices used to exclude site atoms.
    site_chain : str, optional
        Placement chain ID used with ``site_residue``.
    chain_ids : list of str, optional
        Chain IDs aligned with atom names.
    mask_backbone : bool
        Also mask backbone atoms by name.

    Returns
    -------
    numpy.ndarray
        Boolean mask.
    """
    mask = np.zeros(len(atom_names), dtype=bool)
    if residue_indices is None:
        residue_indices = [-1] * len(atom_names)
    wanted = set(residue_indices) if isinstance(site_residue, list) else {site_residue}
    for i, (name, residue_index) in enumerate(zip(atom_names, residue_indices)):
        frame_chain = chain_ids[i] if chain_ids is not None and i < len(chain_ids) else None
        same_chain = site_chain is None or frame_chain is None or str(frame_chain).upper() == str(site_chain).upper()
        if site_residue is not None and same_chain and residue_index in wanted:
            mask[i] = True
            continue
        if mask_backbone and name.upper() in {"CA", "C", "N", "O"}:
            mask[i] = True
    return mask


def _hydrogen_mask(atom_names: list[str]) -> np.ndarray:
    """Return a mask for hydrogen atoms.

    Parameters
    ----------
    atom_names : list[str]
        Atom names.

    Returns
    -------
    numpy.ndarray
        Boolean mask.
    """
    return np.array([_atom_type(name) == "H" for name in atom_names], dtype=bool)


def _protein_charge_mask(atom_names: list[str], resnames: list[str]) -> np.ndarray:
    """Return protein charge codes for Debye-Huckel scoring.

    Parameters
    ----------
    atom_names : list[str]
        Atom names.
    resnames : list[str]
        Residue names.

    Returns
    -------
    numpy.ndarray
        Charge code array.
    """
    charges = np.zeros(len(atom_names), dtype=float)
    for i, (atom, resname) in enumerate(zip(atom_names, resnames)):
        name = atom.upper()
        res = resname.upper()
        if (name == "CZ" and res == "ARG") or (name == "NZ" and res == "LYS"):
            charges[i] = 1.0
        elif (name == "CG" and res == "ASP") or (name == "CD" and res == "GLU"):
            charges[i] = -1.0
        elif name in {"ND1", "NE2"} and res == "HIS":
            charges[i] = 0.25
    return charges


def _rotamer_charge_mask(
    atom_names: list[str],
    positive: list[str] | None = None,
    negative: list[str] | None = None,
    resnames: list[str] | None = None,
) -> np.ndarray:
    """Return rotamer charge codes for Debye-Huckel scoring.

    Parameters
    ----------
    atom_names : list[str]
        Atom names.
    positive : list of str, optional
        Positively charged atom selectors.
    negative : list of str, optional
        Negatively charged atom selectors.
    resnames : list of str, optional
        Residue names aligned with atom names.

    Returns
    -------
    numpy.ndarray
        Charge code array.
    """
    positive = [] if positive is None else list(positive)
    negative = [] if negative is None else list(negative)
    if len(_selector_resnames(positive)) > 1:
        positive = []
    if len(_selector_resnames(negative)) > 1:
        negative = []
    _selector_atom_names(positive, atom_names=atom_names, resnames=resnames)
    _selector_atom_names(negative, atom_names=atom_names, resnames=resnames)
    charges = np.zeros(len(atom_names), dtype=float)
    for i, name in enumerate(atom_names):
        resname = None if resnames is None or i >= len(resnames) else resnames[i]
        is_negative = any(_selector_matches(selector, name, resname) for selector in negative)
        is_positive = any(_selector_matches(selector, name, resname) for selector in positive)
        if is_negative:
            charges[i] = -1.0
        elif is_positive:
            charges[i] = 0.5
    return charges


def _scaled_parameters(
    atom_names: list[str],
    sigma_scaling: float,
    epsilon_scaling: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return scaled LJ parameters for atom names.

    Parameters
    ----------
    atom_names : list[str]
        Atom names.
    sigma_scaling : float
        Sigma scaling factor.
    epsilon_scaling : float
        Epsilon scaling factor.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        Scaled ``Rmin2`` and ``epsilon`` arrays.
    """
    eps = np.array([LJ_PARAMETERS.get(_atom_type(name), LJ_PARAMETERS["C"])["eps"] for name in atom_names], dtype=float)
    rmin2 = np.array(
        [LJ_PARAMETERS.get(_atom_type(name), LJ_PARAMETERS["C"])["p_Rmin2"] for name in atom_names],
        dtype=float,
    )
    return rmin2 * sigma_scaling, eps * epsilon_scaling


def compute_rotamer_score(
    rotamer_coords: np.ndarray,
    protein_coords: np.ndarray,
    protein_atom_names: list[str],
    protein_resnames: list[str],
    rotamer_atom_names: list[str],
    rotamer_metadata: dict[str, Any] | None = None,
    rotamer_resnames: list[str] | None = None,
    protein_residue_indices: list[int] | None = None,
    protein_chain_ids: list[str] | None = None,
    site_residue: list[int] | int | None = None,
    site_chain: str | None = None,
    rotamer_weights: np.ndarray | None = None,
    temperature: float = 300.0,
    ignore_h: bool = True,
    electrostatic: bool = False,
    potential: str = "lj",
    sigma_scaling: float = 0.5,
    epsilon_scaling: float = 1.0,
) -> RotamerScoreResult:
    """Compute Boltzmann weights for a rotamer ensemble.

    Parameters
    ----------
    rotamer_coords : numpy.ndarray
        Rotamer coordinates with shape ``(n_rotamers, n_atoms, 3)`` in Å.
    protein_coords : numpy.ndarray
        Protein coordinates with shape ``(n_atoms, 3)`` in Å.
    protein_atom_names : list of str
        Protein atom names.
    protein_resnames : list of str
        Protein residue names.
    rotamer_atom_names : list of str
        Rotamer atom names.
    rotamer_metadata : dict, optional
        Rotamer metadata containing charged atom selectors.
    rotamer_resnames : list of str, optional
        Rotamer residue names aligned with atom names.
    protein_residue_indices : list of int, optional
        Protein residue numbers aligned with atom names.
    protein_chain_ids : list of str, optional
        Protein chain IDs aligned with atom names.
    rotamer_weights : numpy.ndarray, optional
        External rotamer population weights.
    temperature : float
        Temperature in K.
    ignore_h : bool
        Ignore hydrogen atoms.
    electrostatic : bool
        Include Debye-Huckel electrostatics.
    potential : {'lj', 'gauss'}
        Potential type.
    sigma_scaling : float
        LJ sigma scaling.
    epsilon_scaling : float
        LJ epsilon scaling.
    site_residue : list of int, int, or None
        Placement residue indices used to exclude site atoms.
    site_chain : str, optional
        Placement chain ID used with ``site_residue``.

    Returns
    -------
    RotamerScoreResult
        Normalized weights, partition function, and raw energies.
    """
    metadata = rotamer_metadata or {}
    protein_coords = np.asarray(protein_coords, dtype=np.float64)
    rotamer_coords = np.asarray(rotamer_coords, dtype=np.float64)

    protein_h = _hydrogen_mask(protein_atom_names)
    protein_mask = ~(protein_h if ignore_h else np.zeros(len(protein_atom_names), dtype=bool))
    protein_mask &= ~_site_mask(
        protein_atom_names,
        protein_residue_indices,
        site_residue=site_residue,
        site_chain=site_chain,
        chain_ids=protein_chain_ids,
        mask_backbone=False,
    )

    rotamer_backbone = _site_mask(rotamer_atom_names)
    rotamer_h = _hydrogen_mask(rotamer_atom_names)
    rotamer_mask = ~(rotamer_backbone | (rotamer_h if ignore_h else np.zeros(len(rotamer_atom_names), dtype=bool)))

    protein_idx = np.flatnonzero(protein_mask)
    rotamer_idx = np.flatnonzero(rotamer_mask)
    if protein_idx.size == 0 or rotamer_idx.size == 0:
        weights = np.ones(rotamer_coords.shape[0], dtype=np.float64)
        weights /= weights.sum()
        return RotamerScoreResult(weights=weights, partition=1.0, energies=np.zeros(rotamer_coords.shape[0]))

    protein_types = [_atom_type(name) for name in np.asarray(protein_atom_names)[protein_idx]]
    rotamer_types = [_atom_type(name) for name in np.asarray(rotamer_atom_names)[rotamer_idx]]
    protein_rmin2, protein_eps = _scaled_parameters(protein_types, sigma_scaling, epsilon_scaling)
    rotamer_rmin2, rotamer_eps = _scaled_parameters(rotamer_types, sigma_scaling, epsilon_scaling)

    eps_ij = np.sqrt(np.multiply.outer(rotamer_eps, protein_eps))
    rmin_ij = np.add.outer(rotamer_rmin2, protein_rmin2)

    protein_positions = protein_coords[protein_idx]
    q_rotamer = np.zeros(rotamer_idx.size, dtype=float)
    q_protein = np.zeros(protein_idx.size, dtype=float)

    if electrostatic:
        q_rotamer = _rotamer_charge_mask(
            [rotamer_atom_names[i] for i in rotamer_idx],
            metadata.get("positive"),
            metadata.get("negative"),
            None if rotamer_resnames is None else [rotamer_resnames[i] for i in rotamer_idx],
        )
        q_protein = _protein_charge_mask(
            [protein_atom_names[i] for i in protein_idx],
            [protein_resnames[i] for i in protein_idx],
        )

    pot_energy = np.zeros(rotamer_coords.shape[0], dtype=np.float64)
    dh_energy = np.zeros(rotamer_coords.shape[0], dtype=np.float64)

    for i, coords in enumerate(rotamer_coords):
        distances = np.linalg.norm(coords[rotamer_idx, None, :] - protein_positions[None, :, :], axis=2)
        if electrostatic:
            mask = (q_rotamer[:, None] * q_protein[None, :] != 0.0) & (distances < 20.0)
            if np.any(mask):
                ri, pj = np.nonzero(mask)
                dh_energy[i] += np.sum(q_rotamer[ri] * q_protein[pj] * 7.0 / distances[ri, pj] * np.exp(-distances[ri, pj] / 10.0))

        if potential == "lj":
            mask = distances < 10.0
            if np.any(mask):
                ratio = np.power(rmin_ij[mask] / distances[mask], 6)
                pot_energy[i] += np.sum(eps_ij[mask] * (ratio * ratio - 2.0 * ratio))
        elif potential == "gauss":
            mask = distances < 10.0
            if np.any(mask):
                ratio = np.power(distances[mask] / rmin_ij[mask], 2)
                pot_energy[i] += np.sum(eps_ij[mask] * np.exp(-0.5 * ratio))
        else:
            raise ValueError(f"Unknown potential {potential!r}")

    boltzmann = np.exp(-pot_energy / (_GAS_CONSTANT * temperature) - dh_energy)
    library_weights = np.asarray(rotamer_weights if rotamer_weights is not None else metadata.get("weights", np.ones(rotamer_coords.shape[0])), dtype=np.float64)
    if library_weights.shape != boltzmann.shape:
        library_weights = np.ones_like(boltzmann)
    boltzmann *= library_weights
    boltzmann = np.nan_to_num(boltzmann, nan=0.0, posinf=0.0, neginf=0.0)
    partition = float(np.sum(boltzmann))
    if partition <= 0.0:
        weights = np.ones(rotamer_coords.shape[0], dtype=np.float64) / rotamer_coords.shape[0]
        return RotamerScoreResult(weights=weights, partition=0.0, energies=pot_energy + dh_energy)
    weights = boltzmann / partition
    return RotamerScoreResult(weights=weights, partition=partition, energies=pot_energy + dh_energy)


def kappa2_from_vectors(mu_donor: np.ndarray, mu_acceptor: np.ndarray, r_vectors: np.ndarray) -> np.ndarray:
    """Compute orientation factors from dipole and distance vectors.

    Parameters
    ----------
    mu_donor : numpy.ndarray
        Donor transition-dipole vectors with shape ``(n_donor, 3)``.
    mu_acceptor : numpy.ndarray
        Acceptor transition-dipole vectors with shape ``(n_acceptor, 3)``.
    r_vectors : numpy.ndarray
        Donor-to-acceptor vectors with shape ``(n_donor, n_acceptor, 3)``.

    Returns
    -------
    numpy.ndarray
        ``kappa^2`` matrix.
    """
    r_vectors = np.asarray(r_vectors, dtype=np.float64)
    r_norm = np.linalg.norm(r_vectors, axis=2, keepdims=True)
    r_unit = np.divide(r_vectors, r_norm, out=np.zeros_like(r_vectors), where=r_norm > 0.0)
    cos_da = np.einsum("ik,jk->ij", mu_donor, mu_acceptor)
    cos_dr = np.einsum("ik,ijk->ij", mu_donor, r_unit)
    cos_ar = np.einsum("jk,ijk->ij", mu_acceptor, r_unit)
    return np.power(cos_da - 3.0 * cos_dr * cos_ar, 2)
