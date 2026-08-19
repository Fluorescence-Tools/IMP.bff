"""``RotamerEnsemble``: a screened rotamer library at a site, usable as an AV.

The fps.json position type ``R1`` (PRD-108): a FRETpredict-style 1:1 rotamer
library transformed into the residue's backbone frame and Boltzmann-screened
against the protein. Per rotamer the chromophore centre, the transition
dipole and the weight are kept -- and every atom, so nothing downstream has
to re-derive them -- which is what a FRET *rate distribution* over a pair of
labels needs (R_ij, κ²_ij, w_i·w_j), not just a mean position.

It subclasses :class:`IMP.bff.representation.av.structure.AccessibleVolume` with ``points`` =
(N, 4) centre + weight, so every AV helper in ``fret`` (``av_pair_statistics``,
``histogram_rda``, ``mean_fret_distance``, ...) accepts it unchanged; the pair
kernels live in :mod:`IMP.bff.representation.distance` (``fret_pair_geometry``,
``fret_pair_efficiencies``).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import numpy as np

from IMP.bff.representation.states import States
from IMP.bff.representation.distance import fret_pair_efficiencies, fret_pair_geometry
from IMP.bff.dye.spectra import forster_radius_from_spectra
from IMP.bff.representation.rotamer.io import load_protein_frames, load_rotamer_library
from IMP.bff.scoring.rotamer import compute_rotamer_score

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


def selector_atom_indices(atom_names: Sequence[str], selector) -> list[int]:
    """Indices of the atoms named by a FRETpredict selector (``'C7 and resname A48'``)."""
    names = [selector] if isinstance(selector, str) else list(selector)
    indices: list[int] = []
    for item in names:
        wanted = str(item).split(" and ")[0].strip().upper()
        for i, name in enumerate(atom_names):
            if str(name).upper() == wanted:
                indices.append(i)
                break
    if not indices:
        raise ValueError(f"Atom selector {selector!r} not found in rotamer library")
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

@dataclass(kw_only=True)
class RotamerEnsemble(States):
    """A rotamer library placed and screened at one labelling site (fps ``R1``).

    A **sibling** of :class:`~IMP.bff.representation.AccessibleVolume`, not a
    subclass of it. It used to inherit from the concrete AV, which is why
    distance code worked for rotamers -- by inheritance rather than by design --
    and it then had to carry a grid it does not have, filling
    ``density``, ``grid_step`` and ``grid_shape`` with empty placeholders. No
    consumer ever read them: they all used ``points``, ``mean_position``,
    ``n_points`` or ``has_volume``, which is the
    :class:`~IMP.bff.representation.States` surface both representations share.

    From ``States``: ``points`` (N, 4) chromophore centre + weight,
    ``attachment_point`` (CA), ``orientations`` (the per-rotamer transition
    dipoles), ``position_name``, ``params`` (carries ``simulation_type='R1'``,
    library, temperature, ...). Adds ``atoms`` (N, n_atoms, 3) in the protein
    frame, ``atom_names``, ``resnames``, ``energies``, ``partition`` (Z) and the
    ``library`` name.

    ``mu`` remains as the name the rotamer code uses for the dipoles; it is the
    same array as ``States.orientations``, which is what representation-agnostic
    code asks for.
    """

    mu: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))

    atoms: np.ndarray = field(default_factory=lambda: np.zeros((0, 0, 3)))
    atom_names: tuple = ()
    resnames: tuple = ()
    energies: np.ndarray = field(default_factory=lambda: np.zeros(0))
    partition: float = 0.0
    library: str = ""
    chain: str = ""
    residue: int = 0

    def __post_init__(self):
        # `mu` is what the rotamer code calls the transition dipoles;
        # `States.orientations` is what representation-agnostic code asks for.
        # They are one array, so a caller cannot set one and read a stale other.
        if self.orientations is None and self.mu is not None and len(self.mu):
            object.__setattr__(self, "orientations", self.mu)

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
        centre_idx = selector_atom_indices(names, metadata.get("r", []))[0]
        mu_idx = selector_atom_indices(names, metadata.get("mu", []))
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
    from IMP.bff.representation.rotamer.fps import RotamerPosition

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


__all__ = [
    "SIMULATION_TYPE_R1",
    "RotamerEnsemble",
    "rotamer_ensembles_from_fps",
    "resolve_backbone_site",
    "backbone_rotation",
    "transform_library_to_site",
    "selector_atom_indices",
]
