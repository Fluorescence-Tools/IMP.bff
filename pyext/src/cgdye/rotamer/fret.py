"""Rotamer-based FRET prediction compatible with FRETpredict workflows."""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from IMP.bff.cgdye.rotamer.io import load_protein_frames, load_rotamer_library
from IMP.bff.cgdye.rotamer.r0 import calculate_r0
from IMP.bff.cgdye.rotamer.scoring import compute_rotamer_score, kappa2_from_vectors


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


def _resolve_backbone_site(
    frame: dict[str, Any],
    chain: str | None,
    residue: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Resolve CA, N, and C coordinates from a protein frame.

    Parameters
    ----------
    frame : dict
        Frame dictionary from ``load_protein_frames``.
    chain : str or None
        Chain ID.
    residue : int
        Residue number.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]
        CA, N, and C coordinates.
    """
    coords = np.asarray(frame["coords"], dtype=np.float64)
    atom_names = [str(name).upper() for name in frame["atom_names"]]
    chain_id = (chain or "").upper()
    chain_ids = [str(value).upper() for value in frame.get("chain_ids", [""] * len(atom_names))]
    residue_indices = [int(value) for value in frame.get("residue_indices", [-1] * len(atom_names))]
    matches: dict[str, np.ndarray] = {}
    for coord, atom_name, frame_chain, frame_residue in zip(coords, atom_names, chain_ids, residue_indices):
        if atom_name not in {"CA", "N", "C"} or atom_name in matches:
            continue
        if chain_id and frame_chain and frame_chain != chain_id:
            continue
        if frame_residue != -1 and frame_residue != residue:
            continue
        matches[atom_name] = coord
    missing = [name for name in ["CA", "N", "C"] if name not in matches]
    if missing:
        raise ValueError(f"Missing backbone atoms {missing} for chain {chain or 'A'} residue {residue}")
    return matches["CA"], matches["N"], matches["C"]


def _to_vector(value: Any):
    """Convert a coordinate-like value to an IMP vector.

    Parameters
    ----------
    value : object
        Coordinate-like object.

    Returns
    -------
    object
        IMP vector-compatible value.
    """
    return value


def _transform_coords(coords: np.ndarray, ca: Any, n: Any, c: Any) -> np.ndarray:
    """Transform rotamer coordinates into a protein backbone frame.

    Parameters
    ----------
    coords : numpy.ndarray
        Rotamer coordinates with shape ``(n_rotamers, n_atoms, 3)``.
    ca, n, c
        Backbone CA, N, and C coordinates.

    Returns
    -------
    numpy.ndarray
        Transformed coordinates.
    """
    ca_v = np.asarray([float(ca[0]), float(ca[1]), float(ca[2])], dtype=np.float64)
    n_v = np.asarray([float(n[0]), float(n[1]), float(n[2])], dtype=np.float64)
    c_v = np.asarray([float(c[0]), float(c[1]), float(c[2])], dtype=np.float64)

    x_vector = n_v - ca_v
    x_vector /= np.linalg.norm(x_vector)
    yt_vector = c_v - ca_v
    yt_vector /= np.linalg.norm(yt_vector)
    z_vector = np.cross(x_vector, yt_vector)
    z_vector /= np.linalg.norm(z_vector)
    y_vector = np.cross(z_vector, x_vector)
    rotation = np.vstack([x_vector, y_vector, z_vector])
    return np.tensordot(coords, rotation, axes=([2], [0])) + ca_v


def _selector_atom_indices(atom_names: list[str], selector: str | list[str]) -> list[int]:
    """Find rotamer atom indices from a FRETpredict selector.

    Parameters
    ----------
    atom_names : list of str
        Rotamer atom names.
    selector : str or list of str
        FRETpredict atom selector.

    Returns
    -------
    list of int
        Atom indices.
    """
    if isinstance(selector, str):
        names = [selector]
    else:
        names = list(selector)
    indices: list[int] = []
    for item in names:
        wanted = str(item).split(" and ")[0].strip().upper()
        for i, name in enumerate(atom_names):
            if name.upper() == wanted:
                indices.append(i)
                break
    if not indices:
        raise ValueError(f"Atom selector {selector!r} not found in rotamer library")
    return indices


def _selector_atom_index(atom_names: list[str], selector: str | list[str]) -> int:
    """Find a rotamer atom index from a FRETpredict selector.

    Parameters
    ----------
    atom_names : list of str
        Rotamer atom names.
    selector : str or list of str
        FRETpredict atom selector.

    Returns
    -------
    int
        Atom index.
    """
    return _selector_atom_indices(atom_names, selector)[0]


def _weighted_average_sd_se(values: np.ndarray, weights: np.ndarray) -> tuple[float, float, float]:
    """Compute weighted average, standard deviation, and standard error.

    Parameters
    ----------
    values : numpy.ndarray
        Values to average.
    weights : numpy.ndarray
        Weights.

    Returns
    -------
    tuple[float, float, float]
        Average, standard deviation, and standard error.
    """
    finite = np.isfinite(values)
    values = values[finite]
    weights = weights[finite]
    if values.size == 0:
        return (float("nan"), float("nan"), float("nan"))
    weights = weights / np.sum(weights)
    avg = float(np.average(values, weights=weights))
    variance = float(np.average((values - avg) ** 2, weights=weights))
    return avg, math.sqrt(variance), math.sqrt(variance / values.size)


def _calculate_ws(z_values: np.ndarray) -> np.ndarray:
    """Calculate per-frame weights from partition functions.

    Parameters
    ----------
    z_values : numpy.ndarray
        Array with shape ``(n_frames, 2)``.

    Returns
    -------
    numpy.ndarray
        Per-frame weights.
    """
    z_values = np.asarray(z_values, dtype=np.float64)
    if z_values.shape == (2,):
        return np.array([1.0], dtype=np.float64)
    if z_values.ndim != 2 or z_values.shape[1] != 2:
        raise ValueError(f"Expected Z array with shape (n_frames, 2), got {z_values.shape}")
    z_s = z_values[:, 0] * z_values[:, 1]
    total = np.sum(z_s)
    if total == 0:
        return np.ones(z_values.shape[0], dtype=np.float64) / z_values.shape[0]
    return z_s / total


def _effective_fraction(weights: np.ndarray) -> float:
    """Compute the effective fraction of contributing frames.

    Parameters
    ----------
    weights : numpy.ndarray
        Per-frame weights.

    Returns
    -------
    float
        Effective fraction.
    """
    weights = np.asarray(weights, dtype=np.float64)
    weights = weights[weights != 0]
    if weights.size == 0:
        return 0.0
    uniform = np.ones_like(weights) / weights.size
    entropy = -np.sum(weights * np.log(weights / uniform))
    return float(np.exp(entropy))


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

        if self.verbose:
            logging.basicConfig(filename=kwargs.get("log_file", "log"), level=logging.DEBUG)
            logging.root.setLevel(logging.DEBUG)
        else:
            logging.basicConfig(filename=kwargs.get("log_file", "log"), level=logging.INFO)
            logging.root.setLevel(logging.INFO)

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
        """Resolve a protein labeling site in one frame.

        Parameters
        ----------
        frame : dict
            Frame dictionary from ``load_protein_frames``.
        chain : str or None
            Chain ID.
        residue : int
            Residue number.

        Returns
        -------
        dict
            Mapping with keys ``CA``, ``N``, and ``C``.
        """
        ca, n, c = _resolve_backbone_site(frame, chain, residue)
        return {"CA": ca, "N": n, "C": c}

    def _transform_library(self, library: dict[str, Any], frame: dict[str, Any], chain: str | None, residue: int) -> np.ndarray:
        """Transform a rotamer library onto a protein site.

        Parameters
        ----------
        library : dict
            Rotamer library.
        frame : dict
            Protein frame.
        chain : str or None
            Chain ID.
        residue : int
            Residue number.

        Returns
        -------
        numpy.ndarray
            Transformed rotamer coordinates.
        """
        site = self._resolve_site(frame, chain, residue)
        ca, n, c = [site[key] for key in ["CA", "N", "C"]]
        return _transform_coords(library["coords"], ca, n, c)

    def _frame_fret(
        self,
        frame: dict[str, Any],
        rotamers_1: np.ndarray,
        rotamers_2: np.ndarray,
    ) -> FRETFrameResult:
        """Calculate FRET quantities for one protein frame.

        Parameters
        ----------
        frame : dict
            Protein frame.
        rotamers_1 : numpy.ndarray
            Donor rotamer coordinates.
        rotamers_2 : numpy.ndarray
            Acceptor rotamer coordinates.

        Returns
        -------
        FRETFrameResult
            Per-frame FRET result.
        """
        score_1 = compute_rotamer_score(
            rotamers_1,
            frame["coords"],
            frame["atom_names"],
            frame["resnames"],
            self.lib_1["atom_names"],
            self.lib_1.get("metadata", {}),
            self.lib_1.get("resnames"),
            protein_residue_indices=frame.get("residue_indices"),
            protein_chain_ids=frame.get("chain_ids"),
            site_residue=self.residues[0],
            site_chain=self.chains[0],
            rotamer_weights=self.lib_1.get("weights"),
            temperature=self.temperature,
            ignore_h=self.ignore_h,
            electrostatic=self.electrostatic,
            potential=self.potential,
            sigma_scaling=self.sigma_scaling,
            epsilon_scaling=self.epsilon_scaling,
        )
        score_2 = compute_rotamer_score(
            rotamers_2,
            frame["coords"],
            frame["atom_names"],
            frame["resnames"],
            self.lib_2["atom_names"],
            self.lib_2.get("metadata", {}),
            self.lib_2.get("resnames"),
            protein_residue_indices=frame.get("residue_indices"),
            protein_chain_ids=frame.get("chain_ids"),
            site_residue=self.residues[1],
            site_chain=self.chains[1],
            rotamer_weights=self.lib_2.get("weights"),
            temperature=self.temperature,
            ignore_h=self.ignore_h,
            electrostatic=self.electrostatic,
            potential=self.potential,
            sigma_scaling=self.sigma_scaling,
            epsilon_scaling=self.epsilon_scaling,
        )
        weights_1 = score_1.weights
        weights_2 = score_2.weights
        combined_matrix = np.outer(weights_1, weights_2)

        metadata_1 = dict(self.lib_1.get("metadata", {}) or {})
        metadata_2 = dict(self.lib_2.get("metadata", {}) or {})
        metadata_1.setdefault("weights", self.lib_1.get("weights", self.lib_1.get("weight")))
        metadata_2.setdefault("weights", self.lib_2.get("weights", self.lib_2.get("weight")))
        mu_1_sel = metadata_1.get("mu", [])
        mu_2_sel = metadata_2.get("mu", [])
        center_1_sel = metadata_1.get("r", [])
        center_2_sel = metadata_2.get("r", [])
        mu_1_indices = _selector_atom_indices(self.lib_1["atom_names"], mu_1_sel)
        mu_2_indices = _selector_atom_indices(self.lib_2["atom_names"], mu_2_sel)
        center_1_idx = _selector_atom_index(self.lib_1["atom_names"], center_1_sel)
        center_2_idx = _selector_atom_index(self.lib_2["atom_names"], center_2_sel)

        centers_1 = rotamers_1[:, center_1_idx, :]
        centers_2 = rotamers_2[:, center_2_idx, :]
        r_vectors = centers_1[:, None, :] - centers_2[None, :, :]
        distances_angstrom = np.linalg.norm(r_vectors, axis=2)
        distances_nm = distances_angstrom / 10.0

        mu_1 = rotamers_1[:, mu_1_indices[1], :] - rotamers_1[:, mu_1_indices[0], :] if len(mu_1_indices) >= 2 else rotamers_1[:, 1, :] - rotamers_1[:, 0, :]
        mu_2 = rotamers_2[:, mu_2_indices[0], :] - rotamers_2[:, mu_2_indices[1], :] if len(mu_2_indices) >= 2 else rotamers_2[:, 0, :] - rotamers_2[:, 1, :]
        mu_1 = mu_1 / np.linalg.norm(mu_1, axis=1, keepdims=True)
        mu_2 = mu_2 / np.linalg.norm(mu_2, axis=1, keepdims=True)
        k2 = kappa2_from_vectors(mu_1, mu_2, r_vectors)
        k2_avg = float(np.sum(k2 * combined_matrix))

        if not self.fixed_R0:
            self.r0 = calculate_r0(self.donor, self.acceptor, k2_avg, r0_dir=self.r0lib)
            if self.r0 == 0:
                return FRETFrameResult((float(score_1.partition), float(score_2.partition)), float("nan"), float("nan"), float("nan"), float("nan"))

        ratio6 = np.power(distances_nm / self.r0, 6)
        estatic = 1.0 / (1.0 + 2.0 / 3.0 * np.divide(ratio6, k2))
        edynamic1 = 1.0 / (1.0 + 2.0 / 3.0 / k2_avg * ratio6)
        a_avg = np.sum(3.0 / 2.0 * k2 / ratio6 * combined_matrix)
        edynamic2 = a_avg / (a_avg + 1.0)

        return FRETFrameResult(
            (float(score_1.partition), float(score_2.partition)),
            k2_avg,
            float(np.sum(estatic * combined_matrix)),
            float(np.sum(edynamic1 * combined_matrix)),
            float(edynamic2),
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
            rotamers_1 = self._transform_library(self.lib_1, frame, self.chains[0], self.residues[0])
            rotamers_2 = self._transform_library(self.lib_2, frame, self.chains[1], self.residues[1])
            result = self._frame_fret(frame, rotamers_1, rotamers_2)
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
        logging.debug("Done")
