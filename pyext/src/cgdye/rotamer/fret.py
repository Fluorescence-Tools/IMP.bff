"""Rotamer-based FRET prediction compatible with FRETpredict workflows."""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from IMP.bff.cgdye.rotamer.io import load_protein_frames, load_rotamer_library
from IMP.bff.dye.spectra import forster_radius_from_spectra
from IMP.bff.cgdye.rotamer.ensemble import RotamerEnsemble, resolve_backbone_site, transform_library_to_site
from IMP.bff.representation.distance import fret_pair_efficiencies

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
            self.r0 = forster_radius_from_spectra(self.donor, self.acceptor, k2_avg, r0_dir=self.r0lib)
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
