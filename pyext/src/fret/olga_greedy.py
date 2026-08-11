"""Olga-style greedy informative FRET pair selection (Python port).

This module is a headless implementation of the "informative pair selection" / "experiment planning"
algorithm used by Olga.

Original Olga sources (ChiSurf repo copy):

- `playground/Olga/src/best_dist.h`
  - `greedySelection(...)`
  - `bestPair(...)`
  - `precisionDecay(...)`
  - `rmsdMeanMean(...)`, `rmsdMeanMeanAdd(...)`
  - `chiSquared(...)`

- `playground/Olga/src/chisqdist.hpp`
  - `chisqRTcdf(...)` and gamma-function helpers

- Call flow in GUI:
  `playground/Olga/src/gui/GetInformativePairsDialog.cpp`
  - prepares the efficiency matrix and RMSD matrix
  - performs NaN thresholding/filling
  - calls `greedySelection(...)` then `precisionDecay(...)`

Notes / differences vs Olga:

- This module assumes `effs` contains finite values already (Olga has explicit NaN handling).
- The weighting uses the chi-squared *right-tail* CDF (`chisqRTcdf`). The hand-ported
  expansion was wrong for every **odd** number of degrees of freedom; it is fixed and
  vectorised here, and checked against `scipy.special.gammaincc` -- see `_chisq_rt_cdf`.
- Inputs are accepted as `float32` (the C++/Eigen path's types), but the chi-squared
  accumulator and the weights are `float64`: the weights are a difference of tail
  probabilities and are what the greedy choice turns on.
"""

from __future__ import annotations

import math
from typing import Tuple

import numpy as np
from scipy.special import erfc, gammaincc


def _chisq_rt_cdf(chisq: np.ndarray, ndof: int) -> np.ndarray:
    r"""Chi-squared right-tail CDF, used as the weight.

    Evaluates :math:`Q(\nu/2,\ \chi^2/2)`, the regularized upper incomplete
    gamma, via the closed forms that exist because :math:`\nu/2` is always an
    integer or a half-integer -- elementwise over an array of any shape.

    Parameters
    ----------
    chisq : numpy.ndarray
        Chi-squared values, any shape.
    ndof : int
        Degrees of freedom.

    Returns
    -------
    numpy.ndarray
        The right-tail probability, elementwise, ``float64``.

    Notes
    -----
    **The half-integer branch was wrong for every odd** ``ndof``. Olga takes the
    expansion from Boost, whose half-integer branch loops
    ``for (n = 2; n < a; ++n)`` with ``a`` a half-integer; the port wrote
    ``range(2, int(a))``, and ``int(2.5)`` is ``2`` -- so it ran one iteration
    short every time, zero where Boost runs one. At ``ndof = 5``,
    ``chisq = 3.008`` that returned ``0.3903934`` for a true ``0.6987524`` --
    not a rounding error but very nearly half. ``ndof`` is
    the number of pairs chosen so far, so it is odd on every other greedy step,
    and these weights are exactly what decides which pair looks most
    informative.

    :func:`scipy.special.gammaincc` is the same function and settles the
    direction of that error (the corrected series agrees with it to ``4e-14``),
    but it does not carry the common path: it is general-purpose, and on the
    ``(candidates, n, n)`` arrays this is called with it measured **18x slower**
    end to end than these few-term series. Correctness came from comparing
    against it; speed came from not calling it where the closed form applies.

    ``chisq = 0`` -- the whole diagonal, on every call -- would divide by
    ``sqrt(pi * x)``. ``Q(a, 0) = 1`` exactly, so those entries are filled
    directly rather than computed.
    """
    # 0.5 * chisq allocates, deliberately. Scaling in place would halve the
    # caller's chi-squared accumulator -- which is the same defect this project
    # filed against a library's CDF sampler, and it was written here by hand
    # while trying to save exactly this one allocation.
    x = 0.5 * np.asarray(chisq, dtype=np.float64)
    a = 0.5 * ndof

    if a > 100.0:
        # Olga approximates this tail with a normal, which is off by up to
        # 1.3e-2 near the median. The series would need ~a terms here, but this
        # branch is reached only past 200 selected pairs, so the general
        # function is affordable exactly where the cheap one stops being cheap.
        return gammaincc(a, x)

    if ndof % 2 == 0:
        # a is an integer: Q(a, x) = exp(-x) * sum_{n=0}^{a-1} x**n / n!
        term = np.exp(-x)
        total = term.copy()
        for n in range(1, int(a)):
            term = term * x / n
            total += term
        return total

    # a is a half-integer: Q(a, x) = erfc(sqrt(x)) + exp(-x)/sqrt(pi x) * series
    total = erfc(np.sqrt(x))
    if a > 1.0:
        positive = x > 0.0
        xp = np.where(positive, x, 1.0)
        term = np.exp(-xp) / np.sqrt(np.pi * xp) * xp / 0.5
        series = term.copy()
        n = 2
        while n < a:
            term = term / (n - 0.5) * xp
            series += term
            n += 1
        total = total + np.where(positive, series, 0.0)
    return total


def _weighted_column_means(
        weights: np.ndarray,
        rmsds: np.ndarray,
        diag_weight: float
) -> float:
    """Mean over columns of the weight-averaged RMSD.

    Parameters
    ----------
    weights : numpy.ndarray
        Chi-squared tail weights, shape ``(..., n, n)``; reduced over the row axis.
    rmsds : numpy.ndarray
        Pairwise RMSDs, shape ``(n, n)``.
    diag_weight : float
        Olga's diagonal correction (default 0.99), applied to the denominator.

    Returns
    -------
    float or numpy.ndarray
        The mean over columns; an array if ``weights`` carries leading axes.

    Notes
    -----
    The weighted sum goes through :func:`numpy.einsum` rather than
    ``(weights * rmsds).sum(...)`` because the latter materialises a second
    array the size of ``weights`` -- which at ``(64, n, n)`` is the largest
    allocation in the loop, and it exists only to be summed away.
    """
    sum_w = weights.sum(axis=-2)
    if weights.ndim == 2:
        prod = np.einsum('ij,ij->j', weights, rmsds)
    else:
        prod = np.einsum('kij,ij->kj', weights, rmsds)
    return (prod / (sum_w - 1.0 + diag_weight)).mean(axis=-1)


def _rmsd_mean_mean(
        rmsds: np.ndarray,
        chi2: np.ndarray,
        ndof: int,
        diag_weight: float
) -> float:
    """Mean of the column-wise expected RMSD.

    Port of Olga's ``rmsdMeanMean(rmsds, chi2, ndof, diagWeight)`` in
    ``best_dist.h``.

    Parameters
    ----------
    rmsds : numpy.ndarray
        Pairwise RMSDs, shape ``(n, n)``.
    chi2 : numpy.ndarray
        Accumulated chi-squared, shape ``(n, n)``.
    ndof : int
        Degrees of freedom.
    diag_weight : float
        Diagonal correction.

    Returns
    -------
    float
        Expected mean RMSD.
    """
    return float(_weighted_column_means(_chisq_rt_cdf(chi2, ndof), rmsds, diag_weight))


def _rmsd_mean_mean_add(
        rmsds: np.ndarray,
        chi2: np.ndarray,
        e_add: np.ndarray,
        inv_err_sq: float,
        ndof: int,
        diag_weight: float
) -> np.ndarray:
    """Expected mean RMSD after adding one candidate pair, for many candidates.

    Port of Olga's ``rmsdMeanMeanAdd(...)`` in ``best_dist.h``, evaluated for a
    batch of candidates at once.

    Parameters
    ----------
    rmsds : numpy.ndarray
        Pairwise RMSDs, shape ``(n, n)``.
    chi2 : numpy.ndarray
        Accumulated chi-squared, shape ``(n, n)``.
    e_add : numpy.ndarray
        Per-frame efficiencies of the candidates, shape ``(n,)`` or ``(m, n)``.
    inv_err_sq : float
        ``1 / err**2``; Olga's update is ``d(chi2) = (dE)**2 / err**2``.
    ndof : int
        Degrees of freedom.
    diag_weight : float
        Diagonal correction.

    Returns
    -------
    numpy.ndarray
        Expected mean RMSD, one entry per candidate.
    """
    e = np.atleast_2d(np.asarray(e_add, dtype=np.float64))
    # d[k, row, col] = e[k, row] - e[k, col]; built once and then updated in
    # place into chi2_new, because each of these is the largest array in the
    # loop and a fresh temporary per step is what made the vectorised form
    # memory-bound.
    chi2_new = e[:, :, None] - e[:, None, :]
    np.square(chi2_new, out=chi2_new)
    chi2_new *= inv_err_sq
    chi2_new += chi2
    return _weighted_column_means(_chisq_rt_cdf(chi2_new, ndof), rmsds, diag_weight)


def _add_pair_to_chi2(chi2: np.ndarray, e_pair: np.ndarray, inv_err_sq: float) -> None:
    """Accumulate one selected pair's contribution into chi-squared, in place.

    Matches the ``chiSquared(...)`` accumulation used by Olga's
    ``greedySelection(...)``.

    Parameters
    ----------
    chi2 : numpy.ndarray
        Accumulated chi-squared, shape ``(n, n)``; modified in place.
    e_pair : numpy.ndarray
        Per-frame efficiencies of the selected pair, shape ``(n,)``.
    inv_err_sq : float
        ``1 / err**2``.

    Returns
    -------
    None
    """
    e = np.asarray(e_pair, dtype=np.float64)
    d = e[:, None] - e[None, :]
    chi2 += (d * d) * inv_err_sq


def _best_pair_all_candidates(
        effs: np.ndarray,
        rmsds: np.ndarray,
        chi2: np.ndarray,
        inv_err_sq: float,
        ndof: int,
        diag_weight: float,
        selected_mask: np.ndarray,
        unique_only: bool,
        chunk: int = 64
) -> np.ndarray:
    """Score every candidate pair by the expected mean RMSD it would leave.

    Corresponds to Olga's ``bestPair(...)`` inside ``greedySelection(...)``.

    Parameters
    ----------
    effs : numpy.ndarray
        Efficiencies, shape ``(n_frames, n_pairs)``.
    rmsds : numpy.ndarray
        Pairwise RMSDs, shape ``(n_frames, n_frames)``.
    chi2 : numpy.ndarray
        Accumulated chi-squared, shape ``(n_frames, n_frames)``.
    inv_err_sq : float
        ``1 / err**2``.
    ndof : int
        Degrees of freedom.
    diag_weight : float
        Diagonal correction.
    selected_mask : numpy.ndarray
        Non-zero where a candidate has already been selected.
    unique_only : bool
        If true, an already-selected candidate scores ``+inf``.
    chunk : int
        Candidates evaluated per batch. The intermediate is
        ``(chunk, n_frames, n_frames)``, so this bounds peak memory rather than
        changing any result.

    Returns
    -------
    numpy.ndarray
        Score per candidate, ``float32``.
    """
    m = effs.shape[1]
    out = np.empty(m, dtype=np.float32)
    for start in range(0, m, chunk):
        stop = min(start + chunk, m)
        block = _rmsd_mean_mean_add(
            rmsds, chi2, effs[:, start:stop].T, inv_err_sq, ndof, diag_weight
        )
        out[start:stop] = block.astype(np.float32)
    if unique_only:
        out[np.asarray(selected_mask) != 0] = np.float32(np.inf)
    return out


def _precision_decay(
        selected_pairs: np.ndarray,
        effs: np.ndarray,
        rmsds: np.ndarray,
        inv_err_sq: float,
        diag_weight: float
) -> np.ndarray:
    """Expected mean RMSD after each greedy step.

    Port of Olga's ``precisionDecay(...)`` in ``best_dist.h``.

    Parameters
    ----------
    selected_pairs : numpy.ndarray
        Indices of the selected pairs, in selection order.
    effs : numpy.ndarray
        Efficiencies, shape ``(n_frames, n_pairs)``.
    rmsds : numpy.ndarray
        Pairwise RMSDs, shape ``(n_frames, n_frames)``.
    inv_err_sq : float
        ``1 / err**2``.
    diag_weight : float
        Diagonal correction.

    Returns
    -------
    numpy.ndarray
        Expected mean RMSD after each step, ``float32``.
    """
    n = rmsds.shape[0]
    chi2 = np.zeros((n, n), dtype=np.float64)
    decay = np.empty(selected_pairs.shape[0], dtype=np.float32)
    for i in range(selected_pairs.shape[0]):
        _add_pair_to_chi2(chi2, effs[:, selected_pairs[i]], inv_err_sq)
        decay[i] = np.float32(_rmsd_mean_mean(rmsds, chi2, i + 1, diag_weight))
    return decay


def select_informative_pairs(
    effs: np.ndarray,
    rmsds: np.ndarray,
    err: float,
    max_pairs: int,
    unique_only: bool = True,
    diag_weight: float = 0.99,
) -> Tuple[np.ndarray, np.ndarray]:
    """Olga-style greedy informative pair selection.

    Mapping to Olga:

    - This function implements the control flow of `greedySelection(...)` and
      `precisionDecay(...)` from `playground/Olga/src/best_dist.h`.
    - The chi-squared right-tail weights come from `chisqRTcdf(...)` in
      `playground/Olga/src/chisqdist.hpp`.

    Notes:

    - Olga's GUI performs NaN filtering/filling before calling the selector. This
      function expects finite `effs` values.
    - `ndof` is increased with each added pair; for the very first step we clamp
      to 1 to avoid invalid degrees of freedom.

    Parameters
    ----------
    effs:
        Array of shape (n_frames, n_pairs) with FRET efficiencies per frame.
    rmsds:
        Array of shape (n_frames, n_frames) with pairwise RMSDs (Angstrom).
    err:
        Expected absolute error in FRET efficiency.
    max_pairs:
        Number of pairs to select (will be capped to n_pairs).
    unique_only:
        If True, each candidate pair can be selected at most once.
    diag_weight:
        Weight used in Olga for the diagonal correction (default 0.99).

    Returns
    -------
    selected_pair_indices:
        Indices into the *pair* dimension of `effs`.
    precision_decay:
        Vector of length n_selected with the expected mean RMSD after adding each pair.
    """
    effs_f = np.ascontiguousarray(effs, dtype=np.float32)
    rmsds_f = np.ascontiguousarray(rmsds, dtype=np.float32)

    if effs_f.ndim != 2:
        raise ValueError("effs must be 2D (n_frames, n_pairs)")
    if rmsds_f.ndim != 2 or rmsds_f.shape[0] != rmsds_f.shape[1]:
        raise ValueError("rmsds must be 2D square (n_frames, n_frames)")
    if rmsds_f.shape[0] != effs_f.shape[0]:
        raise ValueError("rmsds size must match effs number of frames")

    n_pairs = effs_f.shape[1]
    max_pairs = int(min(max_pairs, n_pairs))
    if max_pairs <= 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.float32)

    inv_err_sq = float(1.0 / (err * err))

    n = rmsds_f.shape[0]
    chi2 = np.zeros((n, n), dtype=np.float64)
    selected_mask = np.zeros(n_pairs, dtype=np.uint8)
    selected = np.empty(max_pairs, dtype=np.int64)

    for step in range(max_pairs):
        ndof = max(step - 1, 1)
        scores = _best_pair_all_candidates(
            effs_f, rmsds_f, chi2,
            inv_err_sq, ndof, float(diag_weight),
            selected_mask, bool(unique_only),
        )
        best = int(np.argmin(scores))
        selected[step] = best
        selected_mask[best] = 1
        _add_pair_to_chi2(chi2, effs_f[:, best], inv_err_sq)

    decay = _precision_decay(selected, effs_f, rmsds_f, inv_err_sq, float(diag_weight))
    return selected, decay
