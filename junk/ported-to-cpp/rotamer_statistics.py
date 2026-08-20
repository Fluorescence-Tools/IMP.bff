"""The Python that `IMP/bff/RotamerStatistics.h` replaced.

Not importable, not built, not tested -- the record of what the C++ has to
reproduce, kept because these three are the FRET averaging and the exact
handling of the edge cases is the part worth being able to re-read:

* `_calculate_ws` returns *uniform* weights when every partition-function
  product is zero, not undefined ones;
* `_weighted_average_sd_se` drops non-finite values **and their weights**, then
  renormalises what is left, and divides the variance by the surviving frame
  count rather than the effective count;
* `_effective_fraction` drops zero weights before taking the entropy.

Gated on 400 randomised cases including all-zero Z and NaN frames: identical.
"""

import math

import numpy as np


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

