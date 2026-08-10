"""FRET distance metrics and chi-squared scoring.

Pure Numba/numpy utility functions for converting between FRET observables
and model distances, and for scoring model distances against experimental
restraints with asymmetric error bars.

These functions are intended as **non-breaking enhancements** to
``IMP.bff`` and may later migrate into the upstream IMP.bff C++ layer.
"""

from __future__ import annotations

import math
import numpy as np
from IMP.bff._jit import njit as _njit, jit as _jit


@_njit
def chi2_score(
    model_distance: float,
    exp_distance: float,
    error_neg: float,
    error_pos: float,
) -> float:
    r"""Asymmetric chi-squared contribution for a single distance.

    When the model distance is below the experimental value the negative
    error is used; otherwise the positive error is used:

    .. math::

       \chi_i^2 = \begin{cases}
           \left(\frac{d_{\text{mod}} - d_{\text{exp}}}{\sigma_{-}}\right)^2
           & d_{\text{mod}} < d_{\text{exp}} \\
           \left(\frac{d_{\text{mod}} - d_{\text{exp}}}{\sigma_{+}}\right)^2
           & d_{\text{mod}} \geq d_{\text{exp}}
       \end{cases}

    Parameters
    ----------
    model_distance : float
        Model-predicted distance (Å).
    exp_distance : float
        Experimental distance (Å).
    error_neg : float
        Negative (lower) experimental error (Å).
    error_pos : float
        Positive (upper) experimental error (Å).

    Returns
    -------
    float
        Chi-squared contribution.

    Examples
    --------
    >>> chi2_score(50.0, 55.0, 3.0, 5.0)
    1.777777...
    >>> chi2_score(60.0, 55.0, 3.0, 5.0)
    1.0
    """
    d = model_distance - exp_distance
    err = error_neg if d < 0 else error_pos
    if err == 0.0:
        return 0.0 if d == 0.0 else math.inf
    return (d / err) ** 2


@_njit
def fret_efficiency(distance: float, forster_radius: float = 52.0) -> float:
    r"""FRET efficiency *E* for a single donor–acceptor distance.

    .. math::

       E = \frac{1}{1 + \left(\frac{R}{R_0}\right)^6}

    Parameters
    ----------
    distance : float
        Donor–acceptor distance *R* (Å).
    forster_radius : float
        Förster radius *R_0* (Å).  Default 52.0.

    Returns
    -------
    float
        Efficiency *E* ∈ [0, 1].
    """
    return 1.0 / (1.0 + (distance / forster_radius) ** 6.0)


@_njit
def distance_from_fret_efficiency(efficiency: float, forster_radius: float = 52.0) -> float:
    r"""Convert FRET efficiency to distance.

    .. math::

       R = R_0 \left(\frac{1}{E} - 1\right)^{1/6}

    Parameters
    ----------
    efficiency : float
        FRET efficiency *E*.
    forster_radius : float
        Förster radius *R_0* (Å).  Default 52.0.

    Returns
    -------
    float
        Distance *R* (Å).
    """
    if efficiency <= 0.0 or efficiency >= 1.0:
        return 0.0 if efficiency >= 1.0 else math.inf
    return forster_radius * ((1.0 / efficiency) - 1.0) ** (1.0 / 6.0)


@_njit
def gaussian_rmp_to_rda_mean(rmp: float, sigma: float = 6.0) -> float:
    r"""Correct an Rmp (mean-position) distance to approximate <R_DA>.

    For two 3-D Gaussian dye distributions with equal width *σ*, the mean
    inter-dye distance is larger than the distance between the mean
    positions because of the convolution of the distributions:

    .. math::

       \langle R_{DA} \rangle \approx R_{mp}
       + \frac{\sigma^2}{R_{mp}}

    Parameters
    ----------
    rmp : float
        Distance between AV mean positions (Å).
    sigma : float
        Width of the dye distributions (Å).  Default 6.0.

    Returns
    -------
    float
        Sigma-corrected mean distance (Å).
    """
    if rmp <= 0.0:
        return 0.0
    return rmp + sigma ** 2 / rmp


@_njit
def polynomial_transfer(rmp: float, coefficients: np.ndarray) -> float:
    r"""Polynomial transfer function: Rmp -> RDAMean.

    Evaluates a polynomial *P* such that
    :math:`R_{DA,mean} = P(R_{mp})`.

    Parameters
    ----------
    rmp : float
        Mean-position distance (Å).
    coefficients : (n,) np.ndarray
        Polynomial coefficients ``[c0, c1, c2, …]`` for
        :math:`c_0 + c_1 x + c_2 x^2 + \dots`.

    Returns
    -------
    float
    """
    y = 0.0
    x = 1.0
    for c in coefficients:
        y += c * x
        x *= rmp
    return y


@_njit
def av_pair_statistics(
    distances: np.ndarray,
    weights: np.ndarray,
    forster_radius: float = 52.0,
) -> tuple[float, float, float, float]:
    """Compute all four AV distance metrics from a weighted distance sample.

    Parameters
    ----------
    distances : (n_samples,) np.ndarray
        Sampled distances (Å).
    weights : (n_samples,) np.ndarray
        Weight for each sample (product of dye densities).
    forster_radius : float
        Förster radius (Å).

    Returns
    -------
    r_da_mean : float
        Weighted mean distance <R_DA>.
    r_mp : float
        **Always NaN.** R_mp is the distance between the two clouds' *mean
        positions*, and that cannot be recovered from a sample of pair
        distances: |<a> - <b>| is not a function of the distribution of
        |a - b|. Use :func:`mean_position_distance`, which takes the clouds.

        This slot returned ``r_da_mean`` until 2026-07-28, documented as "same
        as r_da_mean here" -- so a caller reading index 1 silently got a
        different quantity under the right name. Measured on 148l E15/E90:
        47.65 A (true R_mp) against 51.53 A (what was returned), 8 % apart.
        NaN propagates instead of lying.
    r_e : float
        FRET-averaged distance R_E.
    mean_efficiency : float
        Weighted mean FRET efficiency <E>.
    """
    w_sum = weights.sum()
    if w_sum == 0.0:
        return 0.0, math.nan, 0.0, 0.0

    r_da_mean = np.dot(distances, weights) / w_sum

    e_vals = 1.0 / (1.0 + (distances / forster_radius) ** 6.0)
    mean_efficiency = np.dot(e_vals, weights) / w_sum

    if mean_efficiency <= 0.0 or mean_efficiency >= 1.0:
        r_e = 0.0 if mean_efficiency >= 1.0 else math.inf
    else:
        r_e = forster_radius * ((1.0 / mean_efficiency) - 1.0) ** (1.0 / 6.0)

    return r_da_mean, math.nan, r_e, mean_efficiency


def mean_position_distance(
    points_a: np.ndarray,
    points_b: np.ndarray,
    weights_a: np.ndarray | None = None,
    weights_b: np.ndarray | None = None,
) -> float:
    """Distance between the mean positions of two point clouds, R_mp.

    This is a different quantity from the mean of the pair distances, and it
    needs the clouds themselves -- ``|<a> - <b>|`` cannot be recovered from the
    distribution of ``|a - b|``, which is why :func:`av_pair_statistics` returns
    NaN for it.

    Contributed from QuEst, where it is ``AV.dRmp``.

    Parameters
    ----------
    points_a, points_b : (n, 3) np.ndarray
        Accessible-volume point clouds in Angstrom.
    weights_a, weights_b : (n,) np.ndarray, optional
        Per-point weights. Uniform if omitted.

    Returns
    -------
    float
        Distance in Angstrom.
    """
    a = np.asarray(points_a, dtype=np.float64).reshape((-1, 3))
    b = np.asarray(points_b, dtype=np.float64).reshape((-1, 3))
    if a.shape[0] == 0 or b.shape[0] == 0:
        raise ValueError("Both point clouds must be non-empty to define R_mp.")
    mean_a = a.mean(axis=0) if weights_a is None else np.average(a, axis=0, weights=weights_a)
    mean_b = b.mean(axis=0) if weights_b is None else np.average(b, axis=0, weights=weights_b)
    return float(np.sqrt(((mean_a - mean_b) ** 2).sum()))
