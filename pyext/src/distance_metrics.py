"""Distance helpers, and two that disagree with their namesakes.

**Read this before using ``polynomial_transfer`` or
``gaussian_rmp_to_rda_mean`` from here.** Both exist under the same names in
:mod:`IMP.bff.representation.distance`, with the same signatures, and give
*different answers*:

===========================  ====================  ==========================
function                     here                  representation.distance
===========================  ====================  ==========================
``gaussian_rmp_to_rda_mean`` ``Rmp + s^2/Rmp``     ``Rmp + s^2/(2 Rmp)``
``polynomial_transfer``      ascending ``c0+c1x``  descending (``np.polyfit``)
===========================  ====================  ==========================

Measured at ``rmp=45, sigma=6``: **45.8 against 45.4** -- a factor of two in the
correction term. And at ``rmp=45, coeffs=[0, 1, 0.02]``: **85.5 against 45.02**,
the polynomial evaluated with its coefficients reversed.

Each is internally consistent and each documents its own convention, so neither
is wrong in isolation. What is wrong is that they share a name in one package.
The versions in ``representation.distance`` are the live ones -- ``fret/engine.py``
calls them, and its ``fit_transfer_polynomial`` produces ``np.polyfit``
coefficients, which only the descending evaluator reads correctly. The versions
here have **no consumers**; only ``chi2_score`` is imported from this module
(by ``restraints/``), and that one is verified identical.

Which ``gaussian_rmp_to_rda_mean`` is right is a question about what ``sigma``
means -- per-component width of an isotropic 3-D cloud, or width of the distance
distribution -- and is left for the owner rather than guessed at. PRD-113 stage
4b removed only what was verified identical: ``fret_efficiency`` and
``distance_from_fret_efficiency``, which are now re-exported from the canonical
module.

What stage 4c then did, with the numba port: the *arithmetic* of the two
disagreeing pairs is now single-sourced in C++, and the difference between them
is stated as an argument rather than duplicated as a body.
:func:`polynomial_transfer` here reverses its coefficients and calls the same
Horner evaluator; :func:`av_pair_statistics` calls the same reduction. The one
thing not resolved is the one that needs a physicist: the factor of two in
:func:`gaussian_rmp_to_rda_mean`.

Also note that ``av_pair_statistics`` is itself two functions under one name --
this one reduces a *sample of distances*, the canonical one takes two *accessible
volumes* -- with different signatures and different return tuples. Only the
canonical one has callers.
"""

from __future__ import annotations

import math
import numpy as np

import IMP.bff
from IMP.bff.representation.distance import (  # noqa: F401
    chi2_score,
    distance_from_fret_efficiency,
    fret_efficiency,
)


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
    # Reversed, then the same Horner evaluator the descending convention uses.
    # The two orders are the whole difference between this and the canonical
    # ``polynomial_transfer``, and saying so in one line beats a second loop.
    return IMP.bff.polynomial_transfer(
        float(rmp), np.asarray(coefficients, dtype=np.float64).ravel()[::-1].copy())


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
    r_da_mean, r_e, mean_efficiency, _sigma = IMP.bff.distance_sample_statistics(
        np.ascontiguousarray(np.asarray(distances, dtype=np.float64).ravel()),
        np.ascontiguousarray(np.asarray(weights, dtype=np.float64).ravel()),
        float(forster_radius))
    if np.asarray(weights).sum() == 0.0:
        return 0.0, math.nan, 0.0, 0.0
    return float(r_da_mean), math.nan, float(r_e), float(mean_efficiency)


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
