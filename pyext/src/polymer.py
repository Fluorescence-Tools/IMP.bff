"""Polymer Radial Distribution Functions (RDF) for dye-linker modelling.

These functions compute the end-to-end distance distribution of ideal and
worm-like polymer chains, which is useful for modelling flexible dye linkers
in FRET experiments.

All functions are Numba-accelerated.
"""

from __future__ import annotations

import math
import numpy as np
from IMP.bff._jit import njit as _njit, jit as _jit

from IMP.bff.distributions import normal_distribution


@_njit
def gaussian_chain_ree(segment_length: float, number_of_segments: int) -> float:
    """Root-mean-square end-to-end distance of a Gaussian chain.

    Parameters
    ----------
    segment_length : float
        Length of one Kuhn segment (Å).
    number_of_segments : int
        Number of Kuhn segments.

    Returns
    -------
    float
        RMS end-to-end distance (Å).
    """
    return segment_length * math.sqrt(float(number_of_segments))


@_njit
def gaussian_chain(
    distances: np.ndarray,
    segment_length: float,
    number_of_segments: int,
) -> np.ndarray:
    r"""Radial distribution function of a Gaussian chain in 3-D.

    .. math::

       P(r) = 4\pi r^2 \left(\frac{3}{2\pi \langle r^2 \rangle}\right)^{3/2}
              \exp\!\left(-\frac{3r^2}{2\langle r^2 \rangle}\right)

    Parameters
    ----------
    distances : np.ndarray
        Values of *r* at which to evaluate the PDF.
    segment_length : float
        Kuhn segment length (Å).
    number_of_segments : int
        Number of Kuhn segments.

    Returns
    -------
    np.ndarray
        PDF values (not normalised to 1 unless the grid covers all
        significant probability mass).
    """
    r2_mean = gaussian_chain_ree(segment_length, number_of_segments) ** 2
    return (
        4.0 * np.pi * distances ** 2
        / (2.0 / 3.0 * np.pi * r2_mean) ** 1.5
        * np.exp(-1.5 * distances ** 2 / r2_mean)
    )


# TODO: needs docstring
@_njit
def _qd(r, kappa):
    return (
        (3.0 / (4.0 * math.pi * kappa)) ** 1.5
        * math.exp(-0.75 * r * r / kappa)
        * (1.0 - 1.25 * kappa + 2.0 * r * r - 0.4125 * r * r * r * r / kappa)
    )


@_njit
def worm_like_chain(
    distances: np.ndarray,
    kappa: float,
    chain_length: float = 0.0,
    normalize: bool = True,
    distance: bool = True,
) -> np.ndarray:
    r"""Radial distribution function of a worm-like chain.

    Uses the multi-piece analytical solution from Becker, Rosa & Everaers
    (Eur Phys J E, 32:53–69, 2010).  The parameter *κ* is the
    dimensionless persistence-length ratio.

    Parameters
    ----------
    distances : np.ndarray
        Values of *r* at which to evaluate the PDF.
    kappa : float
        Stiffness parameter (persistence length / contour length).
    chain_length : float
        Contour length *L*.  If 0, inferred from ``max(distances)``.
    normalize : bool
        If True the returned array is normalised so that ``sum(pr) == 1``.
    distance : bool
        If True the PDF is multiplied by :math:`4\pi r^2` (radial
        distribution).

    Returns
    -------
    np.ndarray
        PDF values.
    """
    if chain_length == 0.0:
        chain_length = float(np.max(distances))

    a = 14.054
    b = 0.473
    c = 1.0 - (1.0 + (0.38 * kappa ** (-0.95)) ** (-5.0)) ** (-0.2)
    pr = np.zeros_like(distances, dtype=np.float64)

    if kappa < 0.125:
        d = kappa + 1.0
    else:
        d = 1.0 - 1.0 / (0.177 / (kappa - 0.111) + 6.4 * math.exp(0.783 * math.log(kappa - 0.111)))

    for i in range(len(distances)):
        r = distances[i]
        if r >= chain_length:
            continue
        r_n = r / chain_length

        pri = ((1.0 - c * r_n ** 2.0) / (1.0 - r_n ** 2.0)) ** 2.5
        pri *= math.exp(
            -d * kappa * a * b * (1.0 + b)
            / (1.0 - (b * r_n) ** 2.0)
            * r_n ** 2.0
        )

        g = (
            ((-0.75) / kappa - 0.5) * r_n ** 2.0
            + ((-0.359375) / kappa + 1.0625) * r_n ** 4.0
            + ((-0.109375) / kappa - 0.5625) * r_n ** 6.0
        )
        pri *= math.exp(g / (1.0 - r_n ** 2.0))
        pri *= math.exp(-d * kappa * a * (1.0 + b) * r_n / (1.0 - (b * r_n) ** 2.0))

        pr[i] = pri

    if distance:
        pr = pr * distances ** 2

    if normalize and pr.sum() > 0.0:
        pr /= pr.sum()

    return pr


@_njit
def worm_like_chain_linker(
    distances: np.ndarray,
    kappa: float,
    chain_length: float = 0.0,
    sigma: float = 6.0,
    normalize: bool = True,
) -> np.ndarray:
    r"""Worm-like chain PDF convolved with dye-linker broadening.

    The WLC end-to-end distribution is broadened by a Gaussian kernel
    that represents the finite size / flexibility of the dye linkers
    attached to the chain ends.

    Parameters
    ----------
    distances : np.ndarray
        *r* axis.
    kappa : float
        WLC stiffness parameter.
    chain_length : float
        Contour length (Å).
    sigma : float
        Linker broadening width (Å).
    normalize : bool
        Normalise the final PDF to sum 1.

    Returns
    -------
    np.ndarray
        Broadened PDF.
    """
    pr = worm_like_chain(
        distances=distances,
        kappa=kappa,
        chain_length=chain_length,
        normalize=False,
        distance=False,
    )
    pn = np.zeros_like(pr)
    for i in range(len(distances)):
        r = distances[i]
        if pr[i] == 0.0:
            continue
        # Gaussian broadening centred on r with width sigma
        broad = normal_distribution(distances, loc=r, scale=sigma, norm=False)
        pn += pr[i] * broad

    if normalize and pn.sum() > 0.0:
        pn /= pn.sum()
    return pn
