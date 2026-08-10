"""Mathematical distributions for FRET/dye-modelling.

Numba-accelerated implementations of normal, skew-normal, Poisson and
distance-between-Gaussians PDFs.  Designed to be upstreamed into
``IMP.bff``.
"""

from __future__ import annotations

import math
import numpy as np
from IMP.bff._jit import njit as _njit, jit as _jit


@_jit(nopython=True, nogil=True)
def poisson_0toN(lam: float, N: int) -> np.ndarray:
    r"""Poisson PMF for :math:`k = 0 \dots N-1` via the recursive relation.

    .. math::

       p_0 = e^{-\lambda}, \qquad
       p_k = p_{k-1} \frac{\lambda}{k}

    Parameters
    ----------
    lam : float
        Rate parameter :math:`\lambda`.
    N : int
        Number of terms to compute.

    Returns
    -------
    np.ndarray
        Shape ``(N,)`` array of probabilities.

    Examples
    --------
    >>> import numpy as np
    >>> p = poisson_0toN(0.2, 5)
    >>> np.round(p, 6)
    array([0.818731, 0.163746, 0.016375, 0.001091, 0.000055])
    """
    p = np.empty(N, dtype=np.float64)
    p[0] = math.exp(-lam)
    for i in range(1, N):
        p[i] = p[i - 1] * lam / i
    return p


@_jit(nopython=True, nogil=True)
def normal_distribution(
    x: np.ndarray,
    loc: float = 0.0,
    scale: float = 1.0,
    norm: bool = True,
) -> np.ndarray:
    r"""Normal (Gaussian) PDF.

    .. math::

       f(x) = \frac{1}{\sqrt{2\pi}\sigma}
              \exp\!\left(-\frac{(x-\mu)^2}{2\sigma^2}\right)

    Parameters
    ----------
    x : np.ndarray
        Evaluation points.
    loc : float
        Mean :math:`\mu`.
    scale : float
        Standard deviation :math:`\sigma`.
    norm : bool
        If True, divide the result by its sum so that ``y.sum() == 1``.
        Useful for discretised distributions.

    Returns
    -------
    np.ndarray
    """
    y = 1.0 / (np.sqrt(2.0 * math.pi) * scale) * np.exp(-((x - loc) ** 2) / (2.0 * scale ** 2))
    if norm and y.sum() > 0.0:
        y /= y.sum()
    return y


@_jit(nopython=True, nogil=True)
def generalized_normal_distribution(
    x: np.ndarray,
    loc: float = 0.0,
    scale: float = 1.0,
    shape: float = 0.0,
    norm: bool = True,
) -> np.ndarray:
    r"""Generalised normal distribution with a shape (skew) parameter.

    When *shape* = 0 the distribution reduces to the standard normal.
    For positive *shape* the distribution is left-skewed; for negative
    *shape* it is right-skewed.

    Parameters
    ----------
    x : np.ndarray
    loc : float
    scale : float
    shape : float
        Skewness parameter.
    norm : bool
        Normalise the result so that ``y.sum() == 1``.

    Returns
    -------
    np.ndarray
    """
    if shape == 0.0:
        z = (x - loc) / scale
    else:
        if scale == 0.0:
            return np.zeros_like(x)
        t = 1.0 - shape * (x - loc) / scale
        t[t < 0.0] = np.spacing(1)
        z = -1.0 / shape * np.log(t)

    y = normal_distribution(z, norm=False)
    if norm and y.sum() > 0.0:
        y /= y.sum()
    return y


@_jit(nopython=True)
def distance_between_gaussian(
    distances: np.ndarray,
    separation_distance: float,
    sigma: float,
    normalize: bool = False,
) -> np.ndarray:
    r"""Distance distribution between two separated 3-D Gaussians.

    For two spherical Gaussian probability clouds whose centres are
    separated by *separation_distance* and whose widths are *sigma*,
    this gives the PDF of the inter-centre distance.

    Parameters
    ----------
    distances : np.ndarray
        *r* axis.
    separation_distance : float
        Centre-to-centre separation.
    sigma : float
        Width (standard deviation) of each Gaussian.
    normalize : bool
        Normalise the result.

    Returns
    -------
    np.ndarray
    """
    if separation_distance > 0.0:
        pr = distances / separation_distance * (
            normal_distribution(distances, loc=separation_distance, scale=sigma, norm=False)
            - normal_distribution(distances, loc=-separation_distance, scale=sigma, norm=False)
        )
    else:
        pr = (
            2.0 * distances ** 2 / sigma ** 2
            * normal_distribution(distances, loc=0.0, scale=sigma, norm=False)
        )
    if normalize and pr.sum() > 0.0:
        pr /= pr.sum()
    return pr
