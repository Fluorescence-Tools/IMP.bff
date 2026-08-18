"""Probability distributions used by the dye and linker models.

**The numerics are C++** (:file:`include/IMP/bff/Distributions.h`). numba is a
prototyping tool in this package, not a runtime dependency, so these are thin
wrappers that hand numpy arrays to the compiled kernels and hand numpy arrays
back. Ported under PRD-113; each function was checked against the numba version
it replaces -- seven of nine cases bit-for-bit identical, the other two within
7e-18, which is the last bit and comes from a different summation order in the
normalisation.

The Python signatures are unchanged, so callers do not know the difference.

.. note::
   The flat names ``IMP.bff.normal_distribution`` and friends resolve to the
   **C++** functions -- SWIG binds them into the package namespace directly, so
   they are the public surface and are deliberately absent from ``api.py``. The
   wrappers here return ``numpy`` arrays rather than SWIG vectors, which is what
   the code inside this package wants.

Filed under ``representation`` in the PRD-113 cleanup: these are the shapes a
label distribution takes. They sat at the package root, outside every domain.
"""

from __future__ import annotations

import numpy as np

import IMP.bff

__all__ = [
    "poisson_0toN",
    "normal_distribution",
    "generalized_normal_distribution",
    "distance_between_gaussian",
]


def poisson_0toN(lam: float, N: int) -> np.ndarray:
    r"""Poisson probabilities for :math:`k = 0 \dots N-1`.

    Uses the recursion :math:`p_0 = e^{-\lambda},\ p_k = p_{k-1}\lambda/k`
    rather than a factorial, which overflows long before the probabilities stop
    mattering.

    :param lam: rate parameter.
    :param N: number of terms.
    :returns: ``(N,)``.

    >>> import numpy as np
    >>> np.round(poisson_0toN(0.2, 5), 6)
    array([0.818731, 0.163746, 0.016375, 0.001091, 0.000055])
    """
    return np.asarray(IMP.bff.poisson_0toN(float(lam), int(N)), dtype=np.float64)


def normal_distribution(
    x: np.ndarray, loc: float = 0.0, scale: float = 1.0, norm: bool = True
) -> np.ndarray:
    """Normal density on *x*.

    :param norm: divide by the sum, so a discretised density sums to one.
    """
    return np.asarray(
        IMP.bff.normal_distribution(
            np.ascontiguousarray(x, dtype=np.float64).ravel(),
            float(loc), float(scale), bool(norm)),
        dtype=np.float64)


def generalized_normal_distribution(
    x: np.ndarray,
    loc: float = 0.0,
    scale: float = 1.0,
    shape: float = 0.0,
    norm: bool = True,
) -> np.ndarray:
    r"""Normal density with a **skew**, applied by transforming the axis.

    :math:`z = -\log(1 - \kappa (x - \mu)/\sigma)/\kappa`, evaluated against the
    standard normal. ``shape = 0`` is the untransformed normal; positive skews
    left, negative right. Not the exponential-power family, despite the name.
    """
    return np.asarray(
        IMP.bff.generalized_normal_distribution(
            np.ascontiguousarray(x, dtype=np.float64).ravel(),
            float(loc), float(scale), float(shape), bool(norm)),
        dtype=np.float64)


def distance_between_gaussian(
    distances: np.ndarray,
    separation_distance: float,
    sigma: float,
    normalize: bool = False,
) -> np.ndarray:
    r"""Distance distribution between two isotropic 3-D Gaussians.

    :math:`p(r) = (r/d)[N(r; d, \sigma) - N(r; -d, \sigma)]`, degenerating at
    :math:`d = 0` to the Maxwell form :math:`2r^2/\sigma^2 \cdot N(r; 0, \sigma)`
    -- which is what the separate branch exists to avoid dividing by.
    """
    return np.asarray(
        IMP.bff.distance_between_gaussian(
            np.ascontiguousarray(distances, dtype=np.float64).ravel(),
            float(separation_distance), float(sigma), bool(normalize)),
        dtype=np.float64)
