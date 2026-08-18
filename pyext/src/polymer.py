"""End-to-end distance distributions of ideal and worm-like chains.

The linker between an attachment point and a dye is a short polymer, and its
end-to-end distribution is what an accessible volume approximates
geometrically. These give it analytically.

**The numerics are C++** (:file:`include/IMP/bff/PolymerChain.h`); these are
wrappers that keep the Python signatures and return numpy arrays. Ported under
PRD-113 -- numba is a prototyping tool in this package, not a runtime
dependency.

The port found a live breakage. ``worm_like_chain_linker`` was numba-jitted and
called ``normal_distribution``, which had just become a C++ delegation that
numba cannot type -- so the function raised ``TypingError`` on any call, and
**nothing in the suite noticed**, because it has no test. The C++ version was
checked against an independent numpy convolution instead: agreement to 7e-18.
"""

from __future__ import annotations

import numpy as np

import IMP.bff

__all__ = [
    "gaussian_chain_ree",
    "gaussian_chain",
    "worm_like_chain",
    "worm_like_chain_linker",
]


def _axis(x) -> np.ndarray:
    return np.ascontiguousarray(x, dtype=np.float64).ravel()


def gaussian_chain_ree(segment_length: float, number_of_segments: int) -> float:
    r"""RMS end-to-end distance of an ideal chain, :math:`b\sqrt{N}`."""
    return float(IMP.bff.gaussian_chain_ree(
        float(segment_length), int(number_of_segments)))


def gaussian_chain(
    distances: np.ndarray, segment_length: float, number_of_segments: int
) -> np.ndarray:
    r"""Radial distribution of an ideal chain.

    :math:`P(r) = 4\pi r^2 (3/2\pi\langle r^2\rangle)^{3/2}
    \exp(-3r^2/2\langle r^2\rangle)`. Not normalised on the given axis unless
    that axis covers the probability mass.
    """
    return np.asarray(IMP.bff.gaussian_chain(
        _axis(distances), float(segment_length), int(number_of_segments)),
        dtype=np.float64)


def worm_like_chain(
    distances: np.ndarray,
    kappa: float,
    chain_length: float = 0.0,
    normalize: bool = True,
    distance: bool = True,
) -> np.ndarray:
    r"""Radial distribution of a worm-like chain.

    The multi-piece analytical solution of Becker, Rosa & Everaers (Eur Phys J E
    32:53-69, 2010); :math:`\kappa` is the dimensionless persistence-length
    ratio and the expression branches at :math:`\kappa = 0.125`.

    :param chain_length: contour length; 0 takes the largest ``r`` on the axis.
    :param distance: multiply by :math:`r^2`, giving a distance distribution
        rather than a density in space.

    Values at or beyond the contour length stay zero -- a chain cannot be longer
    than itself, and the closed form diverges there.
    """
    return np.asarray(IMP.bff.worm_like_chain(
        _axis(distances), float(kappa), float(chain_length),
        bool(normalize), bool(distance)), dtype=np.float64)


def worm_like_chain_linker(
    distances: np.ndarray,
    kappa: float,
    chain_length: float = 0.0,
    sigma: float = 6.0,
    normalize: bool = True,
) -> np.ndarray:
    r"""Worm-like chain broadened by the dye linkers at each end.

    Convolves :func:`worm_like_chain` with a Gaussian of width *sigma*: the
    chain distribution is between the *attachment points*, and what a FRET
    experiment measures is between the *dyes*.
    """
    return np.asarray(IMP.bff.worm_like_chain_linker(
        _axis(distances), float(kappa), float(chain_length), float(sigma),
        bool(normalize)), dtype=np.float64)
