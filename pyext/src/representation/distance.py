"""Distances between two labels, whatever represents them.

Takes anything with :class:`~IMP.bff.representation.States` -- an accessible
volume, a rotamer library, a Gaussian, a coarse-grained ensemble -- so one
implementation serves them all. That is the point: this was the third of **six**
implementations of dye-pair distances in the package.

The six, and what became of them (PRD-113 stage 4b):

* ``fret/distance.py`` -- this module, moved here as the canonical one. It was
  the most complete: pair statistics, the RDA histogram, the transfer
  polynomial, the orientation-resolved pair geometry.
* ``distance_metrics.py`` -- six identically-named functions; now delegates.
* ``av/_kernels.py`` with ``av/basic.py`` -- the numba pair, **still separate**,
  because routing them changes an RNG stream and that is a behaviour change.
* ``representation/label_distribution.py`` -- ``dRmp``/``dRDA``/``dRDAE``/``pRDA`` on
  each of two concrete classes; still separate, same reason.
* the **C++** ``av_distance`` family in ``AV.h``, which takes ``AV`` decorators
  rather than point arrays and is the implementation that should ultimately
  win. Reaching it needs SWIG typemaps for a *second* 2-D input array.

They were checked against each other before any of this: on T4L A132 x A65,
2180 x 3087 points, 2e5 samples, ``<R_DA>`` came out 51.945 against 51.940 A and
``<R_DA>_E`` 51.696 against 51.692 A, with ``Rmp`` identical. **They agree** --
the spread is Monte-Carlo sampling noise -- so this is de-duplication with no
defect hiding in it, unlike the transposed density found in stage 3a.
"""
from __future__ import annotations

from typing import Optional, Tuple, Union

import math

import numpy as np

import IMP.bff

from .types import AccessibleVolume
from IMP.bff.photophysics.kappa2 import kappa2_from_dipoles

N_DISTANCE_SAMPLES: int = 50000

# The pre-move implementation jitted these two kernels with numba. numba is
# not allowed as a dependency of IMP.bff, and vectorised numpy makes the loop
# unnecessary anyway — same Monte-Carlo estimator, no compiler.


def _random_distances(
    p1: np.ndarray,
    p2: np.ndarray,
    n_samples: int,
    seed: int = 0,
) -> np.ndarray:
    """Draw random distance-weight pairs from two AV point clouds.

    Returns (n_samples, 2) array: col 0 = distance, col 1 = weight product.
    """
    rng = np.random.RandomState(seed)
    i1 = rng.randint(0, p1.shape[0], size=n_samples)
    i2 = rng.randint(0, p2.shape[0], size=n_samples)
    diff = p1[i1, :3] - p2[i2, :3]
    result = np.empty((n_samples, 2), dtype=np.float64)
    result[:, 0] = np.sqrt(np.einsum("ij,ij->i", diff, diff))
    result[:, 1] = p1[i1, 3] * p2[i2, 3]
    return result


def _sample_vectors(
    p1: np.ndarray,
    p2: np.ndarray,
    n_samples: int,
    seed: int = 0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Draw random vector-weight pairs from two AV point clouds.

    Parameters
    ----------
    p1 : np.ndarray
        Points in the first accessible volume (N, 4).
    p2 : np.ndarray
        Points in the second accessible volume (M, 4).
    n_samples : int
        Number of random samples to draw.
    seed : int, optional
        Random seed. Default is 0.

    Returns
    -------
    v : np.ndarray
        Sampled vectors from av1 to av2 (n_samples, 3).
    w : np.ndarray
        Sampled weight products (n_samples,).
    """
    rng = np.random.RandomState(seed)
    i1 = rng.randint(0, p1.shape[0], size=n_samples)
    i2 = rng.randint(0, p2.shape[0], size=n_samples)
    v = p2[i2, :3] - p1[i1, :3]
    w = p1[i1, 3] * p2[i2, 3]
    return v, w


def _sample_av_distance(
    av1: AccessibleVolume,
    av2: AccessibleVolume,
    n_samples: int = N_DISTANCE_SAMPLES,
) -> np.ndarray:
    """Sample distances between two AVs."""
    if av1.n_points == 0 or av2.n_points == 0:
        raise ValueError("Cannot sample distance: one or both AVs have no points")
    return _random_distances(av1.points, av2.points, n_samples)


def average_distance(
    av1: AccessibleVolume,
    av2: AccessibleVolume,
    n_samples: int = N_DISTANCE_SAMPLES,
) -> float:
    """Calculate the mean inter-dye distance <R_DA>.

    Uses weighted random sampling from both AV point clouds.
    """
    d = _sample_av_distance(av1, av2, n_samples)
    return float(np.dot(d[:, 0], d[:, 1]) / d[:, 1].sum())


def mean_fret_distance(
    av1: AccessibleVolume,
    av2: AccessibleVolume,
    forster_radius: float = 52.0,
    n_samples: int = N_DISTANCE_SAMPLES,
) -> float:
    """Calculate the FRET-averaged distance R_E.

    R_E = R0 * (1/<E> - 1)^(1/6)
    where <E> is the weighted mean FRET efficiency over the AV distribution.
    """
    d = _sample_av_distance(av1, av2, n_samples)
    r = d[:, 0]
    w = d[:, 1]
    e = 1.0 / (1.0 + (r / forster_radius) ** 6.0)
    mean_e = np.dot(w, e) / w.sum()
    return float(forster_radius * (1.0 / mean_e - 1.0) ** (1.0 / 6.0))


def distance_between_mean_positions(
    av1: AccessibleVolume,
    av2: AccessibleVolume,
) -> float:
    """Calculate the distance between the mean positions (Rmp)."""
    return float(np.sqrt(((av1.mean_position - av2.mean_position) ** 2).sum()))


def standard_deviation_of_distances(
    av1: AccessibleVolume,
    av2: AccessibleVolume,
    n_samples: int = N_DISTANCE_SAMPLES,
) -> float:
    """Calculate the standard deviation of the inter-AV distance distribution."""
    d = _sample_av_distance(av1, av2, n_samples)
    return float(IMP.bff.distance_sample_statistics(
        np.ascontiguousarray(d[:, 0]), np.ascontiguousarray(d[:, 1]))[3])


def av_pair_statistics(
    av1: AccessibleVolume,
    av2: AccessibleVolume,
    forster_radius: float = 52.0,
    n_samples: int = N_DISTANCE_SAMPLES,
) -> Tuple[float, float, float, float]:
    """Calculate distance statistics between two accessible volumes.

    Parameters
    ----------
    av1 : AccessibleVolume
        The first accessible volume.
    av2 : AccessibleVolume
        The second accessible volume.
    forster_radius : float, optional
        Förster radius for FRET-averaged distance calculation. Default is 52.0.
    n_samples : int, optional
        Number of random distance samples. Default is 50000.

    Returns
    -------
    rmp : float
        Distance between the mean positions (Rmp).
    rda_mean : float
        Mean inter-dye distance <R_DA>.
    rda_mean_e : float
        FRET-averaged distance R_E.
    sigma_r : float
        Standard deviation of the inter-AV distance distribution.
    """
    rmp = distance_between_mean_positions(av1, av2)
    if av1.n_points == 0 or av2.n_points == 0:
        return rmp, rmp, rmp, 0.0

    d = _sample_av_distance(av1, av2, n_samples)
    rda_mean, rda_mean_e, _mean_e, sigma_r = IMP.bff.distance_sample_statistics(
        np.ascontiguousarray(d[:, 0]), np.ascontiguousarray(d[:, 1]), forster_radius)
    if d[:, 1].sum() <= 0:
        return rmp, rmp, rmp, 0.0
    return rmp, float(rda_mean), float(rda_mean_e), float(sigma_r)


def fit_transfer_polynomial(
    av1: AccessibleVolume,
    av2: AccessibleVolume,
    distance_type: str,
    forster_radius: float = 52.0,
    degree: int = 3,
    n_samples: int = 10000,
) -> np.ndarray:
    """Fit a polynomial relating Rmp to RDAMean or RDAMeanE by translating the AVs.

    Parameters
    ----------
    av1 : AccessibleVolume
        The first accessible volume.
    av2 : AccessibleVolume
        The second accessible volume.
    distance_type : str
        The target distance type ('RDAMean' or 'RDAMeanE').
    forster_radius : float, optional
        Förster radius for FRET-averaged distance calculation. Default is 52.0.
    degree : int, optional
        Degree of the polynomial. Default is 3.
    n_samples : int, optional
        Number of samples for fit. Default is 10000.

    Returns
    -------
    coeffs : np.ndarray
        Polynomial coefficients [c_0, c_1, c_2, ...].
    """
    rmp_initial = distance_between_mean_positions(av1, av2)
    if av1.n_points == 0 or av2.n_points == 0 or rmp_initial <= 1e-6:
        # Return identity polynomial: y = x
        coeffs = np.zeros(degree + 1)
        coeffs[-2] = 1.0  # slope = 1.0
        return coeffs

    # Draw point cloud relative vectors
    v, w = _sample_vectors(av1.points, av2.points, n_samples)
    w_sum = w.sum()
    if w_sum <= 0:
        coeffs = np.zeros(degree + 1)
        coeffs[-2] = 1.0
        return coeffs

    # Unit vector along the line connecting mean positions
    u = (av2.mean_position - av1.mean_position) / rmp_initial

    # Generate test offsets
    t_min = -min(rmp_initial - 5.0, 15.0)
    t_max = 20.0
    t_values = np.linspace(t_min, t_max, 7)

    rmp_values = rmp_initial + t_values
    r_eff_values = np.empty(len(t_values))

    # Pre-calculate dot product for each vector in sample with u
    v_dot_u = v @ u
    # Pre-calculate squared norm of each vector
    v_sq = np.sum(v ** 2, axis=1)

    for i, t in enumerate(t_values):
        # New distance for each pair is: sqrt(v_sq + 2*t*(v_dot_u) + t^2)
        d_new = np.sqrt(np.maximum(v_sq + 2.0 * t * v_dot_u + t**2, 1e-10))
        if distance_type == "RDAMeanE":
            e = 1.0 / (1.0 + (d_new / forster_radius) ** 6.0)
            mean_e = np.dot(w, e) / w_sum
            if mean_e <= 0:
                r_eff = float(np.dot(d_new, w) / w_sum)
            elif mean_e >= 1:
                r_eff = 0.0
            else:
                r_eff = float(forster_radius * (1.0 / mean_e - 1.0) ** (1.0 / 6.0))
        else:  # RDAMean
            r_eff = float(np.dot(d_new, w) / w_sum)
        r_eff_values[i] = r_eff

    # Fit polynomial
    coeffs = np.polyfit(rmp_values, r_eff_values, degree)
    return coeffs


def polynomial_transfer(
    rmp: Union[float, np.ndarray], coeffs: np.ndarray
) -> Union[float, np.ndarray]:
    """Evaluate polynomial transfer function using Horner's method.

    Parameters
    ----------
    rmp : float or np.ndarray
        Distance(s) to evaluate.
    coeffs : np.ndarray
        Polynomial coefficients (highest power first, i.e., from np.polyfit).

    Returns
    -------
    float or np.ndarray
        Evaluated value(s).
    """
    coeffs = np.asarray(coeffs, dtype=np.float64).ravel()
    if np.ndim(rmp) == 0:
        return IMP.bff.polynomial_transfer(float(rmp), coeffs)
    x = np.asarray(rmp, dtype=np.float64)
    y = np.asarray(IMP.bff.polynomial_transfer_vector(
        np.ascontiguousarray(x.ravel()), coeffs), dtype=np.float64)
    return y.reshape(x.shape)


def gaussian_rmp_to_rda_mean(
    rmp: Union[float, np.ndarray], sigma: float
) -> Union[float, np.ndarray]:
    r"""Correct a mean-position distance to the mean inter-dye distance.

    .. math::
       \langle R_{DA} \rangle \approx R_{mp} + \frac{\sigma^2}{R_{mp}}

    **sigma is the per-component width of the separation vector**, not of one
    dye cloud and not of the distance distribution. That is the convention
    settled on 2026-08-18, and it is the one the expansion gives: writing the
    separation as :math:`d = \Delta + \varepsilon` with
    :math:`\varepsilon \sim N(0, \sigma^2 I_3)`,

    .. math::
       |d| \approx |\Delta| + \varepsilon\cdot\hat\Delta
                 + \frac{|\varepsilon_\perp|^2}{2|\Delta|}

    and since :math:`\langle \varepsilon\cdot\hat\Delta \rangle = 0` while the
    **two** transverse components give
    :math:`\langle|\varepsilon_\perp|^2\rangle = 2\sigma^2`, the correction is
    :math:`\sigma^2/R_{mp}`.

    .. note::
       This returned :math:`R_{mp} + \sigma^2/(2 R_{mp})` until 2026-08-18 --
       half the correction, 45.4 against 45.8 at ``rmp=45, sigma=6``. The other
       form lived under the same name in ``IMP.bff.distance_metrics``; the two
       are now one function. Anything that fitted ``sigma_rda`` *through* this
       code absorbed the factor of two and its fitted sigma will be
       :math:`\sqrt2` too large.

    Parameters
    ----------
    rmp : float or np.ndarray
        Distance between the clouds' mean positions, Angstrom.
    sigma : float
        Per-component width of the separation vector, Angstrom.

    Returns
    -------
    float or np.ndarray
        Corrected distance(s), Angstrom. **Zero where ``rmp`` is zero or
        negative**: the expansion is in :math:`\sigma/R_{mp}` and says nothing
        at coincident mean positions. Clamping the denominator instead -- which
        the canonical version used to do -- returns 3.6e11 A at ``rmp = 0``,
        a number that then propagates as if it meant something.
    """
    r = np.asarray(rmp, dtype=np.float64)
    out = np.where(r > 0.0, r + (sigma ** 2) / np.where(r > 0.0, r, 1.0), 0.0)
    return float(out) if np.ndim(rmp) == 0 else out


def histogram_rda(
    av1: AccessibleVolume,
    av2: AccessibleVolume,
    rda_axis: Optional[np.ndarray] = None,
    rda_min: float = 1.0,
    rda_max: float = 200.0,
    n_rda_bins: int = 100,
    use_log: bool = False,
    n_samples: int = N_DISTANCE_SAMPLES,
    normalize: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute a histogram of the inter-AV distance distribution.

    Returns (histogram, bin_edges).
    """
    if rda_axis is None:
        if use_log:
            rda_axis = np.logspace(
                np.log10(rda_min), np.log10(rda_max), n_rda_bins, dtype=np.float64
            )
        else:
            rda_axis = np.linspace(rda_min, rda_max, n_rda_bins, dtype=np.float64)

    d = _sample_av_distance(av1, av2, n_samples)
    p = np.histogram(d[:, 0], bins=rda_axis, weights=d[:, 1])[0]
    if normalize:
        total = p.sum()
        if total > 0:
            p = p / total
    return p, rda_axis


def model_distance(
    av1: AccessibleVolume,
    av2: AccessibleVolume,
    distance_type: str,
    forster_radius: float = 52.0,
    n_samples: int = N_DISTANCE_SAMPLES,
) -> float:
    """Compute a model distance for the given distance type.

    Supported types: ``Rmp``, ``RDAMean``, ``RDAMeanE``.
    """
    if distance_type == "Rmp":
        return distance_between_mean_positions(av1, av2)
    elif distance_type == "RDAMean":
        return average_distance(av1, av2, n_samples=n_samples)
    elif distance_type == "RDAMeanE":
        return mean_fret_distance(av1, av2, forster_radius=forster_radius, n_samples=n_samples)
    else:
        raise ValueError(f"Unknown distance type: {distance_type}")


def chi2_score(
    model_distance: float,
    experimental_distance: float,
    error_neg: float,
    error_pos: float,
) -> float:
    """Asymmetric chi-squared contribution for one distance restraint.

    ``chi2 = (d_m - d_e)^2 / error^2`` with asymmetric errors.
    """
    delta = model_distance - experimental_distance
    err = error_neg if delta < 0 else error_pos
    if err <= 0:
        return 0.0
    return (delta / err) ** 2


def fret_efficiency(distance: float, forster_radius: float = 52.0) -> float:
    """Single-pair FRET efficiency."""
    return 1.0 / (1.0 + (distance / forster_radius) ** 6.0)


def distance_from_fret_efficiency(
    efficiency: float, forster_radius: float = 52.0
) -> float:
    """Convert FRET efficiency back to distance."""
    return forster_radius * (1.0 / efficiency - 1.0) ** (1.0 / 6.0)


# ---------------------------------------------------------------------------
# Pair distributions over two weighted point sets with dipoles (rotamer
# ensembles, AV clouds): the full N1 x N2 matrix, no sampling.
# ---------------------------------------------------------------------------

def fret_pair_geometry(
    points1: np.ndarray,
    weights1: np.ndarray,
    points2: np.ndarray,
    weights2: np.ndarray,
    mu1: Optional[np.ndarray] = None,
    mu2: Optional[np.ndarray] = None,
) -> dict:
    """Distances, κ² and pair weights over all (i, j) of two weighted point sets.

    ``points`` are (N, 3) centres (Å), ``weights`` (N,) normalised or not,
    ``mu`` optional (N, 3) transition-dipole vectors. Without dipoles κ² is
    the isotropic 2/3 everywhere (an AV cloud). Returns a dict with ``R``
    (N1, N2), ``kappa2`` (N1, N2), ``weight`` (N1, N2, normalised to 1) and
    ``kappa2_avg``.

    The separation vectors used to come back as ``r_vectors`` (N1, N2, 3).
    **Nothing ever read them** -- they were the largest allocation in the call,
    built for no consumer -- and the kernel now forms each one and discards it.
    """
    p1 = np.ascontiguousarray(np.asarray(points1, dtype=np.float64)[:, :3])
    p2 = np.ascontiguousarray(np.asarray(points2, dtype=np.float64)[:, :3])
    w1 = np.asarray(weights1, dtype=np.float64).ravel()
    w2 = np.asarray(weights2, dtype=np.float64).ravel()
    n1, n2 = p1.shape[0], p2.shape[0]

    weight = np.outer(w1, w2)
    total = weight.sum()
    if total > 0:
        weight = weight / total

    oriented = mu1 is not None and mu2 is not None
    empty = np.empty(0)
    # The kernel hands back a numpy view over its own buffer -- no conversion,
    # no copy. Returning a std::vector instead would make SWIG build one Python
    # float per element and numpy walk them back, 35-40 ns each: on a 400x350 pair
    # matrix that was 38 ms against 0.22 ms, for about 1 ms of arithmetic.
    packed = IMP.bff.fret_pair_matrices(
        p1.ravel(), p2.ravel(),
        np.ascontiguousarray(np.asarray(mu1, dtype=np.float64)).ravel() if oriented else empty,
        np.ascontiguousarray(np.asarray(mu2, dtype=np.float64)).ravel() if oriented else empty,
        int(n1), int(n2))
    r = packed[: n1 * n2].reshape(n1, n2)
    kappa2 = packed[n1 * n2:].reshape(n1, n2)

    return {
        "R": r,
        "kappa2": kappa2,
        "weight": weight,
        "kappa2_avg": float(np.sum(kappa2 * weight)),
    }


def fret_pair_efficiencies(
    geometry: dict,
    forster_radius: float,
    tau0: Optional[float] = None,
) -> dict:
    """FRET efficiencies of a pair geometry (see :func:`fret_pair_geometry`).

    ``forster_radius`` is R0 for κ² = 2/3 in the units of the geometry (Å);
    the κ² dependence is applied per pair (rate ∝ (3/2)κ² (R0/R)^6). Returns
    ``static`` (⟨E_ij⟩), ``dynamic1`` (E with ⟨κ²⟩, then averaged),
    ``dynamic2`` (rate-averaged), ``kappa2_avg``, the per-pair ``E`` and
    ``rate_ratio`` (k_FRET/k_rad) matrices and, with ``tau0`` (ns), the
    per-pair FRET rates ``k_fret`` (1/ns).
    """
    r = geometry["R"]
    kappa2 = geometry["kappa2"]
    weight = geometry["weight"]
    kappa2_avg = geometry["kappa2_avg"]
    n1, n2 = r.shape
    packed = np.asarray(IMP.bff.fret_pair_efficiency_matrices(
        np.ascontiguousarray(r, dtype=np.float64).ravel(),
        np.ascontiguousarray(kappa2, dtype=np.float64).ravel(),
        float(forster_radius)), dtype=np.float64)
    e_pair = packed[: n1 * n2].reshape(n1, n2)
    rate_ratio = packed[n1 * n2:].reshape(n1, n2)

    # dynamic1 uses the ensemble-averaged kappa2 rather than the per-pair one,
    # so it is a different expression and stays here -- one pass, no temporaries
    # worth moving.
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio6 = np.power(r / float(forster_radius), 6)
        e_dyn1 = 1.0 / (1.0 + 2.0 / 3.0 / kappa2_avg * ratio6)
    e_dyn1 = np.nan_to_num(e_dyn1, nan=1.0, posinf=1.0)
    rate_avg = float(np.sum(np.nan_to_num(rate_ratio, posinf=0.0) * weight))
    out = {
        "static": float(np.sum(e_pair * weight)),
        "dynamic1": float(np.sum(e_dyn1 * weight)),
        "dynamic2": rate_avg / (rate_avg + 1.0),
        "kappa2_avg": kappa2_avg,
        "E": e_pair,
        "rate_ratio": rate_ratio,
        "weight": weight,
        "R": r,
    }
    if tau0 is not None:
        out["k_fret"] = rate_ratio / float(tau0)
    return out


def fret_pair_distribution(
    points1: np.ndarray,
    weights1: np.ndarray,
    points2: np.ndarray,
    weights2: np.ndarray,
    *,
    forster_radius: float,
    mu1: Optional[np.ndarray] = None,
    mu2: Optional[np.ndarray] = None,
    tau0: Optional[float] = None,
) -> dict:
    """:func:`fret_pair_geometry` followed by :func:`fret_pair_efficiencies` (one call).

    The result carries the flattened per-pair ``R``, ``kappa2``, ``weight``,
    ``E`` (and ``k_fret`` with ``tau0``) plus the ``static``/``dynamic1``/
    ``dynamic2`` averages -- the FRET rate distribution of the pair.
    """
    geometry = fret_pair_geometry(points1, weights1, points2, weights2, mu1, mu2)
    out = fret_pair_efficiencies(geometry, forster_radius, tau0)
    out["kappa2"] = geometry["kappa2"]
    return out


# ---------------------------------------------------------------------------
# The other conventions. Folded in from ``IMP.bff.distance_metrics`` when that
# module was retired (PRD-113 cleanup): each of these is a genuinely different
# question from its neighbour above, and each used to answer to the *same name*
# in the other module. Naming them apart is what makes both usable.
# ---------------------------------------------------------------------------

def polynomial_transfer_ascending(
    rmp: Union[float, np.ndarray], coeffs: np.ndarray
) -> Union[float, np.ndarray]:
    r"""Evaluate a transfer polynomial whose coefficients run **lowest power first**.

    :math:`c_0 + c_1 x + c_2 x^2 + \dots`, the order a human writes by hand.
    :func:`polynomial_transfer` takes the opposite order -- highest power first,
    which is what ``np.polyfit`` returns and therefore what a *fitted*
    calibration carries.

    The two are not interchangeable and the difference is silent: at
    ``rmp = 45`` with ``[0, 1, 0.02]`` this gives **85.5** and
    :func:`polynomial_transfer` gives **45.02**. Both were called
    ``polynomial_transfer``, in two modules, for as long as both existed.

    One evaluator underneath: this reverses and delegates.
    """
    c = np.asarray(coeffs, dtype=np.float64).ravel()[::-1].copy()
    return polynomial_transfer(rmp, c)


def distance_sample_statistics(
    distances: np.ndarray,
    weights: np.ndarray,
    forster_radius: float = 52.0,
) -> Tuple[float, float, float, float]:
    """Distance statistics from a **sample of pair distances**.

    The sibling of :func:`av_pair_statistics`, which takes the two *clouds*.
    This one takes distances that have already been drawn -- from
    :func:`IMP.bff.representation.av._kernels.random_distances`, from a trajectory, or from a
    measurement -- and is the right entry point when the clouds are not at hand.

    :returns: ``(r_da_mean, r_mp, r_e, mean_efficiency)``.

    ``r_mp`` is **always NaN**, and deliberately. The distance between two
    clouds' mean positions is not a function of the distribution of pair
    distances: :math:`|\\langle a\\rangle - \\langle b\\rangle|` cannot be
    recovered from :math:`|a - b|`. Use :func:`mean_position_distance`, which
    takes the clouds.

    That slot returned ``r_da_mean`` until 2026-07-28, documented as "same as
    r_da_mean here" -- so a caller reading index 1 got a different quantity
    under the right name. Measured on 148l E15/E90: 47.65 A against the 51.53 A
    it returned, 8 % apart. NaN propagates instead of lying.
    """
    w = np.asarray(weights, dtype=np.float64).ravel()
    if w.sum() == 0.0:
        return 0.0, math.nan, 0.0, 0.0
    r_da_mean, r_e, mean_efficiency, _sigma = IMP.bff.distance_sample_statistics(
        np.ascontiguousarray(np.asarray(distances, dtype=np.float64).ravel()),
        np.ascontiguousarray(w), float(forster_radius))
    return float(r_da_mean), math.nan, float(r_e), float(mean_efficiency)


def mean_position_distance(
    points_a: np.ndarray,
    points_b: np.ndarray,
    weights_a: Optional[np.ndarray] = None,
    weights_b: Optional[np.ndarray] = None,
) -> float:
    """:math:`R_{mp}` between two **point clouds**, rather than two AV objects.

    The sibling of :func:`distance_between_mean_positions`, which takes
    accessible volumes. Written over arrays so a rotamer library, a
    coarse-grained ensemble or an MD frame reaches it too.

    Contributed from QuEst, where it is ``AV.dRmp``.

    :param points_a: ``(n, 3)`` cloud in Angstrom.
    :param points_b: ``(m, 3)`` cloud in Angstrom.
    :param weights_a: per-point weights; uniform if omitted.
    :param weights_b: per-point weights; uniform if omitted.
    """
    a = np.asarray(points_a, dtype=np.float64).reshape((-1, 3))
    b = np.asarray(points_b, dtype=np.float64).reshape((-1, 3))
    if a.shape[0] == 0 or b.shape[0] == 0:
        raise ValueError("Both point clouds must be non-empty to define R_mp.")
    mean_a = a.mean(axis=0) if weights_a is None else np.average(a, axis=0, weights=weights_a)
    mean_b = b.mean(axis=0) if weights_b is None else np.average(b, axis=0, weights=weights_b)
    return float(np.sqrt(((mean_a - mean_b) ** 2).sum()))
