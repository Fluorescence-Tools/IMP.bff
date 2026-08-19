"""States: positions, weights and orientations -- what every representation supplies.

An accessible-volume grid point, a rotamer, a coarse-grained conformer and an MD
frame are the same kind of thing: *a state the dye can occupy, with a weight, a
position and possibly an orientation*. What differs is only how the states were
generated. Everything downstream -- distances, kappa^2, the interaction terms --
consumes states, so it is written once and works for all of them.

This is the abstraction that was missing. ``RotamerEnsemble`` inherited from the
concrete :class:`~IMP.bff.representation.AccessibleVolume` instead, which is why
distance code happened to work for rotamers: by inheritance, not by design. The
cost showed up in the fields a rotamer library then had to carry and could not
fill -- ``density=zeros((0, 0, 0))``, ``grid_step=0.0``, ``grid_shape=(0, 0, 0)``
-- placeholders for a grid that does not exist. No consumer ever read them:
every one of them used ``points``, ``mean_position``, ``n_points`` or
``has_volume``, which is exactly this surface.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional, TYPE_CHECKING, Tuple, Union
import abc
import math

import numpy as np

import IMP.bff

__all__ = [
    'AccessibleVolume',
    'States',
    'distance_between_gaussian',
    'gaussian_chain',
    'gaussian_chain_ree',
    'generalized_normal_distribution',
    'normal_distribution',
    'poisson_0toN',
    'worm_like_chain',
    'worm_like_chain_linker',
]

# --------------------------------------------------------------------------
# states
# --------------------------------------------------------------------------
"""States: positions, weights and orientations -- what every representation supplies.

An accessible-volume grid point, a rotamer, a coarse-grained conformer and an MD
frame are the same kind of thing: *a state the dye can occupy, with a weight, a
position and possibly an orientation*. What differs is only how the states were
generated. Everything downstream -- distances, kappa^2, the interaction terms --
consumes states, so it is written once and works for all of them.

This is the abstraction that was missing. ``RotamerEnsemble`` inherited from the
concrete :class:`~IMP.bff.representation.AccessibleVolume` instead, which is why
distance code happened to work for rotamers: by inheritance, not by design. The
cost showed up in the fields a rotamer library then had to carry and could not
fill -- ``density=zeros((0, 0, 0))``, ``grid_step=0.0``, ``grid_shape=(0, 0, 0)``
-- placeholders for a grid that does not exist. No consumer ever read them:
every one of them used ``points``, ``mean_position``, ``n_points`` or
``has_volume``, which is exactly this surface.
"""

@dataclass(kw_only=True)
class States:
    """A weighted set of states of one label.

    :param points: ``(N, 4)`` -- ``x, y, z, weight``. One row per state.
    :param attachment_point: ``(3,)`` where the label is tied to the structure.
    :param orientations: ``(N, 3)`` transition dipoles, when the representation
        resolves them. ``None`` for a positional-only model such as an AV, which
        is why kappa^2 from an AV needs an isotropic assumption and kappa^2 from
        a rotamer library does not.
    :param position_name: human-readable label for the site.
    :param params: how these states were produced -- representation parameters,
        not dye or site properties.
    """

    points: np.ndarray
    attachment_point: np.ndarray
    orientations: Optional[np.ndarray] = None
    position_name: str = ""
    params: Dict = field(default_factory=dict)

    @property
    def positions(self) -> np.ndarray:
        """``(N, 3)`` state coordinates."""
        return self.points[:, :3]

    @property
    def weights(self) -> np.ndarray:
        """``(N,)`` state weights, unnormalised."""
        return self.points[:, 3]

    @property
    def n_points(self) -> int:
        return self.points.shape[0] if self.points.ndim == 2 else 0

    @property
    def has_volume(self) -> bool:
        return self.n_points > 0

    @property
    def has_orientations(self) -> bool:
        return self.orientations is not None and len(self.orientations) > 0

    @property
    def mean_position(self) -> np.ndarray:
        """Weight-averaged position, falling back to the attachment point."""
        if self.n_points == 0:
            return np.asarray(self.attachment_point).copy()
        w = self.points[:, 3]
        if w.sum() == 0:
            return np.asarray(self.attachment_point).copy()
        return np.average(self.points[:, :3], axis=0, weights=w)


# --------------------------------------------------------------------------
# types
# --------------------------------------------------------------------------
"""The accessible volume: a region a dye can reach, as one of several representations."""

@dataclass(kw_only=True)
class AccessibleVolume(States):
    """States enumerated as a voxel grid, plus the grid itself.

    One definition. There were two identical ones -- ``IMP.bff.representation.av`` and
    ``IMP.bff.representation.av``, same seven fields in the same order -- which is how they
    came to disagree about the axis order of ``density`` without anything
    noticing (PRD-113 stage 3a).

    :param density: ``(nx, ny, nz)`` accessible density, in the **coordinate**
        axis order. IMP orders its flat tile values with *x* fastest; a C-order
        reshape into ``(nx, ny, nz)`` returns the volume transposed, and a
        mirrored volume has the right voxel count, bounding box and total volume,
        so only a voxel-by-voxel comparison against ``points`` catches it.
    :param grid_origin: ``(3,)`` coordinate of the first voxel centre.
    :param grid_step: voxel edge in Angstrom.
    :param grid_shape: ``(nx, ny, nz)``.
    """

    density: np.ndarray
    grid_origin: np.ndarray
    grid_step: float
    grid_shape: Tuple[int, int, int]


# --------------------------------------------------------------------------
# distance
# --------------------------------------------------------------------------
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


# --------------------------------------------------------------------------
# label_distribution
# --------------------------------------------------------------------------
"""Label distributions: AV-backed and Gaussian representations of where a dye is.

Renamed from ``distribution.py`` in the 2026-08-18 tidy, because it sat beside
``distributions.py`` -- one letter apart, in the same directory, meaning
different things: this one is *where a label is*, that one was a set of
probability density functions. That file is now ``probability.py``. A reader
should not have to open both to find out which is which.

Two representations, alongside the accessible volume itself and the rotamer
library: :class:`LabelDistributionAV` computes an AV and reduces it lazily, and
:class:`DyeDistributionNormal` replaces the cloud with a Gaussian.

Moved here from ``IMP.bff.label`` by PRD-113 stage 3d. That package now means
the *system* -- which dye is attached where -- and a label *distribution* is a
representation of where it can be, which is a different question. The same word
meant both, which is the kind of collision this restructure exists to remove.

.. note::
   These classes re-implement the
   :class:`~IMP.bff.representation.States` surface (``points``,
   ``mean_position``, ``n_points``) and carry a **fourth** copy of the distance
   layer (``dRmp``/``dRDA``/``dRDAE``/``pRDA``, duplicated across both concrete
   classes, and again in ``BasicAV`` and in ``fret/distance.py``). The
   :attr:`LabelDistribution.states` view below is the bridge; folding the four
   distance layers into one is PRD-113 stage 4, and is a behaviour change that
   does not belong in a move.
"""

if TYPE_CHECKING:  # names for annotations only; see _av_types() below
    from IMP.bff.representation.av import ACV, BasicAV


def _av_types():
    """``(BasicAV, ACV, compute_av)``, imported on first use rather than on import.

    The layering here runs ``representation.types`` -> ``av`` ->
    ``representation.label_distribution``: this module *builds* accessible volumes, so
    it sits above the builder, while the builder needs only the dataclass. Both
    edges are real. Written at module level they close a cycle -- and because
    importing ``IMP.bff.representation.types`` also executes the package
    ``__init__``, which imports this module, narrowing the other side does not
    break it. ``import IMP.bff.representation.av`` as a process's first import raised
    ImportError until this was deferred (PRD-113 stage 3).
    """
    from IMP.bff.representation.av import ACV, BasicAV, compute_av
    return BasicAV, ACV, compute_av

# ---------------------------------------------------------------------------
# Helper: find an atom in a coordinate array
# ---------------------------------------------------------------------------

def _find_atom_index(
    atoms_xyz: np.ndarray,
    atoms_vdw: np.ndarray,
    residue_seq_number: int,
    atom_name: str,
    chain_id: Optional[str] = None,
) -> int:
    """Find the index of an atom matching the given criteria.

    This is a simplified replacement for the chisurf
    ``chisurf.core.fio.structure.coordinates.get_atom_index`` routine.
    It first tries to match by all criteria; if *chain_id* is ``None``
    it matches only by residue number and atom name.

    Parameters
    ----------
    atoms_xyz : (N, 3) float64
    atoms_vdw : (N,) float64
    residue_seq_number : int
    atom_name : str
    chain_id : str, optional
    """
    for i in range(len(atoms_xyz)):
        # The minimal signature: residue_seq_number + atom_name.
        # In a real scenario the structured array would contain columns
        # ``residue_seq_number``, ``atom_name``, ``chain_id``.  Here we
        # simply return *i*; subclasses can override with a more
        # sophisticated lookup.
        return i
    return 0


# ---------------------------------------------------------------------------
# Abstract base
# ---------------------------------------------------------------------------

class LabelDistribution(abc.ABC):
    """Abstract base for a 3-D dye label distribution.

    Subclasses must implement :meth:`_compute_density`.

    Attributes
    ----------
    origin : (3,) ndarray
        Attachment site (reference point) in Å.
    density : (ng, ng, ng) ndarray or None
        Density grid (if computed).
    verbose : bool
    simulation_grid_resolution : float
        Grid spacing (Å).
    position_name : str
        Human-readable label.
    """

    density: Optional[np.ndarray] = None
    origin: Optional[np.ndarray] = None
    verbose: bool = True
    simulation_grid_resolution: float
    position_name: str

    def __init__(
        self,
        simulation_type: str = "AV1",
        origin: Optional[np.ndarray] = None,
        simulation_grid_resolution: float = 0.5,
        position_name: str = "",
        verbose: bool = False,
    ):
        self.simulation_type = simulation_type
        self.origin = origin
        self.simulation_grid_resolution = simulation_grid_resolution
        self.position_name = position_name
        self.verbose = verbose
        self._av: Optional["BasicAV"] = None

    @abc.abstractmethod
    def _compute(self):
        """Compute or recompute the underlying accessible volume."""
        ...

    def get_basic_av(self) -> "BasicAV":
        """Return (or create) the underlying ``BasicAV``.

        Returns
        -------
        BasicAV
        """
        if self._av is None:
            self._compute()
        return self._av  # type: ignore

    # Convenience accessors that delegate to BasicAV
    @property
    def points(self) -> np.ndarray:
        """Point cloud ``(N, 4)`` of the dye distribution."""
        return self.get_basic_av().points

    @property
    def mean_position(self) -> np.ndarray:
        """Weighted mean position ``(3,)``."""
        return self.get_basic_av().mean_position

    @property
    def states(self) -> "States":
        """This distribution as :class:`~IMP.bff.representation.States`.

        The representation-agnostic view: whatever produced the cloud, a
        consumer that wants positions and weights asks for this and works for
        an AV, a rotamer library, a Gaussian or an MD trajectory alike.
        """
        # (was: from .distance import ...) -- now in this module
        return States(points=self.points,
                      attachment_point=np.asarray(self.attachment_point)
                      if getattr(self, "attachment_point", None) is not None
                      else self.mean_position)

    @property
    def n_points(self) -> int:
        """Number of points in the dye distribution."""
        return self.get_basic_av().n_points


# ---------------------------------------------------------------------------
# AV-based label distribution
# ---------------------------------------------------------------------------

class LabelDistributionAV(LabelDistribution):
    """Label distribution computed via an Accessible Volume (AV).

    Parameters
    ----------
    atoms_xyz : (N, 3) float64
        Atomic coordinates of the host structure.
    atoms_vdw : (N,) float64
        Van der Waals radii (Å).
    linker_length : float
        Dye linker length (Å).
    linker_width : float
        Linker width (Å).
    dye_radii : tuple (r1, r2, r3)
        Dye-sphere radii for the AV1/AV3 model (Å).
    residue_seq_number : int
        Attachment residue sequence number.
    atom_name : str
        Attachment atom name (e.g. ``"CB"``).
    chain_id : str, optional
        Chain identifier.
    simulation_grid_resolution : float
        AV grid spacing (Å).
    position_name : str
        Optional human-readable label.
    verbose : bool
    """

    def __init__(
        self,
        atoms_xyz: np.ndarray,
        atoms_vdw: np.ndarray,
        linker_length: float = 20.0,
        linker_width: float = 0.5,
        dye_radii: tuple[float, float, float] = (3.5, 0.0, 0.0),
        residue_seq_number: int = 0,
        atom_name: str = "CB",
        chain_id: Optional[str] = None,
        simulation_grid_resolution: float = 1.5,
        position_name: str = "",
        verbose: bool = False,
    ):
        self.atoms_xyz = np.asarray(atoms_xyz, dtype=np.float64)
        self.atoms_vdw = np.asarray(atoms_vdw, dtype=np.float64)
        self.linker_length = float(linker_length)
        self.linker_width = float(linker_width)
        self.dye_radii = tuple(float(r) for r in dye_radii)
        self.residue_seq_number = residue_seq_number
        self.atom_name = atom_name
        self.chain_id = chain_id
        self._attachment_index = _find_atom_index(
            atoms_xyz, atoms_vdw,
            residue_seq_number, atom_name, chain_id,
        )

        super().__init__(
            simulation_type="AV1" if dye_radii[1] == 0.0 else "AV3",
            simulation_grid_resolution=simulation_grid_resolution,
            position_name=position_name,
            verbose=verbose,
        )
        self.origin = atoms_xyz[self._attachment_index].copy()
        # AV is lazily computed in get_basic_av()

    def _compute(self):
        """Compute the AV for this label."""
        if self._av is not None:
            return
        source_xyz = self.atoms_xyz[self._attachment_index]

        BasicAV, _, compute_av = _av_types()
        av_result = compute_av(
            self.atoms_xyz, self.atoms_vdw, source_xyz,
            linker_length=self.linker_length,
            linker_width=self.linker_width,
            dye_radii=self.dye_radii,
            grid_resolution=self.simulation_grid_resolution,
        )
        self.density = av_result.density
        self._av = BasicAV(
            points=av_result.points,
            density=av_result.density,
            grid_origin=av_result.grid_origin,
            grid_step=av_result.grid_step,
            position_name=self.position_name,
        )

    # Distance methods
    def dRmp(self, other: "LabelDistributionAV") -> float:
        """:math:`R_{\\mathrm{mp}}` distance to another label."""
        return self.get_basic_av().dRmp(other.get_basic_av())

    def dRDA(self, other: "LabelDistributionAV", n_samples: int = 50000) -> float:
        """:math:`\\langle R_{DA}\\rangle` mean distance."""
        return self.get_basic_av().dRDA(other.get_basic_av(), n_samples)

    def dRDAE(self, other: "LabelDistributionAV",
              forster_radius: float = 52.0, n_samples: int = 50000) -> float:
        """:math:`R_E` FRET-averaged distance."""
        return self.get_basic_av().dRDAE(
            other.get_basic_av(), forster_radius, n_samples
        )

    def pRDA(self, other: "LabelDistributionAV",
             axis: Optional[np.ndarray] = None,
             n_samples: int = 50000) -> tuple[np.ndarray, np.ndarray]:
        """Distance distribution :math:`p(R_{DA})`."""
        return self.get_basic_av().pRDA(
            other.get_basic_av(), axis, n_samples
        )


# ---------------------------------------------------------------------------
# Normal (Gaussian) dye distribution — no structure needed
# ---------------------------------------------------------------------------

class DyeDistributionNormal(LabelDistribution):
    """Gaussian (normal) dye distribution around a point.

    This distribution does not require a structure; the dye is modelled
    as a 3-D isotropic Gaussian centred at a given point.

    Parameters
    ----------
    origin : (3,) ndarray
        Mean position (Å).
    width : float
        Standard deviation (Å) in each dimension.
    position_name : str
        Optional human-readable label.
    verbose : bool
    """

    def __init__(
        self,
        origin: np.ndarray,
        width: float = 6.0,
        position_name: str = "",
        verbose: bool = False,
    ):
        self.width = float(width)
        super().__init__(
            origin=np.asarray(origin, dtype=np.float64),
            simulation_grid_resolution=1.0,
            position_name=position_name,
            verbose=verbose,
        )
        # Normal distributions don't use a grid; create a point cloud directly
        self._compute()

    def _compute(self):
        """Draw random points from the 3-D Gaussian."""
        n_pts = 50000
        pts = np.random.randn(n_pts, 4).astype(np.float64)
        pts[:, :3] = pts[:, :3] * self.width + self.origin
        # Weight: Gaussian height relative to the distribution centre
        centered = pts[:, :3] - self.origin
        r2 = np.sum(centered ** 2, axis=1)
        pts[:, 3] = np.exp(-0.5 * r2 / (self.width ** 2))
        pts[:, 3] /= pts[:, 3].sum()
        BasicAV, _, _ = _av_types()
        self._av = BasicAV(
            points=pts,
            position_name=self.position_name,
        )

    def dRmp(self, other: "DyeDistributionNormal") -> float:
        return self.get_basic_av().dRmp(other.get_basic_av())

    def dRDA(self, other: "DyeDistributionNormal",
             n_samples: int = 50000) -> float:
        return self.get_basic_av().dRDA(other.get_basic_av(), n_samples)

    def dRDAE(self, other: "DyeDistributionNormal",
              forster_radius: float = 52.0,
              n_samples: int = 50000) -> float:
        return self.get_basic_av().dRDAE(
            other.get_basic_av(), forster_radius, n_samples
        )

    def pRDA(self, other: "DyeDistributionNormal",
             axis: Optional[np.ndarray] = None,
             n_samples: int = 50000) -> tuple[np.ndarray, np.ndarray]:
        return self.get_basic_av().pRDA(
            other.get_basic_av(), axis, n_samples
        )


# --------------------------------------------------------------------------
# polymer
# --------------------------------------------------------------------------
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

Filed under ``representation`` in the PRD-113 cleanup: a linker's end-to-end
distribution *is* a representation of where the dye can be -- the analytic
counterpart of what an accessible volume computes geometrically. It sat at the
package root, which said nothing about that.
"""

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


# --------------------------------------------------------------------------
# probability
# --------------------------------------------------------------------------
"""Probability distributions used by the dye and linker models.

Renamed from ``distributions.py``: it sat one letter away from
``distribution.py`` in the same directory, and the two are unrelated. That one
is now ``label_distribution.py``.

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
