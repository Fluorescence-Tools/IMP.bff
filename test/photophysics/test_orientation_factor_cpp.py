"""kappa^2: the geometry, the wobbling-in-a-cone distribution, and the 2/3 check.

The orientation factor is the dominant systematic in a FRET distance, and it
has exactly one analytic handle: **the isotropic average is 2/3**. That handle
is what these tests are built on, and it is what caught a real defect during the
C++ port -- ``kappasq_all`` sampled dipole directions with ``np.random.random(3)``,
which fills the positive octant of the unit cube rather than the sphere, and
returned 0.333 where the rigid isotropic limit is 0.667.

Two limits pin the sampler from opposite ends:

* ``sD2 = sA2 = 0`` -- freely rotating dyes have no orientation preference, so
  *every* sample is exactly 2/3 and the spread is zero.
* ``sD2 = sA2 = 1`` -- rigid dyes at random mutual orientations, so individual
  values span [0, 4] and only the *mean* is 2/3.

A sampler that passes one and fails the other is wrong in a way a single mean
would hide.
"""

import numpy as np
import pytest

import IMP.bff
import IMP.bff.photophysics as o


# --- the geometry ------------------------------------------------------------

def test_kappa_reproduces_its_documented_value():
    """Perpendicular dipoles offset along y -- the case in ``kappa``'s docstring."""
    donor = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    acceptor = np.array([[0.0, 0.5, 0.0], [0.0, 0.5, 1.0]])
    d, k = o.kappa(donor, acceptor)
    assert d == pytest.approx(0.8660254037844386)
    assert k == pytest.approx(1.0000000000000002)


def test_collinear_dipoles_give_kappa_minus_two():
    """Head to tail along the separation vector: kappa = 1 - 3 = -2, kappa^2 = 4."""
    d, k = o.kappa_distance(np.array([-0.5, 0.0, 0.0]), np.array([0.5, 0.0, 0.0]),
                            np.array([19.5, 0.0, 0.0]), np.array([20.5, 0.0, 0.0]))
    assert d == pytest.approx(20.0)
    assert k == pytest.approx(-2.0)
    assert k ** 2 == pytest.approx(4.0)


def test_parallel_dipoles_perpendicular_to_r_give_kappa_one():
    d, k = o.kappa_distance(np.array([0.0, -0.5, 0.0]), np.array([0.0, 0.5, 0.0]),
                            np.array([30.0, -0.5, 0.0]), np.array([30.0, 0.5, 0.0]))
    assert d == pytest.approx(30.0)
    assert k == pytest.approx(1.0)


def test_kappa_is_signed_and_squaring_loses_that():
    """Anti-parallel dipoles give -1, not 1. The sign is real; kappa^2 discards it."""
    _, k = o.kappa_distance(np.array([0.0, -0.5, 0.0]), np.array([0.0, 0.5, 0.0]),
                            np.array([30.0, 0.5, 0.0]), np.array([30.0, -0.5, 0.0]))
    assert k == pytest.approx(-1.0)


# --- the wobbling form -------------------------------------------------------

def test_wobbling_kappa2_reduces_to_two_thirds_when_nothing_is_ordered():
    for beta1 in (0.0, 0.7, np.pi / 2):
        for beta2 in (0.0, 1.3, np.pi):
            assert o.kappasq(0.9, 0.0, 0.0, beta1, beta2) == pytest.approx(2 / 3)


def test_wobbling_kappa2_matches_the_bare_geometry_when_fully_rigid():
    """At sD2 = sA2 = 1 the cone collapses and eq. 9 must give the static kappa^2."""
    rng = np.random.default_rng(8)
    for _ in range(50):
        d1 = rng.normal(size=3); d1 /= np.linalg.norm(d1)
        d2 = rng.normal(size=3); d2 /= np.linalg.norm(d2)
        delta = np.arccos(np.clip(d1 @ d2, -1, 1))
        beta1 = np.arccos(np.clip(d1[0], -1, 1))
        beta2 = np.arccos(np.clip(d2[0], -1, 1))
        # R_DA along x, which is the convention the module documents
        static = (d1 @ d2 - 3.0 * d1[0] * d2[0]) ** 2
        assert o.kappasq(delta, 1.0, 1.0, beta1, beta2) == pytest.approx(static, abs=1e-9)


# --- the distributions -------------------------------------------------------

def test_free_dyes_give_exactly_two_thirds_every_sample():
    _, _, k2 = o.kappasq_all(0.0, 0.0, n_samples=20000, seed=2)
    assert np.ptp(k2) == 0.0
    assert k2[0] == pytest.approx(2 / 3)


def test_rigid_dyes_average_two_thirds_and_span_the_full_range():
    """The check that caught the octant-sampling defect. Do not loosen it."""
    _, _, k2 = o.kappasq_all(1.0, 1.0, n_bins=81, n_samples=400000, seed=1)
    assert k2.mean() == pytest.approx(2 / 3, abs=0.01)
    assert k2.min() < 0.01 and k2.max() > 3.9, "the sampler is not reaching the sphere"


def test_samples_stay_within_the_physical_bounds():
    _, _, k2 = o.kappasq_all(0.9, 0.9, n_samples=100000, seed=3)
    assert k2.min() >= 0.0 and k2.max() <= 4.0


def test_distribution_is_reproducible_and_seed_dependent():
    a = o.kappasq_all(0.5, 0.4, n_samples=5000, seed=7)[2]
    np.testing.assert_array_equal(a, o.kappasq_all(0.5, 0.4, n_samples=5000, seed=7)[2])
    assert not np.array_equal(a, o.kappasq_all(0.5, 0.4, n_samples=5000, seed=8)[2])


def test_histogram_and_scale_shapes_are_consistent():
    scale, hist, k2 = o.kappasq_all(0.3, 0.3, n_bins=81, n_samples=1000, seed=0)
    assert scale.shape == (81,) and hist.shape == (80,)
    assert hist.sum() == 1000, "every sample must land in a bin, including the top edge"
    assert k2.shape == (1000,)


def test_delta_distribution_is_solid_angle_weighted():
    """Each beta1 ring contributes sin(beta1); without it the poles are over-counted."""
    scale, hist, k2 = o.kappasq_all_delta(delta=0.2, sD2=0.15, sA2=0.25, step=2.0, n_bins=31)
    assert scale.shape == (31,) and hist.shape == (30,)
    assert k2.shape == (45, 180)
    n_beta, n_phi = k2.shape
    weights = np.sin(np.arange(0.001, np.pi / 2, 2.0 * np.pi / 180.0)[:n_beta])
    assert hist.sum() == pytest.approx(n_phi * weights.sum(), rel=1e-9)


def test_delta_distribution_is_deterministic():
    a = o.kappasq_all_delta(0.2, 0.15, 0.25, 2.0, 31)
    b = o.kappasq_all_delta(0.2, 0.15, 0.25, 2.0, 31)
    for x, y in zip(a, b):
        np.testing.assert_array_equal(x, y)


def test_no_numba_left_in_orientation():
    import ast
    import inspect
    tree = ast.parse(inspect.getsource(o))
    imported = {n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.module}
    imported |= {a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
    assert not any("numba" in m or m.endswith("_jit") for m in imported), imported

    # "no decorated function anywhere" was the original proxy for "@njit is
    # gone". It stopped meaning that when photophysics/ merged into one module
    # and brought `@abc.abstractmethod` with it. What the proxy was reaching
    # for is a compilation decorator, so say that instead -- an abstract method
    # is not evidence of numba.
    compiled = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.FunctionDef):
            continue
        for dec in node.decorator_list:
            text = ast.unparse(dec)
            if any(w in text for w in ("njit", "jit", "vectorize", "guvectorize")):
                compiled.append((node.name, text))
    assert not compiled, compiled


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))


# --- the kappa^2 distance-ratio convolution ---------------------------------

def test_the_convolution_matches_the_outer_product_it_replaced():
    """One path now, and it agrees with the array form to machine precision.

    ``convolve_distance_with_k2_ratio`` had two branches. The fast one called
    ``_fast_convolve_loop``, which **was never defined in this package** -- it
    came across in the kappa-squared migration as a call to something that
    stayed behind, so ``use_fast=True`` raised NameError for every input large
    enough to reach it. Nothing noticed, because the module had no consumers.

    The slow branch built the full outer product: ``n x m`` twice over, once for
    the products and once for the weights, purely to hand them to
    ``np.histogram``. The kernel bins without forming either, so there is
    nothing left for a flag to choose between and the flag is gone.
    """
    rng = np.random.default_rng(1)
    for n, m in ((40, 40), (200, 150)):
        r_da = np.linspace(30.0, 70.0, n)
        amp = rng.random(n); amp /= amp.sum()
        ratio = np.linspace(0.7, 1.4, m)
        w = rng.random(m); w /= w.sum()

        products = (r_da[:, None] * ratio[None, :]).ravel()
        weights = (amp[:, None] * w[None, :]).ravel()
        lo, hi = products.min(), products.max()
        want, edges = np.histogram(products, bins=32, range=(lo, hi), weights=weights)
        keep = want > 1e-10 * want.max()

        centres, hist = o.convolve_distance_with_k2_ratio(r_da, amp, ratio, w, n_bins=32)
        np.testing.assert_allclose(hist, want[keep], rtol=0, atol=1e-15)
        np.testing.assert_allclose(centres, (0.5 * (edges[:-1] + edges[1:]))[keep])


def test_the_convolution_conserves_weight():
    """Every product lies inside [min*min, max*max] by construction, so none may
    fall outside the histogram -- including the one exactly at the top edge."""
    r_da = np.linspace(30.0, 70.0, 60)
    ratio = np.linspace(0.7, 1.4, 60)
    amp = np.full(60, 1 / 60)
    w = np.full(60, 1 / 60)
    _, hist = o.convolve_distance_with_k2_ratio(r_da, amp, ratio, w, n_bins=64)
    assert hist.sum() == pytest.approx(1.0, rel=1e-12)


def test_the_dead_fast_flag_is_gone():
    import inspect
    assert "use_fast" not in inspect.signature(o.convolve_distance_with_k2_ratio).parameters


def test_a_single_point_convolution_is_a_single_bin():
    """One product value is a delta, and its range has zero width.

    numpy widens ``(a, a)`` to ``(a - 0.5, a + 0.5)`` rather than returning
    nothing, and that is the contract here too. The weight lands in one bin;
    its centre can only be within half a bin width of the value, which is what
    a histogram is.
    """
    n_bins = 8
    c, h = o.convolve_distance_with_k2_ratio(
        np.array([50.0]), np.array([1.0]), np.array([1.1]), np.array([1.0]),
        n_bins=n_bins)
    assert h.sum() == pytest.approx(1.0)
    assert h.size == 1, "a delta occupies exactly one bin"
    assert abs(c[0] - 55.0) <= 1.0 / n_bins


def test_empty_input_returns_empty():
    c, h = o.convolve_distance_with_k2_ratio(
        np.empty(0), np.empty(0), np.array([1.0]), np.array([1.0]))
    assert c.size == 0 and h.size == 0
