"""One reduction behind four statistics, and the two disagreements it did not resolve.

``distance_sample_statistics`` replaced three separate reductions -- in
``distance_metrics``, in ``representation.distance``, and inside
``av/_kernels`` -- that computed the same four numbers and disagreed at the
limits. The tests here pin the limits, because that is where they differed.

The two *documented* disagreements are deliberately still here. Both are
questions of convention that need an owner, not a refactor:

* ``gaussian_rmp_to_rda_mean`` is ``Rmp + s^2/Rmp`` in one module and
  ``Rmp + s^2/(2 Rmp)`` in the other -- a factor of two that turns on what
  ``sigma`` names.
* ``polynomial_transfer`` reads its coefficients ascending in one and
  descending in the other.

The second is now single-sourced: one Horner evaluator, and the ascending
caller reverses. The first is not touched at all. A test that asserted them
equal would be asserting an answer nobody has given.
"""

import math

import numpy as np
import pytest

import IMP.bff
import IMP.bff.distance_metrics as dm
from IMP.bff.representation import distance as rd


def _reference(d, w, r0=52.0):
    """The four statistics, computed independently in numpy."""
    ws = w.sum()
    mu = (d * w).sum() / ws
    e = (w / (1.0 + (d / r0) ** 6)).sum() / ws
    return (mu,
            r0 * (1.0 / e - 1.0) ** (1.0 / 6.0),
            e,
            math.sqrt(max((d * d * w).sum() / ws - mu * mu, 0.0)))


@pytest.mark.parametrize("n", [10, 997, 50000])
def test_matches_an_independent_numpy_reduction(n):
    rng = np.random.default_rng(n)
    d = rng.uniform(5.0, 140.0, n)
    w = rng.random(n)
    got = IMP.bff.distance_sample_statistics(d, w, 52.0)
    np.testing.assert_allclose(got, _reference(d, w), rtol=1e-12, atol=1e-9)


def test_zero_total_weight_gives_zeros_not_nan():
    """A weightless sample has no statistics; it must not return NaN silently."""
    got = IMP.bff.distance_sample_statistics(np.array([1.0, 2.0]), np.zeros(2), 52.0)
    assert list(got) == [0.0, 0.0, 0.0, 0.0]


def test_delta_sample_has_zero_width_and_R_E_equals_the_distance():
    got = IMP.bff.distance_sample_statistics(np.full(64, 50.0), np.ones(64), 52.0)
    assert got[0] == pytest.approx(50.0)
    assert got[1] == pytest.approx(50.0)
    assert got[3] == 0.0, "the one-pass variance must be clamped, not left slightly negative"


def test_variance_clamp_survives_a_narrow_distribution():
    """``<r^2> - <r>^2`` loses cancellation on a tight, far-out sample.

    Without the clamp this is a small negative number and the width comes back
    NaN, which then propagates through every caller without an error.
    """
    d = 1000.0 + np.full(4096, 1e-9)
    got = IMP.bff.distance_sample_statistics(d, np.ones(4096), 52.0)
    assert not math.isnan(got[3])
    assert got[3] >= 0.0


def test_fret_averaged_distance_limits():
    close = IMP.bff.distance_sample_statistics(np.zeros(8), np.ones(8), 52.0)
    assert close[1] == 0.0, "mean efficiency 1 means zero separation"
    # A merely distant pair does *not* reach the limit: at 1e9 A the efficiency
    # is 1e-43, small but representable, and R_E comes back at the separation.
    distant = IMP.bff.distance_sample_statistics(np.full(8, 1e9), np.ones(8), 52.0)
    assert distant[1] == pytest.approx(1e9, rel=1e-9)
    # The limit needs (r/R0)^6 to overflow, which is the only way <E> is exactly 0.
    far = IMP.bff.distance_sample_statistics(np.full(8, 1e60), np.ones(8), 52.0)
    assert far[2] == 0.0
    assert math.isinf(far[1]), "mean efficiency 0 means no transfer at all"


def test_R_E_is_never_longer_than_the_mean_distance():
    rng = np.random.default_rng(4)
    d = rng.uniform(20.0, 100.0, 20000)
    w = rng.random(20000)
    mean_r, r_e, _, _ = IMP.bff.distance_sample_statistics(d, w, 52.0)
    assert r_e <= mean_r


# --- the calibration evaluator ----------------------------------------------

def test_horner_reads_coefficients_highest_power_first():
    assert IMP.bff.polynomial_transfer(45.0, np.array([0.0, 1.0, 0.02])) == pytest.approx(45.02)
    assert IMP.bff.polynomial_transfer(2.0, np.array([1.0, 0.0, 0.0])) == pytest.approx(4.0)
    assert IMP.bff.polynomial_transfer(7.0, np.array([])) == 0.0


def test_the_two_coefficient_orders_are_one_evaluator():
    coeffs = np.array([0.0, 1.0, 0.02])
    assert rd.polynomial_transfer(45.0, coeffs) == pytest.approx(45.02)
    assert dm.polynomial_transfer(45.0, coeffs) == pytest.approx(85.5)
    # ...and reversing the input makes them agree, which is the whole difference
    assert dm.polynomial_transfer(45.0, coeffs[::-1]) == pytest.approx(
        rd.polynomial_transfer(45.0, coeffs))


def test_polynomial_transfer_is_vectorised_and_shape_preserving():
    x = np.linspace(0.0, 100.0, 24).reshape(4, 6)
    coeffs = np.array([1.0, -2.0, 0.5])
    y = rd.polynomial_transfer(x, coeffs)
    assert y.shape == x.shape
    np.testing.assert_allclose(
        y.ravel(), [rd.polynomial_transfer(float(v), coeffs) for v in x.ravel()], atol=1e-9)


def test_the_sigma_convention_is_still_an_open_question():
    """Not a bug being tolerated -- a decision that has not been made.

    If this ever fails, someone answered it; update the note in
    ``distance_metrics``' module docstring rather than deleting the test.
    """
    assert dm.gaussian_rmp_to_rda_mean(45.0, 6.0) == pytest.approx(45.8)
    assert rd.gaussian_rmp_to_rda_mean(45.0, 6.0) == pytest.approx(45.4)


def test_chi2_score_is_one_function_now():
    assert dm.chi2_score is rd.chi2_score
    assert dm.chi2_score(45.0, 40.0, 2.0, 5.0) == pytest.approx(1.0)
    assert dm.chi2_score(35.0, 40.0, 5.0, 2.0) == pytest.approx(1.0)
    assert dm.chi2_score(50.0, 50.0, 0.0, 0.0) == 0.0, "a zero error must not divide"


def test_no_numba_left_in_the_distance_layer():
    import ast
    import inspect
    for mod in (dm, rd):
        tree = ast.parse(inspect.getsource(mod))
        imported = {n.module for n in ast.walk(tree)
                    if isinstance(n, ast.ImportFrom) and n.module}
        assert not any("numba" in m or m.endswith("_jit") for m in imported), (mod, imported)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
