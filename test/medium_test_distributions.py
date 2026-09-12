"""Tests for ``IMP.bff.representation.distance`` — mathematical distributions."""

from __future__ import annotations

import numpy as np
import pytest

from IMP.bff import (
    poisson_0toN,
    normal_distribution,
    generalized_normal_distribution,
    distance_between_gaussian,
)


class TestPoisson:
    """Poisson probability mass function."""

    def test_p0_is_exp_minus_lam(self):
        p = poisson_0toN(0.2, 5)
        assert p[0] == pytest.approx(np.exp(-0.2), rel=1e-12)

    def test_normalises(self):
        p = poisson_0toN(2.0, 100)
        assert p.sum() == pytest.approx(1.0, abs=1e-6)


class TestNormalDistribution:
    """Gaussian PDF."""

    def test_peak_at_mean(self):
        x = np.array([0.0])
        y = normal_distribution(x, loc=0.0, scale=1.0, norm=False)
        expected = 1.0 / np.sqrt(2.0 * np.pi)
        assert y[0] == pytest.approx(expected, rel=1e-9)

    def test_normalisation(self):
        x = np.linspace(-5, 5, 101)
        y = normal_distribution(x, loc=0.0, scale=1.0, norm=True)
        assert y.sum() == pytest.approx(1.0, abs=1e-3)

    def test_shifted_mean(self):
        x = np.array([0.0, 5.0])
        y = normal_distribution(x, loc=5.0, scale=1.0, norm=False)
        assert y[0] < y[1]


class TestGeneralizedNormal:
    """Skew-normal (generalised normal) PDF."""

    def test_shape_zero_reduces_to_normal(self):
        x = np.linspace(-3, 3, 7)
        y_normal = normal_distribution(x, loc=0.0, scale=1.0, norm=False)
        y_general = generalized_normal_distribution(
            x, loc=0.0, scale=1.0, shape=0.0, norm=False
        )
        np.testing.assert_allclose(y_general, y_normal, atol=1e-12)

    def test_normalisation(self):
        x = np.linspace(-5, 5, 101)
        y = generalized_normal_distribution(
            x, loc=0.0, scale=2.0, shape=-0.3, norm=True
        )
        assert y.sum() == pytest.approx(1.0, abs=1e-2)


class TestDistanceBetweenGaussian:
    """Distance distribution of two separated 3-D Gaussians."""

    def test_non_negative(self):
        r = np.linspace(0, 20, 101)
        pr = distance_between_gaussian(
            r, separation_distance=10.0, sigma=3.0, normalize=False
        )
        assert np.all(pr >= 0.0)

    def test_peak_near_separation(self):
        r = np.linspace(0, 20, 101)
        pr = distance_between_gaussian(
            r, separation_distance=10.0, sigma=1.0, normalize=True
        )
        peak_idx = np.argmax(pr)
        peak_r = r[peak_idx]
        assert 9.0 < peak_r < 11.0

    def test_normalises(self):
        r = np.linspace(0, 50, 201)
        pr = distance_between_gaussian(
            r, separation_distance=20.0, sigma=2.0, normalize=True
        )
        assert pr.sum() == pytest.approx(1.0, abs=1e-2)


# IMP runs every .py under test/ as a standalone script; a file of bare pytest
# functions would import cleanly and exit 0, reporting success without running
# a single assertion. Hand it to pytest so a failure here fails ctest.
if __name__ == "__main__":
    import sys
    try:
        import pytest
    except ImportError:
        print("pytest not installed; skipping", __file__)
        sys.exit(0)
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))


# --- the Gaussian mixture, in one call -------------------------------------
#
# `Gaussians.distribution` used to loop in Python, calling one of these
# kernels per component and summing the returned arrays in numpy: `k`
# crossings of the language boundary to build one curve. The kernel was
# already C++; the *data* were not staying with it. This is the whole mixture
# in one call, and `GaussianDistances::evaluate` uses the same function, so
# the graph path and the Python path cannot drift apart.

import IMP.bff as _b
from IMP.bff import (
    gaussian_distance_mixture, generalized_normal_distribution,
    normal_distribution, distance_between_gaussian,
)

AXIS = np.linspace(1.0, 200.0, 96)
MEANS = np.array([45.0, 70.0, 110.0])
SIGMAS = np.array([6.0, 9.0, 14.0])
SHAPES = np.array([0.0, 0.3, -0.2])
AMPS = np.array([0.5, 0.3, 0.2])
NONE = np.array([])


def _mix_by_hand(component_of, normalize_components):
    total = np.zeros_like(AXIS)
    for i in range(len(MEANS)):
        c = component_of(i)
        if normalize_components and c.sum() > 0:
            c = c / c.sum()
        total = total + AMPS[i] * c
    return total / total.sum()


def test_the_generalised_normal_branch_is_the_weighted_sum():
    ref = _mix_by_hand(
        lambda i: generalized_normal_distribution(AXIS, MEANS[i], SIGMAS[i],
                                                  SHAPES[i], True), False)
    got = gaussian_distance_mixture(AXIS, MEANS, SIGMAS, SHAPES, AMPS,
                                    _b.GAUSSIAN_MIXTURE_GENERALIZED_NORMAL,
                                    True, True)
    np.testing.assert_allclose(got, ref, rtol=1e-14, atol=0)


def test_the_two_cloud_branch_is_the_weighted_sum():
    ref = _mix_by_hand(
        lambda i: distance_between_gaussian(AXIS, MEANS[i], SIGMAS[i]), False)
    got = gaussian_distance_mixture(AXIS, MEANS, SIGMAS, NONE, AMPS,
                                    _b.GAUSSIAN_MIXTURE_DISTANCE_BETWEEN_GAUSSIANS,
                                    False, True)
    np.testing.assert_allclose(got, ref, rtol=1e-14, atol=0)


def test_the_plain_normal_branch_keeps_the_one_over_sigma():
    """GAUSSIAN_MIXTURE_NORMAL is not the generalised normal at shape 0.

    With `norm=False` the generalised form evaluates the *standard* normal at
    z and so drops the `1/scale` factor. Normalising each component hides
    that; leaving them unnormalised does not, and components of unequal width
    are then weighted wrongly against each other. This pins the difference so
    nobody "simplifies" the two branches into one.
    """
    ref = _mix_by_hand(
        lambda i: normal_distribution(AXIS, MEANS[i], SIGMAS[i], False), False)
    got = gaussian_distance_mixture(AXIS, MEANS, SIGMAS, NONE, AMPS,
                                    _b.GAUSSIAN_MIXTURE_NORMAL, False, True)
    np.testing.assert_allclose(got, ref, rtol=1e-14, atol=0)

    wrong = gaussian_distance_mixture(AXIS, MEANS, SIGMAS, NONE, AMPS,
                                      _b.GAUSSIAN_MIXTURE_GENERALIZED_NORMAL,
                                      False, True)
    assert np.max(np.abs(wrong - ref)) > 1e-3, "the two branches must differ"


def test_normalising_components_changes_the_answer_when_widths_differ():
    """The flag is not cosmetic, which is why the caller has to state it."""
    a = gaussian_distance_mixture(AXIS, MEANS, SIGMAS, SHAPES, AMPS,
                                  _b.GAUSSIAN_MIXTURE_GENERALIZED_NORMAL,
                                  True, True)
    b = gaussian_distance_mixture(AXIS, MEANS, SIGMAS, SHAPES, AMPS,
                                  _b.GAUSSIAN_MIXTURE_GENERALIZED_NORMAL,
                                  False, True)
    assert np.max(np.abs(a - b)) > 1e-3


def test_a_short_shape_array_reads_as_zero_skew():
    got = gaussian_distance_mixture(AXIS, MEANS, SIGMAS, NONE, AMPS,
                                    _b.GAUSSIAN_MIXTURE_GENERALIZED_NORMAL,
                                    True, True)
    ref = gaussian_distance_mixture(AXIS, MEANS, SIGMAS, np.zeros(3), AMPS,
                                    _b.GAUSSIAN_MIXTURE_GENERALIZED_NORMAL,
                                    True, True)
    np.testing.assert_allclose(got, ref, rtol=0, atol=0)


def test_the_node_and_the_free_function_agree():
    """The graph path and the Python path must not be able to drift.

    They share one kernel; this is the test that says so, in both branches.
    """
    for two_cloud in (False, True):
        node = _b.GaussianDistances("distances")
        node.set_number_of_components(len(MEANS))
        node.set_distance_between_gaussians(two_cloud)
        node.add_output_port("distances", _b.GraphPort([0.0], False, True))
        node.set_axis_array(np.ascontiguousarray(AXIS))
        for i in range(len(MEANS)):
            node.get_input_port("mean%d" % i).set_value(MEANS[i])
            node.get_input_port("sigma%d" % i).set_value(SIGMAS[i])
            node.get_input_port("shape%d" % i).set_value(
                0.0 if two_cloud else SHAPES[i])
            node.get_input_port("amplitude%d" % i).set_value(AMPS[i])
        node.evaluate()
        interleaved = np.asarray(node.get_distribution(), dtype=float)
        density = interleaved[0::2]
        free = gaussian_distance_mixture(
            AXIS, MEANS, SIGMAS,
            NONE if two_cloud else SHAPES, AMPS,
            _b.GAUSSIAN_MIXTURE_DISTANCE_BETWEEN_GAUSSIANS if two_cloud
            else _b.GAUSSIAN_MIXTURE_GENERALIZED_NORMAL,
            not two_cloud, True)
        np.testing.assert_allclose(density, free, rtol=0, atol=0)
