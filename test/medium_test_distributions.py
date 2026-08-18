"""Tests for ``IMP.bff.representation.probability`` — mathematical distributions."""

from __future__ import annotations

import numpy as np
import pytest

from IMP.bff.representation.probability import (
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
