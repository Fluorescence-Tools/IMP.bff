"""The distance conventions that used to share names across two modules.

``IMP.bff.distance_metrics`` held six functions with the same names as
``IMP.bff``, three of which answered genuinely
different questions and three of which were verified identical. The module was
retired in the PRD-113 cleanup (2026-08-18) and moved to ``junk/`` for
reference; its three distinct functions were folded into the canonical module
under names that say which question they answer:

===================================  ===============================================
was ``distance_metrics.``            is now ``IMP.bff.``
===================================  ===============================================
``polynomial_transfer`` (ascending)  ``polynomial_transfer_ascending``
``av_pair_statistics`` (the clouds)  ``distance_sample_statistics`` (a sample)
``mean_position_distance`` (clouds)  ``mean_position_distance`` (unchanged)
===================================  ===============================================

These tests followed the functions. They exist to make sure the retirement lost
nothing: every value asserted here is the value the retired module produced.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from IMP.bff import (
    av_pair_statistics as _av_pair_statistics_from_volumes,  # noqa: F401
    chi2_score,
    distance_from_fret_efficiency,
    distance_sample_statistics,
    fret_efficiency,
    gaussian_rmp_to_rda_mean,
    mean_position_distance,
    polynomial_transfer_ascending as polynomial_transfer,
)


class TestChi2Score:
    """Asymmetric chi-squared contribution."""

    def test_model_below_exp(self):
        """Negative deviation uses error_neg."""
        s = chi2_score(50.0, 55.0, 3.0, 5.0)
        expected = (5.0 / 3.0) ** 2
        assert s == pytest.approx(expected, rel=1e-9)

    def test_model_above_exp(self):
        """Positive deviation uses error_pos."""
        s = chi2_score(60.0, 55.0, 3.0, 5.0)
        expected = (5.0 / 5.0) ** 2
        assert s == pytest.approx(expected, rel=1e-9)

    def test_perfect_fit(self):
        s = chi2_score(55.0, 55.0, 3.0, 5.0)
        assert s == pytest.approx(0.0, abs=1e-12)


class TestFretEfficiency:
    """FRET efficiency ↔ distance conversions."""

    def test_efficiency_at_foerster_radius(self):
        e = fret_efficiency(52.0, forster_radius=52.0)
        assert e == pytest.approx(0.5, abs=1e-12)

    def test_efficiency_zero_distance(self):
        e = fret_efficiency(0.0, forster_radius=52.0)
        assert e == pytest.approx(1.0, abs=1e-12)

    def test_round_trip(self):
        r0 = 52.0
        for r in [10.0, 30.0, 52.0, 80.0]:
            e = fret_efficiency(r, r0)
            r_back = distance_from_fret_efficiency(e, r0)
            assert r_back == pytest.approx(r, rel=1e-6)


class TestGaussianRmpToRdaMean:
    """Sigma-corrected mean distance."""

    def test_correction_positive(self):
        r = gaussian_rmp_to_rda_mean(50.0, sigma=6.0)
        assert r > 50.0
        expected = 50.0 + 6.0 ** 2 / 50.0
        assert r == pytest.approx(expected, rel=1e-9)

    def test_zero_rmp(self):
        r = gaussian_rmp_to_rda_mean(0.0, sigma=6.0)
        assert r == pytest.approx(0.0, abs=1e-12)


class TestPolynomialTransfer:
    """Polynomial Rmp → RDAMean conversion."""

    def test_identity(self):
        coeffs = np.array([0.0, 1.0], dtype=np.float64)
        r = polynomial_transfer(50.0, coeffs)
        assert r == pytest.approx(50.0, rel=1e-12)

    def test_quadratic(self):
        coeffs = np.array([0.0, 0.0, 1.0], dtype=np.float64)
        r = polynomial_transfer(3.0, coeffs)
        assert r == pytest.approx(9.0, rel=1e-12)


class TestAvPairStatistics:
    """AV distance metric batch computation."""

    def test_uniform_weights(self):
        dists = np.array([40.0, 50.0, 60.0], dtype=np.float64)
        weights = np.ones(3, dtype=np.float64)
        r_da, r_e, mean_e, sigma = distance_sample_statistics(dists, weights)
        assert r_da == pytest.approx(50.0, abs=1e-9)
        assert 0.0 < mean_e < 1.0
        assert r_e > 0.0 and sigma > 0.0

    def test_r_mp_is_not_in_the_tuple(self):
        """A sample of pair distances cannot produce R_mp, so it is not offered.

        The tuple carried a second slot for it, always NaN, because an earlier
        version had returned ``r_da_mean`` under that name -- and measured on
        148l E15/E90 the two are 47.65 A against 51.53 A, 8 % apart. R_mp is
        the distance between the clouds' *mean positions*, and
        ``|<a> - <b>|`` is not a function of the distribution of ``|a - b|``.
        The slot is gone; :func:`mean_position_distance` takes the clouds.
        """
        dists = np.array([40.0, 50.0, 60.0], dtype=np.float64)
        weights = np.ones(3, dtype=np.float64)
        assert len(distance_sample_statistics(dists, weights)) == 4

    def test_mean_position_distance_takes_the_clouds(self):
        a = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        b = np.array([[10.0, 0.0, 0.0], [12.0, 0.0, 0.0]])
        assert mean_position_distance(a.ravel(), b.ravel()) == pytest.approx(10.0)

    def test_mean_position_distance_is_not_the_mean_pair_distance(self):
        """The two differ whenever the clouds are extended, which is always."""
        rng = np.random.default_rng(0)
        a = rng.normal(0.0, 5.0, size=(500, 3))
        b = rng.normal(0.0, 5.0, size=(500, 3)) + np.array([30.0, 0.0, 0.0])

        r_mp = mean_position_distance(a.ravel(), b.ravel())
        pair = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1).ravel()
        r_da_mean = float(pair.mean())

        assert r_mp == pytest.approx(30.0, abs=1.0)
        assert r_da_mean > r_mp + 0.5, (
            "the mean pair distance exceeds the distance between mean "
            "positions for any pair of extended clouds; if these agree the "
            "test geometry is degenerate"
        )

    def test_zero_weights(self):
        dists = np.array([40.0, 50.0], dtype=np.float64)
        weights = np.zeros(2, dtype=np.float64)
        r_da, r_e, mean_e, sigma = distance_sample_statistics(dists, weights)
        assert r_da == pytest.approx(0.0, abs=1e-12)


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
