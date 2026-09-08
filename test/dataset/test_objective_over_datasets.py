"""PRD-140: an objective that asks the data rather than being told.

`ChiSquared` has always been *told* its noise model by whoever built it. That
is workable for one curve and breaks the moment a joint objective holds
members that disagree -- a Poisson decay and a Gaussian correlation curve
weight differently, and only each member knows how.

So a ChiSquared may instead be given a `Dataset`, which carries its family
with it. The old path is untouched: chisurf sets the noise model per fit and
plots the residuals, so nothing there moves.
"""

import numpy as np
import pytest

import IMP.bff


def _chi(name, dataset):
    """The arithmetic is driven directly, as test/chi2 does: compute_chi2 and
    compute_weighted_residuals take the model and need no wiring."""
    c = IMP.bff.ChiSquared(name)
    c.set_dataset(dataset)
    return c


def _poisson_dataset(y):
    d = IMP.bff.Dataset()
    d.set_values(list(np.asarray(y, dtype=float)))
    d.set_noise_family(IMP.bff.NOISE_FAMILY_POISSON)
    return d


def _gaussian_dataset(y, variance):
    d = IMP.bff.Dataset()
    d.set_values(list(np.asarray(y, dtype=float)))
    d.set_noise_family(IMP.bff.NOISE_FAMILY_STORED)
    d.set_stored_variance(list(np.asarray(variance, dtype=float)))
    return d


def test_the_objective_takes_its_weighting_from_the_dataset():
    y = np.array([5.0, 10.0, 20.0])
    mu = np.array([6.0, 9.0, 25.0])
    c = _chi("counts", _poisson_dataset(y))
    assert c.compute_chi2(list(mu)) == pytest.approx(
        _poisson_dataset(y).objective(list(mu)))


def test_the_same_numbers_under_a_stored_variance():
    y = np.array([100.0, 110.0, 90.0])
    var = np.array([25.0, 25.0, 25.0])
    mu = np.array([102.0, 108.0, 95.0])
    c = _chi("measured", _gaussian_dataset(y, var))
    np.testing.assert_allclose(c.compute_chi2(list(mu)),
                               (((y - mu) ** 2) / var).sum(), rtol=1e-12)


def test_two_objectives_that_disagree_about_noise_each_weight_their_own_way():
    """The property a joint fit needs: one Poisson member and one Gaussian
    member, each weighted its own way, summing to the total.

    It drives the arithmetic directly rather than through `JointChiSquared`,
    so it does not yet prove that the *node* composes them -- wiring that is
    the next step, and this test would be overclaiming if it were named for
    it."""
    y_counts = np.array([5.0, 10.0, 20.0])
    mu_counts = np.array([6.0, 9.0, 25.0])
    y_corr = np.array([1.20, 1.10, 1.05])
    var_corr = np.array([0.01, 0.01, 0.01])
    mu_corr = np.array([1.18, 1.12, 1.06])

    decay = _chi("decay", _poisson_dataset(y_counts))
    corr = _chi("corr", _gaussian_dataset(y_corr, var_corr))

    a = decay.compute_chi2(list(mu_counts))
    b = corr.compute_chi2(list(mu_corr))
    expected = (_poisson_dataset(y_counts).objective(list(mu_counts))
                + _gaussian_dataset(y_corr, var_corr).objective(list(mu_corr)))
    assert a + b == pytest.approx(expected, rel=1e-12)

    # each is weighted its own way: scoring the counts as if they had the
    # correlation curve's stored variance gives a different number entirely
    wrong = _chi("wrong", _gaussian_dataset(y_counts, [25.0, 25.0, 25.0]))
    assert wrong.compute_chi2(list(mu_counts)) != pytest.approx(a, rel=1e-6)


def test_the_residual_kind_defaults_to_the_family_and_can_be_asked_for():
    y = np.array([5.0, 10.0, 20.0])
    mu = np.array([6.0, 9.0, 25.0])
    d = _poisson_dataset(y)
    c = _chi("counts", d)
    default = np.asarray(c.compute_weighted_residuals(list(mu)))
    np.testing.assert_allclose(
        default, np.asarray(d.residuals(list(mu), IMP.bff.RESIDUAL_DEVIANCE)))

    c2 = _chi("counts2", d)
    c2.set_residual_kind(IMP.bff.RESIDUAL_PEARSON)
    np.testing.assert_allclose(
        np.asarray(c2.compute_weighted_residuals(list(mu))),
        np.asarray(d.residuals(list(mu), IMP.bff.RESIDUAL_PEARSON)))


def test_a_dataset_and_an_index_window_together_are_refused():
    """Two masking mechanisms that disagree exclude the wrong points quietly."""
    y = np.array([5.0, 10.0, 20.0])
    c = _chi("counts", _poisson_dataset(y))
    c.set_fit_range(1, 3)
    with pytest.raises(ValueError):
        c.compute_weighted_residuals([6.0, 9.0, 25.0])


def test_the_told_path_is_untouched():
    """Existing callers must see exactly what they saw."""
    y = np.array([5.0, 10.0, 20.0])
    mu = np.array([6.0, 9.0, 25.0])
    c = IMP.bff.ChiSquared("old")
    c.set_data(list(y), [])
    c.set_noise_model_name("poisson")
    assert not c.get_has_dataset()
    expected = np.asarray(IMP.bff.weighted_residuals(
        list(y), [], list(mu), 0, -1, "poisson"))
    np.testing.assert_allclose(np.asarray(c.compute_weighted_residuals(list(mu))),
                               expected, rtol=1e-12, atol=1e-12)
