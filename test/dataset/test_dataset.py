"""PRD-140: measured values of any rank, and the noise family that goes with them.

No curve class: a curve is the rank-1 case, an image the rank-2 one, a density
the rank-3 one, and only the shape differs.

The design decision these pin is that a dataset owns the noise **family** and
not the per-point weights. For counts the weight is not a property of the data
at all -- weighting by the data is Neyman and biased, weighting by the model
is Pearson and is what the Fisher information uses -- so the variance is asked
for *given a model*, and it moves as the fit moves. A dataset that owned
sigmas would be right for a Gaussian measurement and quietly wrong for counts,
and the failure is that the fit converges to the wrong place with a reduced
chi-square near one.
"""

import numpy as np
import pytest

import IMP.bff


def _poisson(y, mask=None):
    d = IMP.bff.Dataset()
    d.set_values(list(np.asarray(y, dtype=float)))
    d.set_noise_family(IMP.bff.NOISE_FAMILY_POISSON)
    if mask is not None:
        d.set_mask(list(np.asarray(mask, dtype=float)))
    return d


def test_a_curve_is_the_rank_one_case():
    d = _poisson([1.0, 2.0, 3.0])
    assert d.get_rank() == 1
    assert list(d.get_shape()) == [3]
    assert d.get_size() == 3


def test_the_same_class_carries_an_image_and_a_density():
    for shape in ([4, 5], [3, 4, 5]):
        n = int(np.prod(shape))
        d = IMP.bff.Dataset()
        d.set_values(list(np.arange(n, dtype=float)), list(shape))
        assert d.get_rank() == len(shape)
        assert list(d.get_shape()) == list(shape)
        assert d.get_size() == n


def test_a_shape_that_does_not_describe_the_values_is_refused():
    d = IMP.bff.Dataset()
    with pytest.raises(ValueError):
        d.set_values([1.0, 2.0, 3.0], [2, 2])
    with pytest.raises(ValueError):
        d.set_values([1.0, 2.0], [2, 0])


def test_poisson_variance_is_the_model_not_the_data():
    """The heart of it. Ask twice with different models and get different
    weights from the same data."""
    d = _poisson([10.0, 10.0, 10.0])
    v1 = np.asarray(d.variance([4.0, 9.0, 16.0]))
    v2 = np.asarray(d.variance([1.0, 1.0, 1.0]))
    np.testing.assert_allclose(v1, [4.0, 9.0, 16.0])
    np.testing.assert_allclose(v2, [1.0, 1.0, 1.0])
    assert not np.allclose(v1, v2), "the weights move with the model"


def test_a_stored_variance_is_the_data_s_own():
    d = IMP.bff.Dataset()
    d.set_values([1.0, 2.0, 3.0])
    d.set_noise_family(IMP.bff.NOISE_FAMILY_STORED)
    d.set_stored_variance([4.0, 4.0, 9.0])
    np.testing.assert_allclose(np.asarray(d.variance([0.0, 0.0, 0.0])),
                               [4.0, 4.0, 9.0])


def test_a_stored_family_with_nothing_stored_refuses():
    """It must not fall back to ones. Unweighted least squares on counts over
    four decades fits the peak and ignores the tail, and says nothing."""
    d = IMP.bff.Dataset()
    d.set_values([1.0, 2.0, 3.0])
    d.set_noise_family(IMP.bff.NOISE_FAMILY_STORED)
    with pytest.raises(ValueError):
        d.variance([1.0, 1.0, 1.0])
    with pytest.raises(ValueError):
        d.objective([1.0, 1.0, 1.0])


def test_one_poisson_dataset_gives_both_residuals():
    """A real analysis reports a deviance and runs its runs test on Pearson
    signs. Choosing one at construction is a choice it does not make."""
    y = np.array([5.0, 10.0, 20.0])
    mu = np.array([6.0, 9.0, 25.0])
    d = _poisson(y)
    pearson = np.asarray(d.residuals(list(mu), IMP.bff.RESIDUAL_PEARSON))
    deviance = np.asarray(d.residuals(list(mu), IMP.bff.RESIDUAL_DEVIANCE))
    np.testing.assert_allclose(pearson, (y - mu) / np.sqrt(mu))
    assert not np.allclose(pearson, deviance)
    # they agree in sign, which is what a runs test uses
    np.testing.assert_array_equal(np.sign(pearson), np.sign(deviance))
    # and the deviance residuals square to the objective
    np.testing.assert_allclose((deviance ** 2).sum(), d.objective(list(mu)),
                               rtol=1e-12)


def test_neyman_is_available_and_named():
    """Offered because a great deal of existing analysis used it, and named so
    that nobody reaches for it without meaning to."""
    y = np.array([4.0, 9.0])
    mu = np.array([5.0, 8.0])
    d = _poisson(y)
    np.testing.assert_allclose(
        np.asarray(d.residuals(list(mu), IMP.bff.RESIDUAL_NEYMAN)),
        (y - mu) / np.sqrt(y))


def test_neyman_and_pearson_disagree_and_that_is_the_point():
    """The bias, shown rather than asserted: over a decade of counts, weighting
    by the data and by the model put the misfit in different places."""
    y = np.array([4.0, 400.0])
    mu = np.array([9.0, 405.0])
    d = _poisson(y)
    neyman = np.asarray(d.residuals(list(mu), IMP.bff.RESIDUAL_NEYMAN))
    pearson = np.asarray(d.residuals(list(mu), IMP.bff.RESIDUAL_PEARSON))
    assert abs(neyman[0]) > abs(pearson[0]), "the low bin is weighted harder by 1/y"


def test_a_mask_excludes_points_from_both_residuals_and_objective():
    d = _poisson([5.0, 10.0, 20.0], mask=[1.0, 0.0, 1.0])
    mu = [6.0, 1000.0, 21.0]
    r = np.asarray(d.residuals(mu, IMP.bff.RESIDUAL_PEARSON))
    assert r[1] == 0.0
    assert d.get_number_of_active_points() == 2
    unmasked = _poisson([5.0, 20.0]).objective([6.0, 21.0])
    np.testing.assert_allclose(d.objective(mu), unmasked, rtol=1e-12)


def test_the_objective_is_zero_for_a_perfect_model():
    d = _poisson([5.0, 10.0, 20.0])
    assert d.objective([5.0, 10.0, 20.0]) == pytest.approx(0.0, abs=1e-12)


def test_a_model_of_the_wrong_length_is_refused():
    d = _poisson([1.0, 2.0])
    with pytest.raises(ValueError):
        d.objective([1.0, 2.0, 3.0])


def test_a_three_dimensional_dataset_scores_like_a_flat_one():
    """The N-D case is not an afterthought: shape changes nothing but shape."""
    rng = np.random.default_rng(0)
    y = rng.integers(1, 50, 3 * 4 * 5).astype(float)
    mu = y + rng.normal(0.0, 1.0, y.size)
    mu = np.abs(mu) + 0.5
    nd = IMP.bff.Dataset()
    nd.set_values(list(y), [3, 4, 5])
    nd.set_noise_family(IMP.bff.NOISE_FAMILY_POISSON)
    flat = _poisson(y)
    assert nd.objective(list(mu)) == pytest.approx(flat.objective(list(mu)))
    assert "shape=[3, 4, 5]" in nd.describe()


def test_the_poisson_objective_is_the_one_the_library_already_computes():
    """The migration claim, checked rather than asserted: a Poisson Dataset
    scores a model exactly as `weighted_residuals(..., "poisson")` does, so
    moving an objective onto the dataset does not move any number."""
    rng = np.random.default_rng(3)
    y = rng.integers(0, 500, 200).astype(float)
    mu = np.abs(y + rng.normal(0.0, 8.0, y.size)) + 0.5
    existing = np.asarray(IMP.bff.weighted_residuals(
        list(y), [], list(mu), 0, -1, "poisson"))
    d = _poisson(y)
    from_dataset = np.asarray(d.residuals(list(mu), IMP.bff.RESIDUAL_DEVIANCE))

    # Same magnitudes, so the objective is identical and no fit moves.
    np.testing.assert_allclose(np.abs(from_dataset), np.abs(existing),
                               rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(d.objective(list(mu)), (existing ** 2).sum(),
                               rtol=1e-12)

    # Opposite signs, and the reason is a defect in the existing class rather
    # than a choice here: `weighted_residuals` returns (y - mu)/e under
    # "default" -- positive where the data exceeds the model -- and the
    # opposite under "poisson". So switching noise model there flips every
    # residual, and a residual plot or a runs test flips with it. Dataset is
    # consistent: Pearson, deviance and Neyman all read positive where the
    # data exceeds the model. Not fixed in ChiSquared, because chisurf plots
    # those signs today; recorded in PRD-140.
    np.testing.assert_allclose(from_dataset, -existing, rtol=1e-12, atol=1e-12)


def test_the_stored_objective_is_the_one_the_library_already_computes():
    """And the other path: a stored variance scores as the default noise model
    does when the caller supplies the same errors."""
    rng = np.random.default_rng(4)
    y = rng.normal(100.0, 10.0, 150)
    ey = np.full(y.size, 10.0)
    mu = y + rng.normal(0.0, 3.0, y.size)
    existing = np.asarray(IMP.bff.weighted_residuals(
        list(y), list(ey), list(mu), 0, -1, "default"))
    d = IMP.bff.Dataset()
    d.set_values(list(y))
    d.set_noise_family(IMP.bff.NOISE_FAMILY_STORED)
    d.set_stored_variance(list(ey ** 2))
    from_dataset = np.asarray(d.residuals(list(mu), IMP.bff.RESIDUAL_PEARSON))
    np.testing.assert_allclose(from_dataset, existing, rtol=1e-12, atol=1e-12)
