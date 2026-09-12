"""Uncertainty carried through arithmetic on datasets, correctly under sharing.

Two channels of counts get combined in fluorescence all the time, and the
combination is not Poisson. The constructed magic-angle decay is
`VV + 2 G VH`, whose variance is `VV + 4 G^2 VH` -- the coefficient squares --
while reading the result as counts would say `VV + 2 G VH`. Anisotropy is
worse: `(VV - G VH)/(VV + 2 G VH)` has the same channels above and below the
line, so numerator and denominator are correlated and treating them as
independent is wrong however carefully each half is done.

So propagation here is over the **independent sources**, with the derivatives
composed as the expression is built. The caller writes the arithmetic; the
correlation is handled because both branches trace back to the same VV and VH.
"""

import numpy as np
import pytest

import IMP.bff

G = 1.15


def _counts(values, name):
    d = IMP.bff.FitDataset()
    d.set_values(list(np.asarray(values, dtype=float)))
    d.set_noise_family(IMP.bff.FIT_NOISE_FAMILY_POISSON)
    d.set_as_source(name)
    return d


@pytest.fixture
def channels():
    vv = np.array([900.0, 400.0, 100.0])
    vh = np.array([300.0, 150.0, 40.0])
    return vv, vh, _counts(vv, "VV"), _counts(vh, "VH")


def test_a_source_carries_its_own_variance(channels):
    vv, _, dvv, _ = channels
    np.testing.assert_allclose(np.asarray(dvv.variance(list(vv))), vv)
    assert list(dvv.get_source_names()) == ["VV"]


def test_the_constructed_magic_angle_is_not_poisson(channels):
    """The coefficient squares in the variance and a Poisson reading of the
    result does not know that."""
    vv, vh, dvv, dvh = channels
    magic = IMP.bff.FitDataset.add(dvv, IMP.bff.FitDataset.affine(dvh, 2.0 * G))
    np.testing.assert_allclose(np.asarray(magic.get_values()), vv + 2 * G * vh)

    got = np.asarray(magic.variance(list(magic.get_values())))
    np.testing.assert_allclose(got, vv + 4 * G * G * vh, rtol=1e-12)

    # what a Poisson reading of the same numbers would have claimed
    wrong = vv + 2 * G * vh
    assert np.all(got > wrong), "the perpendicular channel is understated by 2G"
    np.testing.assert_allclose(got[0] / wrong[0], (900 + 4 * G * G * 300) / (900 + 2 * G * 300),
                               rtol=1e-12)


def test_a_ratio_of_the_same_channels_is_correlated(channels):
    """Anisotropy. Propagating over the operands as if they were independent
    is the standard mistake; tracing to VV and VH is right."""
    vv, vh, dvv, dvh = channels
    num = IMP.bff.FitDataset.subtract(dvv, IMP.bff.FitDataset.affine(dvh, G))
    den = IMP.bff.FitDataset.add(dvv, IMP.bff.FitDataset.affine(dvh, 2.0 * G))
    r = IMP.bff.FitDataset.divide(num, den)

    N, D = vv - G * vh, vv + 2 * G * vh
    np.testing.assert_allclose(np.asarray(r.get_values()), N / D, rtol=1e-12)

    # correct: d r/d VV and d r/d VH, then sum of squares times each variance
    dr_dvv = (D - N) / D ** 2
    dr_dvh = (-G * D - N * 2 * G) / D ** 2
    correct = dr_dvv ** 2 * vv + dr_dvh ** 2 * vh
    np.testing.assert_allclose(np.asarray(r.variance(list(r.get_values()))),
                               correct, rtol=1e-12)

    # the mistake, for contrast: numerator and denominator taken as independent
    var_N = vv + G * G * vh
    var_D = vv + 4 * G * G * vh
    naive = var_N / D ** 2 + N ** 2 * var_D / D ** 4
    assert not np.allclose(naive, correct), \
        "if these agreed there would be nothing to get wrong"


def test_both_sources_are_kept_once_not_twice(channels):
    _, _, dvv, dvh = channels
    r = IMP.bff.FitDataset.divide(IMP.bff.FitDataset.subtract(dvv, dvh),
                               IMP.bff.FitDataset.add(dvv, dvh))
    assert sorted(r.get_source_names()) == ["VH", "VV"], \
        "VV appears on both sides and must stay one source"


def test_a_difference_of_one_channel_with_itself_has_no_variance(channels):
    """The sharpest check that correlation is handled: x - x is exactly zero
    and exactly certain. Independent propagation would give it 2*Var(x)."""
    _, _, dvv, _ = channels
    zero = IMP.bff.FitDataset.subtract(dvv, dvv)
    np.testing.assert_allclose(np.asarray(zero.get_values()), 0.0, atol=1e-12)
    np.testing.assert_allclose(
        np.asarray(zero.variance(list(zero.get_values()))), 0.0, atol=1e-12)


def test_products_and_a_general_transform(channels):
    vv, vh, dvv, dvh = channels
    prod = IMP.bff.FitDataset.multiply(dvv, dvh)
    np.testing.assert_allclose(np.asarray(prod.get_values()), vv * vh)
    np.testing.assert_allclose(
        np.asarray(prod.variance(list(prod.get_values()))),
        vh ** 2 * vv + vv ** 2 * vh, rtol=1e-12)

    # sqrt, through the general escape: values and the derivative per point
    root = IMP.bff.FitDataset.transform(dvv, list(np.sqrt(vv)),
                                     list(0.5 / np.sqrt(vv)))
    np.testing.assert_allclose(
        np.asarray(root.variance(list(root.get_values()))),
        (0.5 / np.sqrt(vv)) ** 2 * vv, rtol=1e-12)


def test_an_untracked_operand_is_refused():
    """Silence is the failure mode: an operand carrying no uncertainty would
    contribute none and the result would look more certain than it is."""
    a = _counts([10.0, 20.0], "A")
    plain = IMP.bff.FitDataset()
    plain.set_values([1.0, 2.0])
    with pytest.raises(ValueError):
        IMP.bff.FitDataset.add(a, plain)


def test_a_propagated_dataset_scores_with_its_propagated_variance(channels):
    """And it reaches the objective, which is the point of carrying it."""
    vv, vh, dvv, dvh = channels
    magic = IMP.bff.FitDataset.add(dvv, IMP.bff.FitDataset.affine(dvh, 2.0 * G))
    y = np.asarray(magic.get_values())
    mu = y * 1.02
    var = vv + 4 * G * G * vh
    np.testing.assert_allclose(magic.objective(list(mu)),
                               (((y - mu) ** 2) / var).sum(), rtol=1e-12)
    assert "propagated" in magic.describe()
