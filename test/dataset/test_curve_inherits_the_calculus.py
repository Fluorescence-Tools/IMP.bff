"""A curve class that keeps its own arrays and inherits the probability calculus.

The shape chisurf's `DataCurve` would take: its own storage, its own metadata,
its own object model -- plus `IMP.bff.Dataset` as a base, so that scoring,
weighting and uncertainty propagation come from one implementation rather than
being written again beside every curve class.

What makes it work is that the calculus needs nothing of the curve's storage.
The curve syncs its arrays in when they change, and everything downstream --
`variance`, `residuals`, `objective`, `ChiSquared.set_dataset` -- sees a
Dataset.
"""

import numpy as np
import pytest

import IMP.bff


class Provenance:
    """Stands in for chisurf's ExperimentalData: metadata, and nothing numeric."""

    def __init__(self, filename=""):
        self.filename = filename
        self.meta_data = {}


class DataCurveLike(Provenance, IMP.bff.Dataset):
    """A curve that stores numpy arrays and inherits the calculus."""

    def __init__(self, x, y, ey=None, filename=""):
        Provenance.__init__(self, filename)
        IMP.bff.Dataset.__init__(self)
        self.x = np.asarray(x, dtype=float)
        self.y = np.asarray(y, dtype=float)
        self.ey = None if ey is None else np.asarray(ey, dtype=float)
        self.sync()

    def sync(self):
        IMP.bff.sync_dataset(self, self.y, self.ey, self.x)


def test_it_is_a_dataset_and_still_itself():
    c = DataCurveLike([0.0, 1.0, 2.0], [5.0, 10.0, 20.0], filename="a.txt")
    assert isinstance(c, IMP.bff.Dataset)
    assert c.filename == "a.txt" and c.meta_data == {}
    assert c.get_size() == 3 and c.get_rank() == 1
    np.testing.assert_allclose(np.asarray(c.get_coordinate(0)), [0.0, 1.0, 2.0])


def test_counts_without_errors_are_poisson_and_score_as_such():
    y = np.array([5.0, 10.0, 20.0])
    mu = [6.0, 9.0, 25.0]
    c = DataCurveLike([0.0, 1.0, 2.0], y)
    assert c.get_noise_family() == IMP.bff.NOISE_FAMILY_POISSON
    np.testing.assert_allclose(np.asarray(c.variance(mu)), mu)
    plain = IMP.bff.Dataset()
    plain.set_values(list(y))
    plain.set_noise_family(IMP.bff.NOISE_FAMILY_POISSON)
    assert c.objective(mu) == pytest.approx(plain.objective(mu))


def test_given_errors_it_stores_their_squares():
    y = np.array([100.0, 110.0, 90.0])
    ey = np.array([5.0, 5.0, 5.0])
    c = DataCurveLike([0.0, 1.0, 2.0], y, ey)
    assert c.get_noise_family() == IMP.bff.NOISE_FAMILY_STORED
    np.testing.assert_allclose(np.asarray(c.variance(list(y))), ey ** 2)


def test_the_objective_accepts_it_directly():
    c = DataCurveLike([0.0, 1.0, 2.0], [5.0, 10.0, 20.0])
    chi = IMP.bff.ChiSquared("chi")
    chi.set_dataset(c)
    assert chi.get_has_dataset()
    assert chi.compute_chi2([6.0, 9.0, 25.0]) == pytest.approx(
        c.objective([6.0, 9.0, 25.0]))


def test_a_write_that_is_not_synced_is_the_one_hazard():
    """Two stores, so they can disagree. The curve owns the sync and the test
    says so, because a silent disagreement here would be a wrong objective
    with no symptom."""
    c = DataCurveLike([0.0, 1.0, 2.0], [5.0, 10.0, 20.0])
    before = c.objective([5.0, 10.0, 20.0])
    assert before == pytest.approx(0.0, abs=1e-12)
    c.y = np.array([50.0, 100.0, 200.0])          # written, not synced
    assert c.objective([5.0, 10.0, 20.0]) == pytest.approx(0.0, abs=1e-12)
    c.sync()
    assert c.objective([5.0, 10.0, 20.0]) > 1.0


def test_the_curve_can_be_combined_with_another_and_carry_the_uncertainty():
    """The magic angle, built from two curve objects rather than two bare
    datasets: the calculus does not care which it was given."""
    vv = DataCurveLike([0.0, 1.0], [900.0, 400.0])
    vh = DataCurveLike([0.0, 1.0], [300.0, 150.0])
    vv.set_as_source("VV")
    vh.set_as_source("VH")
    g = 1.15
    magic = IMP.bff.Dataset.add(vv, IMP.bff.Dataset.affine(vh, 2.0 * g))
    np.testing.assert_allclose(np.asarray(magic.get_values()),
                               np.array([900.0, 400.0]) + 2 * g * np.array([300.0, 150.0]))
    np.testing.assert_allclose(
        np.asarray(magic.variance(list(magic.get_values()))),
        np.array([900.0, 400.0]) + 4 * g * g * np.array([300.0, 150.0]), rtol=1e-12)
