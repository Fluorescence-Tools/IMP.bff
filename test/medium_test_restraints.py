"""Tests for ``IMP.bff``'s labelling restraints — AV network and direct.

Both restraints, their measurements and their sites are C++ values
(``include/IMP/bff/LabelingRestraints.h``); this exercises them from the
Python surface SWIG gives them.
"""

from __future__ import annotations

import numpy as np
import pytest

from IMP.bff import AccessibleVolume
from IMP.bff.restraints import (
    SimpleAVNetworkRestraint,
    AVMeasurement,
    DirectLabelingRestraint,
    LabelingSite,
)


class TestAVMeasurement:
    """AVMeasurement, a value."""

    def test_defaults(self):
        m = AVMeasurement(
            av1_name="donor", av2_name="acceptor",
            distance=50.0, error_neg=3.0, error_pos=5.0,
        )
        assert m.av1_name == "donor"
        assert m.distance == 50.0
        assert m.forster_radius == 52.0
        assert m.distance_type == "RDAMean"


class TestSimpleAVNetworkRestraint:
    """Programmatic AV network restraint."""

    @pytest.fixture
    def two_avs(self):
        rng = np.random.RandomState(0)
        pts1 = rng.randn(300, 4).astype(np.float64)
        pts1[:, 3] = np.abs(pts1[:, 3]) + 0.01
        pts2 = rng.randn(300, 4).astype(np.float64)
        pts2[:, 3] = np.abs(pts2[:, 3]) + 0.01
        pts2[:, 0] += 50.0  # shift → ~50 Å apart
        av1 = AccessibleVolume(points=pts1.ravel(), position_name="donor")
        av2 = AccessibleVolume(points=pts2.ravel(), position_name="acceptor")
        return av1, av2

    def _rmp(self, r, distance=50.0, error=3.0):
        r.add_measurement(AVMeasurement(
            av1_name="d", av2_name="a", distance=distance,
            error_neg=error, error_pos=error, distance_type="Rmp"))

    def test_create_empty(self):
        r = SimpleAVNetworkRestraint()
        assert r.evaluate() == 0.0
        assert r.get_n_avs() == 0

    def test_add_av(self, two_avs):
        av1, av2 = two_avs
        r = SimpleAVNetworkRestraint()
        r.add_av("donor", av1)
        r.add_av("acceptor", av2)
        assert r.has_av("donor") and r.has_av("acceptor")
        assert r.get_n_avs() == 2

    def test_evaluate_no_measurements(self, two_avs):
        av1, av2 = two_avs
        r = SimpleAVNetworkRestraint()
        r.add_av("d", av1)
        r.add_av("a", av2)
        assert r.evaluate() == 0.0

    def test_evaluate_one_measurement(self, two_avs):
        av1, av2 = two_avs
        r = SimpleAVNetworkRestraint()
        r.add_av("d", av1)
        r.add_av("a", av2)
        self._rmp(r)
        score = r.evaluate()
        assert score > 0.0
        assert isinstance(score, float)

    def test_weight(self, two_avs):
        av1, av2 = two_avs
        r = SimpleAVNetworkRestraint()
        r.add_av("d", av1)
        r.add_av("a", av2)
        self._rmp(r)
        score_unweighted = r.evaluate()
        r.set_weight(2.0)
        assert r.get_weight() == 2.0
        assert r.evaluate() == pytest.approx(2.0 * score_unweighted, rel=1e-12)

    def test_get_model_distances(self, two_avs):
        av1, av2 = two_avs
        r = SimpleAVNetworkRestraint()
        r.add_av("d", av1)
        r.add_av("a", av2)
        self._rmp(r)
        md = r.get_model_distances()
        assert "d_a" in md
        assert md["d_a"] > 0.0

    def test_multiple_measurements(self, two_avs):
        av1, av2 = two_avs
        r = SimpleAVNetworkRestraint()
        r.add_av("d", av1)
        r.add_av("a", av2)
        self._rmp(r, distance=50.0, error=3.0)
        self._rmp(r, distance=55.0, error=5.0)
        assert r.get_n_measurements() == 2
        assert r.evaluate() > 0.0

    def test_an_unregistered_volume_is_an_error(self, two_avs):
        av1, _ = two_avs
        r = SimpleAVNetworkRestraint()
        r.add_av("d", av1)
        self._rmp(r)
        with pytest.raises(ValueError):
            r.evaluate()

    def test_an_unknown_distance_type_is_an_error(self, two_avs):
        av1, av2 = two_avs
        r = SimpleAVNetworkRestraint()
        r.add_av("d", av1)
        r.add_av("a", av2)
        r.add_measurement(AVMeasurement(
            av1_name="d", av2_name="a", distance=50.0, error_neg=3.0,
            error_pos=3.0, distance_type="nonsense"))
        with pytest.raises(ValueError):
            r.evaluate()


class TestLabelingSite:
    """LabelingSite, a value."""

    def test_defaults(self):
        s = LabelingSite(residue_seq_number=5)
        assert s.residue_seq_number == 5
        assert s.atom_name == "CB"
        assert s.distance == 0.0


class TestDirectLabelingRestraint:
    """Fast direct-labeling restraint."""

    #: A restraint with no coordinates. The array is `(0, 3)` rather than
    #: absent: the kernel takes numpy's shape, and an empty one still has to
    #: say it is three-dimensional.
    NO_COORDS = np.zeros((0, 3))

    def test_create_empty(self):
        r = DirectLabelingRestraint(self.NO_COORDS)
        assert r.evaluate() == 0.0

    def test_no_xyz(self):
        r = DirectLabelingRestraint(self.NO_COORDS)
        r.add_site(LabelingSite(residue_seq_number=0, distance=50.0))
        r.add_site(LabelingSite(residue_seq_number=1, distance=50.0))
        assert r.evaluate() == 0.0

    def test_two_sites(self):
        xyz = np.array([[0., 0., 0.], [30., 0., 0.]], dtype=np.float64)
        r = DirectLabelingRestraint(xyz)
        r.add_site(LabelingSite(residue_seq_number=0, distance=30.0,
                                error_neg=3.0, error_pos=3.0))
        r.add_site(LabelingSite(residue_seq_number=1, distance=30.0,
                                error_neg=3.0, error_pos=3.0))
        # perfect fit → chi2 = 0
        assert r.evaluate() == pytest.approx(0.0, abs=1e-10)

    def test_two_sites_mismatch(self):
        xyz = np.array([[0., 0., 0.], [30., 0., 0.]], dtype=np.float64)
        r = DirectLabelingRestraint(xyz)
        r.add_site(LabelingSite(residue_seq_number=0, distance=35.0,
                                error_neg=3.0, error_pos=3.0))
        r.add_site(LabelingSite(residue_seq_number=1, distance=35.0,
                                error_neg=3.0, error_pos=3.0))
        # |30 - 35| = 5, / 3 → (5/3)^2 ≈ 2.78
        expected = (5.0 / 3.0) ** 2
        assert r.evaluate() == pytest.approx(expected, rel=1e-9)

    def test_three_sites(self):
        xyz = np.array([[0., 0., 0.],
                        [10., 0., 0.],
                        [20., 0., 0.]], dtype=np.float64)
        r = DirectLabelingRestraint(xyz)
        r.add_site(LabelingSite(residue_seq_number=0, distance=10.0,
                                error_neg=1.0, error_pos=1.0))
        r.add_site(LabelingSite(residue_seq_number=1, distance=10.0,
                                error_neg=1.0, error_pos=1.0))
        r.add_site(LabelingSite(residue_seq_number=2, distance=20.0,
                                error_neg=1.0, error_pos=1.0))
        # Sites 0-1: |10-10|=0 → 0
        # Sites 0-2: |20-20|=0 → 0
        # Sites 1-2: |10-20|=10 → (10/1)^2 = 100
        assert r.evaluate() == pytest.approx(100.0, abs=1e-9)

    def test_weight(self):
        xyz = np.array([[0., 0., 0.], [30., 0., 0.]], dtype=np.float64)
        r = DirectLabelingRestraint(xyz, weight=2.0)
        r.add_site(LabelingSite(residue_seq_number=0, distance=25.0,
                                error_neg=5.0, error_pos=5.0))
        r.add_site(LabelingSite(residue_seq_number=1, distance=25.0,
                                error_neg=5.0, error_pos=5.0))
        expected = 2.0 * (5.0 / 5.0) ** 2
        assert r.evaluate() == pytest.approx(expected, rel=1e-9)

    def test_a_site_outside_the_coordinates_is_an_error(self):
        xyz = np.array([[0., 0., 0.], [30., 0., 0.]], dtype=np.float64)
        r = DirectLabelingRestraint(xyz)
        r.add_site(LabelingSite(residue_seq_number=0, distance=30.0))
        r.add_site(LabelingSite(residue_seq_number=7, distance=30.0))
        with pytest.raises(ValueError):
            r.evaluate()


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
