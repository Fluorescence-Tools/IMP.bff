"""Tests for ``IMP.bff.restraints`` — AV and direct-labeling restraints."""

from __future__ import annotations

import numpy as np
import pytest

from IMP.bff.representation.av import BasicAV
from IMP.bff.restraints import (
    SimpleAVNetworkRestraint,
    AVMeasurement,
    DirectLabelingRestraint,
    LabelingSite,
)


class TestAVMeasurement:
    """AVMeasurement dataclass."""

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
        av1 = BasicAV(points=pts1, position_name="donor")
        av2 = BasicAV(points=pts2, position_name="acceptor")
        return av1, av2

    def test_create_empty(self):
        r = SimpleAVNetworkRestraint()
        assert r.evaluate() == 0.0
        assert len(r._avs) == 0

    def test_add_av_object(self, two_avs):
        av1, av2 = two_avs
        r = SimpleAVNetworkRestraint()
        r.add_av_object("donor", av1)
        r.add_av_object("acceptor", av2)
        assert "donor" in r._avs
        assert "acceptor" in r._avs

    def test_evaluate_no_measurements(self, two_avs):
        av1, av2 = two_avs
        r = SimpleAVNetworkRestraint()
        r.add_av_object("d", av1)
        r.add_av_object("a", av2)
        assert r.evaluate() == 0.0

    def test_evaluate_one_measurement(self, two_avs):
        av1, av2 = two_avs
        r = SimpleAVNetworkRestraint()
        r.add_av_object("d", av1)
        r.add_av_object("a", av2)
        r.add_measurement(
            "d", "a",
            distance=50.0, error_neg=3.0, error_pos=3.0,
            distance_type="Rmp",
        )
        score = r.evaluate()
        assert score > 0.0
        assert isinstance(score, float)

    def test_weight(self, two_avs):
        av1, av2 = two_avs
        r = SimpleAVNetworkRestraint()
        r.add_av_object("d", av1)
        r.add_av_object("a", av2)
        r.add_measurement(
            "d", "a",
            distance=50.0, error_neg=3.0, error_pos=3.0,
            distance_type="Rmp",
        )
        score_unweighted = r.evaluate()
        r.weight = 2.0
        score_weighted = r.evaluate()
        assert score_weighted == pytest.approx(2.0 * score_unweighted, rel=1e-12)

    def test_get_model_distances(self, two_avs):
        av1, av2 = two_avs
        r = SimpleAVNetworkRestraint()
        r.add_av_object("d", av1)
        r.add_av_object("a", av2)
        r.add_measurement(
            "d", "a",
            distance=50.0, error_neg=3.0, error_pos=3.0,
            distance_type="Rmp",
        )
        md = r.get_model_distances()
        assert "d_a" in md
        assert md["d_a"] > 0.0

    def test_multiple_measurements(self, two_avs):
        av1, av2 = two_avs
        r = SimpleAVNetworkRestraint()
        r.add_av_object("d", av1)
        r.add_av_object("a", av2)
        r.add_measurement("d", "a", distance=50.0, error_neg=3.0, error_pos=3.0, distance_type="Rmp")
        r.add_measurement("d", "a", distance=55.0, error_neg=5.0, error_pos=5.0, distance_type="Rmp")
        assert r.evaluate() > 0.0

    def test_add_av_from_coords_raises_without_backend(self, monkeypatch):
        """Without IMP.bff's AV decorator, coordinate-based AV creation fails.

        There is one backend since PRD-112 stage 1, so the absence to simulate
        is a build without the decorator. Two earlier revisions of this test
        held for the wrong reason: first because the LabelLib probe looked for
        `LabelLib.AV`, which current builds do not expose, so a usable LabelLib
        was reported missing; then because both flags had to be cleared.
        """
        xyz = np.array([[0., 0., 0.]])
        vdw = np.array([1.5])
        src = np.array([0., 0., 0.])
        r = SimpleAVNetworkRestraint()
        import IMP.bff.representation.av as _compute

        monkeypatch.setattr(_compute, "_HAS_IMP_BFF", False, raising=False)
        with pytest.raises((ImportError, RuntimeError)):
            r.add_av_from_coords("test", xyz, vdw, src)


class TestLabelingSite:
    """LabelingSite dataclass."""

    def test_defaults(self):
        s = LabelingSite(residue_seq_number=5)
        assert s.residue_seq_number == 5
        assert s.atom_name == "CB"
        assert s.distance == 0.0


class TestDirectLabelingRestraint:
    """Fast direct-labeling restraint."""

    def test_create_empty(self):
        r = DirectLabelingRestraint()
        assert r.evaluate() == 0.0

    def test_no_xyz(self):
        r = DirectLabelingRestraint()
        r.add_site(0, distance=50.0)
        r.add_site(1, distance=50.0)
        assert r.evaluate() == 0.0

    def test_two_sites(self):
        xyz = np.array([[0., 0., 0.], [30., 0., 0.]], dtype=np.float64)
        r = DirectLabelingRestraint(xyz)
        r.add_site(0, distance=30.0, error_neg=3.0, error_pos=3.0)
        r.add_site(1, distance=30.0, error_neg=3.0, error_pos=3.0)
        # perfect fit → chi2 = 0
        assert r.evaluate() == pytest.approx(0.0, abs=1e-10)

    def test_two_sites_mismatch(self):
        xyz = np.array([[0., 0., 0.], [30., 0., 0.]], dtype=np.float64)
        r = DirectLabelingRestraint(xyz)
        r.add_site(0, distance=35.0, error_neg=3.0, error_pos=3.0)
        r.add_site(1, distance=35.0, error_neg=3.0, error_pos=3.0)
        # |30 - 35| = 5, / 3 → (5/3)^2 ≈ 2.78
        expected = (5.0 / 3.0) ** 2
        assert r.evaluate() == pytest.approx(expected, rel=1e-9)

    def test_three_sites(self):
        xyz = np.array([[0., 0., 0.],
                        [10., 0., 0.],
                        [20., 0., 0.]], dtype=np.float64)
        r = DirectLabelingRestraint(xyz)
        r.add_site(0, distance=10.0, error_neg=1.0, error_pos=1.0)
        r.add_site(1, distance=10.0, error_neg=1.0, error_pos=1.0)
        r.add_site(2, distance=20.0, error_neg=1.0, error_pos=1.0)
        # Sites 0-1: |10-10|=0 → 0
        # Sites 0-2: |20-20|=0 → 0
        # Sites 1-2: |10-20|=10 → (10/1)^2 = 100
        assert r.evaluate() == pytest.approx(100.0, abs=1e-9)

    def test_weight(self):
        xyz = np.array([[0., 0., 0.], [30., 0., 0.]], dtype=np.float64)
        r = DirectLabelingRestraint(xyz, weight=2.0)
        r.add_site(0, distance=25.0, error_neg=5.0, error_pos=5.0)
        r.add_site(1, distance=25.0, error_neg=5.0, error_pos=5.0)
        expected = 2.0 * (5.0 / 5.0) ** 2
        assert r.evaluate() == pytest.approx(expected, rel=1e-9)


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
