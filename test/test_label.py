"""Tests for ``IMP.bff.label`` — dye label distribution models."""

from __future__ import annotations

import numpy as np
import pytest

from IMP.bff.representation import LabelDistribution, LabelDistributionAV, DyeDistributionNormal


class TestLabelDistribution:
    """Abstract base class."""

    def test_cannot_instantiate(self):
        """LabelDistribution is abstract."""
        with pytest.raises(TypeError):
            LabelDistribution()


class TestDyeDistributionNormal:
    """Gaussian label distribution."""

    def test_create(self):
        dd = DyeDistributionNormal(
            origin=np.array([0.0, 0.0, 10.0]),
            width=6.0,
            position_name="donor",
        )
        assert dd.position_name == "donor"
        assert dd.n_points > 0

    def test_mean_position(self):
        dd = DyeDistributionNormal(
            origin=np.array([5.0, 0.0, 0.0]),
            width=1.0,
        )
        mp = dd.mean_position
        np.testing.assert_allclose(mp, [5.0, 0.0, 0.0], atol=1.0)

    def test_dRmp(self):
        dd1 = DyeDistributionNormal(
            origin=np.array([0.0, 0.0, 0.0]), width=2.0,
        )
        dd2 = DyeDistributionNormal(
            origin=np.array([3.0, 4.0, 0.0]), width=2.0,
        )
        d = dd1.dRmp(dd2)
        assert d == pytest.approx(5.0, abs=2.0)

    def test_dRDA_positive(self):
        dd1 = DyeDistributionNormal(
            origin=np.array([0.0, 0.0, 0.0]), width=4.0,
        )
        dd2 = DyeDistributionNormal(
            origin=np.array([10.0, 0.0, 0.0]), width=4.0,
        )
        d = dd1.dRDA(dd2, n_samples=10000)
        assert d > 0.0


class TestLabelDistributionAV:
    """AV-based label distribution (lazy backend)."""

    def test_create_instantiate(self):
        """Instantiation should succeed even without AV backends."""
        xyz = np.array([[0.0, 0.0, 0.0],
                        [5.0, 0.0, 0.0]], dtype=np.float64)
        vdw = np.array([1.5, 1.5], dtype=np.float64)
        dd = LabelDistributionAV(
            atoms_xyz=xyz, atoms_vdw=vdw,
            residue_seq_number=0, atom_name="CB",
            linker_length=10.0,
        )
        assert dd.position_name == ""
        np.testing.assert_allclose(dd.origin, [0.0, 0.0, 0.0], atol=1e-6)

    def test_get_basic_av_raises_without_backend(self, monkeypatch):
        """Accessing the AV without a backend should raise."""
        xyz = np.array([[0.0, 0.0, 0.0],
                        [5.0, 0.0, 0.0]], dtype=np.float64)
        vdw = np.array([1.5, 1.5], dtype=np.float64)
        dd = LabelDistributionAV(
            atoms_xyz=xyz, atoms_vdw=vdw,
            residue_seq_number=0, atom_name="CB",
        )
        # Force the condition the test is about: a build without IMP.bff's AV
        # decorator, which since PRD-112 stage 1 is the only backend. This
        # used to pass only because `_HAS_LABELLIB` was wrongly False on
        # machines that *did* have a working LabelLib — the check looked for
        # `LabelLib.AV`, which current builds do not expose.
        import IMP.bff.av.compute as _compute

        monkeypatch.setattr(_compute, "_HAS_IMP_BFF", False, raising=False)
        with pytest.raises(ImportError):
            dd.get_basic_av()

    def test_lazy_computation(self):
        """AV is not computed until accessed."""
        xyz = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        vdw = np.array([1.5], dtype=np.float64)
        dd = LabelDistributionAV(
            atoms_xyz=xyz, atoms_vdw=vdw,
            residue_seq_number=0, atom_name="CB",
        )
        assert dd._av is None  # not computed yet


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
