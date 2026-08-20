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
        import IMP.bff.representation.av as _compute

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


def test_find_atom_matches_the_full_scan_it_replaced():
    """`label._find_atom` against the exhaustive scan, on a real structure.

    The reference is the old implementation kept verbatim, because the risk is
    not "does it find an atom" but *which* atom: the match qualifies an atom by
    `AtomType` **or** by name -- a dye's atoms carry types IMP does not
    recognise and are found by name only -- and the scan returned the first in
    hierarchy order. A selection that reordered or over-selected would still
    find something.

    Includes a name that is in no residue ("ZZZ"), so the not-found path is
    covered too.
    """
    import IMP
    import IMP.atom
    from IMP.bff.tools import get_structure_dir
    from IMP.bff.label import _find_atom, _atom_name, _atom_type_from_name

    def exhaustive(hierarchy, chain_id, resnum, atom_name):
        target = _atom_type_from_name(atom_name)
        for a in IMP.atom.get_by_type(hierarchy, IMP.atom.ATOM_TYPE):
            at = IMP.atom.Atom(a)
            if target != IMP.atom.AtomType("UNK") and at.get_atom_type() == target:
                pass
            elif _atom_name(a).upper() == atom_name.upper():
                pass
            else:
                continue
            res_p = a.get_parent()
            if not IMP.atom.Residue.get_is_setup(res_p):
                continue
            if IMP.atom.Residue(res_p).get_index() != resnum:
                continue
            chain_p = res_p.get_parent()
            if not IMP.atom.Chain.get_is_setup(chain_p):
                continue
            if IMP.atom.Chain(chain_p).get_id() != chain_id:
                continue
            return a
        return None

    model = IMP.Model()
    IMP.set_log_level(IMP.SILENT)
    hierarchy = IMP.atom.read_pdb(
        str(get_structure_dir("1DG3.pdb")), model, IMP.atom.NonWaterPDBSelector())
    chain = IMP.atom.Chain(
        IMP.atom.get_by_type(hierarchy, IMP.atom.CHAIN_TYPE)[0]).get_id()
    residues = [IMP.atom.Residue(r)
                for r in IMP.atom.get_by_type(hierarchy, IMP.atom.RESIDUE_TYPE)][:40]

    found = 0
    for residue in residues:
        for name in ("CA", "N", "C", "CB", "O", "ZZZ"):
            want = exhaustive(hierarchy, chain, residue.get_index(), name)
            got = _find_atom(hierarchy, chain, residue.get_index(), name)
            if want is None:
                assert got is None, (residue.get_index(), name)
                continue
            assert got is not None, (residue.get_index(), name)
            assert (IMP.atom.Hierarchy(got).get_particle()
                    == IMP.atom.Hierarchy(want).get_particle())
            found += 1
    assert found > 100, found
