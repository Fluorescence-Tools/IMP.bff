"""Tests for ``IMP.bff.label`` — dye label distribution models."""

from __future__ import annotations

import numpy as np
import pytest

from IMP.bff import LabelDistribution, LabelDistributionAV, ProbeDistributionNormal


class TestLabelDistribution:
    """The base a label distribution shares."""

    def test_cannot_instantiate(self):
        """It has a pure virtual `do_compute`, so SWIG gives it no constructor."""
        with pytest.raises((TypeError, AttributeError)):
            LabelDistribution()


class TestDyeDistributionNormal:
    """Gaussian label distribution."""

    def test_create(self):
        dd = ProbeDistributionNormal(
            origin=np.array([0.0, 0.0, 10.0]),
            width=6.0,
            position_name="donor",
        )
        assert dd.get_position_name() == "donor"
        assert dd.n_points > 0

    def test_mean_position(self):
        dd = ProbeDistributionNormal(
            origin=np.array([5.0, 0.0, 0.0]),
            width=1.0,
        )
        mp = dd.get_mean_position()
        np.testing.assert_allclose(mp, [5.0, 0.0, 0.0], atol=1.0)

    def test_dRmp(self):
        dd1 = ProbeDistributionNormal(
            origin=np.array([0.0, 0.0, 0.0]), width=2.0,
        )
        dd2 = ProbeDistributionNormal(
            origin=np.array([3.0, 4.0, 0.0]), width=2.0,
        )
        d = dd1.dRmp(dd2)
        assert d == pytest.approx(5.0, abs=2.0)

    def test_dRDA_positive(self):
        dd1 = ProbeDistributionNormal(
            origin=np.array([0.0, 0.0, 0.0]), width=4.0,
        )
        dd2 = ProbeDistributionNormal(
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
        dd = LabelDistributionAV(atoms_xyz=xyz, atoms_vdw=vdw,
                                 linker_length=10.0)
        assert dd.get_position_name() == ""
        # No `source_xyz` means the first obstacle.
        np.testing.assert_allclose(dd.get_origin(), [0.0, 0.0, 0.0], atol=1e-6)


    def test_lazy_computation(self):
        """The volume is not built until it is asked for.

        Building one is the expensive thing a label does, and a caller often
        holds several before asking any of them anything. The check is that
        constructing is fast -- the private flag it used to read is now on the
        C++ side.
        """
        import time

        xyz = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        vdw = np.array([1.5], dtype=np.float64)
        t0 = time.perf_counter()
        dd = LabelDistributionAV(atoms_xyz=xyz, atoms_vdw=vdw)
        assert time.perf_counter() - t0 < 0.05
        assert dd.n_points > 0


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


def test_select_atoms_matches_the_full_scan_it_replaced():
    """`select_atoms` against an exhaustive scan, on a real structure.

    The atom lookup is C++ now (`ProbeAttachment.h`), reached through a strip mask
    rather than a private helper. What this pins is unchanged and is the part
    that is easy to get wrong: not *whether* an atom is found but **which** --
    the first in hierarchy order, in the right chain and residue, matched on
    the name as the structure spells it (a dye's atoms carry types IMP does
    not recognise, so a type-based match would miss them).

    "ZZZ" is in no residue, so the not-found path is covered too.
    """
    import IMP
    import IMP.atom
    from IMP.bff import atom_name, get_structure_dir, select_atoms

    def exhaustive(hierarchy, chain_id, resnum, wanted):
        for a in IMP.atom.get_by_type(hierarchy, IMP.atom.ATOM_TYPE):
            if atom_name(IMP.atom.Atom(a)).upper() != wanted.upper():
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
            got = select_atoms(
                hierarchy,
                f"chain {chain} and resid {residue.get_index()} and name {name}")
            if want is None:
                assert not got, (residue.get_index(), name)
                continue
            assert got, (residue.get_index(), name)
            assert got[0] == IMP.atom.Hierarchy(want).get_particle()
            found += 1
    assert found > 100, found
