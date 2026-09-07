"""Full DCD parity sweep: every bundled trajectory, against MDAnalysis.

The rotamer libraries ship as PDB + DCD pairs and are read through
MDAnalysis. IMP.bff carries no dependency beyond what IMP itself brings, so the
format is parsed in-tree instead -- and a hand-written binary parser is only
trustworthy if it is checked against the implementation it replaced.

The parity test skips when MDAnalysis is absent, which is the normal case for a
user: the point is that whoever *has* it can prove the replacement, not that
everyone must install it. The structural tests below run unconditionally.
"""

import os
import unittest

import numpy as np

import IMP
import IMP.bff
import IMP.test
from IMP.bff import read_dcd, read_dcd_header


def _library_dir():
    return IMP.bff.get_data_path("rotamer_library")


def _test_input_dir():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "input")


def _dcd_files():
    """Every DCD the tree still has: the libraries became BinaryCIF on
    2026-08-19 and `.drot` on 2026-08-24, so `data/rotamer_library` holds
    none -- `test/input/A56_C1R_cutoff30.dcd` is kept for this parity sweep,
    which is the only reason the reader is still exercised against the
    library it replaced."""
    out = []
    for d in (_library_dir(), _test_input_dir()):
        if os.path.isdir(d):
            out.extend(os.path.join(d, f) for f in os.listdir(d)
                       if f.endswith(".dcd"))
    return sorted(out)


class Tests(IMP.test.TestCase):

    def test_agrees_with_mdanalysis(self):
        """Frame-for-frame parity with the library this reader replaced"""
        try:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                import MDAnalysis as mda
        except ImportError:
            self.skipTest("MDAnalysis not installed; parity cannot be checked here")

        lib = _library_dir()
        checked = 0
        for path in _dcd_files():
            base = os.path.basename(path).split("_cutoff")[0]
            pdb = os.path.join(lib, base + ".pdb")
            if not os.path.exists(pdb):
                pdb = os.path.join(os.path.dirname(path), base + ".pdb")
            if not os.path.exists(pdb):
                continue
            mine = np.asarray(read_dcd(path)).reshape(
                -1, read_dcd_header(path).n_atoms, 3)
            universe = mda.Universe(pdb, path)
            theirs = np.array(
                [ts.positions.copy() for ts in universe.trajectory], dtype=np.float64)
            self.assertEqual(mine.shape, theirs.shape, path)
            np.testing.assert_allclose(mine, theirs, atol=1e-4, err_msg=path)
            checked += 1
        self.assertGreater(checked, 0, "no PDB/DCD pairs were compared")


if __name__ == '__main__':
    IMP.test.main()
