"""The in-tree DCD reader must agree with MDAnalysis, frame for frame.

The rotamer libraries ship as PDB + DCD pairs and used to be read through
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
    """`test/input/`, not `data/`.

    The shipped rotamer libraries became BinaryCIF on 2026-08-19 and no `.dcd`
    is installed any more. The reader stays, because a user's own library may
    still be one -- so it needs a fixture, and a fixture belongs in
    `test/input/`, which IMP does not install. That is the whole reason the
    directory exists: `em` keeps 68 MB there without a single user seeing it.
    """
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "input")


def _dcd_files():
    d = _library_dir()
    return sorted(
        os.path.join(d, f) for f in os.listdir(d) if f.endswith(".dcd")
    )


class Tests(IMP.test.TestCase):

    def test_header_is_plausible(self):
        """Every bundled DCD has a readable, self-consistent header"""
        files = _dcd_files()
        self.assertGreater(len(files), 0, "no DCD files in the rotamer library")
        # a sample here; expensive_test_dcd_reader.py sweeps every file
        for path in files[:8]:
            head = read_dcd_header(path)
            self.assertGreater(head.n_frames, 0, path)
            self.assertGreater(head.n_atoms, 0, path)
            self.assertIn(head.endianness, ("<", ">"))

    def test_shapes_and_finiteness(self):
        """Coordinates come back as (frames, atoms, 3) and are finite"""
        for path in _dcd_files()[:5]:
            head = read_dcd_header(path)
            coords = np.asarray(read_dcd(path)).reshape(
                -1, head.n_atoms, 3)
            self.assertEqual(
                coords.shape, (head.n_frames, head.n_atoms, 3), path)
            self.assertTrue(np.isfinite(coords).all(), path)

    def test_max_frames_truncates(self):
        """max_frames stops the read without disturbing the values"""
        path = _dcd_files()[0]
        head = read_dcd_header(path)
        full = np.asarray(read_dcd(path)).reshape(-1, head.n_atoms, 3)
        part = np.asarray(read_dcd(path, max_frames=3)).reshape(
            -1, head.n_atoms, 3)
        self.assertEqual(part.shape[0], min(3, full.shape[0]))
        np.testing.assert_array_equal(part, full[:part.shape[0]])

    def test_rejects_a_non_dcd(self):
        """A file that is not a DCD raises rather than returning nonsense"""
        with IMP.test.temporary_directory() as tmp:
            bogus = os.path.join(tmp, "not.dcd")
            with open(bogus, "wb") as fh:
                fh.write(b"\x00" * 512)
            # `IMP::ValueException` is itself a `ValueError`, and the Python's
            # own `DCDFormatError` was a `ValueError` subclass that nothing but
            # this test ever named.
            self.assertRaises(ValueError, read_dcd, bogus)

if __name__ == '__main__':
    IMP.test.main()
