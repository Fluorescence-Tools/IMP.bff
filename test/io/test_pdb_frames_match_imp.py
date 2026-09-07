"""The core PDB reader answers exactly as IMP's reader did, on every shipped PDB.

`load_protein_frames` and `load_structure` used to read a PDB through
`IMP::atom::read_pdb` / `read_multimodel_pdb` and the connection layer's
`protein_frame_from_hierarchy`. Since PRD-137 step 5c the core reads the file
itself (`internal/PdbFrames.h`), reproducing IMP's rules: the `MODEL` split,
the NonWater selector (blank/`A` alternate locations, no `HOH`/`DOD`), the
`HET:`-prefixed raw type for `HETATM`, `UNK` for empty names, the residue
index, the one-character chain.

Both roads exist in the IMP build, so this test walks every `.pdb` under
examples/ and test/ and compares the two answers field by field. It needs
IMP.atom, so it is skipped in a standalone build.
"""
import glob
import os
import unittest

import numpy as np

import IMP.bff as bff

try:
    import IMP
    import IMP.atom
    _HAVE_ATOM = True
except ImportError:  # the standalone core
    _HAVE_ATOM = False

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_PDBS = sorted(p for p in glob.glob(os.path.join(_ROOT, "examples", "**", "*.pdb"), recursive=True)
               + glob.glob(os.path.join(_ROOT, "test", "**", "*.pdb"), recursive=True)
               if "cgprobe" not in p)


@unittest.skipUnless(_HAVE_ATOM, "IMP.atom is the reference reader; not available")
class TestCoreReaderMatchesIMP(unittest.TestCase):

    def _imp_frames(self, pdb):
        m = IMP.Model()
        hs = IMP.atom.read_multimodel_pdb(pdb, m, IMP.atom.NonWaterPDBSelector())
        return [bff.protein_frame_from_hierarchy(h) for h in hs]

    def test_every_shipped_pdb_reads_the_same(self):
        self.assertTrue(_PDBS, "no PDBs found under examples/ or test/")
        for pdb in _PDBS:
            with self.subTest(pdb=os.path.relpath(pdb, _ROOT)):
                ours = bff.load_protein_frames(pdb)
                ref = self._imp_frames(pdb)
                self.assertEqual(len(ours), len(ref), "frame count")
                for a, b in zip(ours, ref):
                    self.assertEqual(a.get_n_atoms(), b.get_n_atoms())
                    np.testing.assert_array_equal(np.asarray(a.get_coords()), np.asarray(b.get_coords()))
                    self.assertEqual(list(a.residue_indices), list(b.residue_indices))
                    self.assertEqual(list(a.atom_names), list(b.atom_names))
                    self.assertEqual(list(a.atom_types), list(b.atom_types))
                    self.assertEqual(list(a.resnames), list(b.resnames))
                    self.assertEqual(list(a.chain_ids), list(b.chain_ids))

    def test_load_structure_is_the_first_frame_without_water(self):
        for pdb in _PDBS:
            with self.subTest(pdb=os.path.relpath(pdb, _ROOT)):
                ours = np.asarray(bff.load_structure(pdb)).reshape(-1, 3)
                m = IMP.Model()
                h = IMP.atom.read_pdb(pdb, m, IMP.atom.NonWaterPDBSelector())
                ref = np.asarray(bff.structure_coordinates(h)).reshape(-1, 3)
                np.testing.assert_array_equal(ours, ref)


if __name__ == "__main__":
    unittest.main()
