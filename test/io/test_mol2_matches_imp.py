"""The core Tripos reader answers exactly as IMP's `read_mol2` road did.

`read_mol2_component` used to read the file through `IMP::atom::read_mol2`
and walk the hierarchy; since PRD-137 step 5c the core walks the file itself
the way IMP does (the MOLECULE header's third line as residue type, ATOM and
BOND sections to an empty line or the next record, a bond kept only when both
atoms are in the molecule). This rebuilds the old road in Python from IMP.atom
and compares every field on every shipped .mol2. Skipped without IMP.atom.
"""
import glob
import os
import unittest

import IMP.bff as bff

try:
    import IMP
    import IMP.atom
    _HAVE_ATOM = True
except ImportError:
    _HAVE_ATOM = False

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_MOL2 = sorted(set(sum((glob.glob(os.path.join(_ROOT, d, "**", "*.mol2"), recursive=True)
                        for d in ("data", "examples", "test")), [])))


def _imp_road(path):
    m = IMP.Model()
    IMP.set_log_level(IMP.SILENT)
    h = IMP.atom.read_mol2(path, m)
    names = dict(bff.read_mol2_atom_names(path))
    atoms = []
    for a in IMP.atom.get_by_type(h, IMP.atom.ATOM_TYPE):
        at = IMP.atom.Atom(a)
        serial = at.get_input_index()
        name = names.get(serial, "")
        res = IMP.atom.get_residue(at, True)
        resname = res.get_residue_type().get_string() if res else "UNK"
        xyz = IMP.core.XYZ(a).get_coordinates()
        atoms.append((serial, name, resname or "UNK", bff.element_from_atom_name(name),
                      xyz[0], xyz[1], xyz[2]))
    bonds = set()
    for b in IMP.atom.get_internal_bonds(h):
        bd = IMP.atom.Bond(b)
        s1 = IMP.atom.Atom(bd.get_bonded(0).get_particle()).get_input_index()
        s2 = IMP.atom.Atom(bd.get_bonded(1).get_particle()).get_input_index()
        bonds.add((min(s1, s2), max(s1, s2)))
    return atoms, sorted(bonds)


@unittest.skipUnless(_HAVE_ATOM, "IMP.atom is the reference reader; not available")
class TestCoreMol2MatchesIMP(unittest.TestCase):

    def test_every_shipped_mol2_reads_the_same(self):
        self.assertTrue(_MOL2, "no .mol2 files under data/, examples/ or test/")
        for path in sorted(set(_MOL2)):
            with self.subTest(mol2=os.path.relpath(path, _ROOT)):
                c = bff.read_mol2_component(path, "X")
                ours = [(a.serial, a.atom_name, a.resname, a.element, a.x, a.y, a.z) for a in c.atoms]
                ref_atoms, ref_bonds = _imp_road(path)
                self.assertEqual(ours, ref_atoms)
                self.assertEqual([tuple(p) for p in c.bonds], ref_bonds)


if __name__ == "__main__":
    unittest.main()
