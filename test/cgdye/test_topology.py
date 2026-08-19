import unittest
import os
from IMP.bff.cgdye.topology import parse_dye_mol2, build_graph, build_angles, build_dihedrals
from IMP.bff.tools import get_structure_dir

class TestTopology(unittest.TestCase):
    def setUp(self):
        self.cx4_mol2 = str(get_structure_dir("cx4.mol2"))

    def test_parse_mol2(self):
        if not os.path.exists(self.cx4_mol2):
            self.skipTest("cx4.mol2 not found")
        
        atoms, bonds = parse_dye_mol2(self.cx4_mol2, "CX4")
        self.assertGreater(len(atoms), 0)
        self.assertGreater(len(bonds), 0)
        
        # Check an atom structure
        first_atom = next(iter(atoms.values()))
        self.assertIn("atom_name", first_atom)
        self.assertIn("element", first_atom)

    def test_graph_and_angles(self):
        if not os.path.exists(self.cx4_mol2):
            self.skipTest("cx4.mol2 not found")
            
        atoms, bonds = parse_dye_mol2(self.cx4_mol2, "CX4")
        graph = build_graph(bonds)
        self.assertEqual(len(graph), len(atoms))
        
        angles = build_angles(graph)
        self.assertGreater(len(angles), 0)
        
        dihedrals = build_dihedrals(graph)
        self.assertGreater(len(dihedrals), 0)

if __name__ == "__main__":
    unittest.main()
