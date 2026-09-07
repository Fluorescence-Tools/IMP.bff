import unittest
import os
import IMP.bff
from IMP.bff import read_mol2_component
from IMP.bff import get_structure_dir

class TestTopology(unittest.TestCase):
    def setUp(self):
        self.cx4_mol2 = str(get_structure_dir("cx4.mol2"))

    def test_parse_mol2(self):
        if not os.path.exists(self.cx4_mol2):
            self.skipTest("cx4.mol2 not found")
        
        component = read_mol2_component(self.cx4_mol2, "CX4")
        self.assertGreater(len(component.atoms), 0)
        self.assertGreater(len(component.bonds), 0)

        # the typed atom, not a dict of the same fields
        first_atom = component.atoms[0]
        self.assertTrue(first_atom.atom_name)
        self.assertTrue(first_atom.element)

    def test_graph_and_angles(self):
        if not os.path.exists(self.cx4_mol2):
            self.skipTest("cx4.mol2 not found")
            
        component = read_mol2_component(self.cx4_mol2, "CX4")
        graph = IMP.bff.MolecularGraph([tuple(b) for b in component.bonds])
        self.assertEqual(len(graph.get_nodes()), len(component.atoms))

        self.assertGreater(len(graph.get_angles()), 0)
        self.assertGreater(len(graph.get_dihedrals()), 0)

if __name__ == "__main__":
    unittest.main()
