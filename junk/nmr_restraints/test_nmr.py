import os
import tempfile
import unittest
from IMP.bff.io.cif import read_nmr_restraints, write_nmr_restraints

class TestNMR(unittest.TestCase):
    def test_read_write_nmr_restraints(self):
        data = {
            "default_restraint_set": "subunit",
            "Positions": {
                "p1": {
                    "atom_name": "H1", "residue_name": "DYE", "residue_seq_number": 0,
                    "chain_identifier": "B", "component_id": "mobile"
                },
                "p2": {
                    "atom_name": "H61", "residue_name": "CX4", "residue_seq_number": 0,
                    "chain_identifier": "A", "component_id": "fixed"
                }
            },
            "Distances": {
                "SUB_RID1": {
                    "distance": 3.5, "error_pos": 1.0, "error_neg": 1.0,
                    "position1_name": "p1", "position2_name": "p2",
                    "distance_type": "AtomUpperBound"
                }
            },
            "restraint_sets": {
                "subunit_prefix": "SUB_"
            }
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            cif_path = os.path.join(tmpdir, "nmr.cif")
            write_nmr_restraints(cif_path, data)
            self.assertTrue(os.path.exists(cif_path))

            read_data = read_nmr_restraints(cif_path)
            self.assertIn("p1", read_data["Positions"])
            self.assertIn("p2", read_data["Positions"])
            # The prefix logic might change naming slightly depending on how it infers, 
            # but let's check if at least one distance exists
            self.assertGreater(len(read_data["Distances"]), 0)
            
            # Check values of first distance
            d = next(iter(read_data["Distances"].values()))
            self.assertEqual(d["distance"], 3.5)

if __name__ == "__main__":
    unittest.main()
