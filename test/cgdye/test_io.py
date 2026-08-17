import os
import tempfile
import unittest
from IMP.bff.cgdye.io.cif import read_dye_forcefield_cif, write_dye_forcefield_cif

class TestIO(unittest.TestCase):
    def test_read_write_ff_system(self):
        # Create a minimal system dictionary
        system = {
            "name": "test_system",
            "components": {
                "comp1": {"mol2": "dummy.mol2", "role": "fixed"}
            },
            "sites": [
                {"id": "comp1/A1", "component": "comp1", "atom_name": "A1", "site_serial": 1, "radius": 1.7, "mass": 12.0}
            ],
            "groups": {"comp1_all": ["comp1/A1"]},
            "rb_groups": {},
            "md_fixed_groups": {},
            "fixed_groups": ["comp1_all"],
            "bond_types": {"B1": {"k": 100.0}},
            "angle_types": {},
            "torsion_types": {},
            "improper_types": {},
            "bonds": [["comp1/A1", "comp1/A1", 1.0, "B1"]], # Dummy self-bond
            "angles": [],
            "dihedrals": [],
            "impropers": [],
            "nonbonded": {"enabled": True, "k": 5.0, "cutoff_A": 6.0},
            "sampling": {
                "temperature_K": 300.0,
                "friction_ps": 10.0,
                "timestep_fs": 0.25,
                "n_steps": 1000,
                "write_every": 100,
                "minimize_steps": 10,
            },
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            cif_path = os.path.join(tmpdir, "test.cif")
            
            # Note: write_dye_forcefield_cif expects mol2 files to exist if it's going to write _atom_site
            # But we can test the basic round-trip of topology categories.
            # We'll mock the _parse_struct_atom_site_rows to avoid file dependency
            import IMP.bff.cgdye.io.cif as cif_mod
            original_parse = cif_mod._parse_struct_atom_site_rows
            cif_mod._parse_struct_atom_site_rows = lambda path, asym, ent: []
            
            try:
                write_dye_forcefield_cif(cif_path, system)
                self.assertTrue(os.path.exists(cif_path))
                
                read_system = read_dye_forcefield_cif(cif_path)
                self.assertEqual(read_system["name"], system["name"])
                self.assertEqual(len(read_system["sites"]), len(system["sites"]))
                self.assertEqual(read_system["sampling"]["n_steps"], system["sampling"]["n_steps"])
            finally:
                cif_mod._parse_struct_atom_site_rows = original_parse

if __name__ == "__main__":
    unittest.main()
