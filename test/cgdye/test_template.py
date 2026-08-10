import os
import tempfile
import unittest
from IMP.bff.cgdye.io.template_cif import read_cgdye_template, write_cgdye_template

class TestTemplate(unittest.TestCase):
    def test_read_write_template(self):
        template = {
            "name": "test_template",
            "features": {
                "f1": {
                    "rb": True,
                    "md_fixed": False,
                    "feature_type": "dof",
                    "region_color": "red",
                    "atoms": [{"name": "C1", "occurrence": 1}]
                }
            },
            "impropers": [{"center_atom": "C1", "type": "ring"}]
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            cif_path = os.path.join(tmpdir, "template.cif")
            write_cgdye_template(cif_path, template)
            self.assertTrue(os.path.exists(cif_path))

            read_template = read_cgdye_template(cif_path)
            self.assertEqual(read_template["name"], template["name"])
            self.assertIn("f1", read_template["features"])
            self.assertTrue(read_template["features"]["f1"]["rb"])
            self.assertEqual(len(read_template["impropers"]), 1)

if __name__ == "__main__":
    unittest.main()
