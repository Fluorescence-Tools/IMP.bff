import os
import tempfile
import unittest
from IMP.bff import ComponentTemplate, Feature, FeatureAtom, Improper
from IMP.bff import read_component_template_cif, write_component_template_cif

class TestTemplate(unittest.TestCase):
    def test_read_write_template(self):
        template = ComponentTemplate()
        template.name = "test_template"
        f1 = Feature()
        f1.rb = True
        f1.md_fixed = False
        f1.feature_type = "dof"
        f1.region_color = "red"
        atom = FeatureAtom()
        atom.name = "C1"
        atom.occurrence = 1
        f1.atoms.append(atom)
        template.features["f1"] = f1
        improper = Improper()
        improper.center_atom = "C1"
        improper.type = "ring"
        template.impropers.append(improper)

        with tempfile.TemporaryDirectory() as tmpdir:
            cif_path = os.path.join(tmpdir, "template.cif")
            write_component_template_cif(cif_path, template)
            self.assertTrue(os.path.exists(cif_path))

            read_template = read_component_template_cif(cif_path)
            self.assertEqual(read_template.name, template.name)
            self.assertIn("f1", read_template.features)
            self.assertTrue(read_template.features["f1"].rb)
            self.assertEqual(read_template.features["f1"].region_color, "red")
            self.assertEqual(read_template.features["f1"].atoms[0].name, "C1")
            self.assertEqual(len(read_template.impropers), 1)
            self.assertEqual(read_template.impropers[0].center_atom, "C1")
            self.assertEqual(read_template.impropers[0].type, "ring")

if __name__ == "__main__":
    unittest.main()
