"""The `_cgprobe_*` component template, read from the files that ship it."""

import unittest

import pytest

import IMP.bff
from IMP.bff import read_component_template_cif


SHIPPED = ["atto655.template.cif", "cx4.template.cif"]


class TestTemplate(unittest.TestCase):
    def test_shipped_templates_parse(self):
        """Every shipped template reads into the fields the topology builder
        asks of it: named features, each with atoms, and impropers naming a
        centre."""
        for name in SHIPPED:
            template = read_component_template_cif(
                str(IMP.bff.get_template_dir(name)))
            self.assertTrue(template.name, name)
            self.assertTrue(template.features, f"{name}: no features")
            for fid, feature in template.features.items():
                self.assertTrue(feature.feature_type, f"{name}/{fid}: no type")
                for atom in feature.atoms:
                    self.assertTrue(atom.name, f"{name}/{fid}: unnamed atom")
            for improper in template.impropers:
                self.assertTrue(improper.center_atom, f"{name}: no centre")

    def test_region_features_are_the_coloured_ones(self):
        """`region_features` selects the features that declare a colour."""
        template = read_component_template_cif(
            str(IMP.bff.get_template_dir("atto655.template.cif")))
        regions = IMP.bff.region_features(template)
        coloured = {fid for fid, f in template.features.items() if f.region_color}
        self.assertEqual(set(regions.keys()), coloured)


def test_a_template_in_the_old_spelling_fails_loudly(tmp_path):
    """`_cgdye_*` became `_cgprobe_*` on 2026-08-27.

    An old template matches no category, so every handler stays silent and the
    reader hands back an empty value -- not even obviously empty, because an
    unnamed template takes the file's stem for a name. It says what happened
    and how to migrate instead.
    """
    shipped = IMP.bff.get_template_dir("atto655.template.cif")
    legacy = tmp_path / "legacy.template.cif"
    legacy.write_text(open(str(shipped)).read().replace("_cgprobe_", "_cgdye_"))

    with pytest.raises(IOError) as caught:
        read_component_template_cif(str(legacy))
    assert "_cgdye_" in str(caught.value)


if __name__ == "__main__":
    unittest.main()
