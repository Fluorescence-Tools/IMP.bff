"""What a system CIF records for its MOL2 files, and whether it reads back.

A system records each component's structure file *relative to the CIF*, so a
directory of systems moves as a unit. The relative path was computed by
counting the **separators** of the base directory rather than its segments,
which is one short whenever the base has no trailing separator -- always. A
system written anywhere that is not a sibling of its MOL2 files therefore
recorded a path one level too shallow, and every command that read it back
failed with a file that "does not exist":

    base   /var/folders/cl/xxxx/T/tmp1        (6 segments)
    mol2   /Users/me/dye.mol2
    wrote  ../../../../../Users/me/dye.mol2   (5 ups)
    read   /var/Users/me/dye.mol2

`imp_bff simulate` was dead for any system built outside the data tree.
"""

import os

import pytest

import IMP.bff

#: `create_forcefield_system` is overloaded (typed specs, or the JSON string),
#: and SWIG turns off keyword arguments for an overloaded function -- so every
#: parameter up to `relative_to` has to be given, in order. These are its own
#: defaults; `bin/imp_bff` passes them the same way.
DEFAULTS = (2000.0, 400.0, 12.0, 1.5, 40.0, 180.0, 120.0, 220.0,
            20000, 100, 1.7, 12.0, 5.0, 6.0, 200)


def _build(components, relative_to=""):
    return IMP.bff.create_forcefield_system(components, *DEFAULTS, relative_to)


def _components():
    return [
        IMP.bff.FFComponentSpec(
            "CX4", str(IMP.bff.get_structure_dir("cx4.mol2")),
            str(IMP.bff.get_template_dir("cx4.template.cif")), "fixed"),
        IMP.bff.FFComponentSpec(
            "atto655", str(IMP.bff.get_structure_dir("atto655.mol2")),
            str(IMP.bff.get_template_dir("atto655.template.cif")), "mobile"),
    ]


def test_a_recorded_path_resolves_from_the_cifs_directory(tmp_path):
    """The test's own `tmp_path` is the case that used to fail: deep, and on
    the other side of the filesystem from the shipped structures."""
    system = _build(_components(), str(tmp_path))
    out = tmp_path / "system.cif"
    IMP.bff.write_probe_forcefield_cif(str(out), system)

    back = IMP.bff.read_forcefield_cif(str(out))
    assert set(back.components) == {"CX4", "atto655"}
    for name, spec in back.components.items():
        assert not os.path.isabs(spec.mol2_path), name
        resolved = os.path.normpath(os.path.join(str(tmp_path),
                                                 spec.mol2_path))
        assert os.path.exists(resolved), f"{name}: {spec.mol2_path}"


def test_a_path_under_the_base_needs_no_climbing(tmp_path):
    """The ordinary case: a system beside its own structures."""
    mol2 = tmp_path / "structures" / "dye.mol2"
    mol2.parent.mkdir()
    mol2.write_bytes(
        open(IMP.bff.get_structure_dir("atto655.mol2"), "rb").read())

    system = _build([IMP.bff.FFComponentSpec("dye", str(mol2), "", "mobile")],
                    str(tmp_path))
    out = tmp_path / "system.cif"
    IMP.bff.write_probe_forcefield_cif(str(out), system)

    spec = IMP.bff.read_forcefield_cif(str(out)).components["dye"]
    assert spec.mol2_path == "structures/dye.mol2"


def test_without_a_base_the_path_is_left_alone(tmp_path):
    system = _build(_components())
    out = tmp_path / "system.cif"
    IMP.bff.write_probe_forcefield_cif(str(out), system)
    for spec in IMP.bff.read_forcefield_cif(str(out)).components.values():
        assert os.path.isabs(spec.mol2_path)
        assert os.path.exists(spec.mol2_path)
