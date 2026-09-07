"""The two `bin/` programs that had no test, run.

`imp_bff_traj2bcif`, `imp_bff_traj2drot` and `imp_bff_potentials2pto` are
exercised by the tests of the formats they write. These two were not, and both
were broken: `imp_bff_probe_pdb2cif` passed `--probe-id`'s `None` default into a
`const std::string&`, so its documented default path raised a `TypeError`
before it read anything, and `imp_bff_labelizer --show` handed a structure to a
container reader, which complained about EBML headers.

They are argparse programs with a `main(argv)`, so a test calls that.
"""

import importlib.machinery
import importlib.util
import sys
from pathlib import Path

import pytest

import IMP.bff

BIN = Path(__file__).resolve().parent.parent.parent / "bin"
PDB = IMP.bff.get_example_path("structure/T4L/3GUN.pdb")


def _program(name):
    path = BIN / name
    loader = importlib.machinery.SourceFileLoader(name, str(path))
    spec = importlib.util.spec_from_file_location(name, str(path),
                                                  loader=loader)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_probe_pdb2cif_converts_with_its_default_id(tmp_path):
    """`--probe-id` is empty by default, meaning "take it from the file". It was
    `None`, which is not a string, so the default invocation never ran."""
    out = tmp_path / "dye.cif"
    program = _program("imp_bff_probe_pdb2cif")
    assert program.main([str(IMP.bff.get_structure_dir("alexa488_r48.pdb")),
                         str(out)]) == 0
    assert out.is_file()
    text = out.read_text()
    assert "_atom_site" in text
    assert text.count("ATOM") + text.count("HETATM") > 10


def test_probe_pdb2cif_takes_an_explicit_id(tmp_path):
    out = tmp_path / "dye.cif"
    program = _program("imp_bff_probe_pdb2cif")
    assert program.main([str(IMP.bff.get_structure_dir("alexa488_r48.pdb")),
                         str(out), "--probe-id", "A48"]) == 0
    assert "A48" in out.read_text()


def test_labelizer_scores_a_structure(tmp_path):
    out = tmp_path / "t4l.mmfdb.pto"
    program = _program("imp_bff_labelizer")
    assert program.main([PDB, "--no-conservation", "-o", str(out)]) == 0
    assert out.is_file()
    assert len(IMP.bff.ll_read_pto_scores(str(out))) > 0


def test_labelizer_show_says_what_it_wants(tmp_path, capsys):
    """`--show` reads a container. Given a structure it used to fail with
    "does not begin with an EBML header", which is true and useless."""
    program = _program("imp_bff_labelizer")
    assert program.main([PDB, "--show"]) == 2
    assert "container" in capsys.readouterr().err


def test_labelizer_show_reads_what_it_wrote(tmp_path, capsys):
    out = tmp_path / "t4l.mmfdb.pto"
    program = _program("imp_bff_labelizer")
    assert program.main([PDB, "--no-conservation", "-o", str(out)]) == 0
    capsys.readouterr()
    assert program.main([str(out), "--show"]) == 0
    assert "arithmetic" in capsys.readouterr().out
