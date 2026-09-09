"""The compiled command line: ``IMP.bff.bin_main`` and its subs.

The bin/ scripts whose subject is the core became C++ subcommands
(src/Bin.cpp), so they ship with the wheel and run with no IMP. These tests
drive the dispatcher the way the console script does -- a list of words, an
exit code -- and check the one property the subs exist for: the round trip
through the container happens inside the library, with no Python in it.
"""

import numpy
import pytest

import IMP.bff


@pytest.fixture
def ensemble():
    """A small ensemble of one connected molecule: six atoms in a chain,
    each 1.4 A from the last, so bond perception finds a spanning tree."""
    rng = numpy.random.default_rng(7)
    base = numpy.array([[i * 1.4, 0.0, 0.0] for i in range(6)])
    xyz = base[None, :, :] + rng.normal(0, 0.05, (40, 6, 3))
    return numpy.ascontiguousarray(xyz, dtype=numpy.float64)


def _write_library(tmp_path, xyz, name):
    names = ["N", "CA", "C", "O", "CB", "SG"]
    path = tmp_path / name
    weights = numpy.ones(len(xyz), dtype=numpy.float64)
    IMP.bff.write_drot(str(path), numpy.ascontiguousarray(xyz).reshape(-1),
                       names, ["N", "C", "C", "O", "C", "S"], ["LYS"] * 6,
                       weights)
    return path


def test_help_lists_the_subs(capfd):
    # capfd, not capsys: the dispatcher writes from C++, at the fd
    assert IMP.bff.bin_main(["help"]) == 0
    out = capfd.readouterr().out
    assert "pdb2cif" in out
    assert "traj2drot" in out


def test_unknown_sub_is_a_usage_error():
    assert IMP.bff.bin_main(["no-such-sub"]) != 0


def test_pdb2cif_writes_an_atom_site_table(tmp_path, ensemble):
    """A PDB written from the library comes back as mmCIF."""
    coords = [c for frame in ensemble[:1] for c in frame.reshape(-1)]
    pdb = tmp_path / "probe.pdb"
    IMP.bff.write_pdb(coords, str(pdb), "A", "LYS")
    cif = tmp_path / "probe.cif"
    assert IMP.bff.bin_main(["pdb2cif", str(pdb), str(cif),
                             "--probe-id", "LYP"]) == 0
    text = cif.read_text()
    assert "_atom_site" in text
    assert "LYP" in text


def test_traj2drot_round_trips_through_the_container(tmp_path, ensemble):
    """The library the sub reads and the one it writes agree to the lossless
    tolerance, with the verification the sub runs itself."""
    src = _write_library(tmp_path, ensemble, "src.drot.pto")
    dst = tmp_path / "out.drot.pto"
    assert IMP.bff.bin_main(["traj2drot", str(src), str(dst)]) == 0
    lib = IMP.bff.read_drot(str(dst))
    assert lib.n_rotamers == 40
    assert lib.n_atoms == 6


def test_traj2drot_bundle_reads_back_by_library_name(tmp_path, ensemble):
    one = _write_library(tmp_path, ensemble, "one.drot.pto")
    two = _write_library(tmp_path, ensemble[:7], "two.drot.pto")
    family = tmp_path / "family.pto"
    assert IMP.bff.bin_main(
        ["traj2drot", "--bundle", str(family), str(one), str(two)]) == 0
    assert IMP.bff.drot_catalog(str(family)) == ("one", "two")
    assert IMP.bff.read_drot(str(family), "two").n_rotamers == 7


def test_bad_usage_exits_nonzero(tmp_path, ensemble):
    src = _write_library(tmp_path, ensemble, "src.drot.pto")
    # src without dst: the sub's own check, exit 1
    assert IMP.bff.bin_main(["traj2drot", str(src)]) == 1
