"""IMP.bff.io — the fps.json reader, the legacy C# formats, and structure IO.

The three used to be one module (``fret/io.py``) and are three now; the tests
reach them through the package, which is the surface a caller sees.
"""

import json
from pathlib import Path

import pytest

import IMP.bff.io as fio

REPO = Path(__file__).resolve().parents[2]


def test_fps_json_roundtrip(tmp_path):
    path = tmp_path / "labels.fps.json"
    positions = {"p1": {"residue_seq_number": 132, "atom_name": "CB",
                        "linker_length": 20.0}}
    distances = {"d1": {"position1_name": "p1", "position2_name": "p1",
                        "distance": 45.0, "error_neg": 2.0, "error_pos": 3.0,
                        "distance_type": "RDAMeanE",
                        "Forster_radius": 52.0}}
    score_sets = {"all": {"distances": ["d1"]}}
    extra = {"FlexFit": {"anything": 1}}
    fio.write_fps_json(path, positions, distances, score_sets, extra,
                       validate=True)
    p, d, s, e = fio.read_fps_json(path, validate=True)
    assert p == positions
    assert d == distances
    assert s == score_sets
    assert e == extra


def test_write_fps_json_validate_refuses_nonconforming(tmp_path):
    path = tmp_path / "bad.fps.json"
    with pytest.raises(ValueError, match="non-conforming"):
        fio.write_fps_json(
            path, {"p": {"linker_length": "long"}}, {}, validate=True)
    assert not path.exists()


def test_read_fps_json_reads_shipped_example():
    p, d, s, extra = fio.read_fps_json(
        REPO / "examples" / "structure" / "T4L" / "fret.fps.json",
        validate=True)
    assert len(p) == 17
    assert len(d) == 99
    assert len(s) == 7
    assert extra == {}


def test_evaluators_roundtrip_with_and_without_factory(tmp_path):
    path = tmp_path / "labels.fps.json"
    fio.write_fps_json(path, {}, {})
    fio.write_evaluators_json(
        path, [{"type": "distance", "cutoff": 30.0}])
    raw = fio.read_evaluators_json(path)
    assert raw == [{"type": "distance", "cutoff": 30.0}]

    class Ev:
        def __init__(self, d):
            self.d = d

    made = fio.read_evaluators_json(path, factory=Ev)
    assert len(made) == 1 and made[0].d["cutoff"] == 30.0

    def broken(d):
        raise RuntimeError("no")

    assert fio.read_evaluators_json(path, factory=broken) == []


def test_read_old_lps_txt_av1_av3_xyz(tmp_path):
    pdb = tmp_path / "mol.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A  12      11.000  12.000  13.000  1.00  0.00\n"
        "ATOM      2  CB  ALA A  12      14.000  15.000  16.000  1.00  0.00\n"
        "END\n")
    lps = tmp_path / "LPs.txt"
    lps.write_text(
        "# comment\n"
        "site1 mol Alexa488 AV1 20.0 4.5 3.5 2\n"
        "site2 mol Alexa647 AV3 21.0 4.0 9.0 3.0 1.5 1\n"
        "fixed mol none XYZ 1.0 2.0 3.0\n")
    (tmp_path / "Distances.txt").write_text(
        "RDAMeanE\n"
        "site1 site2 45.7 4.5 4.6 52.0\n")

    positions, distances, score_sets, extra = fio.read_fps_json(
        lps, pdb_paths=[str(pdb)])
    assert positions["site1"]["atom_name"] == "CB"
    assert positions["site1"]["residue_seq_number"] == 12
    assert positions["site1"]["simulation_type"] == "AV1"
    assert positions["site2"]["simulation_type"] == "AV3"
    assert positions["site2"]["radius2"] == 3.0
    assert positions["fixed"] == {
        "simulation_type": "XYZ", "x": 1.0, "y": 2.0, "z": 3.0, "body_id": 0}
    d = distances["site1_site2"]
    assert d["distance"] == 45.7
    assert d["distance_type"] == "RDAMeanE"
    assert d["Forster_radius"] == 52.0
    assert score_sets == {} and extra == {}
    # the converted payload conforms to the schema
    import IMP.bff.io.fps as fps_schema
    errors, _ = fps_schema.fps_schema_validate(
        {"Positions": positions, "Distances": distances})
    assert not errors, errors


def test_write_pdb_and_rmsd(tmp_path):
    import numpy as np
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
    out = tmp_path / "m.pdb"
    fio.write_pdb(coords, out, transform=np.array([1.0, 1.0, 1.0]))
    text = out.read_text()
    assert text.count("ATOM") == 3
    assert "2.000" in text  # translated x of atom 2

    shifted = coords + 5.0
    assert fio.compute_rmsd(coords, shifted) == pytest.approx(
        np.sqrt(3 * 25.0))
    assert fio.compute_rmsd(coords, shifted, superpose=True) == pytest.approx(
        0.0, abs=1e-9)


def test_write_rmf_is_self_contained(tmp_path):
    import numpy as np
    import RMF
    out = tmp_path / "m.rmf3"
    coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
    fio.write_rmf(coords, out, model_name="probe",
                  metadata={"origin": "test"})
    fh = RMF.open_rmf_file_read_only(str(out))
    assert fh.get_number_of_frames() == 1
    assert "origin" in fh.get_description()


def test_the_forcefield_cif_reader_is_cpp_and_matches_the_python():
    """`IMP.bff.read_forcefield_cif` against the 470-line Python it replaced.

    The Python is gone, so the reference here is a *round trip*: build a system,
    write it, read it back, and check every field of every category survives.
    That is what the old reader was checked against too, and it catches the one
    thing a count comparison does not -- the compact site-number scheme, where a
    term may name a number whose `_ff_site` row appears later in the file.
    """
    import IMP.bff
    from IMP.bff.io.cif import write_dye_forcefield_cif
    from IMP.bff import get_template_dir, get_structure_dir
    from IMP.bff.cgdye.topology import build_dye_protein_system

    system = build_dye_protein_system(
        str(get_structure_dir("cx4.mol2")), str(get_structure_dir("atto655.mol2")),
        "CX4", "atto655",
        protein_template=str(get_template_dir("cx4.template.cif")),
        dye_template=str(get_template_dir("atto655.template.cif")),
    )
    import tempfile, os
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "system.cif")
        write_dye_forcefield_cif(path, system)
        back = IMP.bff.read_forcefield_cif(path)

    assert back.name == system.name
    assert len(back.sites) == len(system.sites) == 138
    assert len(back.bonds) == len(system.bonds) == 146
    assert len(back.angles) == len(system.angles) == 262
    assert len(back.dihedrals) == len(system.dihedrals) == 207
    assert len(back.impropers) == len(system.impropers) == 81

    for a, b in zip(system.sites, back.sites):
        assert (a.id, a.component, a.atom_name, a.site_serial) == \
               (b.id, b.component, b.atom_name, b.site_serial)
        # the writer omits these when uniform, so absent must mean the default
        assert abs(a.radius - b.radius) < 1e-12 and b.radius > 0.0
        assert abs(a.mass - b.mass) < 1e-12 and b.mass > 0.0
    # 5e-4, not 1e-9: `ihm`'s writer emits three decimals ("1.635"), so a round
    # trip is exact only to half of the last place it wrote. That is the file's
    # precision, not the reader's -- reading the same file twice is identical.
    for a, b in zip(system.bonds, back.bonds):
        assert (a.site_a, a.site_b, a.type_id) == (b.site_a, b.site_b, b.type_id)
        assert abs(a.length - b.length) < 5e-4
    for a, b in zip(system.angles, back.angles):
        assert (a.site_a, a.site_b, a.site_c) == (b.site_a, b.site_b, b.site_c)
        assert abs(a.theta - b.theta) < 5e-4
    for name in ("dihedrals", "impropers"):
        for a, b in zip(getattr(system, name), getattr(back, name)):
            assert (a.site_a, a.site_b, a.site_c, a.site_d, a.type_id) == \
                   (b.site_a, b.site_b, b.site_c, b.site_d, b.type_id)
    assert sorted(back.groups) == sorted(system.groups)
    assert sorted(back.lj_types) == sorted(system.lj_types)
    assert sorted(back.improper_types) == sorted(system.improper_types)
    assert back.nonbonded.enabled == system.nonbonded.enabled
    assert abs(back.nonbonded.cutoff - system.nonbonded.cutoff) < 1e-12
    assert back.sampling.n_steps == system.sampling.n_steps
