"""Positions that are a *point* rather than a volume (PRD-121 G3).

FPS's `AVSimlationType.None`. Two spellings reach it and they differ only in
where the point comes from: `XYZ` is a mean dye position measured once, in some
other structure's frame, and `ATOM` is an atom of this structure. Either way
there is no cloud, so `FilterEngine.cs:305-307` scores **every** distance
touching one as `R_mp`, whatever the file's `distance_type` says.

Before this, a network containing an `XYZ` position could not be built at all:
`search_labeling_site` was asked to resolve a position that names no chain, no
residue and no atom, matched the whole structure, and threw *ambiguous
labelling site*. A transported coordinate was computed and reported by the
screening diagnostics and could not enter chi-square.

The reference-frame cases are synthetic on purpose, and for the reason the
sibling file states: a frame has an exact answer only when it was *made* from
the structure by a transform we chose. Here that is pushed one step further --
the fixed point is placed on a known atom of the structure and then written in
the measured frame, so the transported coordinate has to come back to that
atom, and the chi-square has to be the one that separation gives.
"""

import json
import math
from pathlib import Path

import pytest

import IMP
import IMP.algebra
import IMP.atom
import IMP.core
import IMP.bff

REPO = Path(__file__).resolve().parents[2]
T4L = REPO / "examples" / "structure" / "T4L"
PDB = str(T4L / "3GUN.pdb")
FPS_JSON = T4L / "fret.fps.json"

#: A rigid transform with no special structure -- an off-axis rotation and a
#: translation of tens of angstrom, so an identity fit cannot pass by accident.
ROTATION = IMP.algebra.get_rotation_about_axis(
    IMP.algebra.Vector3D(0.3, -0.5, 0.81), 0.7)
TRANSFORM = IMP.algebra.Transformation3D(
    ROTATION, IMP.algebra.Vector3D(12.0, -3.0, 5.5))


@pytest.fixture(scope="module")
def structure():
    IMP.set_log_level(IMP.SILENT)
    m = IMP.Model()
    h = IMP.atom.read_pdb(PDB, m, IMP.atom.NonWaterNonHydrogenPDBSelector())
    return m, h


def _ca(hierarchy, residue, chain="A"):
    ps = IMP.atom.Selection(hierarchy, chain_id=chain, residue_index=residue,
                            atom_type=IMP.atom.AT_CA).get_selected_particles()
    return IMP.core.XYZ(ps[0]).get_coordinates() if ps else None


def _frame_json(hierarchy, residues, transform, chain="A"):
    """`reference_atoms` for a frame that is the structure moved by `transform`."""
    rows = []
    for r in residues:
        c = _ca(hierarchy, r, chain)
        if c is None:
            continue
        v = transform * c
        rows.append({"chain_identifier": chain, "residue_seq_number": r,
                     "atom_name": "CA", "x": v[0], "y": v[1], "z": v[2]})
    return rows


def _one_av_name(payload, distance="19-119_C2"):
    return payload["Distances"][distance]["position1_name"]


def _write(tmp_path, payload, name="points.fps.json"):
    path = tmp_path / name
    path.write_text(json.dumps(payload))
    return str(path)


def _fixed_on_atom(tmp_path, structure, residue=60, distance_exp=30.0,
                   error_neg=2.0, error_pos=3.0, distance_type="RDAMean",
                   name="points.fps.json"):
    """An fps.json whose XYZ position lands, after its fit, on a known atom.

    Returns (path, av position name, the atom's coordinate on the structure).
    """
    _, h = structure
    payload = json.loads(FPS_JSON.read_text())
    av_name = _one_av_name(payload)
    target = _ca(h, residue)
    measured = TRANSFORM * target
    out = {
        "Positions": {
            av_name: payload["Positions"][av_name],
            "fixed": {"simulation_type": "XYZ", "body_id": 0,
                      "x": measured[0], "y": measured[1], "z": measured[2],
                      "reference_atoms": _frame_json(h, range(30, 45),
                                                     TRANSFORM)},
        },
        "Distances": {
            "dd": {"position1_name": av_name, "position2_name": "fixed",
                   "distance": distance_exp, "error_neg": error_neg,
                   "error_pos": error_pos, "distance_type": distance_type,
                   "Forster_radius": 52.0},
        },
        "χ²": {"s": {"distances": ["dd"]}},
    }
    return _write(tmp_path, out, name), av_name, target


# --------------------------------------------------------------------------
# XYZ: the transported coordinate, and the chi-square it produces
# --------------------------------------------------------------------------
def test_an_xyz_position_lands_where_its_frame_says(tmp_path, structure):
    """The fit must undo the transform the frame was written with, exactly.

    The point was chosen to be residue 60's CA *on this structure*, so the
    answer is an atom we can look up rather than a plausible coordinate.
    """
    path, _, target = _fixed_on_atom(tmp_path, structure)
    assembly = IMP.bff.create_docking_assembly([PDB], path, "s")
    network = assembly.get_network()

    assert list(network.get_point_position_names()) == ["fixed"]
    assert list(network.get_atom_position_names()) == []   # XYZ, not ATOM
    got = IMP.core.XYZ(assembly.get_model(),
                       network.get_position_particle_index("fixed")
                       ).get_coordinates()
    # 2.2e-14 when measured; the tolerance is arithmetic noise through an SVD
    # on coordinates of tens of angstrom, not a modelling allowance.
    assert (got - target).get_magnitude() < 1e-9


def test_a_distance_touching_a_fixed_position_enters_chi2(tmp_path, structure):
    """End to end, against a chi-square computed here from the geometry.

    The model distance must be |mean position of the volume - the transported
    point|, and the chi-square must be that separation against the measurement
    with the error bar its sign selects. Nothing here reads the answer off the
    thing being tested.
    """
    path, av_name, target = _fixed_on_atom(tmp_path, structure,
                                           distance_exp=30.0, error_neg=2.0,
                                           error_pos=3.0)
    result = IMP.bff.score_structures([PDB], path, "s")
    assert result.n_distances == 1
    pair = result.pairs[0]

    # the two ends, computed independently. `evaluate()` first: a network with
    # a shared occupancy registry has to be evaluated once before an AV handle
    # can be resampled on its own (the registry's coordinate snapshot is taken
    # there), and reaching into one before that segfaults -- see the report.
    assembly = IMP.bff.create_docking_assembly([PDB], path, "s")
    assembly.evaluate()
    av = assembly.get_network().get_used_av(av_name)
    av.resample()
    expected_r = (av.get_mean_position() - target).get_magnitude()

    assert pair.distance_model == pytest.approx(expected_r, rel=1e-9)
    # the model is short of 30 A, so the *negative* error bar applies
    assert pair.residual < 0.0
    expected_chi2 = ((expected_r - 30.0) / 2.0) ** 2
    assert pair.chi2 == pytest.approx(expected_chi2, rel=1e-9)
    # measured 2026-08-31: R_mp 26.7921 A, chi2 2.5727 -- with the accessible
    # contact volume the position asks for still ignored. Re-measured
    # 2026-09-01 once it was honoured (PRD-121 G9): the cloud is weighted
    # towards the surface, its mean moves, and R_mp shortens to 22.8273 A.
    # Olga's radii were the default for part of that day and gave 13.2917;
    # IMP's own are the default again (AV::set_radii_source) and 12.8620
    # returns exactly.
    assert pair.chi2 == pytest.approx(12.861951, abs=1e-3)
    # the score is half a chi-square per distance plus a clash term, and one
    # rigid body has no clashes
    assert result.score == pytest.approx(0.5 * expected_chi2, rel=1e-6)


def test_the_error_bar_is_still_chosen_by_the_sign(tmp_path, structure):
    """A point end changes the distance convention, not the chi-square."""
    short, _, _ = _fixed_on_atom(tmp_path, structure, distance_exp=40.0,
                                 error_neg=2.0, error_pos=8.0, name="a.json")
    long_, _, _ = _fixed_on_atom(tmp_path, structure, distance_exp=20.0,
                                 error_neg=2.0, error_pos=8.0, name="b.json")
    # bound, not chained: `score_structures(...).pairs` reads a member of a
    # temporary and comes back empty (see the note in pyext/IMP_bff.docking.i)
    ra = IMP.bff.score_structures([PDB], short, "s")
    rb = IMP.bff.score_structures([PDB], long_, "s")
    a, b = ra.pairs[0], rb.pairs[0]
    assert a.residual < 0 and b.residual > 0
    assert a.chi2 == pytest.approx((a.residual / 2.0) ** 2, rel=1e-9)
    assert b.chi2 == pytest.approx((b.residual / 8.0) ** 2, rel=1e-9)


def test_a_point_end_forces_rmp_whatever_the_file_asks_for(tmp_path, structure):
    """`FilterEngine.cs:305-307`. The three conventions must give one number.

    A cloud-to-cloud <R_DA> is several angstrom from R_mp, so if the rule were
    not applied the three files below would disagree -- which is exactly what
    makes this worth asserting rather than assuming.
    """
    models = {}
    for i, dtype in enumerate(("RDAMean", "RDAMeanE", "Rmp")):
        path, _, _ = _fixed_on_atom(tmp_path, structure, distance_type=dtype,
                                    name="t%d.json" % i)
        result = IMP.bff.score_structures([PDB], path, "s")
        pair = result.pairs[0]
        models[dtype] = pair.distance_model
        # and the table says what was scored, not what was asked for
        assert pair.distance_type == "Rmp"
    assert models["RDAMean"] == pytest.approx(models["RDAMeanE"], rel=1e-12)
    assert models["RDAMean"] == pytest.approx(models["Rmp"], rel=1e-12)

    # the same two positions with two *volumes* do not agree like that -- so
    # the equality above is the rule, not a coincidence of this geometry
    payload = json.loads(FPS_JSON.read_text())
    name = "19-119_C2"
    d = dict(payload["Distances"][name])
    keep = {d["position1_name"], d["position2_name"]}
    both = {}
    for i, dtype in enumerate(("RDAMean", "Rmp")):
        out = {"Positions": {k: payload["Positions"][k] for k in keep},
               "Distances": {name: dict(d, distance_type=dtype)},
               "χ²": {"s": {"distances": [name]}}}
        p = _write(tmp_path, out, "vol%d.json" % i)
        r = IMP.bff.score_structures([PDB], p, "s")
        both[dtype] = r.pairs[0].distance_model
    assert abs(both["RDAMean"] - both["Rmp"]) > 1.0


def test_an_xyz_position_without_a_frame_is_taken_as_written(tmp_path,
                                                             structure):
    """No `reference_atoms` means the coordinate is used as it stands.

    That is right when the file was written against this very structure and
    silently wrong otherwise, which is why the schema documents the frame on
    exactly this position type. What is asserted is that nothing is invented.
    """
    payload = json.loads(FPS_JSON.read_text())
    av_name = _one_av_name(payload)
    here = IMP.algebra.Vector3D(1.0, 2.0, 3.0)
    out = {"Positions": {av_name: payload["Positions"][av_name],
                         "fixed": {"simulation_type": "XYZ", "body_id": 0,
                                   "x": 1.0, "y": 2.0, "z": 3.0}},
           "Distances": {"dd": {"position1_name": av_name,
                                "position2_name": "fixed", "distance": 30.0,
                                "error_neg": 2.0, "error_pos": 2.0,
                                "distance_type": "Rmp",
                                "Forster_radius": 52.0}},
           "χ²": {"s": {"distances": ["dd"]}}}
    path = _write(tmp_path, out, "bare.json")
    assembly = IMP.bff.create_docking_assembly([PDB], path, "s")
    got = IMP.core.XYZ(assembly.get_model(),
                       assembly.get_network().get_position_particle_index(
                           "fixed")).get_coordinates()
    assert (got - here).get_magnitude() == 0.0


def test_screening_a_library_with_a_fixed_position_now_runs(tmp_path,
                                                            structure):
    """This is the failure the gap was named after.

    Every structure came back NaN because the network could not be built at
    all, so a screen containing one fixed position looked like a library of
    unscorable models. The `ref_rmsd` column and the score now come from the
    same fit.
    """
    path, _, _ = _fixed_on_atom(tmp_path, structure)
    ranked = list(IMP.bff.screen_structures([PDB], path, "s"))
    assert len(ranked) == 1
    assert math.isfinite(ranked[0].score)
    assert ranked[0].invalid_r == 0
    # the frame fits this structure exactly, because it was made from it
    assert ranked[0].ref_rmsd == pytest.approx(0.0, abs=1e-6)


# --------------------------------------------------------------------------
# ATOM: a point that is an atom of the structure
# --------------------------------------------------------------------------
def test_an_atom_position_is_the_atom_itself(tmp_path, structure):
    """No particle is created and none has to be moved: it *is* the atom."""
    _, h = structure
    payload = json.loads(FPS_JSON.read_text())
    av_name = _one_av_name(payload)
    out = {"Positions": {av_name: payload["Positions"][av_name],
                         "anchor": {"simulation_type": "ATOM", "body_id": 0,
                                    "chain_identifier": "A",
                                    "residue_seq_number": 60,
                                    "atom_name": "CA"}},
           "Distances": {"dd": {"position1_name": av_name,
                                "position2_name": "anchor", "distance": 30.0,
                                "error_neg": 2.0, "error_pos": 2.0,
                                "distance_type": "RDAMean",
                                "Forster_radius": 52.0}},
           "χ²": {"s": {"distances": ["dd"]}}}
    path = _write(tmp_path, out, "atom.json")
    assembly = IMP.bff.create_docking_assembly([PDB], path, "s")
    network = assembly.get_network()
    assert list(network.get_atom_position_names()) == ["anchor"]
    pi = network.get_position_particle_index("anchor")
    got = IMP.core.XYZ(assembly.get_model(), pi).get_coordinates()
    # not exactly zero: the assembly puts its atoms into a rigid-body frame
    # and back, which is a round trip through a quaternion
    assert (got - _ca(h, 60)).get_magnitude() < 1e-9
    # and it is an atom of the assembly, not a copy of one
    assert IMP.atom.Atom.get_is_setup(assembly.get_model(), pi)


def test_a_network_of_two_points_has_no_volumes_at_all(tmp_path, structure):
    """The degenerate case the evaluation pipeline has to survive.

    With both ends points there is no accessible volume to raster, search or
    coarsen, so the threaded evaluation has no AV tasks -- and a pair task that
    still claimed an AV slot would index an empty array.
    """
    _, h = structure
    a, b = _ca(h, 20), _ca(h, 60)
    out = {"Positions": {"p1": {"simulation_type": "ATOM", "body_id": 0,
                                "chain_identifier": "A",
                                "residue_seq_number": 20, "atom_name": "CA"},
                         "p2": {"simulation_type": "ATOM", "body_id": 0,
                                "chain_identifier": "A",
                                "residue_seq_number": 60, "atom_name": "CA"}},
           "Distances": {"dd": {"position1_name": "p1", "position2_name": "p2",
                                "distance": 10.0, "error_neg": 1.0,
                                "error_pos": 1.0, "distance_type": "RDAMean",
                                "Forster_radius": 52.0}},
           "χ²": {"s": {"distances": ["dd"]}}}
    path = _write(tmp_path, out, "two_points.json")
    result = IMP.bff.score_structures([PDB], path, "s")
    assert result.n_avs == 0
    assert result.n_distances == 1
    assert result.pairs[0].distance_model == pytest.approx(
        (a - b).get_magnitude(), rel=1e-12)


# --------------------------------------------------------------------------
# the schema and the legacy reader
# --------------------------------------------------------------------------
def test_the_schema_knows_the_atom_type():
    assert "ATOM" in list(IMP.bff.fps_simulation_types())
    assert "ATOM" in list(IMP.bff.fps_av_simulation_types())
    good = {"simulation_type": "ATOM", "chain_identifier": "A",
            "residue_seq_number": 60, "atom_name": "CA"}
    report = IMP.bff.validate_position(json.dumps(good), "p")
    assert not report.errors, report.errors
    assert not report.warnings, report.warnings

    # it names an atom rather than carrying a coordinate, so it needs a residue
    bad = {"simulation_type": "ATOM", "atom_name": "CA"}
    assert any("residue_seq_number" in e
               for e in IMP.bff.validate_position(json.dumps(bad), "p").errors)
    # and an AV parameter on it is a warning: there is no volume to apply it to
    odd = dict(good, linker_length=20.0)
    assert any("no volume" in w
               for w in IMP.bff.validate_position(json.dumps(odd),
                                                  "p").warnings)


def test_the_legacy_lps_reader_keeps_atom_lines(tmp_path):
    """`Name Molecule D/A ATOM <serial>` was being dropped as an unknown line.

    A legacy file's crosslinks disappeared silently, which is worse than
    refusing the file: the run went on with fewer restraints than it was given.
    """
    lps = tmp_path / "LabelingPositions.txt"
    lps.write_text("site1\tmol\tD\tAV1\t20.0\t4.5\t3.5\t100\n"
                   "anchor1\tmol\t-\tATOM\t250\n"
                   "anchor2\tmol\t-\tATOM\t260\n")
    document = IMP.bff.read_old_lps_txt(str(lps), [])
    positions = json.loads(document.positions)
    assert set(positions) == {"site1", "anchor1", "anchor2"}
    assert positions["anchor1"]["simulation_type"] == "ATOM"
    # no PDB was given, so the serial is carried as the residue key -- and no
    # zeroed AV parameters are invented, which would validate as an AV1 of
    # radius 0
    assert positions["anchor1"]["residue_seq_number"] == 250
    assert "radius1" not in positions["anchor1"]
    assert positions["site1"]["simulation_type"] == "AV1"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
