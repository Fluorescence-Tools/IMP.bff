"""FPS's export files, byte by byte (PRD-121 G5).

`okf/references/fps-export-formats.md` specifies what `SaveForm` writes: a
PyMOL script per result, an Overlay, an OverlayStates, an R table, a chi2 table.
This checks the bytes -- a format test that only asserts a file appeared is not
a format test.

Where this port diverges from FPS it is on purpose, and each divergence is
asserted here so nobody has to diff two implementations to find out:

* **all N transforms are written**, not `1..N-1`. FPS skips molecule 0 because
  its own engine normalises every result so molecule 0 is the identity; a pose
  from anywhere else then lands wrong.
* **`camera=0`** on every `rotate`/`translate`. PyMOL reads both in *camera*
  space by default and FPS writes model-frame vectors, which works only in a
  session whose view has not been touched.
* **`_tmp.pdb` is absolute and removed.** FPS writes the bare name, so it lands
  in PyMOL's working directory and stays there.
* **the chi2 table carries `chi2` and `chi2_r`.** FPS's file writes the raw sum
  under a header that reads as a reduced one while its GUI shows the reduced
  one, so the file and the screen disagree about what `chi2` means.
* **the Filter R table's `Number` is 1-based.** FPS prints the 0-based array
  index there and the 1-based `InternalNumber` everywhere else. Join on `File`.
"""

import json
import math
from pathlib import Path

import pytest

import IMP
import IMP.algebra
import IMP.bff

REPO = Path(__file__).resolve().parents[2]
HIV = REPO / "examples" / "structure" / "HIV_RT"
PDBS = [str(HIV / "protein_1R0A.pdb"), str(HIV / "dna.pdb")]
FPS_JSON = str(HIV / "hiv_rt.fps.json")


@pytest.fixture(scope="module")
def molecules():
    return IMP.bff.fps_molecules(PDBS)


@pytest.fixture(scope="module")
def assembly():
    return IMP.bff.build_docking_assembly(PDBS, FPS_JSON)


def _pair(name, p1, p2, exp, model, err=4.0, dtype="RDAMeanE"):
    pair = IMP.bff.PairDistance()
    pair.name = name
    pair.position1 = p1
    pair.position2 = p2
    pair.distance_exp = exp
    pair.distance_model = model
    pair.error_neg = err
    pair.error_pos = err
    pair.distance_type = dtype
    return pair


def _row(number, poses, chi2, method="Docking", parent=0, converged=True):
    row = IMP.bff.FPSResultRow()
    row.number = number
    row.chi2 = chi2
    row.chi2_bond = 0.5
    row.chi2_clash = 1.25
    row.converged = converged
    row.parent = parent
    row.method = method
    row.poses = poses
    row.pairs = [_pair("d1", "A1", "B1", 50.0, 51.5),
                 _pair("d2", "A2", "B2", 60.0, 58.25)]
    return row


@pytest.fixture(scope="module")
def table(assembly):
    """Two rows: the input pose, and the same with body 1 pushed 10 A along x.

    Synthetic on purpose. The geometry is exactly known, so a wrong number in a
    file is a wrong number and not a converged-to-something-else minimisation.
    """
    base = IMP.bff.capture_poses(assembly)
    moved = json.loads(base)
    moved[1]["t"][0] += 10.0
    t = IMP.bff.FPSResultTable()
    t.n_molecules = 2
    t.reference_number = 1
    t.distance_type = "RDAMeanE"
    t.rows = [_row(1, base, 12.3456789),
              _row(2, json.dumps(moved), 20.5, "ErrorEstimation", 1, False)]
    return t


def _lines(path):
    text = Path(path).read_text()
    assert "\r" not in text, "line endings are LF; FPS's CRLF is Windows only"
    return text.rstrip("\n").split("\n")


# --------------------------------------------------------------------------
# The molecules an export names
# --------------------------------------------------------------------------

def test_a_molecule_is_its_path_its_stem_and_its_frame(molecules):
    names = [m.name for m in molecules]
    assert names == ["protein_1R0A", "dna"]
    assert [m.path for m in molecules] == PDBS
    assert [m.n_atoms for m in molecules] == [8016, 1018]


def test_the_input_frame_is_not_the_identity(molecules, assembly):
    """The trap that makes a naive export silently wrong.

    `IMP.atom.create_rigid_body` starts a body in its principal-axis frame, so
    the pose `capture_poses` records is *not* the transform to apply to the
    coordinates a `load` line reads. A script built from the pose alone rotates
    every structure by that principal-axis rotation -- and still looks like an
    overlay, because every structure is rotated by the same amount.
    """
    poses = json.loads(IMP.bff.capture_poses(assembly))
    angle = IMP.algebra.get_axis_and_angle(
        IMP.algebra.Rotation3D(*poses[0]["q"]))[1]
    assert angle > 0.5, "expected a non-identity principal-axis frame"
    # The centre and the frame's translation are the same point, so the
    # `origin=` of a `rotate` line and the frame the pose was taken in agree.
    assert list(molecules[0].center) == pytest.approx(list(poses[0]["t"]))


def test_the_untransformed_pose_writes_a_zero_transform(tmp_path, table,
                                                        molecules):
    """Row 1 is the input pose, so its rotate/translate must be a no-op."""
    path = IMP.bff.write_fps_pymol_script(str(tmp_path), table.rows[0],
                                          molecules)
    rotates = [l for l in _lines(path) if l.startswith("rotate ")]
    translates = [l for l in _lines(path) if l.startswith("translate ")]
    assert len(rotates) == 2 and len(translates) == 2
    for line in rotates:
        assert ", 0.000, " in line
    for line in translates:
        assert line.startswith("translate [0.000, 0.000, 0.000], ")


# --------------------------------------------------------------------------
# The per-result PyMOL script
# --------------------------------------------------------------------------

def test_the_pymol_script_has_fpss_shape(tmp_path, table, molecules):
    path = IMP.bff.write_fps_pymol_script(str(tmp_path), table.rows[1],
                                          molecules)
    assert Path(path).name == "structure2.pml"
    lines = _lines(path)

    # 1. the energy comment, F8, with FPS's exact parenthesis
    assert lines[0] == "# Energy = 20.50000000 (not converged)"
    # 2. one `load` per molecule, in order, unquoted
    assert lines[1] == "load " + PDBS[0]
    assert lines[2] == "load " + PDBS[1]
    # 3. every molecule's transform, molecule 0 included
    assert lines[3].startswith("rotate [")
    assert lines[3].endswith(", protein_1R0A, origin=[%s], camera=0"
                             % ", ".join("%.3f" % c
                                         for c in molecules[0].center))
    assert lines[4] == "translate [0.000, 0.000, 0.000], protein_1R0A, camera=0"
    assert lines[5].startswith("rotate [")
    assert ", dna, origin=[" in lines[5]
    # body 1 was pushed 10 A along x and not turned, so this is exact
    assert lines[6] == "translate [10.000, 0.000, 0.000], dna, camera=0"
    assert lines[-1] == "deselect"


def test_a_vector_is_three_decimals_comma_space(tmp_path, table, molecules):
    """`Vector3.ToString()` is `{0:F3}, {1:F3}, {2:F3}` (`MatrixVector3.cs:62`)."""
    path = IMP.bff.write_fps_pymol_script(str(tmp_path), table.rows[1],
                                          molecules)
    for line in _lines(path):
        if not line.startswith(("rotate ", "translate ")):
            continue
        inside = line[line.index("[") + 1:line.index("]")]
        parts = inside.split(", ")
        assert len(parts) == 3
        for part in parts:
            assert part.count(".") == 1 and len(part.split(".")[1]) == 3


def test_camera_zero_can_be_turned_off_for_fps_compatibility(tmp_path, table,
                                                             molecules):
    options = IMP.bff.FPSExportOptions()
    options.camera_zero = False
    path = IMP.bff.write_fps_pymol_script(str(tmp_path), table.rows[1],
                                          molecules, options)
    assert "camera=0" not in Path(path).read_text()


def test_save_pdb_writes_an_absolute_path_after_the_transforms(tmp_path, table,
                                                               molecules):
    options = IMP.bff.FPSExportOptions()
    options.save_pdb = True
    path = IMP.bff.write_fps_pymol_script(str(tmp_path), table.rows[1],
                                          molecules, options)
    lines = _lines(path)
    save = [i for i, l in enumerate(lines) if l.startswith("save ")]
    assert len(save) == 1
    assert lines[save[0]] == "save %s/structure2.pdb, sele" % tmp_path
    assert lines[save[0] - 1] == "select all"
    # after every transform, so the PDB it writes is the assembled complex
    assert save[0] > max(i for i, l in enumerate(lines)
                         if l.startswith("translate "))


def test_labelling_positions_become_pseudoatoms(tmp_path, table, molecules,
                                                assembly):
    labels = IMP.bff.fps_label_positions(assembly)
    assert len(labels) == 11        # HIV-RT's eleven sites
    assert sorted(l.name for l in labels)[0].startswith(("p51", "p66", "p_"))
    options = IMP.bff.FPSExportOptions()
    options.labels = labels
    path = IMP.bff.write_fps_pymol_script(str(tmp_path), table.rows[0],
                                          molecules, options)
    text = Path(path).read_text()
    for label in labels:
        assert "pseudoatom %s, pos=[" % label.name in text
        assert 'label %s, "%s"' % (label.name, label.name) in text
        assert "show spheres, %s" % label.name in text
    # Row 1 is the input pose, so a pseudoatom must land on the position
    # itself: the same `R*(lp - CM) + CM + T` the molecules got.
    first = labels[0]
    expected = "pseudoatom %s, pos=[%s]" % (
        first.name, ", ".join("%.3f" % c for c in first.position))
    assert expected in text


def test_the_best_fit_of_a_row_onto_itself_is_the_identity(tmp_path, table,
                                                           molecules,
                                                           assembly):
    """FPS's biggest export trap, asserted as the absence of a 180 degrees.

    FPS reads a `BestFitRotation` no `Save*` method computes; with the best-fit
    box unticked it is the default all-zero matrix, whose axis-angle is
    `atan2(0, -1)` = pi about z. Here the fit is computed against the stated
    reference, and row 1 *is* the reference.
    """
    fitted = IMP.bff.add_best_fit(table, assembly)
    options = IMP.bff.FPSExportOptions()
    options.best_fit = True
    path = IMP.bff.write_fps_pymol_script(str(tmp_path), fitted.rows[0],
                                          molecules, options)
    lines = _lines(path)
    block = lines[lines.index("select all") + 1:]
    assert block[0] == "translate [0.000, 0.000, 0.000], sele, camera=0"
    assert block[1].endswith(", sele, origin=[0, 0, 0], camera=0")
    assert ", 0.000, sele, origin=[0, 0, 0]" in block[1]
    # translate first, then rotate about the origin -- FPS's ordering
    assert block[0].startswith("translate ")
    assert block[1].startswith("rotate ")


# --------------------------------------------------------------------------
# Overlay and OverlayStates
# --------------------------------------------------------------------------

def test_overlay_is_one_object_per_result(tmp_path, table, molecules):
    path = IMP.bff.write_fps_overlay(str(tmp_path), table, molecules)
    assert Path(path).name == "structureOverlay.pml"
    lines = _lines(path)
    assert lines[0] == "# Overlay of structures 1, 2"
    assert lines[1] == "load " + PDBS[0]
    assert lines.count("select _tmp*") == 2
    assert lines.count("create _tmpjoin, sele") == 2
    assert lines.count("delete _tmp*") == 2
    assert "copy structure1, _tmpjoin" in lines
    assert "copy structure2, _tmpjoin" in lines
    # a `copy _tmp<i>` per molecule per result
    assert lines.count("copy _tmp0, protein_1R0A") == 2
    assert lines.count("copy _tmp1, dna") == 2
    # the pristine originals are removed at the end
    assert lines[-2:] == ["delete protein_1R0A", "delete dna"]
    # `select _tmp*`, not FPS's `select object _tmp*`: `object` is not a
    # selection operator in PyMOL.
    assert "select object _tmp*" not in lines


def test_overlay_states_is_one_object_with_named_states(tmp_path, table,
                                                        molecules):
    path = IMP.bff.write_fps_overlay_states(str(tmp_path), table, molecules)
    assert Path(path).name == "structureOverlayStates.pml"
    lines = _lines(path)
    assert lines[0] == "# Overlay of structures 1, 2"
    tmp_pdb = "%s/_tmp.pdb" % tmp_path
    # The scratch file is absolute -- FPS writes the bare `_tmp.pdb`, which
    # lands in PyMOL's working directory and is left behind.
    assert lines.count("save %s, sele" % tmp_pdb) == 2
    assert lines.count("load %s, structure" % tmp_pdb) == 2
    assert "system rm -f %s" % tmp_pdb in lines
    # and each state keeps the InternalNumber FPS drops
    assert 'set_title structure, 1, "structure1"' in lines
    assert 'set_title structure, 2, "structure2"' in lines
    assert "create _tmpjoin, sele" not in lines


def test_the_two_overlays_share_their_head_and_their_copies(tmp_path, table,
                                                            molecules):
    """`SaveForm.cs:314-332` is byte-identical to `:263-287`; so is this."""
    a = _lines(IMP.bff.write_fps_overlay(str(tmp_path / "a"), table, molecules))
    b = _lines(IMP.bff.write_fps_overlay_states(str(tmp_path / "b"), table,
                                                molecules))
    head = 1 + len(molecules) + 3 * len(molecules)
    assert a[:head] == b[:head]


# --------------------------------------------------------------------------
# The tables
# --------------------------------------------------------------------------

def test_the_r_table_is_tab_separated_with_fpss_column_names(tmp_path, table):
    path = IMP.bff.write_fps_r_table(str(tmp_path), table)
    assert Path(path).name == "structureRtable_RDAMeanE.txt"
    lines = _lines(path)
    assert lines[0].split("\t") == ["Structure", "Number", "A1_B1", "A2_B2"]
    assert lines[1].split("\t") == ["structure1", "1", "51.500", "58.250"]
    assert lines[2].split("\t") == ["structure2", "2", "51.500", "58.250"]


def test_a_missing_model_distance_is_written_as_nan(tmp_path, assembly):
    """`Double.NaN.ToString("F3")` is the literal `NaN`; a blank would hide it."""
    t = IMP.bff.FPSResultTable()
    row = _row(1, IMP.bff.capture_poses(assembly), 1.0)
    row.pairs = [_pair("d1", "A1", "B1", 50.0, float("nan"))]
    t.rows = [row]
    path = IMP.bff.write_fps_r_table(str(tmp_path), t)
    assert _lines(path)[1].split("\t")[2] == "NaN"


def test_the_chi2_table_names_both_chi_squares(tmp_path, table, assembly):
    filled = IMP.bff.add_rmsd_columns(table, assembly)
    path = IMP.bff.write_fps_chi2_table(str(tmp_path), filled)
    assert Path(path).name == "structurechi2table.txt"
    lines = _lines(path)
    assert lines[0].split("\t") == [
        "File", "chi2", "chi2_r", "chi2_bond", "chi2_clash", "Converged",
        "Method", "Parent", "RMSD vs previous", "RMSD vs 1"]

    first = lines[1].split("\t")
    assert first[0] == "structure1"
    assert first[1] == "12.34567890"          # F8, raw, not reduced
    assert first[3] == "0.5000" and first[4] == "1.2500"   # F4
    assert first[5] == "True"
    assert first[6] == "Docking"
    assert first[7] == "0"                    # no parent
    assert first[8] == "---"                  # FPS's literal for row 1
    assert first[9] == "0.000000"             # the reference is itself

    second = lines[2].split("\t")
    assert second[0] == "structure2"
    assert second[5] == "False"
    assert second[6] == "ErrorEstimation"
    assert second[7] == "1"
    # body 1 moved 10 A: sqrt(1018/9034) * 10
    expected = 10.0 * math.sqrt(1018.0 / 9034.0)
    assert float(second[8]) == pytest.approx(expected, abs=5e-7)
    assert float(second[9]) == pytest.approx(expected, abs=5e-7)


def test_chi2_r_is_chi2_over_fpss_own_dof(tmp_path, table):
    """`dof = max(Ndist - 6(Nmol - 1), 1)` (`MainForm.cs:683`).

    Two distances and two molecules give `2 - 6 = -4`, floored at 1 -- so the
    two columns coincide here, and that is the point: FPS's dof can go
    non-positive and its floor is what stops the division.
    """
    assert table.get_dof() == 1
    path = IMP.bff.write_fps_chi2_table(str(tmp_path), table)
    row = _lines(path)[1].split("\t")
    assert row[1] == row[2]

    wide = IMP.bff.FPSResultTable()
    wide.n_molecules = 1
    r = IMP.bff.FPSResultRow()
    r.number = 1
    r.chi2 = 20.0
    r.pairs = [_pair("d%d" % i, "A%d" % i, "B%d" % i, 50.0, 50.0)
               for i in range(10)]
    wide.rows = [r]
    assert wide.get_dof() == 10
    path = IMP.bff.write_fps_chi2_table(str(tmp_path / "wide"), wide)
    row = _lines(path)[1].split("\t")
    assert row[1] == "20.00000000" and row[2] == "2.00000000"


def test_e_bond_is_a_subset_of_chi2_not_an_addition(table):
    """`SpringEngine.cs:445` accumulates it inside the total's own loop."""
    row = table.rows[0]
    assert row.chi2_bond <= row.chi2
    payload = json.loads(table.get_json())
    assert payload["rows"][0]["chi2_bond"] == row.chi2_bond
    assert payload["rows"][0]["chi2"] == row.chi2


def test_write_fps_exports_writes_all_five(tmp_path, table, molecules):
    written = list(IMP.bff.write_fps_exports(str(tmp_path), table, molecules))
    assert [Path(p).name for p in written] == [
        "structure1.pml", "structure2.pml", "structureOverlay.pml",
        "structureOverlayStates.pml", "structureRtable_RDAMeanE.txt",
        "structurechi2table.txt"]
    for path in written:
        assert Path(path).stat().st_size > 0


def test_the_table_round_trips_through_json(table):
    payload = json.loads(table.get_json())
    assert payload["n_molecules"] == 2
    assert payload["reference_number"] == 1
    assert payload["distance_type"] == "RDAMeanE"
    assert [r["number"] for r in payload["rows"]] == [1, 2]
    assert payload["rows"][1]["method"] == "ErrorEstimation"
    assert payload["rows"][1]["parent"] == 1
    assert payload["rows"][1]["converged"] is False
    # the poses come back as JSON, not as a string of JSON
    assert isinstance(payload["rows"][0]["poses"], list)
    assert len(payload["rows"][0]["pairs"]) == 2


# --------------------------------------------------------------------------
# Filter mode
# --------------------------------------------------------------------------

def _screened(path, chi2_r, invalid=0, ref=0.0, sigmas=(0, 0, 0)):
    s = IMP.bff.ScreenedStructure()
    s.path = path
    s.chi2_r = chi2_r
    s.invalid_r = invalid
    s.ref_rmsd = ref
    s.sigma1, s.sigma2, s.sigma3 = sigmas
    s.pairs = [_pair("d1", "A1", "B1", 50.0, 49.0)]
    return s


def test_the_filter_chi2_table_is_fpss_seven_columns(tmp_path):
    structures = [_screened("/lib/a.pdb", 1.234, 2, 0.5, (3, 1, 0)),
                  _screened("/lib/b.pdb", 5.678)]
    options = IMP.bff.FPSExportOptions()
    options.prefix = "screening_"
    path = IMP.bff.write_fps_screening_chi2_table(str(tmp_path), structures,
                                                 options)
    assert Path(path).name == "screening_chi2table.txt"
    lines = _lines(path)
    assert lines[0].split("\t") == [
        "File", "Chi2r", "NaNs", "RefRMSD", ">1sigma", ">2sigma", ">3sigma"]
    assert lines[1].split("\t") == ["a.pdb", "1.234", "2", "0.500", "3", "1",
                                    "0"]
    assert lines[2].split("\t") == ["b.pdb", "5.678", "0", "0.000", "0", "0",
                                    "0"]


def test_the_filter_r_tables_number_is_one_based_by_default(tmp_path):
    structures = [_screened("/lib/a.pdb", 1.0), _screened("/lib/b.pdb", 2.0)]
    options = IMP.bff.FPSExportOptions()
    options.prefix = "screening_"
    path = IMP.bff.write_fps_screening_r_table(str(tmp_path), structures,
                                              options)
    lines = _lines(path)
    assert lines[0].split("\t") == ["File", "Number", "A1_B1"]
    assert [l.split("\t")[:2] for l in lines[1:]] == [["a.pdb", "1"],
                                                      ["b.pdb", "2"]]

    # FPS's own numbering, one less than every other Number column it writes.
    options.fps_filter_number = True
    path = IMP.bff.write_fps_screening_r_table(str(tmp_path / "fps"),
                                              structures, options)
    lines = _lines(path)
    assert [l.split("\t")[:2] for l in lines[1:]] == [["a.pdb", "0"],
                                                      ["b.pdb", "1"]]


def test_file_is_the_only_stable_join_key(tmp_path):
    """The `Number` column of a Filter table numbers *the export*, not the run.

    Export a selection and it renumbers from zero, so joining a Filter R table
    to a Filter chi2 table on `Number` silently pairs the wrong rows.
    """
    everything = [_screened("/lib/%s.pdb" % n, 1.0) for n in ("a", "b", "c")]
    selection = [everything[2]]

    options = IMP.bff.FPSExportOptions()
    a = _lines(IMP.bff.write_fps_screening_r_table(str(tmp_path / "all"),
                                                   everything, options))
    b = _lines(IMP.bff.write_fps_screening_r_table(str(tmp_path / "sel"),
                                                   selection, options))
    assert a[3].split("\t")[:2] == ["c.pdb", "3"]
    assert b[1].split("\t")[:2] == ["c.pdb", "1"]   # same structure, new number
