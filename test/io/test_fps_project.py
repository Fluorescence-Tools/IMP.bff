"""The project document: both doors, field by field (PRD-121 G1).

A project is what makes a run repeatable, so what is checked here is that
nothing is lost on the way through either door -- not that "a file appeared".
Two things get their own tests because they are where the legacy format is
easiest to get wrong:

* a **struct** value (`DockParameters = Fps.FPSParameters` and an indented
  ` Field = value` block), and
* an **array** value (`MoleculesPaths = System.String[]` with its elements on
  the *following* line).

The legacy fixture is written here in FPS's exact export shape rather than
taken from a file, so the format under test is `OptionsManager.cs::Export` and
not this module's own writer -- there is no legacy writer, deliberately.
"""

import importlib.machinery
import importlib.util
import json
import os
import sys
from pathlib import Path

import pytest

import IMP
import IMP.bff

REPO = Path(__file__).resolve().parents[2]
HIV = REPO / "examples" / "structure" / "HIV_RT"
BIN = REPO / "bin"


def _program(name):
    path = BIN / name
    loader = importlib.machinery.SourceFileLoader(name, str(path))
    spec = importlib.util.spec_from_file_location(name, str(path),
                                                  loader=loader)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


#: FPS's export of a two-body project, in the exact shape `Export` writes:
#: a banner, a date, `# Category` headers, `Key = value`, an indented field
#: block for a struct, and an array's elements on the line *after* its key.
#: The values are deliberately *not* FPS's defaults, so a field that is
#: silently dropped shows up as a default rather than hiding behind one.
LEGACY = """# FPS v. 1.1.0.0
# Monday, May 19, 2014 11:22:33 AM

# Data
ProjectFPSMode = Dock
MoleculesPaths = System.String[]
 protein_1R0A.pdb dna.pdb
LabelingPositionsPath = LabelingPositions.txt
DistancesPath = Distances.txt

# Search parameters
DockParameters = Fps.FPSParameters
 ViscosityFactor = 0.9
 TimeStepFactor = 0.8
 MaxIterations = 123456
 MaxForce = 401.5
 ClashTolerance = 1.25
 rkT = 11
 ETolerance = 99
 KTolerance = 0.002
 FTolerance = 0.003
 TTolerance = 0.04
 OptimizeSelected = SelectedThenAll

# Refine parameters
RefineParameters = Fps.FPSParameters
 ViscosityFactor = 0.7
 TimeStepFactor = 0.5
 MaxIterations = 500000
 MaxForce = 10000
 ClashTolerance = 0.5
 rkT = 10
 ETolerance = 100
 KTolerance = 0.0005
 FTolerance = 0.0005
 TTolerance = 0.01
 OptimizeSelected = All

# Conversion function
ConversionParameters = Fps.ConversionParameters
 R0 = 49.5
 PolynomOrder = 4

# AV parameters
AVGlobalParameters = Fps.AVGlobalParameters
 GridSize = 0.25
 MinGridSize = 0.45
 LinkerInitialSphere = 0.55
 LinkSearchNodes = 4
 ESamples = 300000

# Selected distances
SelectedDistances = System.Boolean[]
 True False True True
"""


def _legacy(tmp_path, text=LEGACY, name="MyProject.bin.txt"):
    path = tmp_path / name
    path.write_text(text)
    return path


# --------------------------------------------------------------------------
# door one: FPS's legacy project text export
# --------------------------------------------------------------------------
def test_legacy_text_export_reads_every_field(tmp_path):
    """Every key of `OptionsManager.Export`, including a struct and an array."""
    p = IMP.bff.read_fps_project_txt(str(_legacy(tmp_path)))

    # the banner, kept as provenance
    assert p.source == "FPS v. 1.1.0.0"
    assert p.mode == "Dock"

    # the ARRAY value: elements on the line after the key
    assert list(p.structures) == ["protein_1R0A.pdb", "dna.pdb"]
    assert list(p.get_structure_paths()) == [
        str(tmp_path / "protein_1R0A.pdb"), str(tmp_path / "dna.pdb")]

    # the two scalar paths, stored as written and resolved on demand
    assert p.positions_path == "LabelingPositions.txt"
    assert p.distances_path == "Distances.txt"
    assert p.get_labelling_path() == str(tmp_path / "LabelingPositions.txt")
    assert p.get_distances_path() == str(tmp_path / "Distances.txt")

    # the STRUCT value: all eleven fields of FPS's FPSParameters
    d = p.dock
    assert d.viscosity_factor == 0.9
    assert d.time_step_factor == 0.8
    assert d.max_iterations == 123456
    assert d.max_force == 401.5
    assert d.clash_tolerance == 1.25
    assert d.rkt == 11.0
    assert d.e_tolerance == 99.0
    assert d.k_tolerance == 0.002
    assert d.f_tolerance == 0.003
    assert d.t_tolerance == 0.04
    assert d.optimize_selected == "SelectedThenAll"

    r = p.refine
    assert (r.viscosity_factor, r.time_step_factor) == (0.7, 0.5)
    assert (r.max_iterations, r.max_force, r.clash_tolerance) == (
        500000, 10000.0, 0.5)
    assert (r.k_tolerance, r.f_tolerance, r.t_tolerance) == (
        0.0005, 0.0005, 0.01)
    assert r.optimize_selected == "All"

    # the two other structs
    assert p.conversion.forster_radius == 49.5
    assert p.conversion.polynomial_order == 4
    assert p.av.grid_size == 0.25
    assert p.av.min_grid_size == 0.45
    assert p.av.linker_initial_sphere == 0.55
    assert p.av.link_search_nodes == 4
    assert p.av.e_samples == 300000

    # the boolean array
    assert list(p.selected_flags) == [1, 0, 1, 1]
    assert list(p.resolve_selected_distances(["a", "b", "c", "d"])) == [
        "a", "c", "d"]


def test_legacy_blocks_it_does_not_carry_keep_fps_defaults(tmp_path):
    """The fixture has no Sample/Screening/Error block; those stay FPS's."""
    p = IMP.bff.read_fps_project_txt(str(_legacy(tmp_path)))
    assert p.sample.max_iterations == 8000            # ProjectData.cs:123
    assert p.error_estimation.max_force == 10000.0    # ProjectData.cs:106
    assert p.error_estimation.optimize_selected == "All"
    assert p.screening.clash_tolerance == 1.0         # ProjectData.cs:141


def test_legacy_empty_array_writes_a_blank_line_and_reads_as_empty(tmp_path):
    """`Export` writes the element line even for a zero-length array.

    The blank line that follows an array key is its (empty) element list, not
    the category separator -- so it has to be consumed unconditionally or the
    *next* key is eaten instead.
    """
    text = ("# FPS v. 1.1.0.0\n\n# Data\n"
            "ProjectFPSMode = Filter\n"
            "MoleculesPaths = System.String[]\n"
            "\n"
            "LabelingPositionsPath = LPs.txt\n")
    p = IMP.bff.read_fps_project_txt(str(_legacy(tmp_path, text)))
    assert list(p.structures) == []
    assert p.mode == "Filter"
    assert p.positions_path == "LPs.txt"


def test_legacy_reads_a_comma_decimal_separator(tmp_path):
    """`Double.ToString()` follows the writer's culture; a German FPS wrote
    `0,0005`, and `strtod` reads that as 0 and stops."""
    text = LEGACY.replace(" KTolerance = 0.002", " KTolerance = 0,002")
    p = IMP.bff.read_fps_project_txt(str(_legacy(tmp_path, text)))
    assert p.dock.k_tolerance == 0.002


def test_legacy_unknown_keys_and_categories_are_skipped(tmp_path):
    """Other builds of the GUI carried other options; a project that mostly
    parses is worth more than an exception."""
    text = LEGACY.replace("# Data\n",
                          "# Something else\nWhatever = 12\n\n# Data\n")
    p = IMP.bff.read_fps_project_txt(str(_legacy(tmp_path, text)))
    assert p.mode == "Dock"
    assert list(p.structures) == ["protein_1R0A.pdb", "dna.pdb"]


def test_legacy_windows_paths_are_not_glued_to_the_project_directory(tmp_path):
    """`C:\\...` is absolute. Prefixing a POSIX directory makes it unreadable."""
    text = LEGACY.replace("LabelingPositionsPath = LabelingPositions.txt",
                          "LabelingPositionsPath = C:\\data\\LPs.txt")
    p = IMP.bff.read_fps_project_txt(str(_legacy(tmp_path, text)))
    assert p.get_labelling_path() == "C:\\data\\LPs.txt"


# --------------------------------------------------------------------------
# door two: the extended fps.json
# --------------------------------------------------------------------------
def _assert_same_project(a, b):
    """Field by field, including all five parameter blocks."""
    assert a.mode == b.mode
    assert a.source == b.source
    assert list(a.structures) == list(b.structures)
    assert a.labelling_json == b.labelling_json
    assert a.positions_path == b.positions_path
    assert a.distances_path == b.distances_path
    assert a.score_set == b.score_set
    assert list(a.selected_distances) == list(b.selected_distances)
    assert list(a.selected_flags) == list(b.selected_flags)
    assert a.fixed_body == b.fixed_body
    assert a.ev_weight == b.ev_weight
    assert a.sigma_da == b.sigma_da
    assert a.shuffle == b.shuffle
    assert a.refine_av_cycles == b.refine_av_cycles
    # The clash term's own radii (schema 1.7). All three travel together
    # because `coarse_clash` decides whether the other two can be honoured.
    assert a.coarse_clash == b.coarse_clash
    assert a.clash_radii_source == b.clash_radii_source
    assert a.clash_radii_scale == b.clash_radii_scale
    assert json.loads(a.poses) == json.loads(b.poses)
    assert a.conversion.forster_radius == b.conversion.forster_radius
    assert a.conversion.polynomial_order == b.conversion.polynomial_order
    for name in ("grid_size", "min_grid_size", "linker_initial_sphere",
                 "link_search_nodes", "e_samples"):
        assert getattr(a.av, name) == getattr(b.av, name), name
    for mode in IMP.bff.fps_mode_names():
        x, y = a.get_parameters(mode), b.get_parameters(mode)
        for name in ("viscosity_factor", "time_step_factor", "max_iterations",
                     "max_force", "clash_tolerance", "rkt", "e_tolerance",
                     "k_tolerance", "f_tolerance", "t_tolerance",
                     "optimize_selected", "iterations"):
            assert getattr(x, name) == getattr(y, name), (mode, name)


def test_fps_json_round_trip_is_field_for_field(tmp_path):
    """Write, read, compare every field -- including the five blocks."""
    p = IMP.bff.read_fps_project_txt(str(_legacy(tmp_path)))
    p.labelling_json = "labels.fps.json"
    p.positions_path = ""
    p.distances_path = ""
    p.score_set = "resolved"
    p.selected_distances = ["a", "c", "d"]
    p.fixed_body = 1
    p.ev_weight = 0.25
    p.sigma_da = 4.5
    p.shuffle = 3.0
    p.refine_av_cycles = 2
    p.coarse_clash = False
    p.clash_radii_source = "olga"
    p.clash_radii_scale = 0.9
    p.poses = json.dumps([
        {"body_id": 0, "t": [0.0, 0.0, 0.0], "q": [1.0, 0.0, 0.0, 0.0]},
        {"body_id": 1, "t": [1.5, -2.5, 3.5], "q": [0.5, 0.5, 0.5, 0.5]}])
    block = p.get_parameters("Screening")
    block.optimize_selected = "All"
    block.iterations = 7
    p.set_parameters("Screening", block)

    out = tmp_path / "run.fps.json"
    IMP.bff.write_fps_project(str(out), p)
    back = IMP.bff.read_fps_project(str(out))
    _assert_same_project(p, back)

    # and the file is an fps.json, not a private format
    payload = json.loads(out.read_text())
    assert set(payload) >= {"Positions", "Distances", "Project"}
    assert payload["Project"]["schema_version"] == IMP.bff.fps_schema_version()


def test_fps_json_project_paths_stay_relative(tmp_path):
    """A stored path is written as given, so a project directory can be moved.

    Resolving at read time and writing the resolved path back would pin the
    file to one machine, silently.
    """
    p = IMP.bff.FPSProject()
    p.structures = ["a.pdb", "b.pdb"]
    p.labelling_json = "labels.fps.json"
    out = tmp_path / "run.fps.json"
    IMP.bff.write_fps_project(str(out), p)
    assert json.loads(out.read_text())["Project"]["structures"] == [
        "a.pdb", "b.pdb"]
    back = IMP.bff.read_fps_project(str(out))
    assert list(back.structures) == ["a.pdb", "b.pdb"]
    assert list(back.get_structure_paths()) == [
        str(tmp_path / "a.pdb"), str(tmp_path / "b.pdb")]
    assert back.get_labelling_path() == str(tmp_path / "labels.fps.json")


def test_plain_fps_json_without_a_project_section_is_not_an_error(tmp_path):
    """"Run this fps.json" is a thing to want; the defaults point at it."""
    p = IMP.bff.read_fps_project(str(HIV / "hiv_rt.fps.json"))
    assert p.mode == "None"
    assert list(p.structures) == []
    assert p.get_labelling_path() == str(HIV / "hiv_rt.fps.json")


def test_defaults_are_fps_shipped_per_mode_parameters():
    """ProjectData.cs:61-149 -- the settings behind the published numbers."""
    table = {
        # mode: (visc, dt, max_iter, max_force, clash, kTol, tTol, optimize)
        "Dock": (1.0, 1.0, 200000, 400.0, 1.0, 0.001, 0.02, "Selected"),
        "Refine": (0.7, 0.5, 500000, 10000.0, 0.5, 0.0005, 0.01, "All"),
        "Error estimation": (0.7, 1.0, 100000, 10000.0, 0.5, 0.001, 0.02,
                             "All"),
        "Sample": (1.0, 1.0, 8000, 400.0, 1.0, 0.001, 0.02, "Selected"),
        "Screening": (1.0, 1.0, 200000, 400.0, 1.0, 0.001, 0.02, "Selected"),
    }
    assert list(IMP.bff.fps_mode_names()) == list(table)
    project = IMP.bff.FPSProject()
    for mode, expected in table.items():
        for block in (IMP.bff.fps_mode_parameters(mode),
                      project.get_parameters(mode)):
            got = (block.viscosity_factor, block.time_step_factor,
                   block.max_iterations, block.max_force,
                   block.clash_tolerance, block.k_tolerance,
                   block.t_tolerance, block.optimize_selected)
            assert got == expected, mode
            # rkT and ETolerance are the same in all five
            assert (block.rkt, block.e_tolerance) == (10.0, 100.0)
    assert project.conversion.forster_radius == 52.0
    assert project.conversion.polynomial_order == 3
    assert project.av.grid_size == 0.2
    assert project.av.e_samples == 200000


def test_docking_parameters_carry_only_what_this_port_honours():
    """FPS's integrator settings have nowhere to go here (PRD-121 D1)."""
    p = IMP.bff.FPSProject()
    p.score_set = "resolved"
    p.fixed_body = 1
    p.ev_weight = 0.5
    p.sigma_da = 4.0
    p.shuffle = 2.0
    d = p.get_docking_parameters("Refine")
    assert d.max_force == 10000.0        # the Refine block, not Dock's 400
    assert d.clash_tolerance == 0.5
    assert d.optimize_selected == "All"
    assert d.score_set == "resolved"
    assert (d.fixed_body, d.ev_weight, d.sigma_da) == (1, 0.5, 4.0)
    assert d.shuffle_max_translation == 2.0
    # NOT FPS's MaxIterations: 500000 integrator steps is not 500 CG iterations
    assert p.refine.max_iterations == 500000
    assert d.n_frames == p.refine.iterations == 500


# --------------------------------------------------------------------------
# the schema
# --------------------------------------------------------------------------
def test_project_section_validates_and_names_a_missing_distance():
    payload = {
        "Positions": {"a": {"residue_seq_number": 1},
                      "b": {"residue_seq_number": 2}},
        "Distances": {"d": {"position1_name": "a", "position2_name": "b",
                            "distance": 40.0, "error_neg": 2.0,
                            "error_pos": 2.0}},
        "Project": {"structures": ["x.pdb"],
                    "selected_distances": ["d", "ghost"]},
    }
    report = IMP.bff.fps_schema_validate(json.dumps(payload))
    assert any("ghost" in e for e in report.errors), list(report.errors)
    payload["Project"]["selected_distances"] = ["d"]
    report = IMP.bff.fps_schema_validate(json.dumps(payload))
    assert not report.errors, list(report.errors)


def test_project_section_reports_wrong_types_and_unknown_keys():
    report = IMP.bff.validate_project(json.dumps({
        "fixed_body": "first", "mode": "Wobble", "my_note": 1,
        "parameters": {"Dock": {"max_force": "lots"}}}))
    assert any("fixed_body" in e for e in report.errors)
    assert any("mode" in e for e in report.errors)
    assert any("max_force" in e for e in report.errors)
    assert any("my_note" in w for w in report.warnings)


def test_every_project_field_has_a_default():
    """An empty `default_json` is how a field declares it has none, and a
    project is meant to be writable one key at a time."""
    for table in (IMP.bff.fps_project_fields(),
                  IMP.bff.fps_mode_parameter_fields(),
                  IMP.bff.fps_conversion_fields(),
                  IMP.bff.fps_av_global_fields()):
        for spec in table:
            assert spec.default_json, spec.name
            assert not spec.required, spec.name


# --------------------------------------------------------------------------
# the shipped example, and the program
# --------------------------------------------------------------------------
def test_shipped_hiv_rt_project_runs_and_scores_what_the_flags_do():
    """`--project` and the twelve-flag spelling are the same run."""
    project = IMP.bff.read_fps_project(str(HIV / "hiv_rt.project.fps.json"))
    assert not list(project.get_problems()), list(project.get_problems())
    assert project.mode == "Dock"
    assert project.score_set == "resolved"
    assert [os.path.basename(s) for s in project.get_structure_paths()] == [
        "protein_1R0A.pdb", "dna.pdb"]
    # self-contained: the labelling is the project file itself
    assert project.get_labelling_path() == str(HIV / "hiv_rt.project.fps.json")

    from_project = IMP.bff.score_structures(
        list(project.get_structure_paths()),
        project.get_labelling_path(), project.score_set, False,
        project.sigma_da)
    from_flags = IMP.bff.score_structures(
        [str(HIV / "protein_1R0A.pdb"), str(HIV / "dna.pdb")],
        str(HIV / "hiv_rt.fps.json"), "resolved", False, 6.0)
    assert from_project.score == pytest.approx(from_flags.score, abs=1e-9)
    # Olga's radii gave 59.2547 for the hours they were the default on
    # 2026-09-01 (this fixture's DNA carries the old `C1*`/`O1P` spellings,
    # which that name-keyed table does not have, so those 85 atoms took its
    # 1.50 A unknown-name fallback). IMP's own radii are the default again
    # (ProbeAccessibleVolumeDecorator::set_radii_source) and 59.0404 returns exactly.
    assert from_project.score == pytest.approx(59.040389, abs=1e-3)
    assert from_project.n_distances == from_flags.n_distances == 18


def test_shipped_legacy_project_fixture_converts_to_the_resolved_set(tmp_path):
    """The fixture's Boolean[] is the `resolved` score set, by position."""
    program = _program("imp_bff_fps")
    from click.testing import CliRunner
    out = tmp_path / "converted.fps.json"
    result = CliRunner().invoke(program.cli, [
        "project", "convert", str(HIV / "hiv_rt.project.txt"),
        "-o", str(out)])
    assert result.exit_code == 0, result.output
    assert "18 selected of 20 distances" in result.output
    payload = json.loads(out.read_text())
    selected = set(payload["Project"]["selected_distances"])
    resolved = set(json.loads(
        (HIV / "hiv_rt.fps.json").read_text())["χ²"]["resolved"]["distances"])
    assert selected == resolved
    assert len(payload["Positions"]) == 11
    assert len(payload["Distances"]) == 20


def test_positional_selection_is_named_against_the_file_order():
    """FPS's Boolean[] indexes the distances FILE, whose order the fps.json
    object does not keep: its keys come back sorted."""
    order = list(IMP.bff.read_old_distances_order(str(HIV / "Distances.txt")))
    assert order[0] == "p66_Q6C_p_1bp"
    assert order[8:10] == ["p66_K287C_p_1bp", "p66_K287C_p_19bp"]
    document = IMP.bff.read_fps_json(str(HIV / "LabelingPositions.txt"))
    keyed = list(json.loads(document.distances))
    assert sorted(order) == sorted(keyed)
    assert order != keyed          # the trap this function exists for


def test_program_init_show_and_score_through_a_project(tmp_path):
    program = _program("imp_bff_fps")
    from click.testing import CliRunner
    runner = CliRunner()
    out = tmp_path / "run.fps.json"
    result = runner.invoke(program.cli, [
        "project", "init",
        "-p", str(HIV / "protein_1R0A.pdb"), "-p", str(HIV / "dna.pdb"),
        "-j", str(HIV / "hiv_rt.fps.json"), "-c", "resolved",
        "-o", str(out)])
    assert result.exit_code == 0, result.output
    assert "positions: 11   distances: 20" in result.output
    assert "score set: resolved (18 distances)" in result.output

    result = runner.invoke(program.cli, ["project", "show", str(out)])
    assert result.exit_code == 0, result.output
    assert "no problems" in result.output
    assert "Error estimation" in result.output
    assert "score set: resolved" in result.output

    result = runner.invoke(program.cli, [
        "score", "--project", str(out), "--no-pairs"])
    assert result.exit_code == 0, result.output
    # moved with the pin above (back off Olga's radii, 2026-09-01)
    assert "score (chi2): 59.0404" in result.output
    assert "volumes: 10   distances: 18" in result.output


def test_program_reports_a_broken_project_rather_than_crashing(tmp_path):
    """A stale path is a report, not an exception: a run may re-point it."""
    p = IMP.bff.FPSProject()
    p.structures = ["gone.pdb"]
    p.labelling_json = "missing.fps.json"
    out = tmp_path / "broken.fps.json"
    IMP.bff.write_fps_project(str(out), p)
    problems = list(IMP.bff.read_fps_project(str(out)).get_problems())
    assert any("gone.pdb" in x for x in problems)
    assert any("missing.fps.json" in x for x in problems)

    program = _program("imp_bff_fps")
    from click.testing import CliRunner
    result = CliRunner().invoke(program.cli, ["project", "show", str(out)])
    assert result.exit_code == 0, result.output
    assert "MISSING" in result.output


def test_dock_saves_the_pose_and_the_settings_it_actually_used(tmp_path):
    """`--save-project` closes the loop: a run can be stored and continued."""
    program = _program("imp_bff_fps")
    from click.testing import CliRunner
    saved = tmp_path / "after.fps.json"
    result = CliRunner().invoke(program.cli, [
        "dock", "--project", str(HIV / "hiv_rt.project.fps.json"),
        "-o", str(tmp_path / "dock_out"), "-n", "5", "--shuffle", "0",
        "--seed", "1", "--save-project", str(saved)])
    assert result.exit_code == 0, result.output
    after = IMP.bff.read_fps_project(str(saved))
    poses = json.loads(after.poses)
    assert len(poses) == 2
    assert {"body_id", "t", "q"} <= set(poses[0])
    # what the run did, not what the project said
    assert after.shuffle == 0.0
    assert after.dock.iterations == 5
    assert after.score_set == "resolved"
