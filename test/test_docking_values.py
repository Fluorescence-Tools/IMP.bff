"""The values a docking run is told and reports, and the table it writes.

The engine itself (`IMP.bff.dock`, `dock_minimize`, `screen`) has no test in
this tree and did not have one when it was Python either -- it drives IMP.pmi
over a whole assembly, and nothing here builds one. What *can* be pinned is
what it is told and what it reports, which is what moved to C++
(`include/IMP/bff/Docking.h`): the per-pair table, the derived numbers on it,
the fps.json distance vocabulary, and the CSV.

The derived numbers matter more than they look. A residual is `model - data`
(PRD-117, okf/validation/two_chi2_conventions.md) and the asymmetric chi2
picks its error bar from that sign; the engine's own diagnostics table is
where a user reads whether a fit is off, so a sign convention that disagrees
with the restraint being minimised is a table that says the opposite of what
was optimised.
"""

import csv
import json
import math

import pytest

import IMP.bff


def _pair(model=55.0, exp=50.0, neg=1.0, pos=10.0, r0=52.0):
    p = IMP.bff.PairDistance()
    p.name = "d1_a1"
    p.position1 = "d1"
    p.position2 = "a1"
    p.distance_exp = exp
    p.distance_model = model
    p.distance_type = "RDAMean"
    p.forster_radius = r0
    p.error_neg = neg
    p.error_pos = pos
    return p


def test_residual_is_model_minus_data():
    """Too large a model is a *positive* deviation (the 2026-08-24 ruling)."""
    assert _pair(model=55.0, exp=50.0).residual == pytest.approx(5.0)
    assert _pair(model=45.0, exp=50.0).residual == pytest.approx(-5.0)


def test_chi2_is_the_one_kernel_and_picks_the_error_bar_by_sign():
    """The table and the restraint score the same way -- one implementation.

    A model above the experiment is judged against `error_pos`, below against
    `error_neg`. With bars 1 and 10 that is a factor of 100 in chi2, which is
    what makes getting it backwards visible rather than subtle.
    """
    above = _pair(model=60.0, exp=50.0, neg=1.0, pos=10.0)
    below = _pair(model=40.0, exp=50.0, neg=1.0, pos=10.0)
    assert above.chi2 == pytest.approx(1.0)
    assert below.chi2 == pytest.approx(100.0)
    assert above.chi2 == pytest.approx(
        IMP.bff.chi2_score(60.0, 50.0, 1.0, 10.0))


def test_efficiencies_are_the_single_pair_kernel():
    p = _pair(model=52.0, exp=26.0, r0=52.0)
    assert p.efficiency_model == pytest.approx(0.5)          # r = R0
    assert p.efficiency_exp == pytest.approx(
        IMP.bff.fret_efficiency(26.0, 52.0))


def test_the_fps_distance_vocabulary_reads_both_ways():
    """The names an fps.json uses, published once beside the enum.

    They were a private member of the reader, so a table that wanted to report
    which convention it had scored re-derived them from the JSON -- and said
    `RDAMean` for an entry whose absent `distance_type` the reader had
    defaulted to `RDAMeanE`.
    """
    for name in ("RDAMean", "RDAMeanE", "Rmp", "Efficiency", "pRDA"):
        code = IMP.bff.probe_pair_distance_type(name)
        assert code >= 0, name
        assert IMP.bff.probe_pair_distance_type_name(code) == name
    assert IMP.bff.probe_pair_distance_type("not a distance") == -1
    assert IMP.bff.probe_pair_distance_type_name(999) == ""


def test_score_csv_has_a_row_per_pair_and_round_trips(tmp_path):
    pairs = [_pair(model=55.0), _pair(model=45.0)]
    pairs[1].name = "d2_a2"
    out = tmp_path / "scores.csv"
    IMP.bff.write_score_csv(str(out), 12.5, pairs)

    rows = list(csv.reader(out.open()))
    assert rows[0] == ["# total_score", "12.5"]
    assert rows[1][:3] == ["name", "position1", "position2"]
    assert [r[0] for r in rows[2:]] == ["d1_a1", "d2_a2"]
    header = rows[1]
    first = dict(zip(header, rows[2]))
    assert float(first["distance_model"]) == pytest.approx(55.0)
    assert float(first["residual"]) == pytest.approx(5.0)
    assert float(first["chi2"]) == pytest.approx(pairs[0].chi2, abs=1e-3)


def test_a_result_serialises_with_its_derived_numbers(tmp_path):
    """`get_json` is what a caller stores; the derived numbers are in it."""
    result = IMP.bff.DockingResult(
        score=3.25, n_avs=2, n_distances=1, pairs=[_pair()],
        output_dir=str(tmp_path), score_csv="scores.csv",
        extra=json.dumps({"method": "minimize"}),
        poses=json.dumps([{"body_id": 0, "t": [0.0, 0.0, 0.0]}]))
    payload = json.loads(result.get_json())
    assert payload["score"] == pytest.approx(3.25)
    assert payload["extra"]["method"] == "minimize"
    assert payload["poses"][0]["body_id"] == 0
    pair = payload["pairs"][0]
    assert pair["residual"] == pytest.approx(5.0)
    assert pair["chi2"] == pytest.approx(_pair().chi2)
    assert 0.0 < pair["E_model"] < 1.0


def test_default_parameters_are_the_documented_ones():
    """The defaults are the experiment's, so they are worth pinning."""
    p = IMP.bff.DockingParameters()
    assert (p.n_frames, p.mc_steps) == (500, 10)
    assert p.mean_position_restraint is True     # sampling uses the fast one
    assert p.sigma_da == pytest.approx(6.0)
    assert p.coarse_clash is True
    assert p.refine_av_cycles == 0               # opt-in: it re-runs the AVs


def test_pair_distances_at_positions_reads_the_particles_it_is_given(tmp_path):
    """The minimisation path's table comes from the proxies that moved.

    The volumes themselves stay where the run started -- they are not rigid
    body members in that path -- so reading them would report the input pose.
    """
    m = IMP.Model()
    p1 = IMP.Particle(m, "d1")
    p2 = IMP.Particle(m, "a1")
    IMP.core.XYZ.setup_particle(p1, IMP.algebra.Vector3D(0.0, 0.0, 0.0))
    IMP.core.XYZ.setup_particle(p2, IMP.algebra.Vector3D(0.0, 0.0, 40.0))

    e = IMP.bff.AVPairDistanceMeasurement()
    e.position_1, e.position_2 = "d1", "a1"
    e.distance, e.error_neg, e.error_pos = 42.0, 2.0, 2.0
    e.forster_radius = 52.0
    e.distance_type = IMP.bff.PROBE_PAIR_DISTANCE_MP
    distances = IMP.bff.MapStringAVPairDistanceMeasurement()
    distances["d1_a1"] = e

    pairs = list(IMP.bff.pair_distances_at_positions(
        distances, m, ["d1", "a1"], [p1.get_index(), p2.get_index()], 6.0))
    assert len(pairs) == 1
    # Rmp is the separation itself, so the transfer function returns it
    assert pairs[0].distance_model == pytest.approx(40.0, abs=1e-6)
    assert pairs[0].distance_type == "Rmp"
    assert pairs[0].residual == pytest.approx(-2.0, abs=1e-6)

    # a position the caller did not supply cannot be measured, and says so
    pairs = list(IMP.bff.pair_distances_at_positions(
        distances, m, ["d1"], [p1.get_index()], 6.0))
    assert math.isnan(pairs[0].distance_model)


def test_scoring_an_assembly_agrees_with_the_restraint_it_is_built_on(tmp_path):
    """`score()` on the bundled T4L fixture, end to end.

    The number is not arbitrary: it is the same 11.326777201821047 that
    `test_ProbeNetworkRestraint.py` pins for the quadrature score of this
    fixture, so this asserts that assembling a model and scoring it through
    the docking path gives what the restraint gives on its own. It was
    22.69935116408601 until the positions' `strip_mask` began to be honoured
    (PRD-106) and 22.527697710399956 until the accessible contact volume they
    all ask for stopped being ignored (2026-09-01, PRD-121 G9). Olga's
    name-keyed radii were briefly the default the same day and gave
    11.286617190950983; IMP's own are the default again (AV::set_radii_source)
    -- and that both paths moved together each time, and back together, is the
    point of pinning them against each other.
    """
    pdb = IMP.bff.get_example_path("structure/T4L/3GUN.pdb")
    fps = IMP.bff.get_example_path("structure/T4L/fret.fps.json")
    out = tmp_path / "scores.csv"
    result = IMP.bff.score_structures([pdb], fps, score_set="chi2_C2_33p",
                                      output_csv=str(out))

    assert result.score == pytest.approx(11.326777201821047, abs=1e-9)
    assert result.n_avs == 17
    assert result.n_distances == 33
    assert result.score_csv == str(out)

    rows = list(csv.reader(out.open()))
    assert len(rows) == 2 + result.n_distances       # total, header, one each
    names = {r[0] for r in rows[2:]}
    assert "19-132_C2" in names
    # every pair carries a model distance and a chi2 the table can show
    for pair in result.pairs:
        assert math.isfinite(pair.distance_model), pair.name
        assert math.isfinite(pair.chi2), pair.name


def test_minimisation_improves_the_score_and_writes_its_pose(tmp_path):
    """`dock_minimize` end to end, from the same fixture.

    **What this does not test.** The bundled fixture is *one* rigid body, and
    no rigid motion of a single body changes a distance between two of its own
    dyes -- so there is nothing here for a docking optimiser to optimise, and
    the run's score moving at all (16.02 -> 15.48 when measured) is the AV
    resampling between builds, not the minimiser. Nothing in this tree
    exercises multi-body docking; that needs a two-body fixture, which is a
    gap worth filling and not something this test pretends to cover.

    What it does test is that the machinery runs to completion, reports finite
    numbers and the files it promises, and -- the real check -- that the
    docked state it hands back reproduces its own score when applied to a
    fresh assembly. A pose that does not round-trip is a pose that cannot be
    stored and continued from, which is what it exists for.
    """
    pdb = IMP.bff.get_example_path("structure/T4L/3GUN.pdb")
    fps = IMP.bff.get_example_path("structure/T4L/fret.fps.json")
    params = IMP.bff.DockingParameters()
    params.score_set = "chi2_C2_33p"
    params.shuffle_max_translation = 0.0

    result = IMP.bff.dock_minimize([pdb], fps, str(tmp_path), params=params)

    assert math.isfinite(result.score)
    assert result.n_distances == 33
    for name in ("docked.pdb", "scores.csv", "convergence.csv"):
        assert (tmp_path / name).exists(), name

    poses = json.loads(result.poses)
    assert poses and set(poses[0]) == {"body_id", "t", "q"}
    # the state goes back onto an assembly built from the same inputs
    fresh = IMP.bff.build_docking_assembly([pdb], fps,
                                           score_set="chi2_C2_33p")
    IMP.bff.apply_poses(fresh, result.poses)
    assert IMP.bff.score_assembly(fresh).score == pytest.approx(result.score,
                                                               abs=1e-6)


def test_refine_runs_and_writes_a_structure(tmp_path):
    pdb = IMP.bff.get_example_path("structure/T4L/3GUN.pdb")
    fps = IMP.bff.get_example_path("structure/T4L/fret.fps.json")
    result = IMP.bff.refine_docking([pdb], fps, str(tmp_path),
                                    score_set="chi2_C2_33p", steps=25)
    assert math.isfinite(result.score)
    assert (tmp_path / "refined.pdb").exists()
    assert result.score_csv.endswith("scores.csv")


def test_screening_ranks_a_library_and_never_reports_a_silent_nan(tmp_path):
    """`screen` scores each structure and sorts, best first.

    It catches per-structure exceptions and records NaN so one bad file does
    not stop a library -- which also means a broken engine produced a table of
    NaNs and looked like a screening run with unscorable inputs. The NaN
    column is asserted empty here for that reason.
    """
    pdb = IMP.bff.get_example_path("structure/T4L/3GUN.pdb")
    fps = IMP.bff.get_example_path("structure/T4L/fret.fps.json")
    out = tmp_path / "ranked.csv"
    ranked = list(IMP.bff.screen_structures(
        [pdb, pdb], fps, score_set="chi2_C2_33p", output_csv=str(out)))

    assert len(ranked) == 2
    assert all(math.isfinite(entry.score) for entry in ranked)
    assert ranked[0].score == pytest.approx(11.326777201821047, abs=1e-9)
    rows = list(csv.reader(out.open()))
    # The score is followed by FPS's screening diagnostics (PRD-121 G3);
    # test/representation/test_screening_diagnostics.py checks what they say.
    assert rows[0] == ["pdb", "score", "chi2_r", "sigma1", "sigma2", "sigma3",
                       "invalid_r", "ref_rmsd"]
    assert len(rows) == 3


# IMP runs every .py under test/ as a standalone script, so a bare pytest file
# would import cleanly and exit 0 -- reporting success without running.
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))


def test_a_structure_missing_a_labelling_site_scores_nan_and_does_not_crash():
    """Screening a library is exactly where sites go missing.

    `search_labeling_site` read `p_residue[0]` and *then* asked whether
    `p_residue` held one element, and `IMP_USAGE_CHECK` compiles out of a
    release build -- so a structure without the chain a position names indexed
    element zero of an empty vector and the process died. That is the ordinary
    case for a screen: candidate structures differ, and some do not carry every
    site. It must be a NaN row, not a segmentation fault.
    """
    import math

    protein = IMP.bff.get_example_path("structure/HIV_RT/protein_1R0A.pdb")
    dna = IMP.bff.get_example_path("structure/HIV_RT/dna.pdb")
    fps = IMP.bff.get_example_path("structure/HIV_RT/hiv_rt.fps.json")

    # Each path is one candidate, and neither carries all eleven positions:
    # the protein has no chain P, the DNA no chains A or B.
    ranked = list(IMP.bff.screen_structures([protein, dna], fps, "resolved"))
    assert len(ranked) == 2
    assert all(math.isnan(entry.score) for entry in ranked)

    # and the reason is available rather than silent
    with pytest.raises(Exception) as excinfo:
        IMP.bff.build_docking_assembly([dna], fps, "resolved")
    assert "labelling site" in str(excinfo.value)


def test_the_two_bodies_together_do_carry_every_site():
    """The other half of the pair above: assembled, the sites resolve."""
    protein = IMP.bff.get_example_path("structure/HIV_RT/protein_1R0A.pdb")
    dna = IMP.bff.get_example_path("structure/HIV_RT/dna.pdb")
    fps = IMP.bff.get_example_path("structure/HIV_RT/hiv_rt.fps.json")
    result = IMP.bff.score_structures([protein, dna], fps, "resolved")
    assert result.n_avs == 10
    assert result.n_distances == 18
    # Olga's radii gave 59.2547 for the hours they were the default on
    # 2026-09-01: this fixture's DNA uses the old `C1*`/`O1P` spellings, which
    # that name-keyed table does not carry, so those atoms took its 1.50 A
    # unknown-name fallback. IMP's own radii are the default again
    # (AV::set_radii_source) and 59.0404 returns exactly.
    assert result.score == pytest.approx(59.040389, abs=1e-3)


def test_pairs_survive_a_temporary_result():
    """`score_structures(...).pairs[0]` used to raise IndexError.

    SWIG returned a bare `std::vector` **member** of a temporary as a pointer
    into freed memory, and a freed vector reports size 0 -- so the expression
    read as "this structure has no distances", which is the worst possible
    shape for a diagnostics bug: a wrong answer, not an error. Binding the
    result to a name first happened to work, so it looked intermittent.
    Typing the member as `PairDistances` (the IMP_VALUES vector) fixes it.
    """
    pdb = IMP.bff.get_example_path("structure/T4L/3GUN.pdb")
    fps = IMP.bff.get_example_path("structure/T4L/fret.fps.json")
    # the temporary is the point: no intermediate name
    assert len(IMP.bff.score_structures([pdb], fps, "chi2_C2_33p").pairs) == 33
    assert IMP.bff.score_structures([pdb], fps,
                                    "chi2_C2_33p").pairs[0].name == "19-119_C2"
