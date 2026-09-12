"""The two inputs greedy pair selection needs, and the pipeline that feeds it.

``select_probe_pairs`` was reachable but unreachable *from a structure*:
its tests handed it synthetic matrices, and nothing in the tree turned an
ensemble into the ``(n_frames, n_pairs)`` efficiencies and the ``(n, n)`` RMSDs
it takes. The two halves of that gap are covered here --
``ProbeNetworkRestraint.get_pair_efficiencies`` and ``pairwise_rmsd`` -- and
then the whole route end to end through ``imp_bff select-pairs``.
"""

import numpy as np
import pytest

import IMP
import IMP.algebra
import IMP.atom
import IMP.core
import IMP.rmf
import IMP.bff
import RMF


# ---------------------------------------------------------------------------
# pairwise_rmsd
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def ensemble():
    """Four conformers of six atoms: one deformed, one displaced, one rotated.

    Frame 1 is the only one that genuinely differs from frame 0 -- frames 2 and
    3 are rigid motions of it, which is what the superposition has to see
    through.
    """
    rng = np.random.RandomState(7)
    base = rng.normal(0, 5.0, (6, 3))
    deformed = base + rng.normal(0, 0.3, base.shape)
    rotated = base @ np.array([[0.0, -1.0, 0.0],
                               [1.0, 0.0, 0.0],
                               [0.0, 0.0, 1.0]]) + np.array([3.0, -1.0, 2.0])
    return np.stack([base, deformed, base + 4.0, rotated])


def test_pairwise_rmsd_is_a_symmetric_matrix_with_a_zero_diagonal(ensemble):
    got = IMP.bff.pairwise_rmsd(ensemble, False)
    assert got.shape == (len(ensemble), len(ensemble))
    np.testing.assert_allclose(got, got.T, atol=1e-12)
    np.testing.assert_allclose(np.diag(got), 0.0, atol=1e-12)


def test_pairwise_rmsd_agrees_with_compute_rmsd_pair_by_pair(ensemble):
    for superpose in (False, True):
        got = IMP.bff.pairwise_rmsd(ensemble, superpose)
        for i in range(len(ensemble)):
            for j in range(len(ensemble)):
                if i == j:
                    continue
                assert got[i, j] == pytest.approx(
                    IMP.bff.get_rmsd(ensemble[i].ravel(),
                                         ensemble[j].ravel(), [], superpose),
                    abs=1e-9)


def test_superposition_removes_the_rigid_motions_and_nothing_else(ensemble):
    aligned = IMP.bff.pairwise_rmsd(ensemble, True)
    raw = IMP.bff.pairwise_rmsd(ensemble, False)

    # frame 3 is frame 0 rotated and translated: a rigid motion, and the whole
    # point of superposing is that it is not a structural difference
    assert raw[0, 3] > 5.0
    assert aligned[0, 3] == pytest.approx(0.0, abs=1e-8)
    # a translation alone is just as rigid
    assert raw[0, 2] == pytest.approx(4.0 * np.sqrt(3.0), abs=1e-9)
    assert aligned[0, 2] == pytest.approx(0.0, abs=1e-8)
    # what is not a rigid motion survives; superposing can only shrink it
    assert 0.0 < aligned[0, 1] <= raw[0, 1] + 1e-12
    assert aligned[0, 1] > 0.5 * raw[0, 1]


def test_pairwise_rmsd_rejects_coordinates_that_are_not_three_dimensional():
    with pytest.raises(Exception):
        IMP.bff.pairwise_rmsd(np.zeros((3, 4, 2)), False)


def test_an_empty_ensemble_gives_an_empty_matrix():
    assert IMP.bff.pairwise_rmsd(np.zeros((0, 5, 3)), True).shape == (0, 0)


# ---------------------------------------------------------------------------
# the restraint's candidate pairs
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def t4l_restraint():
    IMP.set_log_level(IMP.SILENT)
    m = IMP.Model()
    handle = RMF.open_rmf_file_read_only(
        IMP.bff.get_example_path("structure/T4L/t4l_docking.rmf3"))
    hier = IMP.rmf.create_hierarchies(handle, m)[0]
    IMP.rmf.load_frame(handle, RMF.FrameID(0))
    beads = [IMP.core.XYZ(p) for p in IMP.atom.get_leaves(hier)]
    restraint = IMP.bff.ProbeNetworkRestraint(
        hier, IMP.bff.get_example_path("structure/T4L/fret.fps.json"),
        score_set="chi2_C1_33p")
    # `m` and `hier` are returned rather than dropped: the model is reference
    # counted from Python, and a restraint over a collected model is a fault,
    # not an exception.
    return m, hier, handle, restraint, beads


def test_the_pair_names_are_the_score_sets_distances_in_a_stated_order(
        t4l_restraint):
    _, _, _, restraint, _ = t4l_restraint
    names = list(restraint.get_pair_names())
    used = restraint.get_used_distances()
    assert names == sorted(used.keys()), \
        "the column order is the contract; it has to be the sorted one"
    assert len(names) == 33


def test_the_efficiencies_are_one_per_pair_and_are_efficiencies(t4l_restraint):
    _, _, _, restraint, _ = t4l_restraint
    effs = np.asarray(restraint.get_pair_efficiencies())
    assert effs.shape == (len(restraint.get_pair_names()),)
    assert np.isfinite(effs).all()
    assert ((effs > 0.0) & (effs < 1.0)).all()


def test_an_efficiency_is_what_get_model_distance_reports_for_that_pair(
        t4l_restraint):
    _, _, _, restraint, _ = t4l_restraint
    effs = restraint.get_pair_efficiencies()
    used = restraint.get_used_distances()
    for name, value in zip(restraint.get_pair_names(), effs):
        measurement = used[name]
        assert value == pytest.approx(restraint.get_model_distance(
            measurement.position_1, measurement.position_2,
            measurement.forster_radius, IMP.bff.PROBE_PAIR_EFFICIENCY))


def test_reading_a_frame_does_not_need_it_scored_first(t4l_restraint):
    """`get_pair_efficiencies` re-evaluates, so a freshly loaded frame answers
    for itself rather than for whichever frame was scored last."""
    _, _, handle, restraint, _ = t4l_restraint
    frames = list(handle.get_root_frames())
    IMP.rmf.load_frame(handle, frames[0])
    first = np.asarray(restraint.get_pair_efficiencies())
    IMP.rmf.load_frame(handle, frames[len(frames) // 2])
    middle = np.asarray(restraint.get_pair_efficiencies())
    IMP.rmf.load_frame(handle, frames[0])
    again = np.asarray(restraint.get_pair_efficiencies())

    assert not np.allclose(first, middle), "different frames, same answer"
    np.testing.assert_allclose(first, again, atol=1e-12)


# ---------------------------------------------------------------------------
# ensemble in, ranking out
# ---------------------------------------------------------------------------

def test_the_two_inputs_feed_the_selector(t4l_restraint):
    _, _, handle, restraint, beads = t4l_restraint
    frames = list(handle.get_root_frames())[::10]
    names = list(restraint.get_pair_names())

    effs = np.empty((len(frames), len(names)))
    coords = np.empty((len(frames), len(beads), 3))
    for i, frame in enumerate(frames):
        IMP.rmf.load_frame(handle, frame)
        effs[i] = restraint.get_pair_efficiencies()
        coords[i] = [b.get_coordinates() for b in beads]

    rmsds = IMP.bff.pairwise_rmsd(coords, True)
    selected, decay = IMP.bff.select_probe_pairs(
        effs, rmsds, measurement_error=0.06, max_pairs=5)

    assert len(selected) == len(decay) == 5
    assert len(set(selected)) == 5, "unique_only defaults on"
    assert set(selected) <= set(range(len(names)))

    # The first measurement narrows the answer, and the selection as a whole
    # narrows it further. Step-by-step monotonicity is *not* asserted: the
    # greedy scores each candidate against the chi-squared accumulated so far,
    # and on a small ensemble the best remaining candidate can still leave a
    # larger expectation than the step before it -- observed on a four-member
    # ensemble, where the third selection gives back 0.03 A.
    uninformed = IMP.bff.expected_rmsd(
        rmsds.ravel(), np.zeros(rmsds.size), 1, 0.99, len(frames))
    assert decay[0] < uninformed
    assert decay[-1] < decay[0]
    assert np.isfinite(decay).all() and (decay > 0.0).all()


def test_select_pairs_runs_end_to_end(imp_bff_program, tmp_path):
    from click.testing import CliRunner
    out = tmp_path / "ranking.tsv"
    result = CliRunner().invoke(imp_bff_program.select_pairs, [
        "--fps-json", IMP.bff.get_example_path("structure/T4L/fret.fps.json"),
        "--rmf", IMP.bff.get_example_path("structure/T4L/t4l_docking.rmf3"),
        "--score-set", "chi2_C1_33p",
        "--stride", "10", "--max-pairs", "4",
        "--output", str(out)])
    assert result.exit_code == 0, result.output
    assert out.is_file()

    import pandas as pd
    table = pd.read_csv(out, sep="\t")
    assert list(table.columns) == [
        "rank", "pair", "position_1", "position_2", "forster_radius",
        "expected_rmsd", "gain"]
    assert len(table) == 4
    assert list(table["rank"]) == [1, 2, 3, 4]
    assert table["pair"].is_unique
    assert (table["expected_rmsd"] > 0.0).all()
    assert table["expected_rmsd"].iloc[-1] < table["expected_rmsd"].iloc[0]
    assert table["gain"].iloc[0] > 0.0


def test_select_pairs_takes_a_stack_of_pdbs(imp_bff_program, tmp_path):
    """The other way in: one structure per file rather than one per frame."""
    from click.testing import CliRunner
    source = IMP.bff.get_example_path("structure/T4L/3GUN.pdb")
    rng = np.random.RandomState(3)
    paths = []
    for k in range(3):
        m = IMP.Model()
        h = IMP.atom.read_pdb(source, m,
                              IMP.atom.NonWaterNonHydrogenPDBSelector())
        for p in IMP.atom.get_leaves(h):
            xyz = IMP.core.XYZ(p)
            xyz.set_coordinates(xyz.get_coordinates() + IMP.algebra.Vector3D(
                *rng.normal(0, 0.4 * k, 3)))
        path = tmp_path / f"conf{k}.pdb"
        IMP.atom.write_pdb(h, str(path))
        paths.append(str(path))

    args = ["--fps-json",
            IMP.bff.get_example_path("structure/T4L/fret.fps.json"),
            "--score-set", "chi2_C1_33p", "--max-pairs", "2"]
    for path in paths:
        args += ["--pdb", path]
    result = CliRunner().invoke(imp_bff_program.select_pairs, args)
    assert result.exit_code == 0, result.output
    assert "3 PDB structures" in result.output

    # a file that is not the same molecule is refused, not silently zipped
    short = tmp_path / "short.pdb"
    short.write_text("\n".join(
        open(paths[0]).read().splitlines()[:200]) + "\nEND\n")
    refused = CliRunner().invoke(imp_bff_program.select_pairs,
                                 args + ["--pdb", str(short)])
    assert refused.exit_code != 0
    assert "must correspond" in refused.output


def test_select_pairs_needs_an_ensemble_and_says_so(imp_bff_program):
    from click.testing import CliRunner
    result = CliRunner().invoke(imp_bff_program.select_pairs, [
        "--fps-json", IMP.bff.get_example_path("structure/T4L/fret.fps.json"),
        "--rmf", IMP.bff.get_example_path("structure/T4L/t4l_docking.rmf3"),
        "--stride", "1000"])
    assert result.exit_code != 0
    assert "ensemble" in result.output


def test_select_pairs_wants_exactly_one_source_of_structures(imp_bff_program):
    from click.testing import CliRunner
    fps = IMP.bff.get_example_path("structure/T4L/fret.fps.json")
    both = CliRunner().invoke(imp_bff_program.select_pairs, [
        "--fps-json", fps,
        "--rmf", IMP.bff.get_example_path("structure/T4L/t4l_docking.rmf3"),
        "--pdb", IMP.bff.get_example_path("structure/T4L/3GUN.pdb")])
    neither = CliRunner().invoke(imp_bff_program.select_pairs,
                                 ["--fps-json", fps])
    for result in (both, neither):
        assert result.exit_code != 0
        assert "--rmf" in result.output and "--pdb" in result.output
