"""The second A/B: maltose-binding protein, where the beta sheet lives.

`test_labelizer_ab.py` compares against the reference on 1DDB, which is
**all-helical** -- 96 H and not one strand. Everything DSSP does with bridges,
ladders and bulges is therefore unexercised there, and an implementation that
never assigned `E` at all would pass every test in that file.

This closes it. The paper's Supplementary Data 1 carries a per-residue table
for maltose-binding protein in **two conformations** -- 1OMP (apo, open) and
1ANF (holo, closed) -- with the same CS/SE/CR/SS columns, 370 rows each, plus
the two-state combined score and its delta. MBP is a two-domain alpha/beta
fold, so the reference's own `SS` column contains 72 `E` and 5 `B`.

Two things follow that the 1DDB case could not give:

* the **strand** half of the assignment is compared against real DSSP output,
  not merely exercised;
* the **two-state** layer -- which had no reference output at all before -- can
  be checked against the published combined score.

Precision: the supplement is rounded to five decimals, so comparisons here are
at `1e-5`, not exact. That is a property of the published file, not a tolerance
on the port -- `cr` agrees on all 370 residues to the last digit the supplement
prints.

Source: Gebhardt *et al.*, *Nat. Commun.* **16**, 3305 (2025), Supplementary
Data 1, sheet "MalE - LS". Extracted to CSV in `test/input/labelizer/`.
"""

import csv
import os

import numpy as np
import pytest

import IMP.bff as bff

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "input", "labelizer")

#: The supplement prints five decimals; nothing here can be tighter.
SUPPLEMENT_PRECISION = 1e-5

#: Residues in the MalE table.
N_RESIDUES = 370


def _path(name):
    return os.path.join(DATA, name)


def _reference(tag):
    """The published per-residue scores for one conformation."""
    with open(_path("malE_%s_reference.csv" % tag)) as fh:
        return {int(r["seq_id"]): r for r in csv.DictReader(fh)}


@pytest.fixture(scope="module")
def holo():
    """1ANF -- maltose-bound, the closed conformation."""
    return bff.labelizer_read_structure(_path("1anf.pdb"))


def test_the_reference_table_is_the_shape_the_paper_published():
    reference = _reference("1ANF")
    assert len(reference) == N_RESIDUES
    # Unlike 1DDB's conservation column, every term here is richly populated --
    # this reference is not degenerate and can be compared against.
    for column, least in (("cs", 5), ("se", 5), ("cr", 15), ("ss", 6)):
        distinct = {round(float(r[column]), 5) for r in reference.values()}
        assert len(distinct) >= least, "%s has only %d distinct values" % (
            column, len(distinct))


def test_the_structure_matches_the_table(holo):
    assert len(holo.residues) == N_RESIDUES
    assert {r.chain for r in holo.residues} == {"A"}


def test_cysteine_resemblance_is_exact_on_a_second_protein(holo):
    """A pure lookup, so a second structure is a real check of the mapping."""
    reference = _reference("1ANF")
    model = bff.LabelizerParameterList()
    model.append(bff.LabelizerParameter("cr", "C_CR1_Name", 1))
    got = {r.seq_id: r.value
           for r in bff.labelizer_parameter_scores(holo, model, bff.LabelizerOptions(), {})}

    delta = [abs(got[k] - float(v["cr"])) for k, v in reference.items()
             if k in got]
    assert len(delta) == N_RESIDUES
    assert max(delta) < SUPPLEMENT_PRECISION


def test_the_native_dssp_finds_the_sheet(holo):
    """The assignment the all-helical case cannot reach.

    Recorded 2026-08-24: 359/369 exact (97.3 %). The remaining ten are interior
    strand residues, and they are **not** a threshold artifact -- the published
    -0.5 kcal/mol hydrogen-bond cutoff is optimal here (94.3 % at -0.4, 94.9 %
    at -0.6), so what is left is DSSP's sheet-level bookkeeping rather than its
    energetics. See `okf/validation/labelizer_ab.md`.

    The bar is set below the measurement so drift fails rather than noise.
    """
    reference = _reference("1ANF")
    table = dict(bff.labelizer_load_table("C_SS1_SS").by_key)
    by_value = {round(v, 5): k for k, v in table.items()}

    ours = bff.labelizer_dssp(holo)
    want, got = [], []
    for i, residue in enumerate(holo.residues):
        row = reference.get(residue.seq_id)
        if row is None:
            continue
        letter = by_value.get(round(float(row["ss"]), 5))
        if letter is None:
            continue
        want.append(letter)
        got.append(ours[i])

    assert len(want) > 360
    agreement = sum(a == b for a, b in zip(want, got)) / len(want)
    assert agreement > 0.95, "DSSP agreement regressed to %.3f" % agreement

    # The sheet has to actually be found, not approximated by coils.
    assert want.count("E") > 60, "the reference sees a real sheet here"
    assert got.count("E") > 55, "we must find most of it: %d" % got.count("E")
    # Bulged ladders must be joined: without that, strands fragment into
    # isolated bridges and B outnumbers what the reference sees several times
    # over. This is the assertion that guards the union-find bulge joining.
    assert got.count("B") <= want.count("B") + 5, (
        "too many isolated bridges (%d vs %d) -- ladder joining regressed"
        % (got.count("B"), want.count("B")))


def test_the_three_ten_helix_is_not_a_leftover_stub(holo):
    """A `G` must be three residues, not the remainder after `H` took the rest.

    Assigning the leftover produced one- and two-residue `G` stubs at helix
    C-termini where the reference has `T` -- six of them on this structure.
    """
    reference = _reference("1ANF")
    table = dict(bff.labelizer_load_table("C_SS1_SS").by_key)
    by_value = {round(v, 5): k for k, v in table.items()}
    ours = bff.labelizer_dssp(holo)

    n_ref_g = sum(1 for r in reference.values()
                  if by_value.get(round(float(r["ss"]), 5)) == "G")
    assert ours.count("G") <= n_ref_g + 2, (
        "%d G against the reference's %d -- helix stubs are back"
        % (ours.count("G"), n_ref_g))

    # And no G run shorter than three residues should exist at all.
    runs = [len(r) for r in
            "".join(c if c == "G" else " " for c in ours).split()]
    assert all(n >= 3 for n in runs), "short G runs: %r" % runs


def test_the_two_conformations_differ_where_the_paper_says_they_do(holo):
    """The two-state layer, against published output for the first time.

    Checked here on the published table alone: which terms may differ between
    two conformations of one molecule, and which may not. Conservation and
    residue identity are properties of the *sequence*, so they must be
    identical in both columns; solvent exposure is a property of the *fold*,
    so it must not be. The coordinates themselves are exercised by the pair
    tests below, which use both shipped structures.
    """
    apo, holo_ref = _reference("1OMP"), _reference("1ANF")
    unchanged = [k for k in holo_ref
                 if abs(float(apo[k]["ls"]) - float(holo_ref[k]["ls"])) < 1e-9]
    assert len(unchanged) > 100, "the two conformations share most positions"

    # Every term that does not depend on the coordinates must be identical
    # between the two conformations -- conservation and residue identity are
    # properties of the sequence, not of the fold.
    for k in holo_ref:
        assert apo[k]["cs"] == holo_ref[k]["cs"]
        assert apo[k]["cr"] == holo_ref[k]["cr"]
    # ...and the coordinate-dependent ones must not be.
    assert any(apo[k]["se"] != holo_ref[k]["se"] for k in holo_ref)


def test_the_published_delta_is_the_difference_of_the_two_scores():
    """Pins what the paper's `Delta LS` column means, since the pair layer's
    two-state score is built on the same idea."""
    apo, holo_ref = _reference("1OMP"), _reference("1ANF")
    with open(_path("malE_two_state_reference.csv")) as fh:
        two_state = {int(r["seq_id"]): r for r in csv.DictReader(fh)}

    checked = 0
    for k, row in two_state.items():
        delta = row["delta_ls"]
        if delta in ("", None):
            continue
        expected = float(holo_ref[k]["ls"]) - float(apo[k]["ls"])
        assert abs(float(delta) - expected) < 1e-4, k
        checked += 1
    assert checked > 300


# ---------------------------------------------------------------------------
# The two-state layer, on the real conformational pair
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def apo():
    """1OMP -- maltose-free, the open conformation."""
    return bff.labelizer_read_structure(_path("1OMP.pdb"))


def _coordinate_only_model():
    """A model with no conservation term.

    No ConSurf grades ship for MalE, and under the published model a missing
    term makes the combined score `unavailable` for every position -- correctly,
    but it leaves nothing to pair. These three terms need only the coordinates.
    """
    model = bff.LabelizerParameterList()
    for tag, table in (("se", "N_SE11_MEAN_SURFACE_DIST"),
                       ("cr", "C_CR1_Name"),
                       ("ss", "C_SS1_SS")):
        model.append(bff.LabelizerParameter(tag, table, 1))
    return model


def test_the_apo_structure_scores_as_well_as_the_holo_one(apo):
    """Two independent structures of the same protein, same agreement.

    A DSSP that happened to suit one conformation would show up here.
    """
    reference = _reference("1OMP")
    table = dict(bff.labelizer_load_table("C_SS1_SS").by_key)
    by_value = {round(v, 5): k for k, v in table.items()}
    ours = bff.labelizer_dssp(apo)

    want, got = [], []
    for i, residue in enumerate(apo.residues):
        row = reference.get(residue.seq_id)
        if row is None:
            continue
        letter = by_value.get(round(float(row["ss"]), 5))
        if letter is None:
            continue
        want.append(letter)
        got.append(ours[i])
    agreement = sum(a == b for a, b in zip(want, got)) / len(want)
    assert agreement > 0.95, "apo agreement %.3f" % agreement


def test_the_pair_score_finds_the_hinge_closure():
    """What the two-state score is *for*, on the textbook case.

    Maltose-binding protein closes around its ligand by a hinge bend, so the
    pairs worth measuring are the ones spanning the two domains -- their
    dye-dye distance should *shrink* substantially between apo and holo. If the
    top-ranked pairs did not do that, the score would be ranking noise.
    """
    model = _coordinate_only_model()
    options = bff.LabelizerOptions()
    apo_scores = bff.labelizer_combined_by_key(
        bff.labelizer_score_structure(_path("1OMP.pdb"), model, options, ""))
    holo_scores = bff.labelizer_combined_by_key(
        bff.labelizer_score_structure(_path("1anf.pdb"), model, options, ""))
    assert len(apo_scores) == len(holo_scores) == N_RESIDUES

    fret = bff.LabelizerFRETOptions()
    fret.n_refine = 0
    pairs = list(bff.labelizer_pair_scores_two_states(
        _path("1OMP.pdb"), _path("1anf.pdb"), apo_scores, holo_scores, fret))
    assert len(pairs) > 10000

    top = pairs[:20]
    # Sorted best first, and every top pair must actually move.
    assert all(top[i].value >= top[i + 1].value for i in range(len(top) - 1))
    shifts = [p.distance_2 - p.distance for p in top]
    assert all(abs(d) > 5.0 for d in shifts), \
        "a top two-state pair that barely moves is a ranking failure: %r" % shifts
    # The closure is a contraction, not an expansion.
    assert sum(1 for d in shifts if d < 0) >= 18, \
        "maltose binding closes the cleft; distances should shrink: %r" % shifts

    # Both conformations must be reported, and they must differ.
    for p in top:
        assert p.distance > 0 and p.distance_2 > 0
        assert p.distance != p.distance_2


def test_a_structure_paired_with_itself_has_nothing_to_report():
    """The degenerate case, which must be zero rather than small.

    |E(d) - E(d)| is exactly zero, so every pair scores zero. A non-zero result
    here would mean the two conformations are not being read independently.
    """
    model = _coordinate_only_model()
    scores = bff.labelizer_combined_by_key(
        bff.labelizer_score_structure(_path("1anf.pdb"), model, bff.LabelizerOptions(), ""))
    fret = bff.LabelizerFRETOptions()
    fret.n_refine = 0
    pairs = list(bff.labelizer_pair_scores_two_states(
        _path("1anf.pdb"), _path("1anf.pdb"), scores, scores, fret))
    assert pairs
    assert max(p.value for p in pairs) == 0.0
    assert all(p.distance == p.distance_2 for p in pairs[:100])


def test_the_difference_map_is_antisymmetric_and_finds_the_hinge(apo, holo):
    """The map is `d_holo(i,j) - d_apo(i,j)`, so it is symmetric in (i, j) and
    zero on the diagonal, and its extremes are the domains that move."""
    first, flat = bff.labelizer_cbeta_difference_map(apo, holo)
    matrix = np.asarray(flat)
    n = int(round(np.sqrt(matrix.size)))
    assert n * n == matrix.size
    matrix = matrix.reshape(n, n)

    assert first == 1, "the map starts at the first residue present"
    assert np.allclose(np.diag(matrix), 0.0), "a residue does not move from itself"
    assert np.allclose(matrix, matrix.T), "the map must be symmetric"
    # The hinge is a real, large motion -- several Angstrom, not rounding.
    assert np.abs(matrix).max() > 8.0
    # And it is predominantly a closure.
    off = matrix[~np.eye(n, dtype=bool)]
    assert (off < 0).sum() > (off > 0).sum()


def test_the_accessible_volume_pair_distance_is_not_the_distance_of_the_means():
    """`RDAMeanE` is a property of the two clouds, not of their two means.

    This was wrong: the accessible-volume path built real clouds and then threw
    them away, reporting the distance between their mean positions. The
    reference's `effDistance` inverts the mean *efficiency*
    (`label_lib_functions.py:138`), which is a different number -- measured
    here at **2.2 to 3.0 A** apart on real sites, and in the same direction
    every time. Since the pair score peaks sharply at R = R0, that is a large
    error in the thing being ranked.
    """
    model = _coordinate_only_model()
    scores = bff.labelizer_combined_by_key(
        bff.labelizer_score_structure(_path("1anf.pdb"), model, bff.LabelizerOptions(), ""))

    measured = {}
    for distance_type in ("Rmp", "RDAMean", "RDAMeanE"):
        options = bff.LabelizerFRETOptions()
        options.n_refine = 4
        options.distance_type = distance_type
        pairs = [p for p in bff.labelizer_fret_pair_scores(_path("1anf.pdb"), scores, options)
                 if p.probe_model == bff.PROBE_MODEL_ACCESSIBLE_VOLUME]
        assert pairs, "no pair was refined with a real volume"
        measured[distance_type] = {(p.seq_id_1, p.seq_id_2): p.distance
                                   for p in pairs}

    shared = (set(measured["Rmp"]) & set(measured["RDAMean"])
              & set(measured["RDAMeanE"]))
    assert shared, "the three settings refined no pair in common"

    gaps = []
    for key in shared:
        rmp = measured["Rmp"][key]
        rda = measured["RDAMean"][key]
        r_e = measured["RDAMeanE"][key]
        # Jensen: averaging the distance between two distributions cannot be
        # less than the distance between their means. Equality is possible and
        # is not a failure -- a site so buried that its cloud collapses to a
        # point has no spread to average over, and `av_pair_statistics`
        # documents that an empty cloud yields Rmp for all three types.
        assert rda >= rmp - 1e-9, key
        # The 1/R^6 weighting of the efficiency average pulls R_E back below
        # <RDA>, towards the near side of the distribution.
        assert r_e <= rda + 1e-9, key
        gaps.append(rda - rmp)

    # The point of the test: for ordinary sites the gap is large enough that
    # reporting the wrong one is a real error, not a rounding difference.
    real = [g for g in gaps if g > 1e-6]
    assert len(real) >= max(1, len(gaps) // 2), (
        "almost every cloud was degenerate: %r" % gaps)
    assert max(real) > 1.0, "largest gap only %.2f A" % max(real)


def test_the_point_dye_models_ignore_the_distance_type():
    """A point has one distance to another point, whatever is asked for --
    which is what the reference's SIMPLE and GEBHARDT do too."""
    model = _coordinate_only_model()
    scores = bff.labelizer_combined_by_key(
        bff.labelizer_score_structure(_path("1anf.pdb"), model, bff.LabelizerOptions(), ""))

    got = []
    for distance_type in ("Rmp", "RDAMeanE"):
        options = bff.LabelizerFRETOptions()
        options.n_refine = 0
        options.probe_model = bff.PROBE_MODEL_ALPHA_CONE
        options.distance_type = distance_type
        pairs = list(bff.labelizer_fret_pair_scores(_path("1anf.pdb"), scores, options))[:50]
        got.append([(p.seq_id_1, p.seq_id_2, round(p.distance, 9)) for p in pairs])
    assert got[0] == got[1]
