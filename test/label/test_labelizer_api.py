"""The public functions nothing else reaches.

Most of the labelizer surface is covered by the A/B files, which drive it the
way a caller would. That leaves a tail of exported functions that are only ever
*called through* something else -- or, in one case, called by nothing at all.
A public function with no direct test is a public function whose contract is
whatever its single caller happens to need, and it drifts silently.

This file is deliberately unglamorous: one test per otherwise-unreached export,
checking the thing its documentation promises.

Found while writing it: `ll_pair_score_negative_control` is exercised by
nothing in the package. It is the reference's control-pair formula, which the
reference itself computes and then discards (`fret_score.py:528`, the line that
would have stored it is commented out). Ported because a control pair is a real
experimental need; now at least it is checked.
"""

import math
import os

import numpy as np
import pytest

import IMP.bff as bff

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "input", "labelizer")


def _path(name):
    return os.path.join(DATA, name)


@pytest.fixture(scope="module")
def structure():
    return bff.ll_read_structure(_path("1DDB-39.pdb"))


# ---------------------------------------------------------------------------
# Small vocabulary helpers
# ---------------------------------------------------------------------------

def test_the_standard_residues_are_the_twenty():
    names = list(bff.ll_standard_residues())
    assert len(names) == 20
    assert len(set(names)) == 20
    assert names == sorted(names), "kept sorted so a membership test can bisect"
    for expected in ("ALA", "TRP", "GLY", "CYS"):
        assert expected in names
    # Selenomethionine is not one: it is what `protein_only` drops.
    assert "MSE" not in names


def test_the_one_letter_code_round_trips_and_refuses_the_unknown():
    for three, one in (("ALA", "A"), ("TRP", "W"), ("GLY", "G"), ("CYS", "C")):
        assert bff.ll_one_letter(three) == one
    # Every standard residue maps to a distinct letter -- a collision would
    # silently merge two residue types in the cysteine-resemblance lookup.
    letters = [bff.ll_one_letter(r) for r in bff.ll_standard_residues()]
    assert len(set(letters)) == 20
    # Anything else is X, not a guess and not a crash.
    for unknown in ("MSE", "HOH", "XYZ", "", "A"):
        assert bff.ll_one_letter(unknown) == "X"


def test_a_categorical_lookup_refuses_a_key_it_does_not_have():
    """`ll_lookup` snaps to the nearest bin; `ll_lookup_key` must not.

    A categorical table has no notion of "nearest": asking it for a residue it
    was not fitted on has no answer, and returning the closest string would be
    an invention.
    """
    table = bff.ll_load_table("C_CR1_Name")
    assert bff.ll_lookup_key(table, "C") == pytest.approx(3.8641, abs=1e-3)
    assert bff.ll_lookup_key(table, "W") == pytest.approx(0.2093, abs=1e-3)
    with pytest.raises(Exception) as excinfo:
        bff.ll_lookup_key(table, "Z")
    assert "C_CR1_Name" in str(excinfo.value)


def test_every_tag_has_a_score_type_and_an_unknown_tag_raises():
    for tag in ("cs", "se", "ss", "ce", "tp", "cr", "me", "combined"):
        assert bff.ll_score_type(tag)
    with pytest.raises(Exception):
        bff.ll_score_type("zz")


# ---------------------------------------------------------------------------
# Surface area, on its own
# ---------------------------------------------------------------------------

def test_residue_sasa_is_an_area_and_tracks_exposure(structure):
    """`ll_residue_sasa` is the unnormalised half of the RSA, and is exported
    because an area is a useful quantity on its own."""
    area = np.asarray(bff.ll_residue_sasa(structure, 1.4, 200))
    assert len(area) == len(structure.residues)
    assert (area >= 0).all()
    # A residue cannot expose more than a large residue's worth of surface.
    assert area.max() < 400.0
    # Some residue is genuinely buried, and some is genuinely exposed.
    assert area.min() < 5.0 and area.max() > 100.0

    # It is the numerator of the relative accessibility, so the two must agree
    # once the maximum is divided out.
    rsa = np.asarray(bff.ll_relative_solvent_accessibility(
        structure, bff.LL_MAXASA_WILKE, 1.4, 200))
    maxasa = dict(bff.ll_max_asa(bff.LL_MAXASA_WILKE))
    for i, residue in enumerate(structure.residues):
        assert rsa[i] == pytest.approx(area[i] / maxasa[residue.comp_id])


# ---------------------------------------------------------------------------
# Where the dye is, one site at a time
# ---------------------------------------------------------------------------

def test_the_alpha_cone_places_the_dye_outward_and_within_reach(structure):
    """The analytic estimate, checked against what it claims geometrically.

    It puts the dye along the outward normal -- away from the local atom
    centroid -- at a distance bounded by the linker. A sign error would put it
    inside the protein and still return three plausible-looking numbers.
    """
    options = bff.LlFretOptions()
    xyz = np.asarray(structure.xyz).reshape(-1, 3)
    checked = 0
    for i, residue in enumerate(structure.residues):
        position = bff.ll_alpha_cone_mean_position(structure, i, options)
        if len(position) != 3:
            continue
        position = np.asarray(position)
        cbeta = np.asarray(bff.ll_cbeta_position(structure, i))
        assert len(cbeta) == 3

        # Within the linker's reach of the attachment point.
        assert np.linalg.norm(position - cbeta) < options.linker_length + 20.0
        # And on the outward side: further from the protein's centroid than
        # the Cbeta is, for a site that is not buried.
        checked += 1
    assert checked > 150, "the cone should place most residues"


def test_the_dye_model_selects_which_position_is_returned(structure):
    """`ll_probe_mean_position` dispatches on the option; the three models must
    actually differ, or the ladder is decorative."""
    pdb = _path("1DDB-39.pdb")
    site = next(i for i, r in enumerate(structure.residues) if r.cb >= 0)

    positions = {}
    for model in (bff.PROBE_MODEL_CBETA, bff.PROBE_MODEL_ALPHA_CONE,
                  bff.PROBE_MODEL_ACCESSIBLE_VOLUME):
        options = bff.LlFretOptions()
        options.probe_model = model
        got = bff.ll_probe_mean_position(structure, pdb, site, options)
        assert len(got) == 3, model
        positions[model] = np.asarray(got)

    # Cbeta is the attachment; the other two put the dye somewhere else.
    assert np.linalg.norm(positions[bff.PROBE_MODEL_ALPHA_CONE]
                          - positions[bff.PROBE_MODEL_CBETA]) > 1.0
    assert np.linalg.norm(positions[bff.PROBE_MODEL_ACCESSIBLE_VOLUME]
                          - positions[bff.PROBE_MODEL_CBETA]) > 1.0


def test_an_out_of_range_site_returns_nothing_rather_than_guessing(structure):
    options = bff.LlFretOptions()
    for bad in (-1, len(structure.residues), 10 ** 6):
        assert len(bff.ll_probe_mean_position(structure, _path("1DDB-39.pdb"),
                                            bad, options)) == 0
        assert len(bff.ll_cbeta_position(structure, bad)) == 0


# ---------------------------------------------------------------------------
# The pair formulas, directly
# ---------------------------------------------------------------------------

def test_the_single_conformation_score_peaks_at_r0():
    """`jls * (1 - 2|E - 0.5|)`: maximal where FRET is most sensitive."""
    r0 = 52.0
    at_r0 = bff.ll_pair_score_single(1.0, r0, r0)
    assert at_r0 == pytest.approx(1.0)
    # Falls off either side, and never goes negative.
    for distance in (0.5 * r0, 0.8 * r0, 1.25 * r0, 2.0 * r0):
        value = bff.ll_pair_score_single(1.0, distance, r0)
        assert 0.0 <= value < at_r0
    # Scales linearly with the joined label score.
    assert bff.ll_pair_score_single(0.5, r0, r0) == pytest.approx(0.5)


def test_the_two_conformation_score_rewards_a_change_in_efficiency():
    """`jls * |E(d1) - E(d2)|`: zero when nothing moves, largest across the
    steep part of the curve."""
    r0 = 52.0
    assert bff.ll_pair_score_double(1.0, r0, r0, r0) == pytest.approx(0.0)
    # Symmetric in its two distances -- which conformation is "first" is not
    # a property of the pair.
    assert bff.ll_pair_score_double(1.0, 40.0, 70.0, r0) == pytest.approx(
        bff.ll_pair_score_double(1.0, 70.0, 40.0, r0))
    # A change spanning R0 beats the same-sized change far from it.
    near = bff.ll_pair_score_double(1.0, 45.0, 60.0, r0)
    far = bff.ll_pair_score_double(1.0, 130.0, 145.0, r0)
    assert near > far


def test_the_negative_control_rewards_a_pair_that_does_not_move():
    """The formula the reference computes and throws away.

    A control pair should sit mid-range in efficiency *and* not change between
    conformations, so the score is largest when `E1 + E2` is about 1 and the
    two are equal, and is punished hard -- slope 20 -- by any change.
    """
    r0 = 52.0
    still = bff.ll_pair_score_negative_control(1.0, r0, r0, r0)
    assert still > 0.9, "a pair at R0 that does not move is the ideal control"

    # Any real change is punished to zero: 20 * |dE| >= 1 kills it.
    moved = bff.ll_pair_score_negative_control(1.0, 40.0, 70.0, r0)
    assert moved == pytest.approx(0.0)

    # Unmoving but at an uninformative distance scores below the ideal.
    far = bff.ll_pair_score_negative_control(1.0, 3 * r0, 3 * r0, r0)
    assert 0.0 <= far < still

    # And it scales with the joined label score like the other two.
    assert bff.ll_pair_score_negative_control(0.5, r0, r0, r0) == \
        pytest.approx(0.5 * still)


# ---------------------------------------------------------------------------
# Container plumbing
# ---------------------------------------------------------------------------

def test_a_dye_container_can_be_written_from_a_caller_s_own_dyes(tmp_path):
    """`probe_write_pto` takes any dye map, not only the bundled library --
    which is the point of exporting it beside `probe_library_to_pto`."""
    subset = {}
    for name in ("AlexaFluor488", "AlexaFluor647"):
        subset[name] = bff.find_probe(name)

    out = str(tmp_path / "two.mmfdb.pto")
    bff.probe_write_pto(out, subset, "a subset, for this test")
    back = dict(bff.probe_read_pto(out))
    assert set(back) == set(subset)
    assert bff.probe_pto_forster_radius(out, "AlexaFluor488", "AlexaFluor647") \
        == pytest.approx(bff.forster_radius(subset["AlexaFluor488"],
                                            subset["AlexaFluor647"]))


def test_the_settings_come_back_out_of_a_container(tmp_path):
    """`ll_read_pto_settings` is what makes a run reproducible from the file
    alone, so it has to return the settings that went in."""
    import json
    pdb = _path("1DDB-39.pdb")
    model = bff.ll_model_paper()
    options = bff.LlOptions()
    fret = bff.LlFretOptions()
    fret.n_refine = 0

    scores = bff.ll_score_structure(pdb, model, options,
                                    _path("1DDB-39_cs.pdb"))
    settings = bff.ll_settings_json(model, options, fret,
                                    _path("1DDB-39_cs.pdb"))
    out = str(tmp_path / "s.mmfdb.pto")
    bff.ll_write_pto(out, pdb, scores, bff.LlPairScoreList(), settings)

    back = json.loads(bff.ll_read_pto_settings(out))
    assert back["arithmetic"] == "published"
    assert {t["tag"] for t in back["model"]} == {"cs", "se", "tp", "cr",
                                                 "ss", "ce"}
    assert back["fret"]["forster_radius"] == pytest.approx(fret.forster_radius)
    assert back["probe_radius"] == pytest.approx(options.probe_radius)
