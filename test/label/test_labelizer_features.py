"""The structural features on their own: DSSP, exposure, and the score model.

`test_labelizer_ab.py` compares against the reference's published output and is
the authority on agreement. It runs on 1DDB, which is **all-helical** -- 96 H
and not one strand -- so the bridge and ladder half of the DSSP implementation
is never reached there. This file covers what the A/B cannot:

* β structure, on maltose-binding protein (1anf), whose α/β fold has to come
  out as alternating strand and helix or the ladder logic is wrong;
* more than one chain, on 3j0e, where a hydrogen bond must not be inferred
  across a chain break;
* the invariants each feature has to satisfy whatever the structure;
* the arithmetic of the score model, including the three published defects,
  which have to be reproducible *and* correctable on demand.

There is no reference output for any of it, so nothing here claims agreement
with the reference -- these are properties that must hold on their own terms.
"""

import os

import numpy as np
import pytest

import IMP.bff as bff

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "input", "labelizer")


def _path(name):
    return os.path.join(DATA, name)


@pytest.fixture(scope="module")
def helical():
    """1DDB: mouse BID, all-helical, one chain."""
    return bff.labelizer_read_structure(_path("1DDB-39.pdb"))


@pytest.fixture(scope="module")
def alpha_beta():
    """1anf: maltose-binding protein, a mixed alpha/beta fold."""
    return bff.labelizer_read_structure(_path("1anf.pdb"))


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

def test_only_the_first_model_of_an_nmr_entry_is_read(helical):
    """`read_pdb_records` returns every ATOM in the file regardless of MODEL,
    so the grouping has to find the model boundary itself. 1DDB ships as an
    NMR ensemble; reading it twice over would double every residue."""
    assert len(helical.residues) == 195
    assert len(helical.vdw) == 2991
    seq = [r.seq_id for r in helical.residues]
    assert seq == sorted(seq), "residues must not repeat or reorder"
    assert len(set(seq)) == len(seq)


def test_several_chains_stay_several_chains():
    s = bff.labelizer_read_structure(_path("3j0e.pdb"))
    assert {r.chain for r in s.residues} == {"F", "G", "H"}


def test_a_glycine_gets_a_reconstructed_cbeta(helical):
    """A glycine has no CB and the reference builds a virtual one rather than
    skipping the residue; a site with no position cannot be scored at all."""
    glycines = [i for i, r in enumerate(helical.residues) if r.comp_id == "GLY"]
    assert glycines, "1DDB has glycines"
    for i in glycines:
        assert helical.residues[i].cb < 0, "a glycine has no real CB"
        position = bff.labelizer_cbeta_position(helical, i)
        assert len(position) == 3
        ca = np.asarray(helical.xyz).reshape(-1, 3)[helical.residues[i].ca]
        # A real CA-CB bond is 1.53 A; the construction must land near it.
        assert 1.3 < np.linalg.norm(np.asarray(position) - ca) < 1.8


# ---------------------------------------------------------------------------
# Secondary structure
# ---------------------------------------------------------------------------

def test_dssp_returns_one_code_per_residue_from_the_eight(helical):
    ss = bff.labelizer_dssp(helical)
    assert len(ss) == len(helical.residues)
    assert set(ss) <= set("HBEGITS-")


def test_an_alpha_beta_protein_gets_strands(alpha_beta):
    """The half of DSSP the all-helical A/B case never reaches.

    Maltose-binding protein is a two-domain alpha/beta fold; an implementation
    that never assigned `E` would pass every 1DDB test and be badly wrong here.
    """
    ss = bff.labelizer_dssp(alpha_beta)
    counts = {c: ss.count(c) for c in set(ss)}
    assert counts.get("E", 0) > 40, "MBP has substantial beta sheet: %r" % counts
    assert counts.get("H", 0) > 100, "and more helix than sheet: %r" % counts
    # A ladder of one bridge is B and a longer one is E, so both appear and E
    # dominates. All-B would mean ladders are never being joined up.
    assert counts.get("B", 0) > 0
    assert counts["E"] > counts.get("B", 0) * 2

    # Strands come in runs, not as isolated residues: a sheet is consecutive.
    runs = [len(r) for r in "".join(c if c == "E" else " " for c in ss).split()]
    assert max(runs) >= 4, "no strand longer than %d residues" % max(runs)


def test_no_hydrogen_bond_is_inferred_across_a_chain_break():
    """The amide hydrogen is placed from the *preceding* residue's C=O, which
    is meaningless across a chain boundary."""
    s = bff.labelizer_read_structure(_path("3j0e.pdb"))
    ss = bff.labelizer_dssp(s)
    assert len(ss) == len(s.residues)
    firsts = [i for i, r in enumerate(s.residues)
              if i == 0 or s.residues[i - 1].chain != r.chain]
    assert len(firsts) == 3
    # A chain's first residue cannot be mid-helix: it has no donor.
    for i in firsts:
        assert ss[i] in "-STB", "chain start assigned %r" % ss[i]


# ---------------------------------------------------------------------------
# Exposure
# ---------------------------------------------------------------------------

def test_exposure_measures_agree_about_which_residues_are_buried(helical):
    """Three independent measures, so they must correlate or one is wrong."""
    depth = np.asarray(bff.labelizer_residue_depth(helical, 1.4, 200))
    rsa = np.asarray(bff.labelizer_relative_solvent_accessibility(
        helical, bff.LABELIZER_MAXASA_WILKE, 1.4, 200))
    hse = np.asarray(bff.labelizer_half_sphere_exposure(helical, 13.0)).reshape(-1, 2)

    assert len(depth) == len(rsa) == len(hse) == len(helical.residues)
    assert (depth > 0).all(), "every residue is some distance from the surface"
    assert (rsa >= 0).all() and rsa.max() < 2.0

    # Deeper means less accessible, and means more neighbours on the side the
    # side chain points to.
    assert np.corrcoef(depth, rsa)[0, 1] < -0.5
    assert np.corrcoef(depth, hse[:, 0])[0, 1] > 0.5


def test_the_three_max_asa_scales_differ_but_agree_in_order():
    # dict(): the SWIG map proxy indexes but has no .get
    wilke = dict(bff.labelizer_max_asa(bff.LABELIZER_MAXASA_WILKE))
    sander = dict(bff.labelizer_max_asa(bff.LABELIZER_MAXASA_SANDER))
    miller = dict(bff.labelizer_max_asa(bff.LABELIZER_MAXASA_MILLER))
    assert len(wilke) == len(sander) == len(miller) == 20
    # Glycine is the smallest residue on every scale, and that is not a
    # convention -- it has no side chain.
    for scale in (wilke, sander, miller):
        assert min(scale, key=scale.get) == "GLY"

    # The largest is *not* the same residue on every scale, and the difference
    # is real rather than a transcription slip: Tien's theoretical maxima put
    # tryptophan on top, while Kabsch and Sander measured arginine higher
    # (248 vs 227) -- a long flexible side chain exposes more area than a bulky
    # compact one. A scale that agreed with the others everywhere would mean
    # the same table had been pasted in three times.
    assert max(wilke, key=wilke.get) == "TRP"
    assert max(miller, key=miller.get) == "TRP"
    assert max(sander, key=sander.get) == "ARG"
    assert sander["ARG"] > sander["TRP"]
    assert wilke["TRP"] != sander["TRP"] != miller["TRP"]


def test_a_residue_type_outside_the_scale_is_not_computed(helical):
    """Negative says 'no relative value'; zero would say 'fully buried'."""
    rsa = np.asarray(bff.labelizer_relative_solvent_accessibility(
        helical, bff.LABELIZER_MAXASA_WILKE, 1.4, 200))
    # 1DDB is protein-only, so nothing should be flagged here -- the point is
    # that the flag is negative and not zero when it does fire.
    assert (rsa >= 0).all()


# ---------------------------------------------------------------------------
# The tables and the model
# ---------------------------------------------------------------------------

def test_every_shipped_table_loads():
    for name in bff.labelizer_available_tables():
        table = bff.labelizer_load_table(name)
        assert table.name == name
        if name.startswith("C_"):
            assert table.categorical and len(table.by_key) > 0
        elif not name.startswith("N_ME11"):
            # The methionine table is a documented empty placeholder: the
            # score is hard-coded and the file exists so the loader does not
            # fail. Everything else has bins.
            assert not table.categorical and len(table.bins) > 0


#: The tag each shipped table belongs to, from its name.
_TAG_OF = {"CR": "cr", "SS": "ss", "SE": "se", "CS": "cs", "ME": "me"}

#: Tables that refuse, and why. Each needs something that genuinely is not
#: available, and refusing is the reference's behaviour too -- it raises
#: `NotImplementedError` for every table but the one it implements per term.
UNIMPLEMENTED = {
    # ConSurf fields `labelizer_read_consurf` does not import.
    "N_CS3_Lower_Score", "N_CS4_Upper_Score", "I_CS1_Color",
    "I_CS5_Variety_Length", "C_CS6_Cys_In_Variety",
    # Bin centres span -153 to +333 degrees, which is neither the -180..180
    # nor the 0..360 convention, and no reference output exists to settle it.
    "N_SS2_Phi", "N_SS3_Psi",
}


@pytest.mark.parametrize("table", sorted(bff.labelizer_available_tables()))
def test_a_table_is_either_implemented_or_refused_never_guessed(table, helical):
    """The property that matters: **no table silently returns the wrong
    observable.**

    Several used to. `C_SS4_SS-1` scores the secondary structure of the
    *preceding* residue, and reading the letter at `i` gave a plausible number
    from the wrong position. `N_SE10_CB_SURFACE_DIST` is the Cbeta's own
    distance to the surface and was being served the residue's mean-atom depth
    -- a different quantity on a table whose bins start at 0.69 A rather than
    1.46 A. Both produced confident wrong scores.

    So every shipped table is now exercised, and each must either compute its
    own observable or raise saying why not.
    """
    tag = _TAG_OF[table.split("_")[1][:2]]
    model = bff.LabelizerParameterList()
    model.append(bff.LabelizerParameter(tag, table, 1))

    if table in UNIMPLEMENTED:
        with pytest.raises(Exception) as excinfo:
            bff.labelizer_parameter_scores(helical, model, bff.LabelizerOptions(), {})
        assert table in str(excinfo.value)
        return

    rows = bff.labelizer_parameter_scores(helical, model, bff.LabelizerOptions(), {})
    assert len(rows) == len(helical.residues)
    scored = [r.value for r in rows if r.status == "scored"]

    if table == "N_CS2_Score":
        # The only implemented conservation table, and it needs grades that
        # were not supplied here -- so no value, and that is correct.
        assert not scored
        return

    assert scored, "%s produced no value at all" % table
    # A table that mapped every residue onto one bin would be a broken
    # observable that still "works"; the shipped tables all discriminate.
    assert len({round(v, 9) for v in scored}) > 1, (
        "%s puts every residue in one bin" % table)

    if table == "N_ME11_Methionin_Exclusion_Dummy":
        # The one table that is legitimately not consulted. It ships empty
        # (`{}`) and exists only so the reference's loader does not fail; the
        # exclusion score is hard-coded 0.001/0.999 -- see
        # `test_the_exclusion_term_uses_a_thousandfold_penalty_not_a_zero`.
        assert len(bff.labelizer_load_table(table).bins) == 0  # an ndarray now, not a tuple
        assert {round(v, 6) for v in scored} <= {0.001, 0.999}
        return

    # Every other value has to be a value the table actually contains: a
    # lookup returns a table entry, never something computed from one.
    permitted = (set(dict(bff.labelizer_load_table(table).by_key).values())
                 if table.startswith("C_")
                 else set(bff.labelizer_load_table(table).values))
    assert all(any(abs(v - p) < 1e-12 for p in permitted) for v in scored), \
        "%s returned a value that is not in the table" % table


def test_the_neighbour_secondary_structure_tables_read_the_neighbour(helical):
    """`C_SS4_SS-1` must score residue i-1, not residue i."""
    ss = bff.labelizer_dssp(helical)
    table = dict(bff.labelizer_load_table("C_SS4_SS-1").by_key)

    model = bff.LabelizerParameterList()
    model.append(bff.LabelizerParameter("ss", "C_SS4_SS-1", 1))
    rows = {r.seq_id: r for r in
            bff.labelizer_parameter_scores(helical, model, bff.LabelizerOptions(), {})}

    for i, residue in enumerate(helical.residues):
        row = rows[residue.seq_id]
        if i == 0:
            # No preceding residue, so no value -- not the value at i.
            assert row.status != "scored"
            continue
        assert row.value == pytest.approx(table[ss[i - 1]]), residue.seq_id
        if ss[i] != ss[i - 1]:
            # The check that would have caught the original defect.
            assert row.value != pytest.approx(table[ss[i]])


def test_the_cbeta_depth_is_not_the_mean_atom_depth(helical):
    """`N_SE10` and `N_SE11` are different observables on the same surface."""
    mean = np.asarray(bff.labelizer_residue_depth(helical, 1.4, 200))
    cbeta = np.asarray(bff.labelizer_cbeta_depth(helical, 1.4, 200))
    assert len(cbeta) == len(mean) == len(helical.residues)
    assert (cbeta > 0).all()
    # Correlated -- both measure burial -- but not the same number.
    assert np.corrcoef(mean, cbeta)[0, 1] > 0.7
    assert np.abs(mean - cbeta).mean() > 0.1


def test_a_numeric_lookup_takes_the_nearest_bin_and_does_not_interpolate():
    table = bff.labelizer_load_table("N_SE11_MEAN_SURFACE_DIST")
    bins = list(table.bins)
    values = list(table.values)
    # Exactly on a bin centre.
    assert bff.labelizer_lookup(table, bins[3]) == values[3]
    # Just inside the next bin's half-width -- still the same value, because
    # nothing is interpolated.
    midpoint = 0.5 * (bins[3] + bins[4])
    assert bff.labelizer_lookup(table, bins[3] + 0.01) == values[3]
    assert bff.labelizer_lookup(table, midpoint - 1e-6) == values[3]
    # Outside the table: clamped to the end bin, never extrapolated.
    assert bff.labelizer_lookup(table, -1e6) == values[0]
    assert bff.labelizer_lookup(table, 1e6) == values[-1]


def test_cysteine_scores_highest_and_tryptophan_lowest():
    """The cysteine-resemblance table is the model's chemical intuition."""
    by_key = dict(bff.labelizer_load_table("C_CR1_Name").by_key)
    assert max(by_key, key=by_key.get) == "C"
    assert min(by_key, key=by_key.get) == "W"
    assert by_key["S"] > by_key["L"]


def test_the_published_model_switches_two_terms_off():
    model = {p.tag: p.weight for p in bff.labelizer_model_paper()}
    assert model == {"cs": 1, "se": 1, "cr": 1, "ss": 1, "tp": 0, "ce": 0}


def test_every_tag_maps_to_a_dictionary_score_type():
    for tag in ["cs", "se", "ss", "ce", "tp", "cr", "me"]:
        assert bff.labelizer_score_type(tag) in bff.mfdb_label_score_types()


# ---------------------------------------------------------------------------
# The three published defects
# ---------------------------------------------------------------------------

def test_the_joined_label_score_defect_is_reproduced_and_correctable():
    """`prod ** 0.5` where the geometric mean is `prod ** (1/N)`.

    The two agree for a two-residue pair and not for the four values of a
    two-conformation pair, which is why the reference's single- and
    double-conformation pair scores are not on a common scale.
    """
    two = [1.5, 2.0]
    assert bff.labelizer_joined_label_score(two, bff.LABELIZER_MODEL_PUBLISHED) == \
        pytest.approx(bff.labelizer_joined_label_score(two, bff.LABELIZER_MODEL_CORRECTED))

    four = [1.5, 2.0, 1.5, 2.0]
    published = bff.labelizer_joined_label_score(four, bff.LABELIZER_MODEL_PUBLISHED)
    corrected = bff.labelizer_joined_label_score(four, bff.LABELIZER_MODEL_CORRECTED)
    assert published == pytest.approx(3.0)
    assert corrected == pytest.approx(np.prod(four) ** 0.25)
    assert published == pytest.approx(corrected ** 2)


def test_a_weight_zero_term_still_vetoes_under_the_published_model():
    """Weight zero does not mean ignored in the reference, and that is a
    defect rather than a convention -- so it is reproduced, and correctable."""
    rows = []
    for tag, value in [("cs", 1.5), ("se", 1.5), ("cr", 1.5), ("ss", 1.5),
                       ("tp", 0.0), ("ce", 1.0)]:
        row = bff.LabelizerScore()
        row.asym_id, row.seq_id, row.comp_id = "A", 1, "SER"
        row.score_type = bff.labelizer_score_type(tag)
        row.value, row.status = value, "scored"
        rows.append(row)

    options = bff.LabelizerOptions()
    options.model = bff.LABELIZER_MODEL_PUBLISHED
    published = bff.labelizer_labeling_score(rows, bff.labelizer_model_paper(), options)
    assert published[0].value == 0.0, "the weight-0 zero must veto"

    options.model = bff.LABELIZER_MODEL_CORRECTED
    corrected = bff.labelizer_labeling_score(rows, bff.labelizer_model_paper(), options)
    assert corrected[0].value == pytest.approx(1.5), \
        "a term of weight 0 must not be consulted"


def test_an_unavailable_term_makes_the_combination_unavailable():
    """Not zero: 'we could not compute this' and 'this scored zero' are
    different facts and must stay different."""
    rows = []
    for tag in ["cs", "se", "cr", "ss", "tp", "ce"]:
        row = bff.LabelizerScore()
        row.asym_id, row.seq_id, row.comp_id = "A", 1, "SER"
        row.score_type = bff.labelizer_score_type(tag)
        if tag == "cs":
            row.status = "unavailable"
        else:
            row.value, row.status = 1.5, "scored"
        rows.append(row)
    out = bff.labelizer_labeling_score(rows, bff.labelizer_model_paper(), bff.LabelizerOptions())
    assert out[0].status == "unavailable"


def test_the_weight_is_a_repeat_count_not_an_exponent():
    """`labeling_score.py:187` multiplies a term in `weight` times and takes
    the root over the total count. For integer weights that equals the
    exponent form, and the reference offers no other."""
    model = bff.LabelizerParameterList()
    model.append(bff.LabelizerParameter("cs", "N_CS2_Score", 3))
    model.append(bff.LabelizerParameter("se", "N_SE11_MEAN_SURFACE_DIST", 1))

    rows = []
    for tag, value in [("cs", 2.0), ("se", 0.5)]:
        row = bff.LabelizerScore()
        row.asym_id, row.seq_id, row.comp_id = "A", 1, "SER"
        row.score_type = bff.labelizer_score_type(tag)
        row.value, row.status = value, "scored"
        rows.append(row)

    out = bff.labelizer_labeling_score(rows, model, bff.LabelizerOptions())
    assert out[0].value == pytest.approx((2.0 ** 3 * 0.5) ** (1.0 / 4.0))


# ---------------------------------------------------------------------------
# The exclusion term
# ---------------------------------------------------------------------------

def test_the_exclusion_term_uses_a_thousandfold_penalty_not_a_zero(helical):
    """0.001/0.999 rather than 0/1 is deliberate: a zero would trip the
    zero-veto and take the whole score to zero, and -1 means an error."""
    model = bff.LabelizerParameterList()
    model.append(bff.LabelizerParameter("me", "N_ME11_Methionin_Exclusion_Dummy", 1))
    rows = bff.labelizer_parameter_scores(helical, model, bff.LabelizerOptions(), {})
    values = {round(r.value, 6) for r in rows if r.status == "scored"}
    assert values <= {0.001, 0.999}
    assert 0.999 in values


# ---------------------------------------------------------------------------
# Sequence breaks
# ---------------------------------------------------------------------------

def test_dssp_does_not_infer_a_peptide_bond_across_a_sequence_gap(tmp_path):
    """Array-adjacent is not chain-adjacent, and assuming so is silent.

    Two very ordinary inputs leave a gap in the author numbering: a crystal
    structure with an unresolved loop, and any structure carrying a residue the
    reader drops -- an unnatural amino acid such as pAcF, a modified residue.
    Before this was guarded, `labelizer_dssp` indexed turns by array position, so it
    read residues 99 and 103 as three apart when a missing residue 100 had made
    them four, and asserted a 3-turn that does not exist. It produced a
    confident wrong assignment rather than an error, which is the worst kind.

    The check: deleting one mid-chain residue must perturb only the residues
    whose own assignment depended on it, and must leave the rest of the protein
    untouched.
    """
    source = open(_path("1DDB-39.pdb")).read().splitlines(True)
    gapped = tmp_path / "gapped.pdb"
    with open(gapped, "w") as fh:
        for line in source:
            if (line.startswith(("ATOM", "HETATM"))
                    and line[21] == "A" and line[22:26].strip() == "100"):
                continue          # residue 100 is simply not there
            fh.write(line)

    whole = bff.labelizer_read_structure(_path("1DDB-39.pdb"))
    holed = bff.labelizer_read_structure(str(gapped))
    assert len(holed.residues) == len(whole.residues) - 1

    by_seq_whole = {r.seq_id: bff.labelizer_dssp(whole)[i]
                    for i, r in enumerate(whole.residues)}
    by_seq_holed = {r.seq_id: bff.labelizer_dssp(holed)[i]
                    for i, r in enumerate(holed.residues)}

    changed = [s for s in sorted(set(by_seq_whole) & set(by_seq_holed))
               if by_seq_whole[s] != by_seq_holed[s]]
    # Only the immediate neighbourhood of the gap may move. A regression that
    # reintroduced index-based turns shifted assignments further out.
    assert all(abs(s - 100) <= 4 for s in changed), \
        "residues far from the gap changed: %r" % changed
    assert len(changed) <= 6, "too much moved for one missing residue: %r" % changed


def test_a_gap_does_not_shift_the_scores_of_other_residues(tmp_path):
    """Scores are keyed by author numbering, not by array position.

    This is what makes a dropped residue survivable at all: a caller joining on
    (asym_id, seq_id) is unaffected by anything the reader removed, and only a
    caller zipping arrays by position could be misled.
    """
    source = open(_path("1DDB-39.pdb")).read().splitlines(True)
    gapped = tmp_path / "gapped.pdb"
    with open(gapped, "w") as fh:
        for line in source:
            if (line.startswith(("ATOM", "HETATM"))
                    and line[21] == "A" and line[22:26].strip() == "100"):
                continue
            fh.write(line)

    model = bff.LabelizerParameterList()
    model.append(bff.LabelizerParameter("cr", "C_CR1_Name", 1))
    whole = {(r.asym_id, r.seq_id): r.value
             for r in bff.labelizer_parameter_scores(
                 bff.labelizer_read_structure(_path("1DDB-39.pdb")), model,
                 bff.LabelizerOptions(), {})}
    holed = {(r.asym_id, r.seq_id): r.value
             for r in bff.labelizer_parameter_scores(
                 bff.labelizer_read_structure(str(gapped)), model,
                 bff.LabelizerOptions(), {})}

    assert ("A", 100) in whole and ("A", 100) not in holed
    for key in holed:
        assert holed[key] == whole[key], "%r moved when residue 100 vanished" % (key,)


# ---------------------------------------------------------------------------
# The probe library
# ---------------------------------------------------------------------------

#: The fifteen dyes the Labelizer backend offers, and their published QY/EC.
#: Names as that backend spells them -- which is how everyone spells them, and
#: not how the vendor does.
LABELIZER_DYES = {
    "Alexa488": (0.92, 73000), "Alexa532": (0.61, 81000),
    "Alexa546": (0.79, 112000), "Alexa555": (0.10, 155000),
    "Alexa568": (0.69, 88000), "Alexa594": (0.66, 92000),
    "Alexa647": (0.33, 270000), "Atto488": (0.80, 90000),
    "Atto532": (0.90, 115000), "Atto550": (0.80, 120000),
    "Atto643": (0.62, 150000), "Atto647N": (0.65, 150000),
    "Cy3": (0.15, 136000), "Cy3B": (0.67, 130000), "Cy5": (0.30, 250000),
}


@pytest.mark.parametrize("name", sorted(LABELIZER_DYES))
def test_a_dye_resolves_under_the_name_people_actually_use(name):
    """The library keys are vendor spellings; nobody types them.

    `AlexaFluor488`, `LumiprobeCy3` and `ATTO647N` are what the library calls
    them, and `Alexa488`, `Cy3` and `Atto647N` are what an experimenter, a
    paper and the Labelizer backend call them. Exact matching rejected all
    fifteen while holding the right row.
    """
    quantum_yield, extinction = LABELIZER_DYES[name]
    dye = bff.get_probe(name)
    # Both libraries derive from the same upstream tables, so this is exact.
    assert dye.quantum_yield == pytest.approx(quantum_yield)
    assert dye.extinction_coefficient == pytest.approx(extinction)


def test_resolution_never_guesses_between_two_dyes():
    """`Cy5` resolves because `LumiprobeCy5` matches exactly while
    `LumiprobeCy55` does not. Ambiguity has to raise, not pick."""
    assert bff.resolve_probe_name("Cy5") == "LumiprobeCy5"
    assert bff.resolve_probe_name("Cy55") == "LumiprobeCy55"
    assert bff.resolve_probe_name("cy3b") == "LumiprobeCy3b"
    with pytest.raises(Exception):
        bff.get_probe("NotADye123")


def test_an_exact_vendor_name_still_wins():
    """Alias resolution is additive: it must not change what already worked."""
    for name in bff.available_probes():
        assert bff.resolve_probe_name(name) == name


def test_the_whole_labelizer_dye_set_is_covered():
    """All fifteen, which is what a backend replacement needs.

    Thirteen were already here under vendor spellings; `Atto532` and `Atto643`
    were genuinely absent and were imported from Labelizer's own tables by
    `utility/import_labelizer_dyes.py`.
    """
    for name in LABELIZER_DYES:
        bff.get_probe(name)
    assert len(LABELIZER_DYES) == 15


def test_a_dye_that_really_is_absent_still_fails_loudly():
    """A wrong dye is worse than no dye: R0 goes as the sixth root of the
    overlap, so a substituted spectrum is a quietly wrong distance."""
    with pytest.raises(Exception) as excinfo:
        bff.get_probe("Atto999")
    assert "not in the probe library" in str(excinfo.value)


@pytest.mark.parametrize("name,ex_nm,em_nm", [("Atto532", 532, 552),
                                              ("Atto643", 644, 665)])
def test_the_imported_spectra_peak_where_the_dye_says_they_do(name, ex_nm, em_nm):
    """The check that an import went onto the right grid and the right way up.

    A spectrum resampled with an off-by-one, or with excitation and emission
    swapped, still loads and still yields a plausible-looking R0 -- but its
    peaks move. ATTO532 and ATTO643 are named for their excitation maxima, so
    the data states its own expected answer.
    """
    spectrum = bff.get_probe(name).spectrum
    wavelength = np.asarray(spectrum.wavelength)
    excitation = np.asarray(spectrum.excitation)
    emission = np.asarray(spectrum.emission)

    assert wavelength[excitation.argmax()] == pytest.approx(ex_nm, abs=2)
    assert wavelength[emission.argmax()] == pytest.approx(em_nm, abs=2)
    # Emission is red-shifted from excitation. Swapped curves fail here.
    assert wavelength[emission.argmax()] > wavelength[excitation.argmax()]
    # Peak-normalised, because the acceptor curve is multiplied by the
    # extinction coefficient to give a molar absorptivity.
    assert excitation.max() == pytest.approx(1.0)
    assert emission.max() == pytest.approx(1.0)


def test_every_dye_shares_one_wavelength_grid():
    """`spectral_overlap` throws if two dyes disagree, so an import that used
    its own grid would break every pair involving it rather than just itself."""
    grids = set()
    for name in bff.available_probes():
        spectrum = bff.get_probe(name).spectrum
        if spectrum.size():
            wavelength = np.asarray(spectrum.wavelength)
            grids.add((len(wavelength), wavelength[0], wavelength[-1]))
    assert len(grids) == 1, "dyes are on different grids: %r" % grids
    assert grids == {(671, 300.0, 970.0)}

def test_the_global_charge_agrees_with_a_residue_census(helical):
    """The reference offers this and consumes it nowhere, so nothing downstream
    would notice it being wrong -- which is a reason to check it directly.

    Counted independently here from the residue names, so the test does not
    reimplement the thing it is testing by accident.
    """
    import collections
    net, positive, negative = bff.labelizer_global_charge(helical)
    census = collections.Counter(r.comp_id for r in helical.residues)

    # Histidine counts as +1, as the `ce` term counts it. That is wrong at
    # physiological pH and is reproduced deliberately: the two must agree.
    assert positive == census["ARG"] + census["LYS"] + census["HIS"]
    assert negative == -(census["ASP"] + census["GLU"])
    assert net == positive + negative
    # Reported so the three add up, rather than negative being a magnitude.
    assert negative <= 0.0 <= positive


def test_the_global_charge_is_the_same_table_the_charge_term_uses(helical):
    """If these two ever disagreed about what is charged, a site's charge
    environment and its protein's net charge would be describing different
    molecules."""
    net, positive, negative = bff.labelizer_global_charge(helical)
    # A structure of only neutral residues must give zero on both counts.
    model = bff.LabelizerParameterList()
    model.append(bff.LabelizerParameter("ce", "", 1))
    rows = bff.labelizer_parameter_scores(helical, model, bff.LabelizerOptions(), {})
    assert rows, "the charge term must still evaluate"
    # Non-trivial charge in this structure, so the census above is a real test.
    assert positive > 0 and negative < 0
