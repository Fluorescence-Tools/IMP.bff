"""The native Labelizer against the reference implementation's own output.

The reference (`labelizer`, Gebhardt *et al.*, Nat. Commun. 16, 3305, 2025)
ships a worked example -- 1DDB model 39, chain A, 195 residues -- together with
the per-residue parameter scores it produced, at full double precision. Those
CSVs are the only independent check on this port that exists, and this file is
it.

What is pinned, and why each bar sits where it does:

* **`cr` (cysteine resemblance) and `ss` (secondary structure): exact.** Both
  are a table lookup, so anything short of exact would be a transcription
  error. `ss` being exact is the load-bearing result of the port: the reference
  shells out to the DSSP binary and refuses to run on macOS at all, and
  :func:`IMP.bff.ll_dssp` reproduces its eight-state assignment on all 195
  residues with no disagreement.

* **`cs` (conservation): the lookup is pinned exactly; the shipped example is
  not reproducible and that is the reference's defect, not this port's.** The
  example passes ``prot1_cs=[".../1DDB-39_cs.pdb"]`` and labelizer's
  ``_save_pdb`` then writes its *output* to that same path, so re-running the
  example feeds the previous run's scores back in as conservation grades. The
  lookup has a two-cycle -- 1.63 -> 2.3553..., 2.36 -> 1.6322..., which round
  back to 2.36 and 1.63 -- so the shipped CSV has exactly **two** distinct
  values across 195 residues where the table has ten bins. The two-cycle is
  pinned below to sixteen digits, which verifies the reader and the lookup
  without pretending the science behind the file is recoverable.

* **`se` (solvent exposure): a stated tolerance, because it cannot be exact.**
  The published term is MSMS residue depth and the MSMS binaries labelizer
  ships are 32-bit ppc/i386 Mach-O that do not execute. :func:`ll_residue_depth`
  builds the surface natively instead. Measured against the reference on 1DDB:
  no bias (mean +0.015 A), Pearson r = 0.976, and 92% of residues within one
  table bin. The bars below are set from that distribution and are deliberately
  a little looser than the measurement, so ordinary numerical drift does not
  fail the suite while a real regression does.
"""

import csv
import os

import numpy as np
import pytest

import IMP.bff as bff

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "input", "labelizer")

#: The reference's two-letter tag -> the MMFDB `score_type` it is written as.
TAGS = {
    "cr": "cysteine_resemblance",
    "cs": "conservation",
    "se": "solvent_exposure",
    "ss": "secondary_structure",
}

#: Residues in the shipped example.
N_RESIDUES = 195


def _path(name):
    return os.path.join(DATA, name)


def _reference(tag):
    """The reference parameter scores, keyed `"<chain><seq_id>"`."""
    with open(_path("1DDB-39_%s.csv" % tag)) as fh:
        return {row[0]: float(row[1]) for row in list(csv.reader(fh))[1:]}


@pytest.fixture(scope="module")
def structure():
    return bff.ll_read_structure(_path("1DDB-39.pdb"))


@pytest.fixture(scope="module")
def scored():
    """Every parameter score of the example, as `{score_type: {key: value}}`.

    Conservation is read from `1DDB-39_cs.pdb` -- the file the reference run
    actually consumed, feedback loop and all -- because reproducing the
    reference means feeding it what the reference was fed.
    """
    rows = bff.ll_score_structure(
        _path("1DDB-39.pdb"),
        bff.ll_model_paper(),
        bff.LlOptions(),
        _path("1DDB-39_cs.pdb"),
    )
    out = {}
    for r in rows:
        key = bff.ll_residue_key(r.asym_id, r.seq_id)
        out.setdefault(r.score_type, {})[key] = (
            r.value if r.status == "scored" else None
        )
    return out


def _deltas(scored, tag):
    ref = _reference(tag)
    ours, refs = [], []
    for key, value in ref.items():
        got = scored[TAGS[tag]].get(key)
        if got is None:
            continue
        ours.append(got)
        refs.append(value)
    assert len(ours) == N_RESIDUES, "every residue must be scored"
    return np.abs(np.asarray(ours) - np.asarray(refs))


def test_the_structure_is_the_one_the_reference_scored(structure):
    assert len(structure.residues) == N_RESIDUES
    assert {r.chain for r in structure.residues} == {"A"}
    assert structure.residues[0].comp_id == "MET"


@pytest.mark.parametrize("tag", ["cr", "ss"])
def test_the_table_lookup_terms_are_exact(scored, tag):
    """A lookup that is not exact is a transcription error, not a tolerance."""
    delta = _deltas(scored, tag)
    assert delta.max() == 0.0, "%s disagrees with the reference" % tag


def test_the_native_dssp_reproduces_the_reference_assignment(structure):
    """All 195 residues, eight states, no disagreement.

    The reference's `ss` score is a lookup on the DSSP letter, so inverting the
    published score recovers the letter the reference's DSSP binary assigned
    and the two assignments can be compared directly.
    """
    table = bff.ll_load_table("C_SS1_SS")
    by_value = {round(table.by_key[k], 9): k for k in table.by_key.keys()}
    reference = _reference("ss")

    ours = bff.ll_dssp(structure)
    assert len(ours) == N_RESIDUES

    disagreements = []
    for i, residue in enumerate(structure.residues):
        key = bff.ll_residue_key(residue.chain, residue.seq_id)
        want = by_value[round(reference[key], 9)]
        if ours[i] != want:
            disagreements.append((key, want, ours[i]))
    assert not disagreements, "DSSP letters differ: %r" % (disagreements[:10],)


def test_the_conservation_lookup_is_exact_to_the_last_digit():
    """The two-cycle the shipped example is stuck in, pinned both ways.

    This is what verifies the conservation term. The example itself cannot be
    reproduced -- see the module docstring -- but the map that produced its two
    values is exactly this lookup, and reproducing it to sixteen digits leaves
    nowhere for an error in the reader or the binning to hide.
    """
    table = bff.ll_load_table("N_CS2_Score")
    assert bff.ll_lookup(table, 1.63) == 2.3553071957924936
    assert bff.ll_lookup(table, 2.36) == 1.6322095472510827


def test_the_shipped_conservation_reference_is_degenerate():
    """Guards the claim in the module docstring, so it cannot rot silently.

    If someone replaces the fixture with a real ConSurf run this fails, and the
    right response is to delete this test and tighten `cs` to exact -- not to
    widen it.
    """
    values = set(_reference("cs").values())
    assert len(values) == 2, (
        "the shipped 1DDB conservation reference has two distinct values "
        "across 195 residues because the example overwrites its own input; "
        "got %d" % len(values)
    )
    assert len(set(_reference("cr").values())) == 20, "cr is healthy by contrast"
    assert len(set(_reference("se").values())) == 10, "se is healthy by contrast"


def test_residue_depth_agrees_with_msms_without_bias(structure):
    """The native surface against MSMS, as a distribution rather than a value.

    The reference depth is only recoverable to half a bin (the score is binned
    at 0.268 A), so the comparison is against bin centres and the bars are set
    accordingly.
    """
    table = bff.ll_load_table("N_SE11_MEAN_SURFACE_DIST")
    by_value = {round(table.values[i], 9): table.bins[i]
                for i in range(len(table.bins))}
    reference = _reference("se")

    depth = np.asarray(bff.ll_residue_depth(structure, 1.4, 590))
    ours, refs = [], []
    for i, residue in enumerate(structure.residues):
        key = bff.ll_residue_key(residue.chain, residue.seq_id)
        centre = by_value.get(round(reference[key], 9))
        if centre is None:
            continue
        ours.append(min(depth[i], 4.0))
        refs.append(centre)
    ours, refs = np.asarray(ours), np.asarray(refs)
    residual = ours - refs

    assert abs(residual.mean()) < 0.10, "a bias would be correctable; fix it"
    assert np.corrcoef(ours, refs)[0, 1] > 0.95
    assert (np.abs(residual) < 0.268).mean() > 0.85, "within one table bin"


def test_the_solvent_exposure_score_lands_in_the_reference_bin(scored):
    """The consequence of the depth agreement, on the quantity that is used."""
    delta = _deltas(scored, "se")
    exact = float((delta < 1e-9).mean())
    assert exact > 0.65, "bin-exact fraction regressed: %.3f" % exact
    assert delta.mean() < 0.20


def test_the_combined_score_is_produced_for_every_residue(scored):
    combined = scored["combined"]
    assert len(combined) == N_RESIDUES
    values = [v for v in combined.values() if v is not None]
    assert len(values) == N_RESIDUES
    # A likelihood-ratio product is unbounded above but must stay non-negative.
    assert min(values) >= 0.0
    assert max(values) > 1.0
