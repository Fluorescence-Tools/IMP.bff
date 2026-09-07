"""The pair layer on a structure with more than one chain.

Everything else in `test/label/` runs on a single-chain protein, and that is
exactly where the two defects this file guards against could hide.

**The collision.** `ll_pair_scores_two_states` matched sites across the two
conformations on the residue number alone -- `std::map<int, std::size_t>` keyed
on `seq_id`, the reference's convention (`fret_score.py:512`). On a multimer
residue 100 of chains F, G and H are one key, the last one inserted wins, and
an inter-chain pair therefore gets the *same* site on both sides of the second
state: `distance_2` is 0, `E(0)` is 1, and the score is `jls * |E(d1) - 1|` for
every such pair. Silently, with no status and no warning.

**The missing filter.** There was no way to say "donor in this chain, acceptor
in that one". Without it a homodimer cannot be screened at all: a single
cysteine mutation puts the same residue on *both* protomers, so the only
distance an experiment can measure is i(A)-i(B), and those pairs were drowned
in intra-protomer ones -- or, on a file holding two copies of the assembly,
mixed with lattice contacts that are not biology. The reference Python had the
argument (`fret_score.py:479`, filter at `:728`); the port did not.

3j0e carries three chains (F, G, H), which is enough to show both.
"""

import os

import numpy as np
import pytest

import IMP.bff as bff

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "input", "labelizer")
MULTIMER = os.path.join(DATA, "3j0e.pdb")


def _coordinate_only_model():
    """The three terms that need nothing but the coordinates.

    No ConSurf grades ship for 3j0e, and under the published model a missing
    term makes the whole combined score unavailable -- which would leave
    nothing to pair.
    """
    model = bff.LlParameterList()
    for tag, table in (("se", "N_SE11_MEAN_SURFACE_DIST"),
                       ("cr", "C_CR1_Name"),
                       ("ss", "C_SS1_SS")):
        model.append(bff.LlParameter(tag, table, 1))
    return model


@pytest.fixture(scope="module")
def label_scores():
    return bff.ll_combined_by_key(
        bff.ll_score_structure(MULTIMER, _coordinate_only_model(),
                               bff.LlOptions(), ""))


@pytest.fixture(scope="module")
def chains(label_scores):
    return sorted({key[0] for key in label_scores})


def _fret(**kwargs):
    fret = bff.LlFretOptions()
    fret.n_refine = 0            # the analytic cone; this is about bookkeeping
    for name, value in kwargs.items():
        setattr(fret, name, value)
    return fret


def test_the_fixture_really_has_several_chains(chains):
    assert len(chains) >= 2, \
        "this file tests multi-chain behaviour; %r is not multi-chain" % chains


def test_a_multimer_paired_with_itself_is_exactly_zero(label_scores):
    """The regression test for the collision.

    Pairing a structure against itself must give |E(d) - E(d)| = 0 for every
    pair. Before the fix the inter-chain pairs came back with `distance_2 == 0`
    -- their partner in the second state having been overwritten by another
    chain's residue of the same number -- and scored `jls * |E(d1) - 1|`, which
    is large. The single-chain version of this test (1anf, in
    `test_labelizer_ab_male.py`) passed throughout.
    """
    pairs = list(bff.ll_pair_scores_two_states(
        MULTIMER, MULTIMER, label_scores, label_scores, _fret()))
    assert pairs, "nothing was paired at all"

    inter = [p for p in pairs if p.asym_id_1 != p.asym_id_2]
    assert inter, "no inter-chain pair was enumerated, so this proves nothing"

    assert max(p.value for p in pairs) == 0.0, \
        "a structure paired with itself has nothing to report"
    assert all(p.distance == p.distance_2 for p in pairs)
    assert all(p.distance_2 > 0 for p in inter), \
        "distance_2 == 0 is the collision: the site was matched to itself"


def test_the_inter_protomer_pair_of_one_residue_is_enumerated(label_scores,
                                                              chains):
    """i(A)-i(B): the only pair a single-cysteine homodimer can produce.

    The two sites have the same residue number and different chains, so they
    are exactly what the seq_id-keyed map could not represent.
    """
    first, second = chains[0], chains[1]
    pairs = list(bff.ll_pair_scores_two_states(
        MULTIMER, MULTIMER, label_scores, label_scores,
        _fret(donor_chain=first, acceptor_chain=second)))
    same_residue = [p for p in pairs if p.seq_id_1 == p.seq_id_2]
    assert same_residue, "no i(A)-i(B) pair survived"
    for p in same_residue:
        assert p.asym_id_1 == first and p.asym_id_2 == second
        assert p.distance > 0 and p.distance_2 > 0


def test_the_chain_filter_keeps_only_the_requested_protomer_pair(label_scores,
                                                                 chains):
    first, second = chains[0], chains[1]
    fret = _fret(donor_chain=first, acceptor_chain=second)
    filtered = list(bff.ll_pair_scores(MULTIMER, label_scores, fret))
    everything = list(bff.ll_pair_scores(MULTIMER, label_scores, _fret()))

    assert filtered, "the filter removed everything"
    assert len(filtered) < len(everything), "the filter removed nothing"
    for p in filtered:
        assert {p.asym_id_1, p.asym_id_2} == {first, second}, \
            "a pair outside the requested chains survived: %s" % p
    if len(chains) > 2:
        assert any(chains[2] in (p.asym_id_1, p.asym_id_2)
                   for p in everything), \
            "the unfiltered run should reach the third chain"


def test_an_unset_chain_filter_changes_nothing(label_scores):
    """The default has to be what the module did before the fields existed."""
    fret = _fret()
    assert fret.donor_chain == "" and fret.acceptor_chain == ""
    pairs = list(bff.ll_pair_scores(MULTIMER, label_scores, fret))
    assert any(p.asym_id_1 != p.asym_id_2 for p in pairs), \
        "with no filter, inter-chain pairs must still be produced"
    assert any(p.asym_id_1 == p.asym_id_2 for p in pairs), \
        "with no filter, intra-chain pairs must still be produced"


def test_the_chain_map_renames_the_second_conformation(label_scores, chains):
    """Two files that name the same protomer differently still match up.

    Mapping every chain onto itself must be identical to no map at all; mapping
    a chain onto a different one must change the answer.
    """
    first, second = chains[0], chains[1]
    plain = list(bff.ll_pair_scores_two_states(
        MULTIMER, MULTIMER, label_scores, label_scores,
        _fret(donor_chain=first, acceptor_chain=second)))

    identity = bff.MapStringString()
    for chain in chains:
        identity[chain] = chain
    mapped = list(bff.ll_pair_scores_two_states(
        MULTIMER, MULTIMER, label_scores, label_scores,
        _fret(donor_chain=first, acceptor_chain=second, chain_map=identity)))
    assert len(mapped) == len(plain)
    assert [p.distance_2 for p in mapped] == [p.distance_2 for p in plain]

    swap = bff.MapStringString()
    swap[first] = second
    swap[second] = first
    swapped = list(bff.ll_pair_scores_two_states(
        MULTIMER, MULTIMER, label_scores, label_scores,
        _fret(donor_chain=first, acceptor_chain=second, chain_map=swap)))
    assert swapped, "the swap matched nothing"
    assert [p.distance_2 for p in swapped] != [p.distance_2 for p in plain], \
        "reading the second state through swapped chains changed no distance"


def test_the_cbeta_map_is_meaningless_on_a_multimer_without_a_chain(chains):
    """The map is indexed by residue number, so it describes one chain.

    Left to itself on a multimer every chain lands on the same row and the last
    one wins -- the diagonal is still zero (a residue does not move from
    itself) but the map is a mixture of three proteins. Naming the chain is the
    only reading that means anything, and it must not be silently equal to the
    unnamed one.
    """
    structure = bff.ll_read_structure(MULTIMER)
    first_all, flat_all = bff.ll_cbeta_difference_map(structure, structure)
    first_one, flat_one = bff.ll_cbeta_difference_map(structure, structure,
                                                      chains[0])
    for flat in (flat_all, flat_one):
        matrix = np.asarray(flat)
        n = int(round(np.sqrt(matrix.size)))
        assert n * n == matrix.size
        # Self-comparison: every entry is zero either way.
        assert np.allclose(matrix, 0.0)
    assert first_one >= first_all, \
        "one chain cannot start before the union of all of them"


def test_a_two_state_pair_survives_the_container(label_scores, chains, tmp_path):
    """`distance_2` must come back out of a `.mmfdb.pto` unchanged.

    The container wrote the field and never read it, so a two-state pair round
    tripped through one came back with `distance_2 == 0` -- which is not a
    missing value but a *plausible* one: zero distance is efficiency 1, and the
    score built on it looks perfectly reasonable. Nothing caught it because the
    only writer was the single-conformation path, where the field is NaN and
    correctly omitted.

    The single-state case is checked here too: absent must round-trip to NaN,
    not to zero.
    """
    import math

    first, second = chains[0], chains[1]
    fret = _fret(donor_chain=first, acceptor_chain=second)
    two_state = list(bff.ll_pair_scores_two_states(
        MULTIMER, MULTIMER, label_scores, label_scores, fret))[:50]
    assert two_state

    path = str(tmp_path / "two_state.mmfdb.pto")
    bff.ll_write_pto(path, MULTIMER, bff.LlScoreList(), two_state, "{}")
    back = list(bff.ll_read_pto_pairs(path))
    assert len(back) == len(two_state)
    for wrote, read in zip(two_state, back):
        assert read.seq_id_1 == wrote.seq_id_1
        assert read.seq_id_2 == wrote.seq_id_2
        assert read.asym_id_1 == wrote.asym_id_1
        assert read.distance == pytest.approx(wrote.distance, abs=1e-9)
        assert read.distance_2 == pytest.approx(wrote.distance_2, abs=1e-9)

    single = list(bff.ll_pair_scores(MULTIMER, label_scores, fret))[:20]
    path = str(tmp_path / "single.mmfdb.pto")
    bff.ll_write_pto(path, MULTIMER, bff.LlScoreList(), single, "{}")
    for read in bff.ll_read_pto_pairs(path):
        assert math.isnan(read.distance_2), \
            "one conformation has no second distance; that must not read as 0"
