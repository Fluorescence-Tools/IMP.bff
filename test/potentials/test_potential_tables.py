"""`data/potentials.pto`: one container for every parameter table.

The tables arrived as four loose `.npy` files in a ChiSurf directory. They ship
as one PTO container, which is what this module already does for
`dyes.drot.pto` and `dyes.mmfdb.pto`, and it is built by
`bin/imp_bff_potentials2pto`.

What these tests hold to: the container reads, its manifest describes what is
in it, the two text tables are what `IMP.core.StatisticalPairScore` reads
directly, and the conversion's decisions -- the NaNs, the short-range
divergence, the sentinel in the proline map -- are the ones the manifest claims.
"""

import json

import numpy as np
import pytest

import IMP
import IMP.atom
import IMP.container
import IMP.core
import IMP.bff

NAMES = ["mj", "unres", "hbond", "ramachandran"]


@pytest.fixture(scope="module")
def manifest():
    return json.loads(IMP.bff.read_potential_manifest())


def test_the_container_is_where_it_says_it_is():
    path = IMP.bff.get_potential_container_path()
    assert path.endswith("potentials.pto")
    assert sorted(IMP.bff.potential_table_names()) == sorted(NAMES)


def test_every_table_is_described(manifest):
    assert set(manifest["tables"]) == set(NAMES)
    for name, entry in manifest["tables"].items():
        assert entry["kind"] in ("pot.pmf", "pot.grid")
        assert entry["source"].endswith(".npy")
        assert len(entry["source_sha256"]) == 64
        assert entry["note"], f"{name} says nothing about its conversion"


def test_the_two_contact_tables_are_pmf_text():
    """Miyazawa-Jernigan and UNRES are text in IMP's PMF format, so the
    potential *is* the file: nothing in this module loops over them."""
    for name in ("mj", "unres"):
        t = IMP.bff.read_potential_table(name)
        assert t.kind == "pot.pmf"
        assert not len(t.table)
        header = t.text.split("\n")[0].split()
        assert len(header) == 2
        bin_width, n_types = float(header[0]), int(header[1])
        assert bin_width > 0
        # As wide as the *potential*, not as wide as IMP. An
        # `IMP.atom.ResidueType` table grows when a process meets a residue it
        # has not seen -- read one PDB with a ligand in it and the next index
        # is one further out -- and an index past the side of a square table
        # reads outside it. These are indexed by `IMP.bff.ResidueContactType`,
        # whose names are the ones the table names and no others.
        assert n_types == 20
        assert t.text.count("\n") == 1 + n_types * (n_types + 1) // 2


def test_the_two_grids_have_the_shape_the_manifest_claims(manifest):
    for name in ("hbond", "ramachandran"):
        t = IMP.bff.read_potential_table(name)
        assert t.kind == "pot.grid"
        assert list(t.shape) == manifest["tables"][name]["shape"]
        assert t.table.size == t.get_size()
        assert np.isfinite(t.table).all(), "a NaN in a table poisons an energy"


def test_the_hydrogen_bond_divergence_is_held_flat(manifest):
    """`hb.npy` reaches 7e33 below 1.27 A, where no two atoms are. Those bins
    hold the value at 1.3 A instead."""
    t = IMP.bff.read_potential_table("hbond")
    g = t.table.reshape(*t.shape)
    first = int(round(1.3 / manifest["tables"]["hbond"]["bin_width"]))
    for channel in range(g.shape[0]):
        assert len(np.unique(g[channel, :first])) == 1
        assert g[channel, 0] == g[channel, first]
    assert np.abs(g).max() < 1e6


def test_the_ramachandran_maps_are_log_ratios():
    """Zero is the peak of the distribution and nothing is above it -- which
    is what lets `ramachandran_energy` read the map as an energy and negate
    it. The proline channel's sentinel 1.0 is gone."""
    t = IMP.bff.read_potential_table("ramachandran")
    g = t.table.reshape(*t.shape)
    assert g.shape == (3, 360, 360)
    for channel in range(3):
        assert g[channel].max() == 0.0
        assert g[channel].min() < 0.0


def test_asking_for_a_table_that_is_not_there_says_so():
    with pytest.raises(IOError):
        IMP.bff.read_potential_table("nonesuch")


def test_a_container_of_ones_own_round_trips(tmp_path):
    """A potential this module does not ship is a file, not a patch."""
    grid = IMP.bff.PotentialTable("mine", "pot.grid")
    grid.shape = [2, 3]
    grid.values = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    text = IMP.bff.PotentialTable("words", "pot.pmf")
    text.text = "1.0 2\nALA ALA 0.5 0.5\n"

    path = str(tmp_path / "own.pto")
    IMP.bff.write_potential_tables(path, [grid, text],
                                   json.dumps({"tables": {"mine": {
                                       "kind": "pot.grid", "shape": [2, 3]}}}))
    assert sorted(IMP.bff.potential_table_names(path)) == ["mine", "words"]
    back = IMP.bff.read_potential_table("mine", path)
    np.testing.assert_allclose(back.table, [1, 2, 3, 4, 5, 6])
    assert list(back.shape) == [2, 3]
    assert IMP.bff.read_potential_table("words", path).text == text.text


def test_a_payload_that_compresses_hugely_still_reads(tmp_path):
    """The decompressed size is written into the payload, because the one-shot
    decoder reports a buffer that is too small as an *error* -- so a reader
    that guesses the size from the compressed length (`size * 8`) fails on
    anything that compresses better than eight times. The UNRES text is 2.4 MB
    of digits in a fraction of that."""
    t = IMP.bff.PotentialTable("flat", "pot.grid")
    t.shape = [500000]
    t.values = [0.0] * 500000          # compresses by ~1000x
    path = str(tmp_path / "flat.pto")
    IMP.bff.write_potential_tables(path, [t], "{}")
    back = IMP.bff.read_potential_table("flat", path)
    assert back.table.size == 500000
    assert not back.table.any()


# ---------------------------------------------------------------------------
# The shipped tables, in the scores that read them
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def structure():
    m = IMP.Model()
    h = IMP.atom.read_pdb(IMP.bff.get_example_path("structure/T4L/3GUN.pdb"),
                          m, IMP.atom.NonWaterNonHydrogenPDBSelector())
    return m, h


def test_miyazawa_jernigan_scores_a_real_fold(structure):
    """A folded protein's contacts are favourable, so the score is negative
    and large; the default constructor reads the shipped table."""
    m, h = structure
    typed = IMP.bff.add_residue_type_score_data(h, IMP.atom.AT_CB)
    assert len(typed) > 100
    lsc = IMP.container.ListSingletonContainer(
        m, [p.get_index() for p in typed])
    r = IMP.container.PairsRestraint(
        IMP.bff.MiyazawaJerniganPairScore(),
        IMP.container.ClosePairContainer(lsc, 6.5, 0.0))
    score = IMP.core.RestraintsScoringFunction([r]).evaluate(False)
    assert score < -100.0
    # and nothing is scored past the cutoff
    far = IMP.container.PairsRestraint(
        IMP.bff.MiyazawaJerniganPairScore(1e-3),
        IMP.container.ClosePairContainer(lsc, 6.5, 0.0))
    assert IMP.core.RestraintsScoringFunction([far]).evaluate(False) == 0.0


def test_unres_scores_a_real_fold(structure):
    m, h = structure
    typed = IMP.bff.add_residue_type_score_data(h, IMP.atom.AT_CB)
    lsc = IMP.container.ListSingletonContainer(
        m, [p.get_index() for p in typed])
    r = IMP.container.PairsRestraint(
        IMP.bff.UNRESCentroidPairScore(),
        IMP.container.ClosePairContainer(lsc, 15.0, 0.0))
    assert IMP.core.RestraintsScoringFunction([r]).evaluate(False) < 0.0


def test_the_residue_type_space_stays_closed(structure):
    """A residue the tables say nothing about is not typed at all, so no index
    can land outside the table -- however wide IMP's own residue types grow."""
    m, h = structure
    IMP.atom.ResidueType("XYZZY")        # widen IMP's key table
    typed = IMP.bff.add_residue_type_score_data(h, IMP.atom.AT_CB)
    key = IMP.bff.get_residue_type_key()
    assert typed
    assert max(p.get_value(key) for p in typed) < 20


def test_the_unres_table_is_one_value_per_unordered_pair(manifest):
    """The source is filled in its upper triangle only, and the reference indexed
    it by the residues' order in the chain -- so about one contact in two read
    the empty half and scored zero. A pair potential has one value per
    unordered pair, and that is what the file holds.

    See `okf/validation/hbond_ca_cutoff.md`.
    """
    assert manifest["tables"]["unres"]["triangle"] == "upper"
    t = IMP.bff.read_potential_table("unres")
    rows = {}
    for line in t.text.split("\n")[1:]:
        parts = line.split()
        if len(parts) < 3:
            continue
        rows[(parts[0], parts[1])] = parts[2:]
    # one row per unordered pair, never both orders
    for a, b in rows:
        assert (b, a) not in rows or a == b


def test_the_unres_table_repels_at_short_range(manifest):
    """The Python applied a flat penalty below `min_dist` in its loop; here it
    is in the bins, so the behaviour is the data's."""
    t = IMP.bff.read_potential_table("unres")
    entry = manifest["tables"]["unres"]
    below = int(round(entry["min_dist_A"] / entry["bin_width"]))
    for line in t.text.split("\n")[1:]:
        parts = line.split()
        if len(parts) < 3 or parts[0] != "ALA" or parts[1] != "ALA":
            continue
        bins = [float(v) for v in parts[2:]]
        assert all(v == entry["repulsion"] for v in bins[:below])
        assert any(v != entry["repulsion"] for v in bins[below:])
        break
    else:
        pytest.fail("the table has no ALA-ALA row")


def test_the_ramachandran_map_scores_a_real_backbone(structure):
    m, h = structure
    t = IMP.bff.read_potential_table("ramachandran")
    r = IMP.bff.RamachandranRestraint(m, h, t.table, t.shape[0], t.shape[1])
    score = r.unprotected_evaluate(None)
    n_res = len(IMP.atom.get_by_type(h, IMP.atom.RESIDUE_TYPE))
    assert 0.0 < score < 10.0 * n_res, "every residue costs less than an empty cell"


def test_the_hydrogen_bond_term_needs_hydrogens(structure):
    """A crystal structure read without hydrogens has no amide H, so there is
    nothing to donate and the term is zero. That is not a defect: the term
    scores H bonds, and this structure carries none."""
    m, h = structure
    t = IMP.bff.read_potential_table("hbond")
    r = IMP.bff.HydrogenBondRestraint(m, h, t.table, t.shape[1])
    assert r.unprotected_evaluate(None) == 0.0
    assert r.get_n_hbonds() == 0
