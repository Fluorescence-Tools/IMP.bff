"""The side-chain rotamer family: the Dunbrack table inside the container.

The third rotamer family (PRD-118). Probes and spin labels are ensembles --
conformers with weights, which is what `.drot` stores; a side-chain library is
a *distribution over backbone conformation*, a table indexed by residue and by
the (phi, psi) bin. Same container, different `PtoKind`, which is the whole
claim the shared envelope makes: a new domain is a set of kinds.

Two properties are worth a test and both are exact rather than approximate.
**The conversion is byte-reversible** -- `write_dunbrack_bin` reproduces
FASPR's `dun2010bbdep.bin` bit for bit -- which is what lets the vendored
FASPR engine, which seeks in that file by byte offset, keep running against a
shipped container with its parity pin intact. And **packing from the shipped
container gives the same structure** as packing from the binary, byte for
byte.

The binary itself is not redistributed (its terms are academic-use; PRD-118),
so the tests that need it skip without it. The container ships, so everything
that needs only the container always runs.
"""

import hashlib
from pathlib import Path

import numpy as np
import pytest

import IMP.bff

REPO = Path(__file__).resolve().parent.parent.parent
#: The shipped family container -- this is what a user has.
SIDECHAINS = Path(IMP.bff.get_data_path("rotamer_library")) / "sidechains.drot.pto"
#: FASPR's binary, if this checkout happens to carry one.
DUNBRACK_BIN = REPO / "junk" / "FASPR" / "dun2010bbdep.bin"
#: Residues with chi angles; alanine and glycine have none.
ROTAMERIC = "RNDCQEHILKMFPSTWYV"


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.fixture(scope="module")
def container():
    if not SIDECHAINS.exists():
        pytest.skip(f"{SIDECHAINS} not present")
    return str(SIDECHAINS)


def test_the_catalog_lists_every_rotameric_residue(container):
    assert sorted(IMP.bff.drot_catalog(container)) == sorted(ROTAMERIC)


def test_a_bin_holds_the_rotamers_that_backbone_allows(container):
    """A helix-ish backbone gives lysine several rotamers, ordered by weight."""
    got = IMP.bff.read_dunbrack_rotamers(container, "K", -60.0, -45.0)
    assert got.n_chi == 4 and got.n_rotamers > 1
    chi = np.asarray(got.get_chi()).reshape(got.n_rotamers, got.n_chi)
    probability = np.asarray(got.get_probability())
    sigma = np.asarray(got.get_sigma()).reshape(got.n_rotamers, got.n_chi)

    assert np.all(np.diff(probability) <= 0), "the library orders by probability"
    assert probability.sum() <= 1.0 + 1e-6 and probability.min() >= 0.01
    assert np.all(np.abs(chi) <= 180.0)
    assert np.all(sigma > 0.0), "every chi carries a spread"


def test_the_cuts_are_the_callers(container):
    """Taking the bin whole gives at least as much as the default cuts."""
    default = IMP.bff.read_dunbrack_rotamers(container, "F", -60.0, -45.0)
    whole = IMP.bff.read_dunbrack_rotamers(container, "F", -60.0, -45.0, 0.0, 1.0)
    assert whole.n_rotamers >= default.n_rotamers


def test_backbone_actually_selects(container):
    """A different (phi, psi) is a different distribution, not the same one."""
    helix = IMP.bff.read_dunbrack_rotamers(container, "W", -60.0, -45.0, 0.0, 1.0)
    sheet = IMP.bff.read_dunbrack_rotamers(container, "W", -120.0, 130.0, 0.0, 1.0)
    a = np.asarray(helix.get_probability())
    b = np.asarray(sheet.get_probability())
    n = min(a.size, b.size)
    assert n and not np.allclose(a[:n], b[:n]), "the table is backbone-dependent"


def test_alanine_and_glycine_are_refused(container):
    for residue in ("A", "G", "X"):
        with pytest.raises(Exception):
            IMP.bff.read_dunbrack_rotamers(container, residue, -60.0, -45.0)


def test_every_residue_answers_at_every_corner_of_the_bin_grid(container):
    """The bin index wraps; the corners are where an off-by-one shows up."""
    for residue in ROTAMERIC:
        for phi in (-180.0, -175.0, 0.0, 175.0, 180.0):
            for psi in (-180.0, 180.0):
                got = IMP.bff.read_dunbrack_rotamers(container, residue, phi, psi,
                                                     0.0, 1.0)
                assert got.n_rotamers > 0, (residue, phi, psi)
                assert len(got.chi) == got.n_rotamers * got.n_chi


@pytest.mark.skipif(not DUNBRACK_BIN.exists(),
                    reason="dun2010bbdep.bin is not redistributed (PRD-118)")
def test_the_conversion_is_byte_reversible(tmp_path, container):
    """The container holds the binary, not a rendering of it.

    This is the property the FASPR parity pin rests on: the engine seeks in
    `dun2010bbdep.bin` by byte offset, so a container that reproduces it
    exactly can stand in for it and a container that merely reproduces its
    *numbers* cannot.
    """
    back = tmp_path / "dun2010bbdep.bin"
    IMP.bff.write_dunbrack_bin(container, str(back))
    assert _sha256(back) == _sha256(DUNBRACK_BIN)


@pytest.mark.skipif(not DUNBRACK_BIN.exists(),
                    reason="dun2010bbdep.bin is not redistributed (PRD-118)")
def test_writing_the_container_reproduces_the_shipped_one(tmp_path, container):
    """Re-converting is a no-op in a diff: nothing in the file is timestamped."""
    again = tmp_path / "sidechains.drot.pto"
    IMP.bff.write_dunbrack_library(str(DUNBRACK_BIN), str(again))
    assert _sha256(again) == _sha256(container)


@pytest.mark.skipif(not DUNBRACK_BIN.exists(),
                    reason="dun2010bbdep.bin is not redistributed (PRD-118)")
def test_packing_from_the_container_matches_packing_from_the_binary(tmp_path,
                                                                    container):
    """`faspr_pack` takes the shipped container and packs the same structure."""
    import sys
    sys.path.insert(0, str(REPO / "test" / "faspr"))
    from test_faspr_port import _clean_pdb

    source = REPO / "examples" / "structure" / "T4L" / "3GUN.pdb"
    if not source.exists():
        pytest.skip(f"{source} not present")
    clean = tmp_path / "clean.pdb"
    _clean_pdb(source, clean)

    from_bin = tmp_path / "from_bin.pdb"
    from_container = tmp_path / "from_container.pdb"
    IMP.bff.faspr_pack(str(clean), str(from_bin), str(DUNBRACK_BIN))
    IMP.bff.faspr_pack(str(clean), str(from_container), container)
    assert _sha256(from_container) == _sha256(from_bin)


# IMP runs every .py under test/ as a standalone script, and a file of bare
# pytest functions would import cleanly and exit 0 -- reporting success without
# running a single assertion. Hand the file to pytest explicitly so a failure
# here is a failure in ctest.
if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
