"""The FPS strip: obstacles removed around the attachment site.

Pinned here because the strip is what fps.json positions are calibrated
against -- without it a declared ``allowed_sphere_radius`` of a few
Angstrom is walled in by the attachment residue's own atoms and the AV
comes back empty. The strip is expressed as the PyMOL selection dialect
fps documents carry (``chain A and resid 36 and not name N+CA+C+O+CB``);
anything outside the dialect is refused loudly, never silently ignored.
"""

import numpy as np
import pytest

import IMP.bff as av

_PDB = """\
ATOM      1  N   ILE A   3      11.104   6.134  -6.504  1.00  0.00           N
ATOM      2  CA  ILE A   3      10.000   6.873  -6.071  1.00  0.00           C
ATOM      3  C   ILE A   3       9.104   6.208  -5.005  1.00  0.00           C
ATOM      4  O   ILE A   3       9.400   5.105  -4.551  1.00  0.00           O
ATOM      5  CB  ILE A   3      10.522   8.304  -5.757  1.00  0.00           C
ATOM      6  CG1 ILE A   3       9.548   9.247  -6.470  1.00  0.00           C
ATOM      7  CG2 ILE A   3      11.829   8.570  -6.467  1.00  0.00           C
ATOM      8  CD1 ILE A   3       9.977  10.637  -6.154  1.00  0.00           C
ATOM      9  N   LYS A   4       8.100   6.971  -4.512  1.00  0.00           N
ATOM     10  CA  LYS A   4       7.200   6.500  -3.500  1.00  0.00           C
END
"""


@pytest.fixture(name="pdb")
def _pdb(tmp_path):
    path = tmp_path / "mini.pdb"
    path.write_text(_PDB)
    return str(path)


def _kept_names(path, chain="A", resi="3"):
    return [
        line[12:16].strip()
        for line in open(path)
        if line.startswith(("ATOM  ", "HETATM"))
        and line[21] == chain
        and line[22:26].strip() == resi
    ]


def test_the_default_strip_removes_the_side_chain_not_the_backbone(pdb):
    """Backbone + attachment stay; the rest of the side chain goes."""
    assert av.default_strip_mask("A", 3, "CB") == (
        "chain A and resid 3 and not name N+CA+C+O+CB"
    )
    stripped = av.stripped_pdb_for(pdb, "A", 3, "CB")
    assert _kept_names(stripped) == ["N", "CA", "C", "O", "CB"]


def test_a_declared_mask_is_honoured_as_given(pdb):
    stripped = av.stripped_pdb_for(
        pdb, "A", 3, "CB", "chain A and resid 3 and not name CA+CB+C+N+O"
    )
    assert _kept_names(stripped) == ["N", "CA", "C", "O", "CB"]
    # A whole-residue mask leaves only the attachment atom.
    whole = av.stripped_pdb_for(pdb, "A", 3, "CB", "chain A and resid 3")
    assert _kept_names(whole) == ["CB"]


def test_the_attachment_atom_survives_any_mask(pdb):
    """The attachment is resolved from this file; it is never stripped."""
    survived = av.stripped_pdb_for(pdb, "A", 3, "CB", "name CB")
    assert _kept_names(survived) == ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"]


@pytest.mark.parametrize(
    "mask",
    [
        "name CA CB C N O",  # space-separated lists: a parse error in PyMOL
        "chain A or resid 3",
        "within 5 of name CB",
        "chain A and resid 3 and not name N+CA+C+O+CB and name X",
        "chain A and chain B",
    ],
)
def test_a_mask_outside_the_dialect_is_refused_loudly(pdb, mask):
    with pytest.raises(ValueError, match="strip_mask"):
        av.stripped_pdb_for(pdb, "A", 3, "CB", mask)


def test_a_mask_never_silently_selects_nothing(pdb):
    """A well-formed mask naming a absent residue strips nothing -- but the
    attachment atom is still kept and the file is still valid."""
    stripped = av.stripped_pdb_for(
        pdb, "A", 3, "CB", "chain A and resid 99 and name N"
    )
    assert _kept_names(stripped) == ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"]


def test_the_strip_cache_answers_the_same_file_per_site_and_mask(pdb):
    assert av.stripped_pdb_for(pdb, "A", 3, "CB") == av.stripped_pdb_for(
        pdb, "A", 3, "CB"
    )
    other = av.stripped_pdb_for(pdb, "A", 3, "CB", "chain A and resid 3")
    assert other != av.stripped_pdb_for(pdb, "A", 3, "CB")
