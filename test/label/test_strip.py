"""The one strip engine (PRD-106): grammar, PDB lines, obstacle arrays, hierarchies.

``fret.av`` strips PDB text for the AV build and cgdye strips an IMP
hierarchy before attaching an explicit dye; both go through
``IMP.bff.label`` so the same mask selects the same atoms whatever the
representation.
"""

import numpy as np
import pytest

import IMP
import IMP.atom

import IMP.bff.label as strip

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
ATOM     11  N   ILE B   3      21.104   6.134  -6.504  1.00  0.00           N
ATOM     12  CB  ILE B   3      20.522   8.304  -5.757  1.00  0.00           C
END
"""
_LINES = _PDB.splitlines(keepends=True)


def _names(lines, chain, resi):
    return [l[12:16].strip() for l in lines
            if l.startswith("ATOM") and l[21] == chain and int(l[22:26]) == resi]


def test_grammar():
    sel = strip.parse_strip_mask("chain A and resid 3 and not name N+CA+C+O")
    assert (sel.get_chain() == "A" and list(sel.get_resids()) == [3]
            and sel.get_negate())
    assert sel.matches("A", 3, "CB") and not sel.matches("A", 3, "CA")
    assert not sel.matches("B", 3, "CB") and not sel.matches("A", 4, "CB")
    sel = strip.parse_strip_mask("resi 3 and name CB")
    # An unstated chain is empty rather than `None`: a C++ string has no
    # third state, and "any chain" is what both spellings meant.
    assert sel.get_chain() == "" and sel.matches("B", 3, "cb")
    for bad in ("chain A or resid 3", "(resid 3)", "resid 3 and name CA CB", "resname HOH"):
        with pytest.raises(ValueError):
            strip.parse_strip_mask(bad)


def test_masks_are_spelled_in_the_dialect():
    assert strip.default_strip_mask("A", 36, "CB") == "chain A and resid 36 and not name N+CA+C+O+CB"
    assert strip.site_strip_mask("", 481, ("N", "CA", "C", "O", "OXT")) == "resid 481 and not name N+CA+C+O+OXT"
    # both round-trip through the parser
    strip.parse_strip_mask(strip.default_strip_mask("A", 36, "CB"))
    strip.parse_strip_mask(strip.site_strip_mask("A", 481, ("N", "CA", "C", "O", "OXT")))


def test_pdb_lines_and_attachment_protection():
    kept = strip.strip_pdb_lines(_LINES, "chain A and resid 3 and not name N+CA+C+O")
    assert _names(kept, "A", 3) == ["N", "CA", "C", "O"]
    assert _names(kept, "A", 4) == ["N", "CA"] and _names(kept, "B", 3) == ["N", "CB"]
    kept = strip.strip_pdb_lines(_LINES, "chain A and resid 3", "A", 3, "CB")
    assert _names(kept, "A", 3) == ["CB"]
    kept = strip.strip_pdb_lines(_LINES, "resid 3 and name CB")  # any chain
    assert _names(kept, "A", 3) == ["N", "CA", "C", "O", "CG1", "CG2", "CD1"]
    assert _names(kept, "B", 3) == ["N"]


def test_obstacle_arrays():
    records = []
    atoms = []
    for l in _LINES:
        if l.startswith("ATOM"):
            records.append((l[21].strip(), int(l[22:26]), l[12:16].strip(),
                            float(l[30:38]), float(l[38:46]), float(l[46:54]), 1.7))
            atoms.append([float(l[30:38]), float(l[38:46]), float(l[46:54]), 1.7])
    atoms = np.array(atoms)
    out = strip.strip_obstacles(atoms, records, "chain A and resid 3 and not name N+CA+C+O+CB")
    assert out.shape == (atoms.shape[0] - 3, 4)
    with pytest.raises(ValueError):
        strip.strip_obstacles(atoms[:-1], records, "resid 3")


def test_hierarchy_non_destructive_and_in_place(tmp_path):
    pdb = tmp_path / "mini.pdb"
    pdb.write_text(_PDB)
    m = IMP.Model()
    hier = IMP.atom.read_pdb(str(pdb), m, IMP.atom.AllPDBSelector())
    n_before = len(IMP.atom.get_by_type(hier, IMP.atom.ATOM_TYPE))

    sel = strip.select_atoms(hier, "chain A and resid 3 and not name N+CA+C+O")
    assert sorted(a.get_atom_type().get_string() for a in sel) == ["CB", "CD1", "CG1", "CG2"]

    stripped, n = strip.strip_hierarchy(hier, "chain A and resid 3 and not name N+CA+C+O")
    assert n == 4
    assert len(IMP.atom.get_by_type(hier, IMP.atom.ATOM_TYPE)) == n_before  # untouched
    assert len(IMP.atom.get_by_type(stripped, IMP.atom.ATOM_TYPE)) == n_before - 4

    same, n = strip.strip_hierarchy(hier, "resid 3 and name CB", inplace=True)
    assert n == 2 and same == hier
    assert len(IMP.atom.get_by_type(hier, IMP.atom.ATOM_TYPE)) == n_before - 2


def test_cgdye_and_av_defaults_differ_only_in_cb(tmp_path):
    """The two consumers own their keep-sets; the engine is shared (PRD-106)."""
    from IMP.bff.label import SITE_KEEP_ATOM_NAMES, strip_sidechain_at_site
    pdb = tmp_path / "mini.pdb"
    pdb.write_text(_PDB)
    m = IMP.Model()
    hier = IMP.atom.read_pdb(str(pdb), m, IMP.atom.AllPDBSelector())
    assert "CB" not in SITE_KEEP_ATOM_NAMES
    removed = strip_sidechain_at_site(hier, "A", 3)
    assert removed == 4  # CB CG1 CG2 CD1
    keep_cb = strip_sidechain_at_site(
        IMP.atom.read_pdb(str(pdb), m, IMP.atom.AllPDBSelector()), "A", 3,
        keep_atom_names=SITE_KEEP_ATOM_NAMES + ("CB",))
    assert keep_cb == 3
    av_kept = strip.strip_pdb_lines(_LINES, strip.default_strip_mask("A", 3, "CB"))
    assert _names(av_kept, "A", 3) == ["N", "CA", "C", "O", "CB"]


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
