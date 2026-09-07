"""The one strip engine (PRD-106): grammar, PDB lines, obstacle arrays, hierarchies.

``fret.av`` strips PDB text for the AV build and cgprobe strips an IMP
hierarchy before attaching an explicit dye; both go through
``IMP.bff.label`` so the same mask selects the same atoms whatever the
representation.
"""

import numpy as np
import pytest

import IMP
import IMP.atom

import IMP.bff as strip

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


def _atom(chain, resi, name):
    import IMP.bff
    a = IMP.bff.SelectionAtom()
    a.chain, a.resi, a.name = chain, resi, name
    return a


def test_grammar():
    """A mask is a selection expression, and `matches` answers one atom."""
    sel = strip.parse_strip_mask("chain A and resid 3 and not name N+CA+C+O")
    assert sel.matches(_atom("A", 3, "CB"))
    assert not sel.matches(_atom("A", 3, "CA"))
    assert not sel.matches(_atom("B", 3, "CB"))
    assert not sel.matches(_atom("A", 4, "CB"))

    # An unstated chain matches any, and names are matched case-insensitively.
    sel = strip.parse_strip_mask("resi 3 and name CB")
    assert sel.matches(_atom("B", 3, "cb"))

    # What the old four-term dialect refused, and this reads:
    for good in ("chain A or resid 3", "(resid 3)", "resid 3 and name CA CB",
                 "resname HOH", "resid 3 to 8", "same residue as name CB"):
        strip.parse_strip_mask(good)

    # What is refused is what cannot be read or cannot be answered.
    for bad in ("chain A and", "(resid 3", "wobble 7", "color red"):
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
    # `strip_keep_mask` says which rows survive; numpy does the slicing.
    keep = np.asarray(strip.strip_keep_mask(
        [r[0] for r in records], [r[1] for r in records], [r[2] for r in records],
        "chain A and resid 3 and not name N+CA+C+O+CB"), dtype=bool)
    out = atoms[keep]
    assert out.shape == (atoms.shape[0] - 3, 4)
    with pytest.raises(ValueError):
        strip.strip_keep_mask([r[0] for r in records], [r[1] for r in records],
                              [r[2] for r in records][:-1], "resid 3")


def test_hierarchy_non_destructive_and_in_place(tmp_path):
    pdb = tmp_path / "mini.pdb"
    pdb.write_text(_PDB)
    m = IMP.Model()
    hier = IMP.atom.read_pdb(str(pdb), m, IMP.atom.AllPDBSelector())
    n_before = len(IMP.atom.get_by_type(hier, IMP.atom.ATOM_TYPE))

    sel = strip.select_atoms(hier, "chain A and resid 3 and not name N+CA+C+O")
    assert sorted(IMP.atom.Atom(a).get_atom_type().get_string() for a in sel) == ["CB", "CD1", "CG1", "CG2"]

    # Stripping is in place and says how many atoms went. The old default
    # cloned the whole structure to leave the input untouched, which is a cost
    # a caller should choose: clone first if that is what you want.
    clone = IMP.atom.create_clone(hier)
    n = strip.strip_hierarchy(clone, "chain A and resid 3 and not name N+CA+C+O")
    assert n == 4
    assert len(IMP.atom.get_by_type(hier, IMP.atom.ATOM_TYPE)) == n_before  # untouched
    assert len(IMP.atom.get_by_type(clone, IMP.atom.ATOM_TYPE)) == n_before - 4

    n = strip.strip_hierarchy(hier, "resid 3 and name CB")
    assert n == 2
    assert len(IMP.atom.get_by_type(hier, IMP.atom.ATOM_TYPE)) == n_before - 2


def test_cgprobe_and_av_defaults_differ_only_in_cb(tmp_path):
    """The two consumers own their keep-sets; the engine is shared (PRD-106)."""
    from IMP.bff import site_keep_atom_names, strip_sidechain_at_site
    keep_names = list(site_keep_atom_names())
    pdb = tmp_path / "mini.pdb"
    pdb.write_text(_PDB)
    m = IMP.Model()
    hier = IMP.atom.read_pdb(str(pdb), m, IMP.atom.AllPDBSelector())
    assert "CB" not in keep_names
    removed = strip_sidechain_at_site(hier, "A", 3)
    assert removed == 4  # CB CG1 CG2 CD1
    keep_cb = strip_sidechain_at_site(
        IMP.atom.read_pdb(str(pdb), m, IMP.atom.AllPDBSelector()), "A", 3,
        keep_atom_names=keep_names + ["CB"])
    assert keep_cb == 3
    av_kept = strip.strip_pdb_lines(_LINES, strip.default_strip_mask("A", 3, "CB"))
    assert _names(av_kept, "A", 3) == ["N", "CA", "C", "O", "CB"]


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))


def test_one_backbone_list_serves_both_the_language_and_the_site():
    """Requirement 3: the keep-set exists in one place.

    `backbone` in a selection and the atoms a labelling site keeps are the
    same five names, so they are one list -- and its *order* is part of it,
    because a mask written from it should read the way a person writes one.
    """
    import IMP.bff

    assert list(IMP.bff.protein_backbone_atom_names()) == [
        "N", "CA", "C", "O", "OXT"]
    assert (list(IMP.bff.site_keep_atom_names()) ==
            list(IMP.bff.protein_backbone_atom_names()))
    assert (strip.site_strip_mask("A", 481, IMP.bff.site_keep_atom_names()) ==
            "chain A and resid 481 and not name N+CA+C+O+OXT")

    # The AV default is deliberately a different set: it keeps CB, because the
    # attachment atom stays, and drops OXT. That difference is data.
    assert list(strip.backbone_atom_names()) == ["N", "CA", "C", "O"]
    assert (strip.default_strip_mask("A", 36, "CB") ==
            "chain A and resid 36 and not name N+CA+C+O+CB")
