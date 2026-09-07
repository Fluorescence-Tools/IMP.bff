"""The selection language: both spellings, and the masks that ship with fps.json.

The grammar is one parser with two vocabularies. One is what a molecular
viewer takes (`resi`, `resn`, `name CA+CB`, `50-60`, `byres`); the other is
what the pteros/VMD family takes (`resid`, `resname`, `name CA CB`,
`50 to 60`, `same residue as`). Files in this repository carry both, sometimes
in one expression, so both are accepted.
"""

import pytest

import IMP
import IMP.atom
import IMP.bff


def _t4l():
    m = IMP.Model()
    h = IMP.atom.read_pdb(IMP.bff.get_example_path("structure/T4L/3GUN.pdb"), m,
                          IMP.atom.NonWaterNonHydrogenPDBSelector())
    return m, h


def _atoms(chain="A", resi=(1, 2), names=("N", "CA", "C", "O", "CB")):
    out = []
    i = 0
    for r in resi:
        for n in names:
            a = IMP.bff.SelectionAtom()
            a.chain, a.resi, a.name = chain, r, n
            a.resn = "ALA"
            a.elem = n[0]
            i += 1
            a.index, a.id = i, i
            out.append(a)
    return out


# --------------------------------------------------------------------------
# the two spellings
# --------------------------------------------------------------------------

@pytest.mark.parametrize("expression", [
    "chain A and resi 1 and not name CA+CB",       # one vocabulary
    "chain A and resid 1 and not name CA CB",      # the other
    "chain A and resi 1 and not name CA CB",       # mixed
    "chain A and resid 1 and not name CA+CB",
])
def test_the_two_spellings_select_the_same_atoms(expression):
    atoms = _atoms()
    selected = IMP.bff.SelectionExpression(expression).evaluate(atoms)
    got = {(atoms[i].resi, atoms[i].name)
           for i, s in enumerate(selected) if s}
    assert got == {(1, "N"), (1, "C"), (1, "O")}


@pytest.mark.parametrize("expression", [
    "resi 2-3", "resid 2 to 3", "resi 2+3", "resid 2,3",
])
def test_ranges_and_lists_agree(expression):
    atoms = _atoms(resi=(1, 2, 3, 4))
    selected = IMP.bff.SelectionExpression(expression).evaluate(atoms)
    assert {atoms[i].resi for i, s in enumerate(selected) if s} == {2, 3}


def test_wildcards():
    atoms = _atoms(names=("N", "CA", "CB", "CG", "OD1"))
    sel = IMP.bff.SelectionExpression("name C*")
    selected = sel.evaluate(atoms)
    assert {atoms[i].name for i, s in enumerate(selected) if s} == {"CA", "CB", "CG"}


def test_precedence_is_not_then_and_then_or():
    atoms = _atoms(resi=(1, 2))
    # `not` binds tighter than `and`: not(name CA) and resi 1
    sel = IMP.bff.SelectionExpression("not name CA and resi 1")
    got = {(atoms[i].resi, atoms[i].name)
           for i, s in enumerate(sel.evaluate(atoms)) if s}
    assert got == {(1, "N"), (1, "C"), (1, "O"), (1, "CB")}
    # `and` binds tighter than `or`
    sel = IMP.bff.SelectionExpression("resi 1 and name CA or name CB")
    got = {(atoms[i].resi, atoms[i].name)
           for i, s in enumerate(sel.evaluate(atoms)) if s}
    assert got == {(1, "CA"), (1, "CB"), (2, "CB")}


def test_byres_and_same_residue_as_agree():
    atoms = _atoms(resi=(1, 2))
    a = IMP.bff.SelectionExpression("byres name CB").evaluate(atoms)
    b = IMP.bff.SelectionExpression("same residue as name CB").evaluate(atoms)
    assert list(a) == list(b)
    assert sum(a) == len(atoms)  # every atom of both residues


def test_byres_binds_loosest():
    """`byres chain A and name CA` extends the whole conjunction."""
    atoms = _atoms(resi=(1, 2))
    for a in atoms[5:]:
        a.chain = "B"
    sel = IMP.bff.SelectionExpression("byres chain A and name CA")
    got = {(atoms[i].chain, atoms[i].resi)
           for i, s in enumerate(sel.evaluate(atoms)) if s}
    assert got == {("A", 1)}


# --------------------------------------------------------------------------
# the masks that ship
# --------------------------------------------------------------------------

def test_the_shipped_t4l_masks_parse_and_apply():
    """Every `strip_mask` in `examples/structure/T4L/fret.fps.json`.

    They are the reason this parser exists: each is a parenthesised group,
    an `or`, and a space-separated `resname` list.
    """
    import json
    from pathlib import Path

    repo = Path(__file__).resolve().parents[2]
    fps = json.loads((repo / "examples/structure/T4L/fret.fps.json").read_text())
    masks = [p["strip_mask"] for p in fps["Positions"].values()
             if p.get("strip_mask")]
    assert masks, "the shipped file carries no strip_mask"

    m, h = _t4l()
    for mask in masks:
        selection = IMP.bff.selection_from_expression(h, mask)
        # the residue named in the mask loses its side chain, and nothing else
        picked = selection.get_selected_particle_indexes()
        assert len(picked) > 0, mask
        names = {IMP.atom.Atom(m.get_particle(p)).get_atom_type().get_string()
                 for p in picked}
        assert not (names & {"CA", "CB", "C", "N", "O"}), mask


def test_an_olga_mask_parses_and_applies():
    """Olga writes the same intent in the other vocabulary."""
    m, h = _t4l()
    selection = IMP.bff.selection_from_expression(
        h, "chain A and resid 115 and not name CA CB C N O")
    picked = selection.get_selected_particle_indexes()
    assert len(picked) > 0
    for p in picked:
        atom = IMP.atom.Atom(m.get_particle(p))
        residue = IMP.atom.Residue(atom.get_parent())
        assert residue.get_index() == 115
        assert atom.get_atom_type().get_string() not in ("CA", "CB", "C", "N", "O")


# --------------------------------------------------------------------------
# it is an IMP selection
# --------------------------------------------------------------------------

def test_an_expression_is_an_imp_selection():
    """The compile target is `IMP.atom.Selection`, so it composes with IMP."""
    m, h = _t4l()
    expression = IMP.bff.selection_from_expression(h, "resi 10-12 and name CA")
    native = IMP.atom.Selection(h, residue_indexes=[10, 11, 12],
                                atom_type=IMP.atom.AT_CA)
    assert (sorted(expression.get_selected_particle_indexes()) ==
            sorted(native.get_selected_particle_indexes()))


def test_within_measures_a_distance():
    m, h = _t4l()
    near = IMP.bff.selection_from_expression(h, "within 5 of resi 10")
    far = IMP.bff.selection_from_expression(h, "beyond 5 of resi 10")
    n_near = len(near.get_selected_particle_indexes())
    n_far = len(far.get_selected_particle_indexes())
    assert n_near > 0 and n_far > 0
    assert n_near + n_far == len(IMP.atom.get_leaves(h))


# --------------------------------------------------------------------------
# what it refuses
# --------------------------------------------------------------------------

@pytest.mark.parametrize("expression,fragment", [
    ("chain A and", "needs a right-hand side"),
    ("(chain A", "unbalanced"),
    ("chain A) and resi 5", "unexpected"),
    ("wobble 3", "not a selection keyword"),
    ("color red", "cannot answer"),
    ("within 5 pbc of resi 3", "no cell"),
    ("name", "needs a value"),
])
def test_a_bad_expression_says_what_is_wrong(expression, fragment):
    with pytest.raises((ValueError, IMP.ValueException)) as caught:
        IMP.bff.SelectionExpression(expression)
    assert fragment in str(caught.value)


# --------------------------------------------------------------------------
# which spelling an expression is in
# --------------------------------------------------------------------------

@pytest.mark.parametrize("expression,dialect", [
    ("chain A and resi 3 and not name N+CA+C+O", IMP.bff.SELECTION_DIALECT_PLUS_LIST),
    ("byres name CB", IMP.bff.SELECTION_DIALECT_PLUS_LIST),
    ("resn ALA", IMP.bff.SELECTION_DIALECT_PLUS_LIST),
    ("chain A and resid 115 and not name CA CB C N O",
     IMP.bff.SELECTION_DIALECT_SPACE_LIST),
    ("resname HOH SOL WAT", IMP.bff.SELECTION_DIALECT_SPACE_LIST),
    ("same residue as name CB", IMP.bff.SELECTION_DIALECT_SPACE_LIST),
    ("resid 50 to 60", IMP.bff.SELECTION_DIALECT_SPACE_LIST),
    ("chain A and name CA", IMP.bff.SELECTION_DIALECT_UNMARKED),
    ("resid 3", IMP.bff.SELECTION_DIALECT_UNMARKED),
    ("(resid 132 and not name CA+CB+C+N+O) or resname HOH SOL WAT",
     IMP.bff.SELECTION_DIALECT_MIXED),
])
def test_the_spelling_is_detected(expression, dialect):
    assert IMP.bff.selection_dialect(expression) == dialect
    assert IMP.bff.SelectionExpression(expression).get_dialect() == dialect


def test_index_counts_from_one_or_zero_by_spelling():
    """The one place the two spellings disagree about meaning.

    `index 1` is the first atom in one and the second in the other, so the
    expression alone does not say which atom is meant -- the spelling does.
    """
    atoms = _atoms(resi=(1,), names=("N", "CA", "C"))
    # marked as the space-list spelling by `resname`, so index counts from 0
    space = IMP.bff.SelectionExpression("resname ALA XXX and index 0")
    assert space.get_dialect() == IMP.bff.SELECTION_DIALECT_SPACE_LIST
    got = [atoms[i].name for i, s in enumerate(space.evaluate(atoms)) if s]
    assert got == ["N"]

    # marked as the other spelling by `resn`, so index counts from 1
    plus = IMP.bff.SelectionExpression("resn ALA and index 1")
    assert plus.get_dialect() == IMP.bff.SELECTION_DIALECT_PLUS_LIST
    got = [atoms[i].name for i, s in enumerate(plus.evaluate(atoms)) if s]
    assert got == ["N"]


def test_a_caller_can_declare_the_spelling():
    """Detection is a default, not a decision taken away from the caller."""
    atoms = _atoms(resi=(1,), names=("N", "CA", "C"))
    expression = "index 1"
    assert (IMP.bff.selection_dialect(expression) ==
            IMP.bff.SELECTION_DIALECT_UNMARKED)
    declared = IMP.bff.SelectionExpression(
        expression, IMP.bff.SELECTION_DIALECT_SPACE_LIST)
    got = [atoms[i].name for i, s in enumerate(declared.evaluate(atoms)) if s]
    assert got == ["CA"]  # counted from 0, so 1 is the second atom
