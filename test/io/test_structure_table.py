"""``read_structure_table``: IMP's readers behind flat columns.

The point of the door is that an application gets IMP's answer -- its radii,
its masses, its element table, its mmCIF reader -- without importing IMP and
without a SWIG call per atom. So the tests here check three things: that the
columns agree with IMP's own hierarchy walk atom for atom, that mmCIF is read
at all, and that the refusals are refusals rather than empty structures.

The whole file names ``IMP.atom``, so ``conftest.py`` deselects it in the
IMP-free lane, which is correct: the door does not exist there.
"""

import math
from pathlib import Path

import numpy as np
import pytest

import IMP
import IMP.atom
import IMP.core
import IMP.bff

REPO = Path(__file__).resolve().parents[2]
PDB = REPO / "examples" / "structure" / "T4L" / "3GUN.pdb"
CIF = REPO / "examples" / "structure" / "GBP" / "1F5N.cif"


def _imp_reference(path, keep_water=False):
    """The same read, done the long way through IMP's Python."""
    model = IMP.Model()
    selector = (IMP.atom.NonAlternativePDBSelector() if keep_water
                else IMP.atom.NonWaterPDBSelector())
    if str(path).lower().endswith((".pdb", ".ent")):
        h = IMP.atom.read_pdb(str(path), model, selector)
    else:
        h = IMP.atom.read_mmcif(str(path), model, selector)
    rows = []
    scale = 2.0 ** (-1.0 / 6.0)
    for a in IMP.atom.get_by_type(h, IMP.atom.ATOM_TYPE):
        atom = IMP.atom.Atom(a)
        residue = IMP.atom.Residue(atom.get_parent())
        res_name = residue.get_residue_type().get_string()
        if not IMP.bff.is_standard_residue(res_name):
            continue
        xyzr = IMP.core.XYZR(a)
        v = xyzr.get_coordinates()
        rows.append((res_name, residue.get_index(),
                     (v[0], v[1], v[2]),
                     xyzr.get_radius() * scale,
                     IMP.atom.Mass(a).get_mass()))
    return rows


@pytest.mark.skipif(not PDB.exists(), reason="example structure not present")
def test_columns_match_imps_own_walk():
    table = IMP.bff.read_structure_table(str(PDB))
    reference = _imp_reference(PDB)

    assert table.n_atoms == len(reference)
    assert table.xyz.shape == (len(reference), 3)

    for i, (res_name, res_id, xyz, radius, mass) in enumerate(reference):
        assert table.res_name[i] == res_name
        assert table.res_id[i] == res_id
        assert table.xyz[i] == pytest.approx(xyz)
        assert table.radius[i] == pytest.approx(radius)
        assert table.mass[i] == pytest.approx(mass)


@pytest.mark.skipif(not PDB.exists(), reason="example structure not present")
def test_every_column_is_the_same_length():
    table = IMP.bff.read_structure_table(str(PDB))
    n = table.n_atoms
    assert n > 0
    assert (len(table.radius) == len(table.mass) == len(table.bfactor)
            == len(table.atom_id) == len(table.res_id) == n)
    assert (len(table.chain) == len(table.res_name) == len(table.atom_name)
            == len(table.element) == n)
    assert table.xyz.shape == (n, 3)


@pytest.mark.skipif(not PDB.exists(), reason="example structure not present")
def test_rmin_and_the_zero_of_the_potential_differ_by_the_lj_factor():
    """The radius the AV has always been given is 2^(-1/6) Rmin, not Rmin."""
    at_zero = IMP.bff.read_structure_table(str(PDB), radius_no_interaction=True)
    at_min = IMP.bff.read_structure_table(str(PDB), radius_no_interaction=False)
    assert at_zero.n_atoms == at_min.n_atoms
    ratio = np.asarray(at_zero.radius) / np.asarray(at_min.radius)
    assert ratio == pytest.approx(2.0 ** (-1.0 / 6.0))


@pytest.mark.skipif(not PDB.exists(), reason="example structure not present")
def test_water_is_dropped_by_default_and_kept_on_request():
    without = IMP.bff.read_structure_table(str(PDB), keep_water=False)
    with_water = IMP.bff.read_structure_table(str(PDB), keep_water=True)
    # Whatever the file holds, keeping water can only add atoms.
    assert with_water.n_atoms >= without.n_atoms


@pytest.mark.skipif(not PDB.exists(), reason="example structure not present")
def test_non_standard_residues_are_dropped_only_when_asked():
    filtered = IMP.bff.read_structure_table(str(PDB),
                                            only_standard_residues=True)
    unfiltered = IMP.bff.read_structure_table(str(PDB),
                                              only_standard_residues=False)
    assert unfiltered.n_atoms >= filtered.n_atoms
    assert all(IMP.bff.is_standard_residue(r) for r in filtered.res_name)


@pytest.mark.skipif(not CIF.exists(), reason="example mmCIF not present")
def test_mmcif_is_read():
    """The format the core cannot parse at all -- this door is why it can."""
    table = IMP.bff.read_structure_table(str(CIF))
    assert table.n_atoms > 0
    reference = _imp_reference(CIF)
    assert table.n_atoms == len(reference)
    assert table.xyz[0] == pytest.approx(reference[0][2])


def test_an_unreadable_extension_raises_rather_than_returning_nothing():
    """An empty structure is the failure mode this refuses to have.

    A trajectory handed to a coordinate reader used to come back with zero
    atoms, and whatever was built from it was simply blank.
    """
    with pytest.raises(Exception):
        IMP.bff.read_structure_table("/tmp/does_not_matter.dcd")


def test_a_compressed_file_is_refused_by_name():
    """IMP's TextInput does not decompress, and a gzip stream read as text
    yields an empty structure rather than an error -- so it is refused."""
    with pytest.raises(Exception):
        IMP.bff.read_structure_table("/tmp/whatever.pdb.gz")


@pytest.mark.skipif(not PDB.exists(), reason="example structure not present")
def test_atom_names_carry_no_HET_prefix():
    """IMP prints an unrecognised type as ``HET: C3``; a five-character field
    truncates that to ``HET: `` and every ligand atom shares one name."""
    table = IMP.bff.read_structure_table(str(PDB),
                                         only_standard_residues=False)
    assert not any(n.upper().startswith("HET:") for n in table.atom_name)
    assert not any(n.startswith('"') for n in table.atom_name)


def test_is_standard_residue_knows_both_nucleotide_spellings():
    for name in ("ALA", "GLY", "HIS", "DA", "DT", "A", "U"):
        assert IMP.bff.is_standard_residue(name)
    for name in ("HOH", "ATP", "MG", "NAG"):
        assert not IMP.bff.is_standard_residue(name)
    # whitespace and case are the file's business, not the caller's
    assert IMP.bff.is_standard_residue(" ala ")
