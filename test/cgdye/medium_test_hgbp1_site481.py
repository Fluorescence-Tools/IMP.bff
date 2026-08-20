#!/usr/bin/env python3

import sys
from pathlib import Path
from IMP.bff import get_template_dir, get_structure_dir

import IMP
import IMP.atom
import IMP.core


from IMP.bff.label import attach_dyes, resolve_dye_site


def _site_atoms(hier, chain_id, resnum):
    out = []
    for a in IMP.atom.get_by_type(hier, IMP.atom.ATOM_TYPE):
        res_p = a.get_parent()
        if not IMP.atom.Residue.get_is_setup(res_p):
            continue
        res = IMP.atom.Residue(res_p)
        if res.get_index() != resnum:
            continue
        chain_p = res_p.get_parent()
        if not IMP.atom.Chain.get_is_setup(chain_p):
            continue
        chain = IMP.atom.Chain(chain_p)
        if chain.get_id() != chain_id:
            continue
        out.append(a)
    return out


def test_hgbp1_site_481_exists_and_is_attachable():
    model = IMP.Model()
    protein = IMP.atom.read_pdb(
        str(get_structure_dir("1DG3.pdb")), model, IMP.atom.NonWaterPDBSelector()
    )
    site = resolve_dye_site(protein, "A", 481)
    assert "CA" in site and "N" in site and "C" in site


def test_attach_alexa488_r48_to_hgbp1_site_481():
    model = IMP.Model()
    protein = IMP.atom.read_pdb(
        str(get_structure_dir("1DG3.pdb")), model, IMP.atom.NonWaterPDBSelector()
    )
    dye = IMP.atom.read_mol2(str(get_structure_dir("alexa488_r48.mol2")), model)

    before = []
    for a in IMP.atom.get_by_type(dye, IMP.atom.ATOM_TYPE):
        c = IMP.core.XYZ(a).get_coordinates()
        before.append((c[0], c[1], c[2]))

    attached = attach_dyes(protein, [(dye, "A", 481)], strip_site_sidechain=True)
    assert len(attached) == 1
    assert attached[0]["resnum"] == 481

    after = []
    for a in IMP.atom.get_by_type(dye, IMP.atom.ATOM_TYPE):
        c = IMP.core.XYZ(a).get_coordinates()
        after.append((c[0], c[1], c[2]))

    changed = any(
        abs(b[0] - a[0]) > 1e-8 or abs(b[1] - a[1]) > 1e-8 or abs(b[2] - a[2]) > 1e-8
        for b, a in zip(before, after)
    )
    assert changed


def test_attach_strips_sidechain_atoms_at_site_481():
    model = IMP.Model()
    protein = IMP.atom.read_pdb(
        str(get_structure_dir("1DG3.pdb")), model, IMP.atom.NonWaterPDBSelector()
    )
    dye = IMP.atom.read_mol2(str(get_structure_dir("alexa488_r48.mol2")), model)

    backbone = {"N", "CA", "C", "O", "OXT"}
    before_names = {
        IMP.atom.Atom(a).get_atom_type().get_string().upper()
        for a in _site_atoms(protein, "A", 481)
    }
    before_sidechain = {n for n in before_names if n not in backbone}
    assert before_sidechain

    attach_dyes(protein, [(dye, "A", 481)], strip_site_sidechain=True)

    after_names = {
        IMP.atom.Atom(a).get_atom_type().get_string().upper()
        for a in _site_atoms(protein, "A", 481)
    }
    after_sidechain = {n for n in after_names if n not in backbone}
    assert not after_sidechain


# IMP runs every .py under test/ as a standalone script, and a file of bare
# pytest functions would import cleanly and exit 0 -- reporting success without
# running a single assertion. Hand the file to pytest explicitly so a failure
# here is a failure in ctest.
if __name__ == "__main__":
    import sys
    try:
        import pytest
    except ImportError:
        print("pytest not installed; skipping", __file__)
        sys.exit(0)
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
