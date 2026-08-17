"""Test dye attachment logic against reference implementation."""

import os
from pathlib import Path

import IMP
import IMP.atom
import IMP.core
import numpy as np
import pytest

from IMP.bff.cgdye.labeling.attachment import attach_dyes, place_dye_from_coords
from IMP.bff.cgdye.utils import get_structure_dir, get_output_dir

# The FRETpredict reference libraries are IMP.bff module data
# (data/rotamer_library): one <dye>_<linker>.pdb per dye/linker combination,
# authored in the residue's backbone frame (CA at the origin, x along CA->N,
# y in the N-CA-C plane).
def _reference_lib_dir():
    import IMP.bff
    return Path(IMP.bff.get_data_path("rotamer_library"))


def _supported_reference_combinations():
    """Return list of (dye_name, pdb_path) for available reference dyes."""
    ref_dir = _reference_lib_dir()
    if not ref_dir.exists():
        return []
    return sorted((pdb.stem, pdb) for pdb in ref_dir.glob("*.pdb"))


def _reference_transform(ca, n, c, points):
    """Reference-compatible backbone frame transformation (numpy)."""
    origin = ca
    v1 = n - origin
    v1 /= np.linalg.norm(v1)
    v2 = c - origin
    v2 -= np.dot(v2, v1) * v1
    v2 /= np.linalg.norm(v2)
    v3 = np.cross(v1, v2)
    R = np.vstack([v1, v2, v3]).T
    return np.dot(points, R.T) + origin


def _coords(hier):
    return np.array([IMP.core.XYZ(a).get_coordinates()
                     for a in IMP.atom.get_by_type(hier, IMP.atom.ATOM_TYPE)])


@pytest.fixture(scope="module")
def hgbp1_site481():
    model = IMP.Model()
    protein = IMP.atom.read_pdb(str(get_structure_dir("1DG3.pdb")), model, IMP.atom.NonWaterPDBSelector())
    xyz = {}
    for name in ("CA", "N", "C"):
        p = IMP.atom.Selection(protein, chain_id="A", residue_index=481,
                               atom_type=IMP.atom.AtomType(name)).get_selected_particles()[0]
        xyz[name] = np.array(IMP.core.XYZ(p).get_coordinates())
    return model, protein, xyz


@pytest.mark.parametrize("combo_name,pdb_path", _supported_reference_combinations())
def test_attach_all_reference_dyes_to_hgbp1_site481(combo_name, pdb_path, hgbp1_site481):
    """cgdye placement equals the reference backbone-frame transform for every library."""
    model, protein, xyz = hgbp1_site481
    dye = IMP.atom.read_pdb(str(pdb_path), model, IMP.atom.AllPDBSelector())
    base = _coords(dye)
    assert base.shape[0] > 0

    # attach_dyes moves the dye in place; the protein is not modified here
    # (strip_site_sidechain=False) so the module-scoped protein stays reusable.
    attach_dyes(protein, [(dye, "A", 481)], strip_site_sidechain=False)
    placed = _coords(dye)
    expected = _reference_transform(xyz["CA"], xyz["N"], xyz["C"], base)
    np.testing.assert_allclose(placed, expected, atol=1e-6)
    IMP.atom.destroy(dye)


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
