import os
import pytest
import numpy as np
import random
import IMP
import IMP.atom
import IMP.core
from pathlib import Path

from IMP.bff import read_rotamer_library_rmf
from IMP.bff import apply_coordinates, sample_weighted_index
from IMP.bff import linker_geometry_from_mol2, sample_linker
from IMP.bff import get_template_dir, get_structure_dir

ROTAMER_LIB_DIR = Path(get_template_dir("rotamer"))
STRUCTURES_DIR = Path(get_structure_dir())

def get_all_dyes():
    dyes = []
    for f in os.listdir(ROTAMER_LIB_DIR):
        if f.endswith(".rmf3"):
            dyes.append(f.replace(".rmf3", ""))
    return sorted(dyes)

@pytest.mark.parametrize("dye_name", get_all_dyes())
def test_rotamer_sampling_all_dyes(probe_name):
    """Test that every rotamer library can be loaded and sampled."""
    lib_path = os.path.join(ROTAMER_LIB_DIR, f"{probe_name}.rmf3")
    # An `IMP.bff.ProbeRotamerLibrary` -- the same value every other reader
    # returns. This one used to hand back a dict, because it was Python.
    lib = read_rotamer_library_rmf(lib_path)

    weights = lib.weights
    assert len(weights) > 0
    assert lib.coords.shape == (lib.n_rotamers, lib.n_atoms, 3)
    assert len(weights) == lib.n_rotamers

    # Test minimal sampling
    rng = random.Random(42)
    idx = sample_weighted_index(weights, seed=rng.randint(0, 2**31 - 1))
    coords = lib.coords[idx]
    assert coords.shape[0] == len(lib.atom_names)

@pytest.mark.parametrize("dye_mol2", [
    "alexa488_r48.mol2",
    "atto655.mol2",
    "cx4.mol2"
])
def test_linker_sampling_available_mol2(dye_mol2):
    """Metropolis sampling over the linker runs for every bundled MOL2."""
    mol2_path = os.path.join(STRUCTURES_DIR, dye_mol2)
    geometry = linker_geometry_from_mol2(mol2_path)

    # the molecule has something that turns
    assert (geometry.get_number_of_torsions() > 0
            or geometry.get_number_of_angles() > 0)

    result = sample_linker(mol2_path, n_steps=10, write_every=2, seed=42)
    assert result.n_frames == 5                       # 10 / 2
    assert result.coordinates.shape == (5, result.n_atoms, 3)
    assert 0.0 <= result.acceptance <= 1.0


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
