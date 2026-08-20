import os
import pytest
import numpy as np
import random
import IMP
import IMP.atom
import IMP.core
from pathlib import Path

from IMP.bff.io.structure import read_rotamer_library_rmf
from IMP.bff.sampling import apply_rotamer_coordinates, sample_rotamer_index
from IMP.bff.cgdye.sampling import LinkerSampler
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
def test_rotamer_sampling_all_dyes(dye_name):
    """Test that every rotamer library can be loaded and sampled."""
    lib_path = os.path.join(ROTAMER_LIB_DIR, f"{dye_name}.rmf3")
    lib = read_rotamer_library_rmf(lib_path)
    
    assert "coords" in lib
    assert "weight" in lib or "weights" in lib # read_rotamer_library_rmf returns 'weight'
    weights = lib.get("weight") or lib.get("weights")
    assert len(weights) > 0
    assert len(lib["coords"]) == len(weights)
    
    # Test minimal sampling
    rng = random.Random(42)
    idx = sample_rotamer_index(weights, rng=rng)
    coords = lib["coords"][idx + 1]
    assert coords.shape[0] == len(lib["atom_names"])

@pytest.mark.parametrize("dye_mol2", [
    "alexa488_r48.mol2",
    "atto655.mol2",
    "cx4.mol2"
])
def test_linker_sampling_available_mol2(dye_mol2):
    """LinkerSampler (Metropolis over linker torsions) runs for available MOL2 files."""
    mol2_path = os.path.join(STRUCTURES_DIR, dye_mol2)
    sampler = LinkerSampler(mol2_path)
    
    # Verify internal DOFs are found
    assert len(sampler.rot_bonds) > 0 or len(sampler.rot_angles) > 0
    
    # Test minimal sampling (10 steps)
    n_steps = 10
    write_every = 2
    frames = sampler.sample(n_steps=n_steps, write_every=write_every, seed=42)
    
    # 10 / 2 = 5 frames
    assert frames.shape[0] == 5
    assert frames.shape[1] == len(sampler.atoms)


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
