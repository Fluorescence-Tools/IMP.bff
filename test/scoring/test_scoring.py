import numpy as np
from IMP.bff.scoring import BoundingBoxFilter, DyeInternalEnergyEvaluator, lj_score

def test_bounding_box_filter():
    bbf = BoundingBoxFilter(pad=3.5)
    
    # Create rotamer coords: 2 frames, 1 atom each
    coords = np.array([
        [[0.0, 0.0, 0.0]], # Frame 0 at origin
        [[100.0, 100.0, 100.0]] # Frame 1 far away
    ])
    
    # Create reference coords: 1 atom at origin
    ref_coords = np.array([[0.0, 0.0, 0.0]])
    
    surviving, mask = bbf.filter_frames(coords, ref_coords)
    
    # Frame 0 should survive, Frame 1 should be filtered out
    assert mask.shape == (2,)
    assert mask[0] == True
    assert mask[1] == False
    assert surviving.shape == (1, 1, 3)

def test_internal_energy_evaluator_vectorized():
    system = {
        'sites': [
            {'id': 'a', 'atom_name': 'C1'},
            {'id': 'b', 'atom_name': 'C2'}
        ],
        'bonds': []
    }
    evaluator = DyeInternalEnergyEvaluator(system)
    
    # 2 frames, 2 atoms
    coords = np.array([
        [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]], # Frame 0: distance 3.0
        [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]]  # Frame 1: distance 5.0 (should be 0 energy if > rmin)
    ])
    
    energies = evaluator.evaluate_batch(coords)
    assert energies.shape == (2,)
    
    # Manually calculate expected energy for Frame 0
    # Assuming C-C interaction uses default params from CHARMM36_LJ
    expected_energy_0 = evaluator.evaluate(coords[0])
    expected_energy_1 = evaluator.evaluate(coords[1])
    
    np.testing.assert_allclose(energies[0], expected_energy_0)
    np.testing.assert_allclose(energies[1], expected_energy_1)

def test_evaluate_batch_filtered():
    system = {
        'sites': [
            {'id': 'a', 'atom_name': 'C1'},
            {'id': 'b', 'atom_name': 'C2'}
        ],
        'bonds': []
    }
    evaluator = DyeInternalEnergyEvaluator(system)
    
    # Frame 0 is close, Frame 1 is far
    coords = np.array([
        [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
        [[100.0, 0.0, 0.0], [103.0, 0.0, 0.0]]
    ])
    
    ref_coords = np.array([[0.0, 0.0, 0.0]])
    
    energies, mask = evaluator.evaluate_batch_filtered(coords, ref_coords)
    
    assert mask.shape == (2,)
    assert mask[0] == True
    assert mask[1] == False
    assert energies[1] == 0.0 # Energy should be 0 since it didn't overlap
    assert energies[0] != 0.0 # Energy should be evaluated


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
