import numpy as np
from IMP.bff.scoring.mean_field import rotamer_mean_field_weights, rotamer_mean_field_weights_multi_dye

def test_mean_field_weights_single():
    n_clusters = 5
    n_dye_atoms = 8
    n_prot = 20
    
    np.random.seed(42)
    rotamer_coords = np.random.randn(n_clusters, n_dye_atoms, 3) * 5
    initial_weights = np.ones(n_clusters) / n_clusters
    protein_coords = np.random.randn(n_prot, 3) * 3
    
    dye_elems = ['C'] * n_dye_atoms
    prot_elems = ['C'] * n_prot
    
    w = rotamer_mean_field_weights(
        rotamer_coords, initial_weights, protein_coords,
        dye_elems, prot_elems, K=1.0, n_iter=5
    )
    
    assert w.shape == (n_clusters,)
    assert abs(w.sum() - 1.0) < 1e-9

def test_mean_field_weights_multi():
    n_clusters_1 = 5
    n_clusters_2 = 3
    n_dye_atoms_1 = 8
    n_dye_atoms_2 = 6
    n_prot = 20
    
    np.random.seed(42)
    rot_1 = np.random.randn(n_clusters_1, n_dye_atoms_1, 3) * 5
    rot_2 = np.random.randn(n_clusters_2, n_dye_atoms_2, 3) * 5 + 10 # Far away
    
    w1_init = np.ones(n_clusters_1) / n_clusters_1
    w2_init = np.ones(n_clusters_2) / n_clusters_2
    
    protein_coords = np.random.randn(n_prot, 3) * 3
    
    dye1_elems = ['C'] * n_dye_atoms_1
    dye2_elems = ['C'] * n_dye_atoms_2
    prot_elems = ['C'] * n_prot
    
    w_list = rotamer_mean_field_weights_multi_dye(
        [rot_1, rot_2],
        [w1_init, w2_init],
        protein_coords,
        [dye1_elems, dye2_elems],
        prot_elems, K=1.0, n_iter=3
    )
    
    assert len(w_list) == 2
    assert w_list[0].shape == (n_clusters_1,)
    assert w_list[1].shape == (n_clusters_2,)
    
    assert abs(w_list[0].sum() - 1.0) < 1e-9
    assert abs(w_list[1].sum() - 1.0) < 1e-9


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
