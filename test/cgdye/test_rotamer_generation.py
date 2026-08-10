"""Tests for rotamer library generation (Phase 5)."""

import unittest
import numpy as np
import os
import tempfile

from IMP.bff.cgdye.sampling.clustering import cluster_leader, cluster_assignment
from IMP.bff.cgdye.sampling.boltzmann import compute_boltzmann_weights, cluster_weights
from IMP.bff.cgdye.io.rotamer_cif import write_rotamer_library, read_rotamer_library


class TestRotamerGeneration(unittest.TestCase):
    def test_clustering_and_boltzmann(self):
        # Create 100 frames of 2 atoms
        n_frames = 100
        n_atoms = 2
        coords = np.zeros((n_frames, n_atoms, 3))
        
        # Two distinct groups of conformations
        coords[:50, 1, 0] = 1.0  # Group 1
        coords[50:, 1, 0] = 5.0  # Group 2
        
        # Add some noise
        coords += np.random.normal(0, 0.1, coords.shape)
        
        # Clustering
        centers = cluster_leader(coords, threshold=1.0)
        self.assertEqual(len(centers), 2)
        
        assignments = cluster_assignment(coords, centers)
        self.assertEqual(len(np.unique(assignments)), 2)
        
        # Boltzmann scoring
        # Group 1 has lower energy
        energies = np.zeros(n_frames)
        energies[:50] = 0.0
        energies[50:] = 10.0 # High energy group
        
        weights = compute_boltzmann_weights(energies, temperature=298.15)
        self.assertAlmostEqual(np.sum(weights), 1.0)
        self.assertGreater(weights[0], weights[50])
        
        c_weights = cluster_weights(assignments, weights, len(centers))
        self.assertAlmostEqual(np.sum(c_weights), 1.0)
        self.assertGreater(c_weights[0], c_weights[1])

    def test_io_cycle(self):
        library = {
            "weight": [0.7, 0.3],
            "atom_names": ["C1", "C2"],
            "coords": {
                1: np.array([[0, 0, 0], [1, 0, 0]]),
                2: np.array([[0, 0, 0], [0, 1, 0]])
            }
        }
        
        with tempfile.TemporaryDirectory() as tmpdir:
            base = os.path.join(tmpdir, "test_lib")
            write_rotamer_library(base, library)
            
            # Verify files exist
            self.assertTrue(os.path.exists(base + "_coords.npy"))
            self.assertTrue(os.path.exists(base + "_weights.txt"))
            self.assertTrue(os.path.exists(base + "_atoms.txt"))
            
            read_lib = read_rotamer_library(base)
            self.assertEqual(len(read_lib["weight"]), 2)
            self.assertEqual(read_lib["atom_names"], ["C1", "C2"])
            np.testing.assert_array_almost_equal(read_lib["coords"][1], library["coords"][1])


if __name__ == "__main__":
    unittest.main()
