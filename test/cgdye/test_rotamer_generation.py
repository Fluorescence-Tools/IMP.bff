"""Tests for rotamer library generation (Phase 5)."""

import unittest
import numpy as np
import os
import tempfile

from IMP.bff.cgdye.sampling.clustering import cluster_frames_leader, assign_frames_to_clusters
from IMP.bff.scoring.boltzmann import boltzmann_weights, rotamer_cluster_weights
from IMP.bff.io.rotamer_cif import write_rotamer_library, read_rotamer_library


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
        centers = cluster_frames_leader(coords, threshold=1.0)
        self.assertEqual(len(centers), 2)
        
        assignments = assign_frames_to_clusters(coords, centers)
        self.assertEqual(len(np.unique(assignments)), 2)
        
        # Boltzmann scoring
        # Group 1 has lower energy
        energies = np.zeros(n_frames)
        energies[:50] = 0.0
        energies[50:] = 10.0 # High energy group
        
        weights = boltzmann_weights(energies, temperature=298.15)
        self.assertAlmostEqual(np.sum(weights), 1.0)
        self.assertGreater(weights[0], weights[50])
        
        c_weights = rotamer_cluster_weights(assignments, weights, len(centers))
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


class TestLinkerSamplerPhysics(unittest.TestCase):
    """The linker LJ score excludes bonded neighbours; a short run is pinned."""

    def test_bonded_pairs_are_excluded_from_the_linker_score(self):
        from IMP.bff.cgdye.topology.builder import parse_dye_mol2
        from IMP.bff.scoring.dye_lj import (
            DyeInternalEnergyEvaluator, compute_exclusions, dye_internal_system)
        from IMP.bff.tools.paths import get_structure_dir
        atoms, bonds = parse_dye_mol2(str(get_structure_dir("alexa488_r48.mol2")), "dye")
        system = dye_internal_system(atoms, bonds)
        self.assertEqual(len(system["sites"]), len(atoms))
        self.assertEqual(len(system["bonds"]), len(bonds))
        excluded = compute_exclusions(system)
        # every bond, angle end pair and dihedral end pair is excluded
        for a, b, _l, _t in system["bonds"]:
            self.assertIn(frozenset({a, b}), excluded)
        for a, _b, c, _th, _t in system["angles"]:
            self.assertIn(frozenset({a, c}), excluded)
        for a, _b, _c, d, _t in system["dihedrals"]:
            self.assertIn(frozenset({a, d}), excluded)
        with_excl = DyeInternalEnergyEvaluator(system)
        without = DyeInternalEnergyEvaluator({"sites": system["sites"], "bonds": []})
        self.assertLess(len(with_excl.pairs), len(without.pairs))
        # site ids are unique although MOL2 atom names are not (83 atoms, 68 names)
        self.assertEqual(len({s["id"] for s in system["sites"]}), len(atoms))
        # the MOL2 geometry scores a few tens of kcal/mol (repulsive-only 12-6
        # on non-bonded pairs) once bonded pairs are excluded; scored against
        # bonded neighbours (the old behaviour) it is ~1e6
        coords = np.array([[a["x"], a["y"], a["z"]] for a in sorted(atoms.values(), key=lambda x: x["serial"])])
        self.assertLess(with_excl.evaluate(coords), 50.0)
        self.assertGreater(without.evaluate(coords), 1e5)

    def test_generate_rotamers_pins(self):
        import hashlib
        import json
        from pathlib import Path
        from IMP.bff.cgdye.sampling.library_gen import generate_linker_rotamers
        from IMP.bff.tools.paths import get_structure_dir
        pins_path = Path(__file__).resolve().parents[1] / "references" / "cgdye_sampler_pins.json"
        with open(pins_path) as fh:
            pins = json.load(fh)
        lib = generate_linker_rotamers(str(get_structure_dir("alexa488_r48.mol2")),
                                n_steps=300, write_every=10, cluster_threshold=1.0, seed=42)
        weights = np.asarray(lib["weight"])
        self.assertEqual(len(weights), pins["n_clusters"])
        self.assertAlmostEqual(float(weights.sum()), 1.0, places=12)
        self.assertTrue((weights > 0).all())
        np.testing.assert_allclose(weights, pins["weights"], rtol=0, atol=1e-10)
        transitions = np.asarray(lib["transitions"])
        self.assertEqual(transitions.shape, (len(weights), len(weights)))
        self.assertEqual(int(transitions.sum()), int(np.asarray(pins["transitions"]).sum()))
        np.testing.assert_array_equal(transitions, np.asarray(pins["transitions"]))
        coords = np.stack([lib["coords"][i + 1] for i in range(len(weights))])
        self.assertEqual(coords.shape[1], pins["n_atoms"])
        self.assertEqual(hashlib.sha256(np.round(coords, 6).tobytes()).hexdigest(), pins["centre_coords_sha256"])
        self.assertEqual(lib["atom_names"][:5], pins["atom_names_first5"])


if __name__ == "__main__":
    unittest.main()
