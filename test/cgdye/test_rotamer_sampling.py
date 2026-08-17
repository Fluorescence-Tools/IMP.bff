"""Tests for rotamer-library based sampling."""

import os
import sys
import tempfile
import unittest
from pathlib import Path

import IMP
import IMP.atom
import IMP.core
import numpy as np

from IMP.bff.cgdye.sampling.rotamer import (
    apply_rotamer_coords,
    find_reference_rotamer_files,
    load_reference_rotamers,
    sample_rotamer_index,
)


class TestRotamerSampling(unittest.TestCase):
    def test_load_reference_rotamers(self):
        """Verify that we can load rotamers from PDB+DCD."""
        # The reference libraries are IMP.bff module data (data/rotamer_library).
        import IMP.bff
        lib_dir = IMP.bff.get_data_path("rotamer_library")

        pdb, dcd, weights = find_reference_rotamer_files(
            lib_dir, "A48_C1R", cutoff=10
        )
        data = load_reference_rotamers(pdb, dcd, weights, max_frames=8)

        self.assertEqual(data["coords"].shape[0], 8)
        self.assertEqual(len(data["weights"]), 8)
        self.assertAlmostEqual(float(data["weights"].sum()), 1.0)

    def test_apply_rotamer_coords(self):
        """Verify coordinates are correctly applied to an IMP hierarchy."""
        model = IMP.Model()
        p1 = IMP.Particle(model, "at1")
        IMP.core.XYZ.setup_particle(p1)
        IMP.atom.Atom.setup_particle(p1, IMP.atom.AtomType("C"))
        h = IMP.atom.Hierarchy.setup_particle(p1)

        new_coords = np.array([[10.0, 20.0, 30.0]])
        apply_rotamer_coords(h, new_coords)

        xyz = IMP.core.XYZ(p1).get_coordinates()
        self.assertAlmostEqual(xyz[0], 10.0)
        self.assertAlmostEqual(xyz[1], 20.0)
        self.assertAlmostEqual(xyz[2], 30.0)

    def test_sample_rotamer_index(self):
        """Verify sampling from weights."""
        weights = [0.0, 1.0, 0.0]
        idx = sample_rotamer_index(weights)
        self.assertEqual(idx, 1)

        weights = [0.1, 0.1, 0.8]
        samples = [sample_rotamer_index(weights) for _ in range(100)]
        self.assertIn(2, samples)


if __name__ == "__main__":
    unittest.main()
