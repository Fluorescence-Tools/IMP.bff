import unittest
import numpy as np
import IMP
import IMP.core
import IMP.algebra
import tempfile
import os

import pytest

@pytest.mark.usefixtures("imp_bff_program")
class TestAnalysis(unittest.TestCase):
    """The density analysis is a program, so the fixture loads `bin/imp_bff`."""

    @pytest.fixture(autouse=True)
    def _program(self, imp_bff_program):
        self.prog = imp_bff_program

    def test_compute_fixed_axis(self):
        model = IMP.Model()
        # Create a flat ring of particles in the XY plane
        particles = []
        for i in range(4):
            p = IMP.Particle(model)
            xyz = IMP.core.XYZ.setup_particle(p)
            # Coordinates: (1,0,0), (0,1,0), (-1,0,0), (0,-1,0)
            x = 1.0 if i == 0 else (-1.0 if i == 2 else 0.0)
            y = 1.0 if i == 1 else (-1.0 if i == 3 else 0.0)
            xyz.set_coordinates(IMP.algebra.Vector3D(x, y, 0.0))
            particles.append(p)
        
        # The minimal variance axis for a flat XY ring should be Z (0,0,1)
        center, axis = self.prog._compute_fixed_axis(particles)
        
        self.assertAlmostEqual(center[0], 0.0)
        self.assertAlmostEqual(center[1], 0.0)
        self.assertAlmostEqual(center[2], 0.0)
        
        # Axis should be [0,0,1] or [0,0,-1]
        self.assertAlmostEqual(abs(axis[2]), 1.0)
        self.assertAlmostEqual(axis[0], 0.0)
        self.assertAlmostEqual(axis[1], 0.0)

    def test_compute_long_axis(self):
        model = IMP.Model()
        # Create a line of particles along the X axis
        particles = []
        for x in [-5.0, 5.0]:
            p = IMP.Particle(model)
            xyz = IMP.core.XYZ.setup_particle(p)
            xyz.set_coordinates(IMP.algebra.Vector3D(x, 0.0, 0.0))
            particles.append(p)
            
        # Long axis (max variance) should be X (1,0,0)
        axis = self.prog._compute_long_axis(particles)
        self.assertAlmostEqual(abs(axis[0]), 1.0)
        self.assertAlmostEqual(axis[1], 0.0)
        self.assertAlmostEqual(axis[2], 0.0)

    def test_write_radial_histogram(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_path = os.path.join(tmpdir, "hist.csv")
            distances = [1.0, 2.0, 3.0, 1.5]
            self.prog.write_radial_histogram(distances, csv_path, bin_width=1.0)
            
            self.assertTrue(os.path.exists(csv_path))
            with open(csv_path, 'r') as f:
                lines = f.readlines()
                # Header + 3 bins (0-1, 1-2, 2-3) - wait, max is 3.0, so 0-1, 1-2, 2-3.
                # Actually min is 1.0, max is 3.0. 
                # bins: [0,1), [1,2), [2,3]
                self.assertGreater(len(lines), 1)

if __name__ == "__main__":
    unittest.main()
