import utils
import os
import unittest
import sys
import json

import numpy as np
import IMP
import IMP.bff

from constants import *

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)


class Tests(unittest.TestCase):

    def make_AccessibleVolume(self, voxel_size):
        m = IMP.Model()
        ps = [IMP.Particle(m)]
        IMP.core.XYZR.setup_particle(ps[0], IMP.algebra.Sphere3D(IMP.algebra.Vector3D(1, 2, 3), 4))
        IMP.atom.Mass.setup_particle(ps[0], 1.0)
        av = IMP.bff.AccessibleVolume(ps, 10.0, voxel_size)
        return m, av

    # def test_AccessibleVolume_init(self):
    #     name_test = "TestName"
    #     av1 = IMP.bff.AccessibleVolume()
    #     av2 = IMP.bff.AccessibleVolume(name_test)
    #     self.assertNotEqual(av1.get_name(), name_test)
    #     self.assertEqual(av2.get_name(), name_test)



if __name__ == '__main__':
    unittest.main()
