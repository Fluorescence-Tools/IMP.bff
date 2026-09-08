"""`IMP.bff.labellib`: LabelLib's API, name for name.

LabelLib is the small library that computes an accessible volume and the
observables read off one. bff computes the same volumes, so a caller should
not have to install both, nor rewrite a working script to move over:
`from IMP.bff import labellib as ll` is meant to stand in for
`import LabelLib as ll` unchanged.

Two things are checked. The **interface**, always: every LabelLib name
exists, takes LabelLib's arguments in LabelLib's order, accepts atoms as
(4, N) and as (N, 4), and returns a `Grid3D` with LabelLib's attributes and
a `points()` of shape (4, n). And the **numbers**, wherever the real
LabelLib is importable: the two agree on the observables to within the
discretisation -- they are different lattice searches over the same physics,
not the same code, so the volumes differ voxel by voxel while the mean
distance and the mean efficiency agree closely.
"""

import unittest

import numpy as np

import IMP.bff
from IMP.bff import labellib as ll

try:
    import LabelLib as _reference
except ImportError:
    _reference = None


# five atoms, not four: a (4, 4) is ambiguous and LabelLib reads it as
# (4, N), so a square array would not test the orientation at all
ATOMS = np.array([[0.0, 0.0, 0.0, 1.5],
                  [5.0, 0.0, 0.0, 1.5],
                  [0.0, 5.0, 0.0, 1.5],
                  [-4.0, -3.0, 2.0, 1.7],
                  [2.0, -6.0, 1.0, 1.6]])
SOURCE = [0.0, 0.0, 3.0]
FAR_SOURCE = [12.0, 0.0, 3.0]
L, W, R, STEP = 10.0, 2.0, 3.5, 1.0


class TestTheNames(unittest.TestCase):

    def test_every_labellib_name_is_here(self):
        for name in ("Grid3D", "dyeDensityAV1", "dyeDensityAV3",
                     "minLinkerLength", "addWeights", "meanDistance",
                     "meanEfficiency", "sampleDistanceDistInv"):
            self.assertTrue(hasattr(ll, name), name)
            self.assertTrue(hasattr(ll, name + "_arr") or name == "Grid3D", name)

    def test_it_is_importable_as_a_module(self):
        import sys
        self.assertIn("IMP.bff.labellib", sys.modules)
        from IMP.bff import labellib
        self.assertIs(labellib, ll)


class TestTheGrid(unittest.TestCase):

    def setUp(self):
        self.g = ll.dyeDensityAV1(ATOMS.T, SOURCE, L, W, R, STEP)

    def test_a_grid_carries_labellibs_attributes(self):
        self.assertEqual(len(self.g.shape), 3)
        self.assertEqual(len(self.g.originXYZ), 3)
        self.assertEqual(self.g.discStep, STEP)
        self.assertEqual(len(self.g.grid),
                         self.g.shape[0] * self.g.shape[1] * self.g.shape[2])

    def test_points_are_four_by_n_and_agree_with_the_grid(self):
        p = self.g.points()
        self.assertEqual(p.shape[0], 4)
        self.assertGreater(p.shape[1], 0)
        # x fastest, as LabelLib's flat grid is: rebuilding the coordinates
        # from the cube must reproduce points() exactly
        cube = np.asarray(self.g.grid).reshape(self.g.shape, order="F")
        idx = np.argwhere(cube > 0)
        xyz = np.asarray(self.g.originXYZ) + idx * self.g.discStep
        self.assertEqual(sorted(map(tuple, np.round(xyz, 6))),
                         sorted(map(tuple, np.round(p[:3].T, 6))))

    def test_the_cloud_clears_the_obstacles(self):
        p = self.g.points()[:3].T
        for x, y, z, r in ATOMS:
            d = np.linalg.norm(p - np.array([x, y, z]), axis=1)
            self.assertGreater(d.min(), r, "a point sits inside an atom")


class TestTheArrayConventions(unittest.TestCase):

    def test_atoms_are_taken_either_way_round(self):
        a = ll.dyeDensityAV1(ATOMS.T, SOURCE, L, W, R, STEP)
        b = ll.dyeDensityAV1(ATOMS, SOURCE, L, W, R, STEP)
        self.assertEqual(a.shape, b.shape)
        np.testing.assert_allclose(a.grid, b.grid)

    def test_a_wrong_shape_is_refused(self):
        with self.assertRaises(ValueError):
            ll.dyeDensityAV1(np.zeros((5, 7)), SOURCE, L, W, R, STEP)

    def test_av3_takes_three_radii(self):
        g = ll.dyeDensityAV3(ATOMS.T, SOURCE, L, W, [3.5, 4.0, 2.0], STEP)
        self.assertGreater((np.asarray(g.grid) > 0).sum(), 0)
        with self.assertRaises(ValueError):
            ll.dyeDensityAV3(ATOMS.T, SOURCE, L, W, [3.5, 4.0], STEP)


class TestTheObservables(unittest.TestCase):

    def setUp(self):
        self.g1 = ll.dyeDensityAV1(ATOMS.T, SOURCE, L, W, R, STEP)
        self.g2 = ll.dyeDensityAV1(ATOMS.T, FAR_SOURCE, L, W, R, STEP)

    def test_mean_distance_is_between_the_extremes(self):
        d = ll.meanDistance(self.g1, self.g2, 200000)
        p1, p2 = self.g1.points()[:3].T, self.g2.points()[:3].T
        alld = np.linalg.norm(p1[:, None, :] - p2[None, :, :], axis=2)
        self.assertGreater(d, alld.min())
        self.assertLess(d, alld.max())

    def test_efficiency_falls_with_the_forster_radius(self):
        e_near = ll.meanEfficiency(self.g1, self.g2, 60.0)
        e_far = ll.meanEfficiency(self.g1, self.g2, 20.0)
        self.assertGreater(e_near, e_far)
        for e in (e_near, e_far):
            self.assertGreaterEqual(e, 0.0)
            self.assertLessEqual(e, 1.0)

    def test_samples_are_distances(self):
        s = np.asarray(ll.sampleDistanceDistInv(self.g1, self.g2, 5000))
        self.assertEqual(s.ndim, 1)
        self.assertEqual(len(s), 5000)
        # the two clouds overlap here, so a zero distance is a real draw
        self.assertGreaterEqual(s.min(), 0.0)
        self.assertGreater(s.max(), 0.0)

    def test_add_weights_raises_the_weights_it_covers(self):
        before = np.asarray(self.g1.grid)
        after = np.asarray(ll.addWeights(
            self.g1, np.array([[0.0, -5.0, -5.0, 6.0, 2.0]]).T).grid)
        changed = after != before
        self.assertGreater(changed.sum(), 0)
        np.testing.assert_allclose(after[changed] - before[changed], 2.0)
        # only occupied voxels take a weight
        self.assertTrue((before[changed] > 0).all())

    def test_min_linker_length_reaches_and_refuses(self):
        m = ll.minLinkerLength(ATOMS.T, SOURCE, L, W, R, STEP)
        v = np.asarray(m.grid)
        self.assertEqual(m.shape, self.g1.shape)
        self.assertLessEqual(v.max(), L + m.discStep)
        self.assertGreater((v > 0).sum(), 0)
        self.assertGreater((v < 0).sum(), 0, "unreachable voxels are negative")
        # every voxel of the volume is a voxel the linker reaches
        cube_av = np.asarray(self.g1.grid).reshape(self.g1.shape, order="F")
        cube_ml = v.reshape(m.shape, order="F")
        self.assertTrue((cube_ml[cube_av > 0] >= 0).all())


@unittest.skipIf(_reference is None, "LabelLib is not installed here")
class TestAgainstLabelLib(unittest.TestCase):
    """The numbers, where the real LabelLib is importable.

    Two different lattice searches over the same physics: the clouds differ
    voxel by voxel (different windows, different stencils), so what is
    compared is what a caller reads off them.
    """

    def setUp(self):
        a = np.ascontiguousarray(ATOMS.T, dtype=np.float32)
        s1 = np.asarray(SOURCE, dtype=np.float32)
        s2 = np.asarray(FAR_SOURCE, dtype=np.float32)
        self.ref1 = _reference.dyeDensityAV1(a, s1, L, W, R, STEP)
        self.ref2 = _reference.dyeDensityAV1(a, s2, L, W, R, STEP)
        self.g1 = ll.dyeDensityAV1(ATOMS.T, SOURCE, L, W, R, STEP)
        self.g2 = ll.dyeDensityAV1(ATOMS.T, FAR_SOURCE, L, W, R, STEP)

    def test_the_grid_attributes_have_labellibs_meaning(self):
        for ours, theirs in ((self.g1, self.ref1),):
            self.assertEqual(len(ours.shape), len(theirs.shape))
            self.assertEqual(ours.discStep, theirs.discStep)
            self.assertEqual(ours.points().shape[0], theirs.points().shape[0])

    def test_the_clouds_occupy_the_same_region(self):
        ours = self.g1.points()[:3].T
        theirs = np.asarray(self.ref1.points())[:3].T
        np.testing.assert_allclose(ours.mean(axis=0), theirs.mean(axis=0), atol=1.0)

    def test_the_mean_distance_agrees(self):
        ours = ll.meanDistance(self.g1, self.g2, 200000)
        theirs = _reference.meanDistance(self.ref1, self.ref2, 200000)
        self.assertLess(abs(ours - theirs), 0.5, (ours, theirs))

    def test_the_mean_efficiency_agrees(self):
        ours = ll.meanEfficiency(self.g1, self.g2, 50.0, 200000)
        theirs = _reference.meanEfficiency(self.ref1, self.ref2, 50.0, 200000)
        self.assertLess(abs(ours - theirs), 0.02, (ours, theirs))


if __name__ == "__main__":
    unittest.main()
