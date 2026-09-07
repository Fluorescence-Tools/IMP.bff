"""`DensityGrid`: the lattice `PathMap` is built on, and the AV output it must not move.

`PathMap` derived from `IMP::em::SampledDensityMap` and is now built on
`IMP.bff.DensityGrid`, a lattice of this module's own. That swap was made
under one rule -- **the numbers do not change** -- and this file is where the
rule is enforced rather than remembered.

Two kinds of test, and they catch different things.

**The golden AV records** pin `compute_av_from_structure` on three structures
at three grid steps, as SHA-256 of the density and of the point cloud, taken
from the `IMP.em`-backed build before the swap. They are what a change to the
sampling, the bounding box or the voxel/location arithmetic shows up in. They
did not, on their own, catch four regressions the swap introduced -- which is
what the second kind is for.

**The unit tests** each pin one of those four, at the level it broke:

* `set_origin` must rebuild the coordinate caches, not only invalidate them.
  They are `unique_ptr`s, so the failure is a null dereference, not a wrong
  number, and it surfaced only in the lattice fast path.
* `PathMap` must keep the methods it inherits. SWIG gives a class the methods
  of a base it has *seen*; when the base was not in the interface `PathMap`
  wrapped, imported, and had no `get_number_of_voxels`.
* The header is float-backed. `IMP::em::DensityHeader` stores spacing and
  origin as `float`; holding them as `double` computed voxel locations *more*
  accurately and moved AV mean positions in the eighth significant figure.
* Voxel-index arithmetic is IMP's, promotion for promotion. Two extra float
  roundings in `get_dim_index_by_location` pick a different voxel at a
  half-voxel boundary, and the legacy lattice search seeds itself there.
"""

import hashlib
import os
import unittest

import numpy as np

import IMP.algebra
import IMP.bff as bff

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _sha(a):
    return hashlib.sha256(
        np.ascontiguousarray(a, dtype=np.float64).tobytes()).hexdigest()[:16]


# (label, pdb, chain, resseq, atom, linker_length, linker_width, radius1)
_SITES = [
    ("T4L-132", "examples/structure/T4L/3GUN.pdb", "", 132, "CB", 22.0, 3.5, 3.5),
    ("T4L-150", "examples/structure/T4L/3GUN.pdb", "", 150, "CB", 22.0, 2.5, 3.5),
    ("GBP-400", "examples/structure/GBP/mGBP2A.pdb", "A", 400, "CB", 22.0, 4.5, 3.5),
]

# Recorded 2026-09-07 from the IMP.em-backed build, before PathMap moved onto
# DensityGrid. key -> (n_density, density_sha, n_points, points_sha, ng).
_GOLDEN = {
    "T4L-132@1.5": (29791, "39742a2572f3a12b", 14356, "269334fbb8e053f2", 31),
    "T4L-132@1.0": (91125, "2536d73005254013", 53120, "90063dc25c2c4327", 45),
    "T4L-132@0.5": (704969, "dd4b08907ac0f84b", 331156, "a3304f75b05374f6", 89),
    "T4L-150@1.5": (29791, "368f4274767d1eca", 10060, "818e27628107e959", 31),
    "T4L-150@1.0": (91125, "f357118e3ab9fb87", 30676, "8b2c027f91058e72", 45),
    "T4L-150@0.5": (704969, "8e539b6d042fc7d4", 248400, "5525fe5619ebf10a", 89),
    "GBP-400@1.5": (29791, "07a95aeb323c606d", 13916, "4a0b177fab0c45ea", 31),
    "GBP-400@1.0": (91125, "071cbd6b2870ed1a", 52448, "74138fbfdd50abb5", 45),
    "GBP-400@0.5": (704969, "268ae0b87ec6b037", 371144, "bb3a458acc320dea", 89),
}


class TestGoldenAV(unittest.TestCase):
    """The AV solver's output, bit for bit, across the structures shipped."""

    def test_av_records_are_unchanged(self):
        for label, rel, chain, resseq, atom, ll, lw, r1 in _SITES:
            pdb = os.path.join(_ROOT, rel)
            if not os.path.exists(pdb):
                self.skipTest("no %s in this checkout" % rel)
            for step in (1.5, 1.0, 0.5):
                key = "%s@%.1f" % (label, step)
                with self.subTest(key=key):
                    av = bff.compute_av_from_structure(
                        pdb, chain, resseq, atom, ll, lw, r1, 0.0, 0.0, step)
                    d = np.asarray(av.get_density())
                    p = np.asarray(av.get_points())
                    n_d, sha_d, n_p, sha_p, ng = _GOLDEN[key]
                    self.assertEqual(int(av.get_ng()), ng)
                    self.assertEqual(d.size, n_d)
                    self.assertEqual(p.size, n_p)
                    self.assertEqual(_sha(d), sha_d, "density moved: " + key)
                    self.assertEqual(_sha(p), sha_p, "points moved: " + key)


def _grid(n=10, spacing=1.0, origin=(0.0, 0.0, 0.0)):
    g = bff.DensityGrid("g")
    h = g.get_header_writable()
    h.set_spacing(spacing)
    h.update_map_dimensions(n, n, n)
    g.resize(g.get_number_of_voxels())
    g.set_origin(IMP.algebra.Vector3D(*origin))
    return g


def _occupied(g):
    return sum(g.get_value(i) for i in range(g.get_number_of_voxels()))


class TestSampling(unittest.TestCase):
    """IMP's `BINARIZED_SPHERE`, reproduced: strict, accumulating."""

    def test_sphere_of_radius_two_marks_the_27_interior_points(self):
        # Lattice points strictly inside r=2 about a lattice point: every
        # (dx,dy,dz) in {-1,0,1}^3 has dx^2+dy^2+dz^2 <= 3 < 4, and any
        # coordinate of 2 gives >= 4. So 27, and nothing else.
        g = _grid()
        g.set_spheres([bff.GridSphere(IMP.algebra.Vector3D(5, 5, 5), 2.0)])
        g.resample()
        self.assertEqual(_occupied(g), 27)

    def test_the_boundary_is_strict(self):
        # r == 1 exactly: the six neighbours sit at distance 1, and 1 < 1 is
        # false, so only the centre is inside. IMP's kernel is `<`, not `<=`.
        g = _grid()
        g.set_spheres([bff.GridSphere(IMP.algebra.Vector3D(5, 5, 5), 1.0)])
        g.resample()
        self.assertEqual(_occupied(g), 1)

    def test_overlap_accumulates(self):
        # Two identical spheres: every voxel inside reads 2, not 1. IMP added
        # (`data[ivox] += value`) rather than assigned, and nothing reading the
        # map cares -- but the threshold that follows is `>`, so a change here
        # to assignment would still pass a naive occupancy count.
        g = _grid()
        s = bff.GridSphere(IMP.algebra.Vector3D(5, 5, 5), 2.0)
        g.set_spheres([s, s])
        g.resample()
        values = {g.get_value(i) for i in range(g.get_number_of_voxels())}
        self.assertEqual(values, {0.0, 2.0})

    def test_zero_radius_blocks_nothing(self):
        g = _grid()
        g.set_spheres([bff.GridSphere(IMP.algebra.Vector3D(5, 5, 5), 0.0)])
        g.resample()
        self.assertEqual(_occupied(g), 0)


class TestIndexing(unittest.TestCase):
    """Voxel <-> index <-> location, with IMP's rounding."""

    def test_index_roundtrip(self):
        g = _grid(n=7)
        for ix, iy, iz in ((0, 0, 0), (6, 6, 6), (3, 1, 5)):
            v = g.xyz_ind2voxel(ix, iy, iz)
            self.assertEqual([g.get_dim_index_by_voxel(v, d) for d in range(3)],
                             [ix, iy, iz])

    def test_location_roundtrip(self):
        g = _grid(n=7, spacing=1.5, origin=(-3.0, 2.0, 0.5))
        for v in (0, 17, 342):
            loc = g.get_location_by_voxel(v)
            self.assertEqual(g.get_voxel_by_location(loc), v)

    def test_half_voxel_rounds_up_as_imp_did(self):
        # IMP: floor(0.5 + (loc - origin) / spacing). Half a voxel past the
        # origin is index 1; just under half is index 0. The subtraction and
        # division are done in double after the float inputs, which is what
        # keeps a coordinate a hair under the boundary from flipping.
        g = _grid(n=7, spacing=2.0, origin=(0.0, 0.0, 0.0))
        self.assertEqual(g.get_dim_index_by_location(1.0, 0), 1)
        self.assertEqual(g.get_dim_index_by_location(0.98, 0), 0)
        self.assertEqual(g.get_dim_index_by_location(2.99, 0), 1)
        self.assertEqual(g.get_dim_index_by_location(3.0, 0), 2)


class TestTheFourRegressions(unittest.TestCase):
    """One test per bug the swap shipped past the golden records."""

    def test_set_origin_rebuilds_the_caches(self):
        # Invalidate-only left `x_loc_` null; reading a location afterwards
        # was a crash. The lattice fast path moves the origin and reads
        # without asking for a recompute in between.
        g = _grid(n=5)
        g.set_origin(IMP.algebra.Vector3D(10.0, 20.0, 30.0))
        loc = g.get_location_by_voxel(0)
        np.testing.assert_array_equal(np.asarray(list(loc)), [10.0, 20.0, 30.0])

    def test_pathmap_keeps_its_inherited_api(self):
        for name in ("get_number_of_voxels", "get_value", "get_spacing",
                     "get_location_by_voxel", "get_dim_index_by_voxel",
                     "get_voxel_by_location", "xyz_ind2voxel", "get_header",
                     "resample"):
            self.assertTrue(hasattr(bff.PathMap, name), name)

    def test_header_is_float_backed(self):
        # 0.1 is not representable; a float-backed header hands back the
        # float32 rounding, a double-backed one would hand back 0.1.
        h = _grid().get_header_writable()
        h.set_spacing(0.1)
        self.assertEqual(h.get_spacing(), float(np.float32(0.1)))
        self.assertNotEqual(h.get_spacing(), 0.1)

    def test_pathmapheader_carries_a_gridheader(self):
        # The em header it used to carry serialized ~40 EM fields; this one
        # serializes eleven. `get_density_header` is the by-pointer accessor
        # IMP's was, so the extent is readable from Python.
        pmh = bff.PathMapHeader()
        h = pmh.get_density_header()
        self.assertIsInstance(h, bff.GridHeader)


if __name__ == "__main__":
    unittest.main()
