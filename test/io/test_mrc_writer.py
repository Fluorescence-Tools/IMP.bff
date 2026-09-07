"""The core writes a correct MRC2014 map -- checked with an independent reader.

`write_mrc` (DensityGrid.h) replaced the IMP::em writers behind
`write_map_feature` (PRD-137). The gate is not "the same bytes as IMP" but
"a map any MRC reader understands and that says what the grid is": mode 2,
cell lengths of `n * spacing`, zero starts with the origin in ORIGIN, real
statistics, the MAP tag, a machine stamp, NVERSION 20140. `mrcfile` is the
independent reader; IMP.em is a second one when present.
"""
import os
import struct
import tempfile
import unittest

import numpy as np

import IMP.algebra
import IMP.bff as bff

try:
    import mrcfile
    _HAVE_MRCFILE = True
except ImportError:
    _HAVE_MRCFILE = False
try:
    import IMP.em
    _HAVE_EM = True
except ImportError:
    _HAVE_EM = False


def _grid(n=(6, 5, 4), spacing=0.5, origin=(-1.5, 2.0, 3.25)):
    g = bff.DensityGrid("g")
    h = g.get_header_writable()
    h.set_spacing(spacing)
    h.update_map_dimensions(*n)
    g.resize(g.get_number_of_voxels())
    g.set_origin(IMP.algebra.Vector3D(*origin))
    # a value that tells x, y and z apart: v = x + 10 y + 100 z
    for v in range(g.get_number_of_voxels()):
        ix, iy, iz = (g.get_dim_index_by_voxel(v, d) for d in range(3))
        g.set_value(v, float(ix + 10 * iy + 100 * iz))
    return g


class TestMrcWriter(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.path = os.path.join(self.tmp, "grid.mrc")
        self.grid = _grid()
        self.grid.write_mrc(self.path)

    def test_the_header_is_mrc2014(self):
        with open(self.path, "rb") as fh:
            raw = fh.read()
        words = struct.unpack("<56i", raw[:224])
        floats = struct.unpack("<56f", raw[:224])
        self.assertEqual(words[0:3], (6, 5, 4))                  # NX NY NZ
        self.assertEqual(words[3], 2)                             # MODE
        self.assertEqual(words[4:7], (0, 0, 0))                   # starts
        self.assertEqual(words[7:10], (6, 5, 4))                  # MX MY MZ
        np.testing.assert_allclose(floats[10:13], (3.0, 2.5, 2.0))  # cell = n * spacing
        np.testing.assert_allclose(floats[13:16], (90.0, 90.0, 90.0))
        self.assertEqual(words[16:19], (1, 2, 3))                 # MAPC MAPR MAPS
        self.assertEqual(words[22], 1)                            # ISPG: a volume
        self.assertEqual(words[27], 20140)                        # NVERSION (word 28)
        np.testing.assert_allclose(floats[49:52], (-1.5, 2.0, 3.25))  # ORIGIN
        self.assertEqual(raw[208:212], b"MAP ")
        self.assertEqual(raw[212:214], b"\x44\x44")                # little-endian stamp
        self.assertEqual(words[55], 1)                            # one label
        self.assertEqual(len(raw), 1024 + 4 * 6 * 5 * 4)
        values = np.array([self.grid.get_value(v) for v in range(self.grid.get_number_of_voxels())])
        np.testing.assert_allclose(floats[19:22], (values.min(), values.max(), values.mean()), rtol=1e-6)
        np.testing.assert_allclose(floats[54], values.std(), rtol=1e-6)

    @unittest.skipUnless(_HAVE_MRCFILE, "mrcfile not installed")
    def test_mrcfile_reads_it_back_with_the_axes_the_right_way_round(self):
        self.assertTrue(mrcfile.validate(self.path))
        with mrcfile.open(self.path) as m:
            data = m.data  # (nz, ny, nx)
            self.assertEqual(data.shape, (4, 5, 6))
            self.assertAlmostEqual(float(m.voxel_size.x), 0.5, places=6)
            self.assertAlmostEqual(float(m.voxel_size.z), 0.5, places=6)
            np.testing.assert_allclose((m.header.origin.x, m.header.origin.y, m.header.origin.z), (-1.5, 2.0, 3.25))
            # v = x + 10 y + 100 z at (ix, iy, iz): a transposed write would fail this
            self.assertEqual(float(data[3, 2, 1]), 1 + 20 + 300)
            self.assertEqual(float(data[0, 4, 5]), 5 + 40)

    @unittest.skipUnless(_HAVE_EM, "IMP.em not available")
    def test_imp_em_reads_it_back(self):
        dm = IMP.em.read_map(self.path, IMP.em.MRCReaderWriter())
        h = dm.get_header()
        self.assertEqual((h.get_nx(), h.get_ny(), h.get_nz()), (6, 5, 4))
        self.assertAlmostEqual(h.get_spacing(), 0.5, places=6)
        np.testing.assert_allclose((h.get_xorigin(), h.get_yorigin(), h.get_zorigin()), (-1.5, 2.0, 3.25))
        for v in (0, 7, 119):
            self.assertAlmostEqual(dm.get_value(v), self.grid.get_value(v), places=5)

    def test_write_map_feature_writes_mrc_only(self):
        # the PathMap door writes MRC through the same writer; other suffixes
        # are the em bridge's business and say so
        pm = bff.PathMap(bff.PathMapHeader(5.0, 1.0))
        with self.assertRaises(Exception):
            bff.write_map_feature(pm, os.path.join(self.tmp, "x.xplor"), bff.PM_TILE_DENSITY, (0.0, 5.0))


if __name__ == "__main__":
    unittest.main()
