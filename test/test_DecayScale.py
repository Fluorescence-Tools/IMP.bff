from __future__ import division
import unittest

import numpy as np
import numpy.testing
import IMP.bff

x = np.arange(0, 16)
y = np.ones_like(x) * 100
y2 = y * 10
data = IMP.bff.DecayCurve(x, y2)
settings = {
    "data": data,
    "constant_background": 0,
    "active": True,
    "start": 0,
    "stop": -1
}


class Tests(unittest.TestCase):

    def test_scale_init(self):
        ds = IMP.bff.DecayScale(**settings)
        self.assertEqual(np.allclose(ds.data.y, data.y), True)
        self.assertEqual(np.allclose(ds.data.x, data.x), True)
        self.assertEqual(ds.constant_background, 0.0)
        self.assertEqual(ds.active, True)

    def test_scale_scale_1(self):
        import IMP.bff
        x = np.arange(0, 16)
        y = np.ones_like(x) * 100
        y2 = y * 10

        model = IMP.bff.DecayCurve(x, y)
        data = IMP.bff.DecayCurve(x, y2)

        settings = {
            "data": data,
            "constant_background": 0,
            "active": True,
            "start": 0,
            "stop": -1
        }
        ds = IMP.bff.DecayScale(**settings)
        ds(model)
        np.testing.assert_allclose(model.y, data.y)

        # if ds.active = False, nothing happend
        ds.active = False
        self.assertEqual(ds.active, False)
        model = IMP.bff.DecayCurve(x, y)
        ds(model)
        self.assertEqual(np.allclose(model.y, data.y), False)
        np.testing.assert_allclose(model.y, y)

    def test_scale_scale_2(self):
        import numpy as np
        import IMP.bff

        x = np.arange(0, 16)
        y = np.ones_like(x) * 100
        y2 = y * 10
        data = IMP.bff.DecayCurve(x, y2)
        settings = {
            "data": data,
            "constant_background": 0,
            "active": True,
            "start": 0,
            "stop": -1
        }

        ds = IMP.bff.DecayScale(**settings)
        model = IMP.bff.DecayCurve(x, y)
        np.testing.assert_allclose(ds.data.y, np.ones_like(ds.data.y) * 1000)
        np.testing.assert_allclose(model.y, np.ones_like(ds.data.y) * 100)

        ds(model)  # equivalent to ds.add(model)
        np.testing.assert_almost_equal(np.sum(model.y), np.sum(ds.data.y))

    def test_scale_scale_3(self):
        ds = IMP.bff.DecayScale(**settings)
        # Scaling considering noise & background of data
        ds.constant_background = 100
        # Set number of photons to negative value for scaling
        model = IMP.bff.DecayCurve(x, y)
        ds.add(model)
        np.testing.assert_allclose(model.y, ds.data.y - 100)
        np.testing.assert_almost_equal(ds.number_of_photons, 16000)
