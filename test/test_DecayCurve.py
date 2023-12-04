from __future__ import division
import unittest

import numpy as np
import numpy.testing
import math
import IMP.bff


def norm_pdf(x, mu, sigma):
    variance = sigma**2
    num = x - mu
    denom = 2*variance
    pdf = ((1/(np.sqrt(2*np.pi)*sigma))*np.exp(-(num**2)/denom))
    return pdf


class Tests(unittest.TestCase):

    def test_shift(self):
        time_axis = np.linspace(0, 12, 25)
        irf_position = 6.0
        irf_width = 1.0
        irf = norm_pdf(time_axis, irf_position, irf_width)
        shifts = [2.1, 1.1, 0.5, 0.0, -0.5, -1.1, -2.1, 26.0, -26.2]
        for shift in shifts:
            shifted_irf_ip = IMP.bff.DecayCurve.shift_array(
                input=irf,
                shift=shift
            )
            shifted_irf_ip = np.array(shifted_irf_ip)
            i0 = math.floor(shift)
            i1 = math.ceil(shift)
            f_shift = shift - i0
            ref = np.roll(irf, -i0) * (1 - f_shift) + f_shift * np.roll(irf, -i1)
            np.testing.assert_allclose(shifted_irf_ip, ref)

    def test_shift_2(self):
        x = np.arange(0, 8, dtype=np.float64)
        y = np.zeros_like(x)
        y[5] = 1.0
        dc = IMP.bff.DecayCurve(x, y)
        self.assertEqual(dc.shift, 0.0)
        dc.shift = 2
        np.testing.assert_allclose(np.roll(y, -2), dc.y)
        dc.shift = 2.5
        np.testing.assert_allclose(np.array([0., 0., 0.5, 0.5, 0., 0., 0., 0.]), dc.y)

    def test_DecayCurve_init(self):
        x = np.linspace(0, 10, 10)
        y = np.sin(x) * 1000
        curve1 = IMP.bff.DecayCurve(x, y)
        curve2 = IMP.bff.DecayCurve(x=x, y=y)
        np.testing.assert_allclose(curve1.x, curve2.x)
        np.testing.assert_allclose(curve1.y, curve2.y)
        # Noise is by default Poisson noise
        ey = np.clip(np.sqrt(np.abs(y)), 1, 1e12)
        curve1 = IMP.bff.DecayCurve(x, y)
        np.testing.assert_allclose(curve1.ey, ey)

    def test_DecayCurve_operator(self):
        x = np.linspace(0, 10, 10)
        y = np.sin(x)
        ey = np.clip(np.sqrt(np.abs(y)), 1, 1e12)

        # Multiplication
        curve1 = IMP.bff.DecayCurve(x, y)
        curve2 = curve1 * curve1
        curve3 = curve1 * 2
        curve1 *= 2
        np.testing.assert_allclose(curve1.y, y * 2)
        np.testing.assert_allclose(curve1.ey, ey * 2)
        np.testing.assert_allclose(curve3.y, curve1.y)
        np.testing.assert_allclose(curve2.x, curve1.x)
        np.testing.assert_allclose(curve3.ey, curve1.ey)
        np.testing.assert_allclose(curve2.y, y * y)

        # Addition
        curve1 = IMP.bff.DecayCurve(x, y)
        curve3 = curve1 + curve1
        curve2 = curve1 + 2
        curve1 += 2
        np.testing.assert_allclose(curve1.y, y + 2)
        np.testing.assert_allclose(curve1.ey, ey)
        np.testing.assert_allclose(curve2.y, curve1.y)
        np.testing.assert_allclose(curve2.x, curve1.x)
        np.testing.assert_allclose(curve2.ey, curve1.ey)
        np.testing.assert_allclose(curve3.y, y + y)

    def test_decay_curve_resize(self):
        d = IMP.bff.DecayCurve()
        d.x = np.array([0, 0.1])
        d.resize(5)
        np.testing.assert_allclose(d.x, np.array([0. , 0.1, 0.2, 0.3, 0.4]))
        d.resize(6)
        np.testing.assert_allclose(d.x, np.array([0. , 0.1, 0.2, 0.3, 0.4, 0.5]))
