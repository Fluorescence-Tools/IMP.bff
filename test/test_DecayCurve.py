from __future__ import division
import unittest

import numpy as np
import numpy.testing
import scipy.stats
import IMP.bff


class Tests(unittest.TestCase):

    def test_shift(self):
        time_axis = np.linspace(0, 12, 25)
        irf_position = 6.0
        irf_width = 1.0
        irf = scipy.stats.norm.pdf(time_axis, loc=irf_position, scale=irf_width)
        # integer shift
        shifted_irf_ip = IMP.bff.DecayCurve.shift_array(
            input=irf,
            shift=1.0
        )
        ref = np.array(
            [0.00000000e+00, 1.48671951e-06, 1.59837411e-05, 1.33830226e-04,
             8.72682695e-04, 4.43184841e-03, 1.75283005e-02, 5.39909665e-02,
             1.29517596e-01, 2.41970725e-01, 3.52065327e-01, 3.98942280e-01,
             3.52065327e-01, 2.41970725e-01, 1.29517596e-01, 5.39909665e-02,
             1.75283005e-02, 4.43184841e-03, 8.72682695e-04, 1.33830226e-04,
             1.59837411e-05, 1.48671951e-06, 1.07697600e-07, 6.07588285e-09,
             6.07588285e-09]
        )        
        np.testing.assert_allclose(shifted_irf_ip, ref)
        shifted_irf_in = IMP.bff.DecayCurve.shift_array(
            input=irf,
            shift=-1.0
        )
        ref = np.array(
            [6.07588285e-09, 6.07588285e-09, 1.07697600e-07, 1.48671951e-06,
             1.59837411e-05, 1.33830226e-04, 8.72682695e-04, 4.43184841e-03,
             1.75283005e-02, 5.39909665e-02, 1.29517596e-01, 2.41970725e-01,
             3.52065327e-01, 3.98942280e-01, 3.52065327e-01, 2.41970725e-01,
             1.29517596e-01, 5.39909665e-02, 1.75283005e-02, 4.43184841e-03,
             8.72682695e-04, 1.33830226e-04, 1.59837411e-05, 1.48671951e-06,
             1.07697600e-07]
        )
        np.testing.assert_allclose(shifted_irf_in, ref)

        # floating shift
        shifted_irf_fp = IMP.bff.DecayCurve.shift_array(
            input=irf,
            shift=1.5
        )
        ref = np.array(
            [0.00000000e+00, 8.73523031e-06, 7.49069834e-05, 5.03256460e-04,
             2.65226555e-03, 1.09800745e-02, 3.57596335e-02, 9.17542811e-02,
             1.85744160e-01, 2.97018026e-01, 3.75503804e-01, 3.75503804e-01,
             2.97018026e-01, 1.85744160e-01, 9.17542811e-02, 3.57596335e-02,
             1.09800745e-02, 2.65226555e-03, 5.03256460e-04, 7.49069834e-05,
             8.73523031e-06, 7.97208558e-07, 5.68867416e-08, 6.07588285e-09,
             5.68867416e-08]
        )
        np.testing.assert_allclose(ref, shifted_irf_fp)

        shifted_irf_fn = IMP.bff.DecayCurve.shift_array(
            input=irf,
            shift=-1.5
        )
        ref = np.array(
            [6.07588285e-09, 5.68867416e-08, 7.97208558e-07, 8.73523031e-06,
             7.49069834e-05, 5.03256460e-04, 2.65226555e-03, 1.09800745e-02,
             3.57596335e-02, 9.17542811e-02, 1.85744160e-01, 2.97018026e-01,
             3.75503804e-01, 3.75503804e-01, 2.97018026e-01, 1.85744160e-01,
             9.17542811e-02, 3.57596335e-02, 1.09800745e-02, 2.65226555e-03,
             5.03256460e-04, 7.49069834e-05, 8.73523031e-06, 7.97208558e-07,
             5.68867416e-08]
        )
        np.testing.assert_allclose(ref, shifted_irf_fn)

        # rollover
        shifted_irf_irp = IMP.bff.DecayCurve.shift_array(
            input=irf,
            shift=26.0
        )
        np.testing.assert_allclose(shifted_irf_irp, shifted_irf_ip)

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
