from __future__ import division
import unittest
import platform

import numpy as np
import numpy.testing
import IMP.bff

x = np.linspace(0, 20, 32)
irf_position = 2.0
irf_width = 0.1


def norm_pdf(x, mu, sigma):
    variance = sigma**2
    num = x - mu
    denom = 2*variance
    pdf = ((1/(np.sqrt(2*np.pi)*sigma))*np.exp(-(num**2)/denom))
    return pdf



irf_y = norm_pdf(x, irf_position, irf_width)


class Tests(unittest.TestCase):

    def test_DecayConvolution_init(self):
        irf = IMP.bff.DecayCurve(x, irf_y)
        lh = IMP.bff.DecayLifetimeHandler([1, 4])
        settings = {
            "lifetime_handler": lh,
            "instrument_response_function": irf,
            "convolution_method": IMP.bff.DecayConvolution.FAST,
            "excitation_period": 100,
            "irf_shift_channels": 0.0,
            "irf_background_counts": 0.0,
            "start": 0, "stop": -1,
            "active": True
        }
        dc = IMP.bff.DecayConvolution(**settings)
        decay = IMP.bff.DecayCurve(x)
        dc.add(decay)
        ref = np.array([1.19203776e-10, 1.59816046e-05, 2.61648339e-02, 4.63676222e-01,
                        8.17250308e-01, 7.51713472e-01, 6.39828454e-01, 5.44522781e-01,
                        4.63413367e-01, 3.94385610e-01, 3.35639884e-01, 2.85644631e-01,
                        2.43096423e-01, 2.06885985e-01, 1.76069274e-01, 1.49842868e-01,
                        1.27523017e-01, 1.08527821e-01, 9.23620548e-02, 7.86042610e-02,
                        6.68957600e-02, 5.69312992e-02, 4.84510951e-02, 4.12340602e-02,
                        3.50920390e-02, 2.98649029e-02, 2.54163750e-02, 2.16304778e-02,
                        1.84085090e-02, 1.56664688e-02, 1.33328693e-02, 1.13468712e-02])
        np.testing.assert_allclose(decay.y, ref)
        ref = np.array([ 0.        ,  0.64516129,  1.29032258,  1.93548387,  2.58064516,
                         3.22580645,  3.87096774,  4.51612903,  5.16129032,  5.80645161,
                         6.4516129 ,  7.09677419,  7.74193548,  8.38709677,  9.03225806,
                         9.67741935, 10.32258065, 10.96774194, 11.61290323, 12.25806452,
                         12.90322581, 13.5483871 , 14.19354839, 14.83870968, 15.48387097,
                         16.12903226, 16.77419355, 17.41935484, 18.06451613, 18.70967742,
                         19.35483871, 20.        ])
        np.testing.assert_allclose(decay.x, ref)

    def test_DecayConvolution_init(self):
        irf = IMP.bff.DecayCurve(x, irf_y)
        lh = IMP.bff.DecayLifetimeHandler([1, 4])
        settings = {
            "instrument_response_function": irf,
            "lifetime_handler": lh,
            "convolution_method": 0,
            "excitation_period": 100,
            "irf_shift_channels": 0,
            "irf_background_counts": 0,
            "start": 0,
            "stop": -1
        }
        ref = [
            [4.70950187e-11, 9.33162000e-01, 2.56789386e-01, 7.06638168e-02],
            [1.78095108e-87, 9.33162000e-01, 2.56789386e-01, 7.06638168e-02],
            [4.49386205e-11, 9.27072464e-01, 2.65610708e-01, 7.60987419e-02],
            [1.72529636e-87, 9.27072464e-01, 2.65610708e-01, 7.60987419e-02],
            [1.72529636e-87, 9.27072464e-01, 2.65610708e-01, 7.60987419e-02],
            [4.49386205e-11, 9.27072464e-01, 2.65610708e-01, 7.60987419e-02]
        ]
        conv_methods = [
            IMP.bff.DecayConvolution.FAST_PERIODIC_TIME,
            IMP.bff.DecayConvolution.FAST_TIME,
            IMP.bff.DecayConvolution.FAST_PERIODIC,
            IMP.bff.DecayConvolution.FAST
        ]
        conv_methods_fast = [
            IMP.bff.DecayConvolution.FAST_AVX,
            IMP.bff.DecayConvolution.FAST_PERIODIC_AVX
        ]
        if "AMD64" in platform.machine():
            if platform.system() == "Linux":
                conv_methods += conv_methods_fast
            elif platform.system() == "Windows":
                conv_methods += conv_methods_fast
        for i in conv_methods:
            settings["convolution_method"] = i
            dc = IMP.bff.DecayConvolution(**settings)
            decay = IMP.bff.DecayCurve(x)
            dc.add(decay)
            np.testing.assert_allclose(decay.y[::8], ref[i])

    def test_irf(self):
        irf = IMP.bff.DecayCurve(x=x, y=irf_y)
        dc = IMP.bff.DecayConvolution()
        dc.irf = irf
        dc.irf_background_counts = 1e-3
        np.testing.assert_allclose(dc.irf.y, irf_y)
        self.assertEqual(dc.irf_background_counts, 1e-3)

        np.allclose(dc.corrected_irf.y, np.clip(dc.irf.y - dc.irf_background_counts, 0, 1e30))
        dc.irf_background_counts = 0.0
        dc.irf_shift_channels = 0.0
        np.testing.assert_allclose(
            dc.corrected_irf.y,
            np.array([5.52094836e-087, 5.51600473e-040, 4.61808726e-011, 3.23985966e+000,
                      1.90465720e-007, 9.38283342e-033, 3.87326980e-076, 1.33982492e-137,
                      3.88369320e-217, 0.00000000e+000, 0.00000000e+000, 0.00000000e+000,
                      0.00000000e+000, 0.00000000e+000, 0.00000000e+000, 0.00000000e+000,
                      0.00000000e+000, 0.00000000e+000, 0.00000000e+000, 0.00000000e+000,
                      0.00000000e+000, 0.00000000e+000, 0.00000000e+000, 0.00000000e+000,
                      0.00000000e+000, 0.00000000e+000, 0.00000000e+000, 0.00000000e+000,
                      0.00000000e+000, 0.00000000e+000, 0.00000000e+000, 0.00000000e+000]),
            atol=5e-9
        )
        dc.irf_shift_channels = 10.0
        self.assertEqual(dc.irf_shift_channels, 10.0)
        np.testing.assert_allclose(
            dc.corrected_irf.y,
            np.array([0.00000000e+000, 0.00000000e+000, 0.00000000e+000, 0.00000000e+000,
                      0.00000000e+000, 0.00000000e+000, 0.00000000e+000, 0.00000000e+000,
                      0.00000000e+000, 0.00000000e+000, 0.00000000e+000, 0.00000000e+000,
                      0.00000000e+000, 0.00000000e+000, 0.00000000e+000, 0.00000000e+000,
                      0.00000000e+000, 0.00000000e+000, 0.00000000e+000, 0.00000000e+000,
                      0.00000000e+000, 0.00000000e+000, 5.52094836e-087, 5.51600473e-040,
                      4.61808726e-011, 3.23985966e+000, 1.90465720e-007, 9.38283342e-033,
                      3.87326980e-076, 1.33982492e-137, 3.88369320e-217, 0.00000000e+000]),
            atol=5e-9
        )

    def test_getter_setter(self):
        irf = IMP.bff.DecayCurve(x, irf_y)
        lh = IMP.bff.DecayLifetimeHandler()
        ls = np.array([1., 10.0])

        dc = IMP.bff.DecayConvolution()
        dc.irf = irf

        self.assertEqual(dc.convolution_method, IMP.bff.DecayConvolution.FAST)
        dc.convolution_method = IMP.bff.DecayConvolution.FAST_PERIODIC
        self.assertEqual(dc.convolution_method, IMP.bff.DecayConvolution.FAST_PERIODIC)

        self.assertEqual(dc.excitation_period, 100.0)
        dc.excitation_period = 32.0
        self.assertEqual(dc.excitation_period, 32.0)

        self.assertEqual(dc.start, 0)
        dc.start = 5
        self.assertEqual(dc.start, 5)

        self.assertEqual(dc.stop, -1)
        dc.stop = 20
        self.assertEqual(dc.stop, 20)

        dc.lifetime_handler = lh
        dc.lifetime_handler.lifetime_spectrum = ls
        np.allclose(dc.lifetime_handler.lifetime_spectrum, ls)

    def test_mean_lifetime(self):
        irf = IMP.bff.DecayCurve(x, irf_y)
        lh = IMP.bff.DecayLifetimeHandler([1, 4])
        settings = {
            "instrument_response_function": irf,
            "lifetime_handler": lh,
            "convolution_method": 0,
            "excitation_period": 100,
            "irf_shift_channels": 0,
            "irf_background_counts": 0,
            "start": 0,
            "stop": -1
        }
        dc = IMP.bff.DecayConvolution(**settings)
        y = np.ones_like(x)
        decay = IMP.bff.DecayCurve(x, y)
        self.assertAlmostEqual(dc.get_mean_lifetime(decay), 7.812499963266247)

