import unittest
import numpy as np
import numpy.testing
import IMP.bff
import IMP.test

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


class Tests(IMP.test.TestCase):

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
        # Only test AVX functionality if it was enabled
        if IMP.bff.IMP_BFF_HAS_AVX:
            conv_methods_fast = [
                IMP.bff.DecayConvolution.FAST_AVX,
                IMP.bff.DecayConvolution.FAST_PERIODIC_AVX
            ]
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


if __name__ == '__main__':
    IMP.test.main()
