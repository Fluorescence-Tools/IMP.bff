from __future__ import division
import unittest

import numpy as np
import IMP.bff


class Tests(unittest.TestCase):

    def test_DecayLifetimeHandler_init(self):
        lt = IMP.bff.DecayLifetimeHandler()
        self.assertEqual(len(lt.lifetime_spectrum), 0)

        settings = {
            'lifetime_spectrum': [1, 4, -0.01, 2],
            'use_amplitude_threshold': False,
            'abs_lifetime_spectrum': False,
            'amplitude_threshold': 0.1
        }
        lt = IMP.bff.DecayLifetimeHandler(**settings)
        self.assertEqual(
            np.allclose(lt.lifetime_spectrum, np.array([1, 4, -0.01, 2])), True
        )

        settings = {
            'lifetime_spectrum': [1, 4, -0.01, 2],
            'use_amplitude_threshold': False,
            'abs_lifetime_spectrum': True,
            'amplitude_threshold': 0.1
        }
        lt = IMP.bff.DecayLifetimeHandler(**settings)
        self.assertEqual(
            np.allclose(lt.lifetime_spectrum, np.array([1, 4, 0.01, 2])), True
        )

        settings = {
            'lifetime_spectrum': [1, 4, -0.01, 2],
            'use_amplitude_threshold': True,
            'abs_lifetime_spectrum': True,
            'amplitude_threshold': 0.1
        }
        lt = IMP.bff.DecayLifetimeHandler(**settings)
        self.assertEqual(
            np.allclose(lt.lifetime_spectrum, np.array([1, 4, 0.0, 2])), True
        )

        lt = IMP.bff.DecayLifetimeHandler([1., 5])
        self.assertEqual(
            np.allclose(lt.lifetime_spectrum, np.array([1, 5])), True
        )

    def test_DecayLifetimeHandler_setter_getter(self):
        import IMP.bff
        import numpy as np
        lt = IMP.bff.DecayLifetimeHandler()

        self.assertEqual(lt.use_amplitude_threshold, False)
        lt.use_amplitude_threshold = True
        self.assertEqual(lt.use_amplitude_threshold, True)

        self.assertAlmostEqual(lt.amplitude_threshold, 2.220446049250313e-16)
        lt.amplitude_threshold = 0.111
        self.assertEqual(lt.amplitude_threshold, 0.111)

        l = np.array([1., 4])
        lt.lifetime_spectrum = l
        self.assertEqual(np.allclose(lt.lifetime_spectrum, l), True)
