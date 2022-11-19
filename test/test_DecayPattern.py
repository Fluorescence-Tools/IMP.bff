from __future__ import division
import unittest

import numpy as np
import numpy.testing
import IMP.bff


class Tests(unittest.TestCase):

    def test_constant_background(self):
        x = np.arange(0, 4, dtype=np.float64)
        y = np.zeros_like(x)
        dc = IMP.bff.DecayCurve(x, y)
        np.testing.assert_allclose(dc.y, y)
        bg = IMP.bff.DecayPattern(5.0)
        bg.add(dc)
        np.testing.assert_allclose(dc.y, y + 5)
        bg.add(dc)
        np.testing.assert_allclose(dc.y, y + 10)
        bg.constant_offset = -100
        self.assertEqual(bg.constant_offset, -100)
        bg(dc)
        np.testing.assert_allclose(dc.y, y - 90)

    def test_pattern(self):
        x = np.linspace(0, 10, 8, dtype=np.float64)
        y = np.sin(x) + 1.0
        dc = IMP.bff.DecayCurve(x, y)
        np.testing.assert_allclose(dc.y, y)

        bg_pattern = x
        bg_pattern_curve = IMP.bff.DecayCurve(x, bg_pattern)
        f = 0.8
        offset = 1.0
        bg = IMP.bff.DecayPattern(
            constant_offset=offset,
            pattern=bg_pattern_curve,
            pattern_fraction=f
        )
        bg.add(dc)
        r = y * (1-f) + bg_pattern * f * y.sum() / np.sum(bg_pattern) + offset
        np.testing.assert_allclose(r, dc.y)
