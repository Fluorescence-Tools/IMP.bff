from __future__ import division
import unittest

import numpy as np
import pandas as pd

import IMP.bff


class Tests(unittest.TestCase):

    def test_lin(self):
        dt = 0.0141
        df = pd.read_csv(
            IMP.bff.get_example_path("spectroscopy/hgbp1/eTCSPC_whitelight.txt"),
            skiprows=6,
            sep='\t'
        )
        lin_data = IMP.bff.DecayCurve(x=df['Chan'] * dt, y=df['Data'])
        model = IMP.bff.DecayCurve(
            x=df['Chan'] * dt,
            y=np.ones_like(lin_data.y)
        )
        dl = IMP.bff.DecayLinearization(
            linearization_table=lin_data,
            start=200, stop=3700, active=True, n_window=20
        )
        lin = dl.get_linearization_table()
        dl.add(model)

        ref_lin = [1., 0.99886486, 1.00740838, 0.99532836, 0.99401942, 1.0006662,
                   0.99965218, 0.99696761, 1.00444274, 0.99765081, 0.99634392, 1.00284555,
                   1.00678725, 1.00884042, 0.9865171, 1.]
        ref_mdl = [1., 0.99886486, 1.00740838, 0.99532836, 0.99401942, 1.0006662,
                   0.99965218, 0.99696761, 1.00444274, 0.99765081, 0.99634392, 1.00284555,
                   1.00678725, 1.00884042, 0.9865171, 1.]
        np.testing.assert_array_almost_equal(ref_lin, lin.y[::256])
        np.testing.assert_array_almost_equal(ref_mdl, model.y[::256])
