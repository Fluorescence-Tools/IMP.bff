"""Numerical contract frozen before the PRD-141 public API rename."""

import numpy as np
import pytest

import IMP.bff


@pytest.mark.parametrize("measurement_kind", ["FRET efficiency", "EPR distance", "PRE rate"])
def test_predicted_measurements_preserve_selection_and_precision(measurement_kind):
    # The selector receives predictions, not a modality-specific forward model.
    # Identical normalized predictions/errors must give identical decisions.
    predicted_measurements = np.array([
        [.1, .4, .3, .2], [.5, .3, .7, .1],
        [.8, .9, .2, .6], [.2, .1, .6, .9],
    ])
    rmsds = np.array([
        [0., 2, 6, 3], [2, 0, 5, 4], [6, 5, 0, 7], [3, 4, 7, 0],
    ])
    pair_sites = np.array([[0, 0], [1, 1], [0, 1], [1, 2]], dtype=np.int32)
    pairs, pair_precision = IMP.bff.select_probe_pairs(
        predicted_measurements=predicted_measurements, rmsds=rmsds,
        measurement_error=.2, max_pairs=4)
    positions, position_precision = IMP.bff.select_probe_positions(
        predicted_measurements=predicted_measurements, rmsds=rmsds,
        pair_sites=pair_sites, measurement_error=.2, max_sites=3)
    np.testing.assert_array_equal(pairs, [1, 0, 3, 2], err_msg=measurement_kind)
    np.testing.assert_array_equal(positions, [1, 0, 2], err_msg=measurement_kind)
    np.testing.assert_allclose(pair_precision, [
        .8821218178384835, .6803200132679036,
        .1810654733262066, .07491362762366367,
    ], rtol=1e-13, atol=1e-15, err_msg=measurement_kind)
    np.testing.assert_allclose(position_precision, [
        .8821218178384835, .6974618288556984, .07491362762366367,
    ], rtol=1e-13, atol=1e-15, err_msg=measurement_kind)
