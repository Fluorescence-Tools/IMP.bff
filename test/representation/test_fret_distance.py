"""IMP.bff distances / fret.engine / greedy Olga — algorithm sanity checks."""

import numpy as np
import pytest

import IMP.bff

from IMP.bff import ProbeAccessibleVolume
import IMP.bff as fdist


def _point_av(xyz):
    pts = np.array([[xyz[0], xyz[1], xyz[2], 1.0]] * 4, dtype=np.float64)
    return ProbeAccessibleVolume(
        points=pts,
        density=np.zeros((1, 1, 1)),
        grid_origin=np.zeros(3),
        grid_step=1.0,
        attachment_point=np.asarray(xyz, dtype=np.float64),
    )


def test_point_cloud_distances_match_geometry():
    av1 = _point_av([0.0, 0.0, 0.0])
    av2 = _point_av([50.0, 0.0, 0.0])
    assert fdist.distance_between_mean_positions(av1, av2) == pytest.approx(50.0)
    assert fdist.states_average_distance(av1, av2, 200) == pytest.approx(50.0)
    r_e = fdist.states_mean_fret_distance(av1, av2, 52.0, 200)
    assert r_e == pytest.approx(50.0, abs=1e-6)
    assert fdist.standard_deviation_of_distances(
        av1, av2, 200) == pytest.approx(0.0, abs=1e-9)


def test_model_distance_types_and_chi2():
    av1 = _point_av([0.0, 0.0, 0.0])
    av2 = _point_av([40.0, 0.0, 0.0])
    for dtype in ("Rmp", "RDAMean", "RDAMeanE"):
        assert fdist.model_distance(
            av1, av2, dtype, n_samples=100) == pytest.approx(40.0, abs=1e-6)
    with pytest.raises(ValueError):
        fdist.model_distance(av1, av2, "nope")
    assert fdist.chi2_score(45.0, 40.0, 2.0, 5.0) == pytest.approx(1.0)
    assert fdist.chi2_score(35.0, 40.0, 5.0, 2.0) == pytest.approx(1.0)


def test_fret_efficiency_roundtrip():
    e = fdist.fret_efficiency(52.0, forster_radius=52.0)
    assert e == pytest.approx(0.5)
    assert fdist.distance_from_fret_efficiency(
        e, forster_radius=52.0) == pytest.approx(52.0)


def test_the_named_transfer_functions_dispatch():
    """The dispatch an fps.json calibration names, in one place.

    It used to be a method on a `DistanceRestraint` dataclass whose only other
    consumer was a rigid-body docking engine that had already been removed.
    """
    assert fdist.effective_distance(40.0, "None", 6.0) == pytest.approx(40.0)
    # sigma is the per-component width of the separation vector, so the
    # correction is s^2/Rmp -- settled 2026-08-18. This read 36/80, i.e. half.
    assert fdist.effective_distance(40.0, "Gaussian", 6.0) == pytest.approx(
        40.0 + 36.0 / 40.0)
    # Coefficients run lowest power first here, as a hand-written calibration
    # does; `polynomial_transfer` is the highest-power-first sibling.
    assert fdist.effective_distance(
        40.0, "Polynomial", 6.0, [0.0, 1.0, 0.02]) == pytest.approx(
            40.0 + 0.02 * 1600.0)
    # A calibration that names a polynomial and carries none has not been
    # fitted; sigma is the parametric stand-in.
    assert fdist.effective_distance(40.0, "Polynomial", 6.0) == pytest.approx(
        40.0 + 36.0 / 40.0)
    # A transfer function nobody implements is no transfer function.
    assert fdist.effective_distance(40.0, "Spline", 6.0) == pytest.approx(40.0)


def test_select_probe_pairs_prefers_discriminating_pair():
    rng = np.random.RandomState(0)
    n_frames = 20
    # pair 0 separates two conformation clusters; pair 1 is pure noise
    effs = np.empty((n_frames, 2), dtype=np.float64)
    labels = np.arange(n_frames) % 2
    effs[:, 0] = 0.2 + 0.6 * labels + rng.normal(0, 0.01, n_frames)
    effs[:, 1] = 0.5 + rng.normal(0, 0.01, n_frames)
    rmsds = np.abs(labels[:, None] - labels[None, :]) * 10.0
    selected, decay = IMP.bff.select_probe_pairs(
        effs, rmsds, measurement_error=0.05, max_pairs=2)
    assert selected[0] == 0
    assert decay.shape == (2,)
    assert decay[0] <= rmsds.mean() + 1e-6
