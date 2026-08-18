"""IMP.bff.representation.distance / fret.engine / olga_greedy — algorithm sanity checks."""

import numpy as np
import pytest

from IMP.bff.fret.av import AccessibleVolume
from IMP.bff.representation import distance as fdist
from IMP.bff.fret.engine import DistanceRestraint, RigidBody
from IMP.bff.fret.olga_greedy import select_informative_pairs


def _point_av(xyz):
    pts = np.array([[xyz[0], xyz[1], xyz[2], 1.0]] * 4, dtype=np.float64)
    return AccessibleVolume(
        points=pts,
        density=np.zeros((1, 1, 1), dtype=np.float32),
        grid_origin=np.zeros(3),
        grid_step=1.0,
        grid_shape=(1, 1, 1),
        attachment_point=np.asarray(xyz, dtype=np.float64),
    )


def test_point_cloud_distances_match_geometry():
    av1 = _point_av([0.0, 0.0, 0.0])
    av2 = _point_av([50.0, 0.0, 0.0])
    assert fdist.distance_between_mean_positions(av1, av2) == pytest.approx(50.0)
    assert fdist.average_distance(av1, av2, n_samples=200) == pytest.approx(50.0)
    r_e = fdist.mean_fret_distance(av1, av2, forster_radius=52.0,
                                   n_samples=200)
    assert r_e == pytest.approx(50.0, abs=1e-6)
    assert fdist.standard_deviation_of_distances(
        av1, av2, n_samples=200) == pytest.approx(0.0, abs=1e-9)


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


def test_rigid_body_global_coords_and_restraint():
    body = RigidBody(
        name="b",
        atoms_local=np.array([[1.0, 0.0, 0.0, 1.7]]),
        com=np.array([10.0, 0.0, 0.0]),
        rotation=np.eye(3),
        translation=np.array([10.0, 0.0, 0.0]),
    )
    assert body.global_coords()[0] == pytest.approx([11.0, 0.0, 0.0])
    r = DistanceRestraint(
        name="d", body_a=0, offset_a=np.zeros(3), body_b=0,
        offset_b=np.zeros(3), distance_exp=40.0, error_neg=2.0,
        error_pos=2.0, transfer_function_type="None")
    assert r.get_effective_distance(40.0) == pytest.approx(40.0)
    r.transfer_function_type = "Gaussian"
    r.sigma_rda = 6.0
    # sigma is the per-component width of the separation vector, so the
    # correction is s^2/Rmp -- settled 2026-08-18. This read 36/80, i.e. half.
    assert r.get_effective_distance(40.0) == pytest.approx(
        40.0 + 36.0 / 40.0)


def test_select_informative_pairs_prefers_discriminating_pair():
    rng = np.random.RandomState(0)
    n_frames = 20
    # pair 0 separates two conformation clusters; pair 1 is pure noise
    effs = np.empty((n_frames, 2), dtype=np.float64)
    labels = np.arange(n_frames) % 2
    effs[:, 0] = 0.2 + 0.6 * labels + rng.normal(0, 0.01, n_frames)
    effs[:, 1] = 0.5 + rng.normal(0, 0.01, n_frames)
    rmsds = np.abs(labels[:, None] - labels[None, :]) * 10.0
    selected, decay = select_informative_pairs(
        effs, rmsds, err=0.05, max_pairs=2)
    assert selected[0] == 0
    assert decay.shape == (2,)
    assert decay[0] <= rmsds.mean() + 1e-6
