"""The fused walk->rate->photons path: identical to the three calls, without the trajectory.

A PET-quenching prediction is three steps, and the middle one used to cross the
Python boundary twice carrying the largest object in the chain. At the default
``t_max`` the trajectory is 20 million doubles, handed out so that a rate map
can be read along it and the trace handed back, to produce a few thousand
photons. Nothing downstream of the decay wants the coordinates.

:meth:`QuenchedDonorDecay.photons_fused` does the same computation with the
intermediates left where they are made. The claim these tests defend is that it
is the *same computation*, not an approximation: given the same seeds it must
reproduce the three-call path photon for photon. That is achievable because both
paths run the same walk core (``internal/RandomWalk.h``) with the same generator
and the same seeds -- so this is a real equality assertion, not a tolerance.

The three-call path stays, and should. A trajectory is what a correlation
function or a visualisation needs, and the rate trace is worth inspecting on its
own.
"""

import time

import numpy as np
import pytest

from IMP.bff.quenching.model import QuenchedDonorDecay
from IMP.bff.quenching.pet import amino_acid_quenching_defaults

ATOM_DTYPE = [
    ("chain", "U4"), ("res_id", "i8"), ("res_name", "U4"),
    ("atom_name", "U4"), ("coord", "f8", 3),
]


class _FakeAV:
    def __init__(self, density, grid_step, attachment_point):
        self.density = density
        self.grid_step = grid_step
        self.attachment_point = attachment_point


def _sphere(ng=41, radius=15):
    i = np.arange(ng) - (ng - 1) // 2
    x, y, z = np.meshgrid(i, i, i, indexing="ij")
    return ((x ** 2 + y ** 2 + z ** 2) < radius ** 2).astype(np.uint8)


def _atoms():
    rows = [
        ("A", 1, "TRP", "CB", [5.0, 0.0, 0.0]),
        ("A", 1, "TRP", "NE1", [7.0, 0.0, 0.0]),
        ("A", 1, "TRP", "CD2", [7.0, 1.0, 0.0]),
        ("A", 2, "ALA", "CB", [-6.0, 0.0, 0.0]),
    ]
    a = np.zeros(len(rows), dtype=ATOM_DTYPE)
    for i, r in enumerate(rows):
        a[i] = r
    return a


def _model(**kw):
    av = _FakeAV(_sphere(), 1.0, np.zeros(3))
    options = dict(
        tau0=4.0, quenching_table=amino_acid_quenching_defaults(),
        critical_distance=7.0, slow_radius=10.0,
        t_max=400.0, t_step=0.004, n_photons=40000,
        n_trajectories=1, random_seed=7,
    )
    options.update(kw)
    return QuenchedDonorDecay(av, _atoms(), **options)


def test_fused_reproduces_the_three_call_path_photon_for_photon():
    """The claim. Not a tolerance -- the same walk, the same draws."""
    slow = _model()
    slow.simulate_diffusion()
    delays_a, emitted_a = slow.simulate_photons()

    fast = _model()
    delays_b, emitted_b = fast.photons_fused()

    np.testing.assert_array_equal(emitted_a, emitted_b)
    np.testing.assert_array_equal(delays_a, delays_b)


def test_fused_reports_the_statistics_the_trajectory_would_have_given():
    """It skips the coordinates, so it has to hand back what they were read for."""
    slow = _model()
    slow.simulate_diffusion()
    walk = slow.diffusion

    fast = _model()
    fast.photons_fused()

    assert fast._fused_stats["n_frames"] == walk.n_frames
    # Exactly, because the fused kernel holds the rate trace as float32 too.
    # It did not at first, and the mean came out 2.6e-8 off -- harmless in
    # itself, but it meant the two paths could in principle disagree about a
    # photon. Matching `sample_grid`'s float32 removed that and halved the
    # memory the photon race walks, which is the bottleneck at long
    # trajectories: two reasons pointing the same way.
    assert fast._fused_stats["mean_k_quench"] == pytest.approx(
        float(walk.k_quench.mean()), rel=1e-12)
    assert fast._fused_stats["collision_fraction"] == pytest.approx(
        walk.collision_fraction, rel=1e-12)
    assert fast.diffusion.n_accepted == walk.n_accepted
    assert fast.diffusion.n_rejected == walk.n_rejected


def test_fused_agrees_across_several_trajectories():
    """Concatenation, not averaging: each photon draws its own start frame."""
    slow = _model(n_trajectories=3)
    slow.simulate_diffusion()
    delays_a, emitted_a = slow.simulate_photons()

    fast = _model(n_trajectories=3)
    delays_b, emitted_b = fast.photons_fused()

    np.testing.assert_array_equal(emitted_a, emitted_b)
    np.testing.assert_array_equal(delays_a, delays_b)


def test_fused_agrees_with_a_per_voxel_mobility_field():
    """The stickiness field is the branch most likely to diverge silently.

    ``_slow_factor_map`` is set directly rather than through ``update_grids``,
    because the point here is the mobility *branch* in each path, not how the
    field was built.
    """
    slow_map = np.where(_sphere(41, 8) > 0, 0.05, 1.0)

    slow = _model()
    slow.update_grids()
    slow._slow_factor_map = slow_map
    slow.simulate_diffusion()
    delays_a, emitted_a = slow.simulate_photons()

    fast = _model()
    fast.update_grids()
    fast._slow_factor_map = slow_map
    delays_b, emitted_b = fast.photons_fused()

    np.testing.assert_array_equal(emitted_a, emitted_b)
    np.testing.assert_array_equal(delays_a, delays_b)


def test_the_quantum_yield_and_lifetime_come_out_the_same():
    slow = _model()
    slow.simulate_diffusion()
    slow.simulate_photons()

    fast = _model()
    fast.photons_fused()

    assert fast.quantum_yield == pytest.approx(slow.quantum_yield, rel=1e-12)
    assert fast.fluorescence_lifetime == pytest.approx(
        slow.fluorescence_lifetime, rel=1e-12)


def test_an_empty_volume_yields_no_photons_either_way():
    empty = np.zeros((21, 21, 21), dtype=np.uint8)
    m = _model()
    m.diffusion.density = empty
    m._quenching_rate_map = np.zeros((21, 21, 21))
    delays, emitted = m.photons_fused()
    assert delays.size == 0 and emitted.size == 0


def test_fused_is_faster():
    """A timing assertion can rot on a loaded machine, so this one is loose.

    Measured best-of-3 with the grids already stamped, so the comparison is of
    the walk-and-race and nothing else:

        1e6 steps    1891 ms -> 813 ms    2.33x     32 MB not crossed
        5e6 steps    8419 ms -> 1582 ms   5.32x    160 MB
        2e7 steps   43065 ms -> 4102 ms  10.50x    640 MB

    The ratio grows with trajectory length because what is removed *is* the
    trajectory. The test asserts only a modest factor, at a size small enough to
    stay quick.

    Two mistakes cost an hour here and are worth remembering. Timing a fresh
    model per repetition put ``update_grids`` inside the measurement, which
    dominates at small sizes and made fusing look 4x *slower*. And
    ``photons_fused`` first read ``self.diffusion``, a property that runs the
    split walk on access -- so the fused path was doing both walks.
    """
    kw = dict(t_max=4000.0, t_step=0.002, n_photons=20000)

    a = _model(**kw)
    a.update_grids()
    t0 = time.perf_counter()
    a.simulate_diffusion()
    a.simulate_photons()
    t_split = time.perf_counter() - t0

    b = _model(**kw)
    b.update_grids()
    t0 = time.perf_counter()
    b.photons_fused()
    t_fused = time.perf_counter() - t0

    np.testing.assert_array_equal(a.photon_trace[0], b.photon_trace[0])
    assert t_fused < t_split, f"fused {t_fused:.3f}s vs split {t_split:.3f}s"


def test_fused_does_not_run_the_walk_twice():
    """``self.diffusion`` runs the split walk on first access. Reading it inside
    the fused path made it do both, and only a benchmark caught that."""
    m = _model()
    assert m._diffusion is None
    m.photons_fused()
    # a simulation object exists, but it was never asked to produce a trajectory
    assert m._diffusion is not None
    assert m._diffusion.trajectory is None


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
