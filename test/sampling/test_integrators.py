"""The particle model in C++: a confined Brownian walk, photons, and the curve.

These three kernels carry a random number generator, so none of them can be
gated bit-for-bit against the numba they replace -- no C++ generator reproduces
numba's stream. They are gated against **analytic results** instead, which is
the stronger check anyway:

* a free walk's mean squared displacement,
* a confined walk's equilibrium distribution -- uniform on the accessible
  region, which is what makes rejection sampling correct,
* the emitted fraction and mean delay time under a constant quenching rate,
  both of which have closed forms,
* the decay curve's rate, which must be ``1/tau0 + kq``.

The walk and the field solver are the same dynamics -- a Langevin trajectory
and its Fokker-Planck density -- so they must agree on what ``D`` means. They
did not until 2026-08-18: the step width was ``sqrt(2*D*3*dt)`` per Cartesian
component, the total three-dimensional MSD used as one component's width, so
the walk diffused at ``3D``. Traced to a 2019 QuEst docstring that transcribes a
Berkeley teaching page's step *magnitude* as a per-component sigma; fixed, and
``D`` left alone, since ``D = 40 A^2/ns`` is free Alexa488 entered correctly
(1 A^2/ns = 10 um^2/s) and no calibration of ``D`` exists anywhere in the stack.
Full record: ``okf/validation/particle_vs_field_diffusion.md``.

The three tests below measure both models and assert they agree, so this cannot
drift apart again.
"""

import numpy as np
import pytest

import IMP.bff
import IMP.bff as dif
import IMP.bff as ph
from IMP.bff import GridDiffusionSolver, diffusion_stability_limit


def _ball(ng, radius):
    i = np.arange(ng)
    x, y, z = np.meshgrid(i, i, i, indexing="ij")
    c = (ng - 1) // 2
    return (((x - c) ** 2 + (y - c) ** 2 + (z - c) ** 2) < radius ** 2).astype(np.uint8)


# --- the walk ----------------------------------------------------------------

def test_free_walk_step_variance_is_2Ddt_per_component():
    """The exact form, measured per step so no wall can touch it.

    A 45-voxel box, not 151. The step width is 0.28 voxels, so even 4000 steps
    stay within a few voxels of the start and the walls are never reached -- and
    the occupancy grid is marshalled once per call, which at 151 cubed is 3.4
    million elements and 112 ms of pure conversion per walk.

    A displacement over many steps is noisy and boundary-sensitive; the
    single-step variance is neither, and it is what the convention actually
    fixes. The 3-D total is then 6 D dt, which is the number the old code put
    into one component.
    """
    D, t_step = 8.0, 0.005
    box = np.ones((45, 45, 45), dtype=np.uint8)
    steps = [np.diff(t.xyz, axis=0) for t in
             (dif.simulate_dye_diffusion(box, dg=1.0, t_max=4000 * t_step,
                                         t_step=t_step, D=D, random_seed=s)
              for s in range(40))
             if t.acceptance_ratio > 0.9999]
    st = np.concatenate(steps)
    assert st.var(axis=0).mean() == pytest.approx(2 * D * t_step, rel=0.02)
    assert (st ** 2).sum(axis=1).mean() == pytest.approx(6 * D * t_step, rel=0.02)


def test_free_walk_diffuses_at_the_stated_D():
    """<dx^2> = 2 D t over many steps, which is the observable claim."""
    D, t_step, n = 8.0, 0.005, 60
    box = np.ones((45, 45, 45), dtype=np.uint8)
    d = []
    for seed in range(300):
        t = dif.simulate_dye_diffusion(box, dg=1.0, t_max=n * t_step,
                                       t_step=t_step, D=D, random_seed=seed)
        if t.acceptance_ratio > 0.999:      # untouched by the walls
            d.append(t.xyz[n - 1] - t.xyz[0])
    per_component = (np.array(d) ** 2).mean(axis=0).mean()
    assert per_component == pytest.approx(2 * D * n * t_step, rel=0.15)


def test_grid_solver_uses_the_standard_convention():
    """The other half of the mismatch: the field model gives 2 D t, exactly."""
    ng, dg, D = 81, 1.0, 8.0
    bounds = np.zeros((ng, ng, ng)); bounds[2:-2, 2:-2, 2:-2] = 1.0
    density = np.zeros((ng, ng, ng)); density[40, 40, 40] = 1.0
    t_step = 0.4 * diffusion_stability_limit(D, dg)
    s = GridDiffusionSolver(diffusion_map=np.full((ng, ng, ng), D), bounds=bounds,
                            density=density, t_step=t_step, dg=dg)
    s.run(300, n_out=1)
    p = np.asarray(s.get_density()).reshape(ng, ng, ng)
    p = p / p.sum()
    axis = (np.arange(ng) - 40) * dg
    msd_x = (p.sum(axis=(1, 2)) * axis ** 2).sum()
    assert msd_x == pytest.approx(2 * D * 300 * t_step, rel=1e-3)


def test_the_two_models_agree_on_D():
    """The walk and the solver, measured side by side at the same D.

    These are one dynamics expressed two ways; a disagreement here is a
    disagreement about the model, not about numerics.
    """
    D, t_step, n = 8.0, 0.005, 60
    box = np.ones((45, 45, 45), dtype=np.uint8)
    d = [t.xyz[n - 1] - t.xyz[0] for t in
         (dif.simulate_dye_diffusion(box, dg=1.0, t_max=n * t_step, t_step=t_step,
                                     D=D, random_seed=s) for s in range(300))
         if t.acceptance_ratio > 0.999]
    walk_msd = (np.array(d) ** 2).mean(axis=0).mean()

    ng, dg = 81, 1.0
    bounds = np.zeros((ng, ng, ng)); bounds[2:-2, 2:-2, 2:-2] = 1.0
    density = np.zeros((ng, ng, ng)); density[40, 40, 40] = 1.0
    solver_step = 0.4 * diffusion_stability_limit(D, dg)
    solver = GridDiffusionSolver(diffusion_map=np.full((ng, ng, ng), D), bounds=bounds,
                                 density=density, t_step=solver_step, dg=dg)
    solver.run(300, n_out=1)
    p = np.asarray(solver.get_density()).reshape(ng, ng, ng)
    p = p / p.sum()
    axis = (np.arange(ng) - 40) * dg
    field_msd = (p.sum(axis=(1, 2)) * axis ** 2).sum()

    # Compare the diffusion coefficient each one implies, not the raw MSDs:
    # the two ran for different lengths of time.
    assert walk_msd / (2 * n * t_step) == pytest.approx(
        field_msd / (2 * 300 * solver_step), rel=0.15)


def test_walk_never_leaves_the_accessible_region():
    ng = 31
    ball = _ball(ng, 12)
    t = dif.simulate_dye_diffusion(ball, dg=1.0, t_max=2000.0, t_step=0.005,
                                   D=8.0, random_seed=1)
    idx = np.floor(t.xyz).astype(int) + (ng - 1) // 2
    inside = ball[np.clip(idx[:, 0], 0, ng - 1),
                  np.clip(idx[:, 1], 0, ng - 1),
                  np.clip(idx[:, 2], 0, ng - 1)]
    assert inside.all()


def test_confined_walk_samples_the_region_uniformly():
    """Uniform in a ball gives p(r) ~ r^2, so <r> = 3R/4.

    This is what rejection sampling buys, and it is why a rejected step still
    emits a frame: dropping them would bias the walk away from the boundary.
    """
    t = dif.simulate_dye_diffusion(_ball(31, 12), dg=1.0, t_max=4000.0,
                                   t_step=0.005, D=8.0, random_seed=1)
    r = np.linalg.norm(t.xyz, axis=1)
    assert r.mean() == pytest.approx(0.75 * 12.0, abs=1.5)


def test_walk_emits_one_frame_per_step_accepted_or_not():
    t = dif.simulate_dye_diffusion(_ball(21, 4), dg=1.0, t_max=100.0,
                                   t_step=0.01, D=40.0, random_seed=2)
    assert t.n_frames == 10000
    assert t.accepted.shape == (10000,)
    assert t.n_accepted + t.n_rejected == 10000
    assert 0.0 < t.acceptance_ratio < 1.0, "a small ball must reject some steps"


def test_walk_is_reproducible_and_seed_dependent():
    kw = dict(dg=1.0, t_max=500.0, t_step=0.005, D=8.0)
    ball = _ball(31, 12)
    a = dif.simulate_dye_diffusion(ball, random_seed=5, **kw)
    b = dif.simulate_dye_diffusion(ball, random_seed=5, **kw)
    np.testing.assert_array_equal(a.xyz, b.xyz)
    assert a.n_accepted == b.n_accepted
    c = dif.simulate_dye_diffusion(ball, random_seed=6, **kw)
    assert not np.array_equal(a.xyz, c.xyz)


def test_mobility_field_and_scalar_mask_agree():
    """The two Python kernels this replaced were the same walk twice over."""
    ng = 31
    ball = _ball(ng, 12)
    i = np.arange(ng); x, y, z = np.meshgrid(i, i, i, indexing="ij")
    outer = (((x - 15) ** 2 + (y - 15) ** 2 + (z - 15) ** 2) > 64)
    kw = dict(dg=1.0, t_max=2000.0, t_step=0.005, D=8.0, random_seed=2)

    field = dif.simulate_dye_diffusion(ball, slow_fact=np.where(outer, 0.1, 1.0), **kw)
    scalar = dif.simulate_dye_diffusion(ball, slow_density=outer.astype(np.uint8),
                                        slow_fact=0.1, **kw)
    free = dif.simulate_dye_diffusion(ball, **kw)

    step = lambda t: np.linalg.norm(np.diff(t.xyz, axis=0), axis=1).mean()
    assert step(field) == pytest.approx(step(scalar), rel=0.05)
    assert step(field) < step(free), "a slow region must shorten the mean step"


def test_no_accessible_voxel_gives_an_empty_trajectory():
    t = dif.simulate_dye_diffusion(np.zeros((11, 11, 11), dtype=np.uint8),
                                   dg=1.0, t_max=10.0, t_step=0.01, D=8.0, random_seed=1)
    assert t.n_accepted == 0 and t.n_rejected == 0
    assert t.acceptance_ratio == 0.0
    assert not t.xyz.any()


# --- photons -----------------------------------------------------------------

@pytest.mark.parametrize("kq_val", [0.0, 0.25, 0.5, 2.0])
def test_constant_rate_reproduces_the_closed_form(kq_val):
    """At a constant quenching rate both observables have exact answers."""
    tau0 = 4.0
    dts, emitted = ph.simulate_photon_trace(
        200000, np.full(2000, kq_val), 0.01, tau0, random_seed=11)
    assert emitted.mean() == pytest.approx(1.0 / (1.0 + kq_val * tau0), abs=0.006)
    assert dts[emitted > 0].mean() == pytest.approx(1.0 / (1.0 / tau0 + kq_val), rel=0.02)


def test_no_photon_is_emitted_before_it_was_excited():
    """The epsilon this code inherited made 2.4e-4 of delay times negative."""
    dts, _ = ph.simulate_photon_trace(500000, np.zeros(10), 0.01, 4.0, random_seed=3)
    assert (dts >= 0).all()


def test_quenched_events_report_zero_delay():
    dts, emitted = ph.simulate_photon_trace(20000, np.full(500, 5.0), 0.01, 4.0, random_seed=1)
    assert (dts[emitted == 0] == 0.0).all()
    assert (dts[emitted > 0] > 0.0).all()


def test_photon_trace_is_reproducible_and_seed_dependent():
    kq = np.full(100, 0.3)
    a = ph.simulate_photon_trace(5000, kq, 0.01, 4.0, random_seed=4)[0]
    np.testing.assert_array_equal(a, ph.simulate_photon_trace(5000, kq, 0.01, 4.0, random_seed=4)[0])
    assert not np.array_equal(a, ph.simulate_photon_trace(5000, kq, 0.01, 4.0, random_seed=5)[0])


def test_unseeded_traces_are_independent_samples():
    """Consecutive unseeded runs returning the same answer reads as precision
    that is not there -- which is what numba's per-thread state used to do."""
    kq = np.full(100, 0.3)
    a = ph.simulate_photon_trace(5000, kq, 0.01, 4.0)[0]
    b = ph.simulate_photon_trace(5000, kq, 0.01, 4.0)[0]
    assert not np.array_equal(a, b)


# --- the decay curve ---------------------------------------------------------

@pytest.mark.parametrize("kq_val", [0.0, 0.5])
def test_decay_falls_at_the_total_rate(kq_val):
    n_bins, dt_tac, tau0 = 512, 0.02, 4.0
    decay = np.zeros(n_bins)
    ph.simulate_quenched_decay(300, decay, dt_tac, np.full(4000, kq_val), 0.01, tau0,
                               random_seed=2)
    t = (np.arange(n_bins) + 0.5) * dt_tac
    keep = decay > decay.max() * 1e-4
    rate = -np.polyfit(t[keep], np.log(decay[keep]), 1)[0]
    assert rate == pytest.approx(1.0 / tau0 + kq_val, rel=0.02)


def test_decay_accumulates_into_the_caller_array():
    decay = np.full(256, 7.0)
    ph.simulate_quenched_decay(10, decay, 0.05, np.full(500, 0.4), 0.01, 4.0, random_seed=1)
    assert (decay >= 7.0).all()
    assert decay.sum() > 7.0 * 256


def test_decay_is_reproducible_and_seed_dependent():
    rng = np.random.default_rng(0)
    kq = np.abs(rng.normal(0.5, 0.4, 500))   # a rate that varies along the trajectory
    out = []
    for seed in (9, 9, 10):
        d = np.zeros(512)
        ph.simulate_quenched_decay(50, d, 0.02, kq, 0.01, 4.0, random_seed=seed)
        out.append(d)
    np.testing.assert_array_equal(out[0], out[1])
    assert not np.array_equal(out[0], out[2])


def test_a_constant_rate_makes_the_start_frame_irrelevant():
    """Not a seeding bug: with kq constant every starting frame is the same."""
    a, b = np.zeros(512), np.zeros(512)
    ph.simulate_quenched_decay(50, a, 0.02, np.full(500, 0.5), 0.01, 4.0, random_seed=9)
    ph.simulate_quenched_decay(50, b, 0.02, np.full(500, 0.5), 0.01, 4.0, random_seed=10)
    np.testing.assert_array_equal(a, b)


def test_decay_accumulation_is_not_racy():
    """The bin index is data-dependent, so a shared accumulator loses updates.

    The numba ancestor did exactly that and dropped up to 2 % of the intensity,
    differently each run. Blocked accumulation makes the total exact and the
    reduction order fixed.
    """
    totals = set()
    for _ in range(6):
        d = np.zeros(256)
        ph.simulate_quenched_decay(200, d, 0.05, np.full(2000, 0.3), 0.01, 4.0, random_seed=1)
        totals.add(round(float(d.sum()), 9))
    assert len(totals) == 1, f"the total varies between identical runs: {totals}"


def test_no_numba_left_in_the_particle_model():
    import ast
    import inspect
    for mod in (dif, ph):
        tree = ast.parse(inspect.getsource(mod))
        imported = {n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.module}
        imported |= {a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
        assert not any("numba" in m or m.endswith("_jit") for m in imported), (mod, imported)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))


def test_the_accept_flag_rides_in_the_returned_array():
    """Not in a SWIG out-parameter, and the reason is measurable.

    A returned ``std::vector`` becomes a Python tuple, which numpy converts at
    C speed. An out-parameter stays a wrapper object that numpy walks one
    ``__getitem__`` at a time -- about 480 ns per element. On a 500 000-step
    walk that cost 170 ms against 55 ms for the entire simulation: three
    quarters of the wall clock spent handing back one bit per step. Packing the
    flag as a fourth column took a 20-site scan from 5.17 s to 2.45 s.

    This test pins the contract that made that possible: xyz and accepted come
    from one array and therefore always have the same length.
    """
    t = dif.simulate_dye_diffusion(_ball(21, 6), dg=1.0, t_max=200.0,
                                   t_step=0.01, D=40.0, random_seed=3)
    assert t.xyz.shape == (20000, 3)
    assert t.accepted.shape == (20000,)
    assert t.accepted.dtype == np.uint8
    assert set(np.unique(t.accepted)) <= {0, 1}
    # the flag agrees with the trajectory: a rejected step does not move
    moved = np.any(np.diff(t.xyz, axis=0) != 0.0, axis=1)
    np.testing.assert_array_equal(moved, t.accepted[1:].astype(bool))
    assert int(t.accepted.sum()) == t.n_accepted


def test_the_grid_goes_through_without_being_copied():
    """The occupancy and mobility grids reach the kernel as numpy's own buffer.

    Converting a numpy array into a ``std::vector`` costs about 34 ns per
    element. On a 101-cubed grid that is 34 ms of pure marshalling per call --
    which for a short walk was the entire wall clock, and made the cost scale
    with the *grid* rather than with the number of steps:

        ng=41    2.31 ms  ->  0.05 ms
        ng=101  33.77 ms  ->  0.34 ms
        ng=151 111.52 ms  ->  0.36 ms

    Zero-copy is easy to get subtly wrong, so the properties that make it safe
    are asserted rather than assumed: the caller's array is not written to, the
    result is unchanged, and an array of the wrong dtype or layout still works
    (the adapter converts it, and then *that* is what passes through).
    """
    ng = 31
    ball = _ball(ng, 12)
    original = ball.copy()
    kw = dict(dg=1.0, t_max=500.0, t_step=0.005, D=8.0, random_seed=5)

    first = dif.simulate_dye_diffusion(ball, **kw)
    np.testing.assert_array_equal(ball, original)

    second = dif.simulate_dye_diffusion(ball, **kw)
    np.testing.assert_array_equal(first.xyz, second.xyz)
    assert first.n_accepted == second.n_accepted

    awkward = np.asfortranarray(ball).astype(np.float32)
    third = dif.simulate_dye_diffusion(awkward, **kw)
    np.testing.assert_array_equal(first.xyz, third.xyz)


def test_a_short_walk_on_a_large_grid_is_no_longer_dominated_by_the_grid():
    """It still scales with the grid, an order of magnitude more weakly.

    The SWIG marshalling is gone, but the adapter still casts the caller's array
    to int32 -- a 1 MB to 4 MB copy at 101 cubed, about 0.3 ms. So the honest
    claim is a bound, not independence: 60 steps on a 101-cubed grid took
    33.8 ms and now take about 0.34 ms. A caller that already holds int32 pays
    nothing at all.

    Loose enough to survive a busy machine; it is guarding an
    order of magnitude, not a stopwatch.
    """
    import time

    box = np.ones((101, 101, 101), dtype=np.uint8)
    kw = dict(dg=1.0, t_max=0.3, t_step=0.005, D=8.0)
    dif.simulate_dye_diffusion(box, random_seed=1, **kw)              # warm
    t0 = time.perf_counter()
    for seed in range(5):
        dif.simulate_dye_diffusion(box, random_seed=seed, **kw)
    per_call = (time.perf_counter() - t0) / 5

    assert per_call < 5e-3, (
        f"{per_call*1000:.2f} ms for a 60-step walk on a 101-cubed grid; "
        "it was 33.8 ms when the grid was marshalled element by element")
