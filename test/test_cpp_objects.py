"""The model objects that moved to C++, and the things that moved with them.

Four objects left Python for C++ over 2026-08-19: ``LifetimeSpectrum``,
``BasicAV``/``ACV``, ``DyeDiffusionSimulation`` and ``GridDiffusionSolver``,
along with the kappa^2 sampler. Each was gated against the Python it replaced
*before* that Python was deleted -- but a gate against code that no longer
exists cannot be re-run, so what it established has to be written down as
invariants instead.

Everything here is either a property that must hold on its own terms, or a
thing that was **actually wrong at some point** during the move and would be
silent if it came back. The per-surface behaviour is already covered by the
suites those objects came with (``test_lifetime_spectrum.py``,
``medium_test_av.py``, ``test_quenching_field.py``, ``test_integrators.py``),
which were not modified by the port and are the real regression net.
"""

import numpy as np
import pytest

import IMP.bff


# --------------------------------------------------------------------------
# GridDiffusionSolver
# --------------------------------------------------------------------------
def _grid(ng=15, r=None):
    ax = np.arange(ng) - (ng - 1) / 2
    x, y, z = np.meshgrid(ax, ax, ax, indexing="ij")
    return np.sqrt(x ** 2 + y ** 2 + z ** 2)


def test_the_stability_bound_does_not_include_the_rate():
    """`k` is carried as ``exp(-k dt)``, so it constrains nothing.

    That is the *exact* solution of ``dp/dt = -k p`` over the step; it is the
    ``1 - k dt`` form that diverges once ``k dt > 1``, and this scheme does not
    use it. Passing ``k_max`` into the bound cut the allowed step by a third on
    a site with ``k_max = 96 1/ns`` and rejected steps that are perfectly
    stable. The check is the two-argument call, and a large rate must not move
    it.
    """
    assert IMP.bff.diffusion_stability_limit(8.0, 1.0) == pytest.approx(1.0 / 48.0)

    ng = 15
    bounds = (_grid(ng) < ng / 3.0).astype(np.float64)
    dmap = 8.0 * bounds
    dt = 0.9 / 48.0                       # inside the diffusion-only bound
    for k in (0.0, 96.0):                 # a rate that would fail the wrong bound
        rate = k * bounds
        solver = IMP.bff.GridDiffusionSolver(
            dmap, bounds, bounds / bounds.sum(), rate, dt, 1.0)
        solver.run(10, 5)                 # must not raise


def test_a_step_past_the_diffusion_bound_still_raises():
    ng = 15
    bounds = (_grid(ng) < ng / 3.0).astype(np.float64)
    solver = IMP.bff.GridDiffusionSolver(
        8.0 * bounds, bounds, bounds / bounds.sum(), None, 1.1 / 48.0, 1.0)
    with pytest.raises(ValueError, match="stability limit"):
        solver.run(1, 1)


def test_the_constructor_takes_the_diffusion_map_first():
    """The parameter order callers pass positionally.

    It is *not* the order the C++ constructor takes -- the shadow bridges them
    -- and getting it wrong produced a solver whose `D` and `k` were swapped,
    which fails as a stability error a long way from the cause.
    """
    ng = 15
    bounds = (_grid(ng) < ng / 3.0).astype(np.float64)
    dmap = 8.0 * bounds
    solver = IMP.bff.GridDiffusionSolver(dmap, bounds, bounds / bounds.sum())
    np.testing.assert_allclose(solver.diffusion_map, dmap)
    np.testing.assert_allclose(solver.bounds, bounds)


def test_a_domain_touching_the_outer_shell_raises():
    """The 7-point stencil cannot be evaluated there, so population would be
    discarded silently."""
    ng = 11
    bounds = np.ones((ng, ng, ng))
    with pytest.raises(ValueError, match="outer shell"):
        IMP.bff.GridDiffusionSolver(
            bounds, bounds, bounds / bounds.sum(), None, 1e-4, 1.0).run(1, 1)


# --------------------------------------------------------------------------
# DyeDiffusionSimulation
# --------------------------------------------------------------------------
def _walk(ng=15):
    r = _grid(ng)
    occ = (r < ng / 3.0).astype(np.uint8)
    return IMP.bff.DyeDiffusionSimulation(
        density=occ, dg=1.0, x0=np.zeros(3),
        quenching_rate_map=np.where(r < ng / 6.0, 0.3, 0.0))


def test_one_seed_gives_one_walk_and_a_different_seed_a_different_one():
    a, b, c = _walk(), _walk(), _walk()
    ta = a.run(t_max=40.0, n_trajectories=1, random_seed=3)
    tb = b.run(t_max=40.0, n_trajectories=1, random_seed=3)
    tc = c.run(t_max=40.0, n_trajectories=1, random_seed=4)
    np.testing.assert_array_equal(ta, tb)
    assert not np.array_equal(ta, tc)


def test_concatenated_walks_do_not_repeat_one_trajectory():
    """Seeds are strided by a large prime, not reused.

    Reusing the base seed would run n identical trajectories and concatenate
    them, which looks like n times the sampling and is none of it.
    """
    w = _walk()
    w.run(t_max=40.0, n_trajectories=4, random_seed=11)
    n = w.n_frames // 4
    first = w.trajectory[:n]
    for k in range(1, 4):
        assert not np.array_equal(first, w.trajectory[k * n:(k + 1) * n])


def test_replacing_the_volume_discards_the_trajectory():
    """It was a walk in the *old* volume; keeping it would let a caller read
    positions the new occupancy forbids. The Python allowed exactly that."""
    w = _walk()
    w.run(t_max=40.0, n_trajectories=1, random_seed=5)
    assert w.n_frames > 0
    w.density = np.zeros((15, 15, 15), dtype=np.uint8)
    assert w.n_frames == 0
    assert w.trajectory is None
    assert w.n_accepted == 0


def test_each_grid_takes_its_shape_from_its_own_length():
    """The fields are stored separately and are only conventionally the same
    size, so swapping the volume for a smaller one must not make the rate map
    unreadable."""
    w = _walk(ng=15)
    w.density = np.zeros((11, 11, 11), dtype=np.uint8)
    assert w.density.shape == (11, 11, 11)
    assert w.quenching_rate_map.shape == (15, 15, 15)


def test_sample_grid_uses_the_integer_centre_and_floors():
    """One convention, and both halves of it were bugs once.

    The offset is the *integer* ``grid_center_index``, which differs from the
    float ``(ng - 1) / 2`` on every even edge length -- and even is the normal
    case. The conversion is ``floor``, not ``trunc``: ``trunc`` maps
    ``[-1, 0)`` to 0, so a position up to one voxel below the grid would read
    voxel 0's value.
    """
    assert IMP.bff.grid_center_index(15) == 7
    assert IMP.bff.grid_center_index(16) == 7          # not 7.5
    w = _walk()
    w.run(t_max=40.0, n_trajectories=1, random_seed=5)
    ramp = np.arange(15 ** 3, dtype=np.float64).reshape(15, 15, 15)
    got = w.sample_grid(ramp)
    centre = IMP.bff.grid_center_index(15)
    idx = np.floor(w.trajectory / w.dg + centre).astype(np.int64)
    inside = np.all((idx >= 0) & (idx < 15), axis=1)
    want = np.zeros(w.n_frames, dtype=np.float64)
    want[inside] = ramp[tuple(idx[inside].T)]
    np.testing.assert_array_equal(got.astype(np.float64), want)


# --------------------------------------------------------------------------
# kappa^2
# --------------------------------------------------------------------------
def test_the_isotropic_limit_is_exactly_two_thirds():
    """Drawing directions with three standard normals is load-bearing.

    Normalising a *uniform* draw fills the cube's positive octant rather than
    the sphere, which halved <kappa^2> to 0.333 here. This limit is the check
    that caught it: with both order parameters zero there is no wobble at all,
    so every sample must be 2/3 exactly.
    """
    out = np.asarray(IMP.bff.sample_kappa2_diffusion_with_traps(
        0.0, 0.0, 0.5, 20000, 31, 0.0, 4.0, 7))
    samples = out[61:]
    np.testing.assert_allclose(samples, 2.0 / 3.0, rtol=1e-12)


def test_the_fully_trapped_limit_spans_the_whole_range():
    out = np.asarray(IMP.bff.sample_kappa2_diffusion_with_traps(
        1.0, 1.0, 0.5, 50000, 31, 0.0, 4.0, 7))
    samples = out[61:]
    assert samples.min() < 0.05 and samples.max() > 3.9
    assert samples.mean() == pytest.approx(2.0 / 3.0, abs=0.05)


def test_the_distance_ratio_transform_is_zero_outside_the_sampled_range():
    """``np.interp(..., left=0, right=0)``, not clamped.

    The transform puts a 5 % margin either side of the transformed points, so
    the first and last bins are outside the range *by construction*. Clamping
    instead planted the extreme weight on a ratio no kappa^2 in the input maps
    to -- a 0.12 error in the weights while the axis and <kappa^2> still
    matched to 2e-16.
    """
    k2_val = np.array([0.5, 1.0, 2.0])
    k2_amp = np.array([1.0, 1.0, 1.0])
    axis, w, mean = np.split(
        np.asarray(IMP.bff.kappa2_distance_ratio_transform(k2_amp, k2_val, 32)),
        [32, 64])
    assert w[0] == 0.0 and w[-1] == 0.0
    assert w.sum() == pytest.approx(1.0)


def test_the_distance_ratio_transform_carries_its_jacobian():
    """Transforming the abscissa and carrying the weights across unchanged
    would be wrong wherever the mapping is non-linear, which is everywhere."""
    k2_val = np.array([0.4, 0.8, 1.6, 3.2])
    k2_amp = np.ones(4)
    out = np.asarray(IMP.bff.kappa2_distance_ratio_transform(k2_amp, k2_val, 64))
    mean = out[128]
    assert mean == pytest.approx(np.mean(k2_val))
    # r = (<k2>/k2)^(1/6) is decreasing in k2, so the smallest k2 gives the
    # largest ratio -- the transform must not silently reverse it.
    axis = out[:64]
    assert axis[0] < (mean / k2_val.max()) ** (1 / 6.0)
    assert axis[-1] > (mean / k2_val.min()) ** (1 / 6.0)


@pytest.mark.parametrize("amp,val,match", [
    (np.ones(3), -np.ones(3), "strictly positive"),
    (np.ones(3), np.ones(4), "same shape"),
    (-np.ones(3), np.ones(3), "non-negative"),
    (np.zeros(3), np.ones(3), "must be positive"),
])
def test_the_transform_refuses_impossible_input(amp, val, match):
    with pytest.raises(ValueError, match=match):
        IMP.bff.kappa2_distance_ratio_transform(amp, val, 8)


# --------------------------------------------------------------------------
# the objects are values
# --------------------------------------------------------------------------
def test_a_spectrum_is_a_value_not_a_handle():
    """`IMP_SWIG_VALUE`, so a copy is independent."""
    s = IMP.bff.LifetimeSpectrum(np.array([0.4, 0.2]), np.array([0.25, 1.0]))
    t = IMP.bff.LifetimeSpectrum(s.amplitudes, s.rate_constants, s.exact)
    assert t.n_species == 2
    np.testing.assert_array_equal(t.amplitudes, s.amplitudes)


def test_an_av_carries_no_grid_when_it_was_built_from_points():
    """``density`` is ``None``, not an empty array: every consumer tests it."""
    pts = np.ascontiguousarray(np.random.default_rng(0).random((40, 4)))
    av = IMP.bff.BasicAV(points=pts, position_name="donor")
    assert av.density is None
    assert av.grid_origin is None
    assert av.n_points == 40
    assert av.points.shape == (40, 4)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))


def test_the_molecular_graph_is_the_systems_own():
    """`get_bonded_neighbors`, `get_exclusions`, `find_rings`, `is_within_bonds`.

    All four read only the bonds, angles and torsions the system already holds,
    so they are methods on it rather than Python functions taking it apart.
    Before, three modules had their own copies -- `cgdye.sim`, `cgdye.topology`
    and `scoring` -- which agreed on the shipped system and had drifted into
    three container shapes for one answer (a set of tuples, a set of frozensets,
    a defaultdict of sets).

    The numbers here are the Python's, recorded when it was replaced.
    """
    import IMP.bff
    from IMP.bff.tools import get_template_dir, get_structure_dir
    from IMP.bff.cgdye.topology import build_dye_protein_system, build_graph, find_cycles

    system = build_dye_protein_system(
        str(get_structure_dir("cx4.mol2")), str(get_structure_dir("atto655.mol2")),
        "CX4", "atto655",
        protein_template=str(get_template_dir("cx4.template.cif")),
        dye_template=str(get_template_dir("atto655.template.cif")),
    )

    # the graph, against the builder that still serves raw bond lists
    from_bonds = {k: set(v) for k, v in
                  build_graph([(b.site_a, b.site_b) for b in system.bonds]).items()}
    from_system = {k: set(v) for k, v in system.get_bonded_neighbors().items()}
    assert from_system == from_bonds

    # rings, against the Python cycle finder over that graph
    assert ({tuple(sorted(c)) for c in system.find_rings(8)}
            == {tuple(sorted(c)) for c in find_cycles(from_bonds, max_len=8)})
    assert len(system.find_rings(8)) == 9

    # exclusions: 1-2 from bonds, 1-3 from angle ends, 1-4 from dihedral ends
    excl = system.exclusions()
    assert len(excl) == 600
    assert all(a <= b for a, b in excl), "each pair is sorted"
    for bond in system.bonds:
        assert tuple(sorted((bond.site_a, bond.site_b))) in excl
    for angle in system.angles:
        assert tuple(sorted((angle.site_a, angle.site_c))) in excl
    for tor in system.dihedrals:
        assert tuple(sorted((tor.site_a, tor.site_d))) in excl

    # impropers are built now (they were dropped by this builder until the
    # three builders were collapsed into one), so `include_impropers` is live:
    # it adds each improper's a-c, a-d and b-d pairs, and every one of those is
    # already 1-2 or 1-3 bonded on these components -- so the flag changes
    # nothing here while no longer being inert by accident.
    assert len(system.impropers) == 81
    assert system.exclusions(False) <= excl
    for tor in system.impropers:
        for x, y in ((tor.site_a, tor.site_c), (tor.site_a, tor.site_d),
                     (tor.site_b, tor.site_d)):
            assert tuple(sorted((x, y))) in excl

    # reachability, bounded by bond count
    first = system.bonds[0]
    assert system.is_within_bonds(first.site_a, first.site_b, 1)
    assert system.is_within_bonds(first.site_a, first.site_a, 0)
    # depth actually bounds: the far end of an angle is two bonds away, so it
    # is reachable at 2 and not at 1
    angle = next(a for a in system.angles if a.site_a != a.site_c)
    assert system.is_within_bonds(angle.site_a, angle.site_c, 2)
    assert not system.is_within_bonds(angle.site_a, angle.site_c, 1)


def test_the_molecular_graph_matches_the_python_it_replaced():
    """`MolecularGraph` against the Python derivations, on three real dyes.

    Adjacency, angles, torsions and rings were Python in three modules. The
    numbers below were produced by that Python and are what the C++ has to
    reproduce; the shipped MOL2 files are the fixture because a hand-built
    graph would not have found the two things that did go wrong -- key order
    (the Python's `defaultdict` is in bond-insertion order, and a seeded
    sampler walks it) and node type (`scoring` builds graphs over site-id
    strings, not MOL2 serials).
    """
    import IMP.bff
    from IMP.bff.tools import get_structure_dir
    from IMP.bff.cgdye.topology import (
        build_angles, build_dihedrals, build_graph, find_cycles, parse_dye_mol2)

    expected = {                      # atoms, angles, dihedrals, rings
        "atto655.mol2": (70, 138, 207, 5),
        "cx4.mol2": (68, 124, 176, 4),
        "alexa488_r48.mol2": (83, 151, 218, 5),
    }
    for mol2, (n_atoms, n_ang, n_dih, n_ring) in expected.items():
        atoms, bonds = parse_dye_mol2(str(get_structure_dir(mol2)), "X")
        assert len(atoms) == n_atoms
        graph = build_graph(bonds)
        assert len(build_angles(graph)) == n_ang
        assert len(build_dihedrals(graph)) == n_dih
        assert len(find_cycles(graph, max_len=8)) == n_ring

        # and the same through the C++ type directly
        cpp = IMP.bff.MolecularGraph(sorted(bonds))
        assert {tuple(a) for a in cpp.get_angles()} == set(build_angles(graph))
        assert {tuple(d) for d in cpp.get_dihedrals()} == set(build_dihedrals(graph))

    # key order is first-appearance, not ascending
    bonds = [(9, 4), (4, 7), (7, 1)]
    assert list(build_graph(bonds)) == [9, 4, 7, 1]

    # nodes need not be integers: scoring builds graphs over site ids
    named = [("dye:C1", "dye:C2"), ("dye:C2", "dye:C3")]
    assert build_angles(build_graph(named)) == [("dye:C1", "dye:C2", "dye:C3")]
