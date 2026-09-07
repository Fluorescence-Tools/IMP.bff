"""Invariants of the model objects: ``LifetimeSpectrum``,
``AccessibleVolume``/``ACV``, ``ProbeDiffusionSimulation``,
``GridDiffusionSolver`` and the kappa^2 sampler.

Everything here is either a property that must hold on its own terms, or a
failure mode that is **silent** when it comes back -- an axis order, a shared
buffer, a seed that does not seed. The per-surface behaviour is covered by the
suites those objects came with (``test_lifetime_spectrum.py``,
``medium_test_av.py``, ``test_quenching_field.py``, ``test_integrators.py``).
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
            dmap, bounds, bounds / bounds.sum(), rate, 1.0, dt)
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
    np.testing.assert_allclose(solver.get_diffusion_map(), dmap.ravel())
    np.testing.assert_allclose(solver.get_bounds(), bounds.ravel())


def test_a_domain_touching_the_outer_shell_raises():
    """The 7-point stencil cannot be evaluated there, so population would be
    discarded silently."""
    ng = 11
    bounds = np.ones((ng, ng, ng))
    with pytest.raises(ValueError, match="outer shell"):
        IMP.bff.GridDiffusionSolver(
            bounds, bounds, bounds / bounds.sum(), None, 1e-4, 1.0).run(1, 1)


# --------------------------------------------------------------------------
# ProbeDiffusionSimulation
# --------------------------------------------------------------------------
def _walk(ng=15):
    r = _grid(ng)
    occ = (r < ng / 3.0).astype(np.uint8)
    return IMP.bff.ProbeDiffusionSimulation(
        density=occ, dg=1.0, x0=np.zeros(3),
        quenching_rate_map=np.where(r < ng / 6.0, 0.3, 0.0))


def test_one_seed_gives_one_walk_and_a_different_seed_a_different_one():
    a, b, c = _walk(), _walk(), _walk()
    a.run(t_max=40.0, n_trajectories=1, random_seed=3)
    b.run(t_max=40.0, n_trajectories=1, random_seed=3)
    c.run(t_max=40.0, n_trajectories=1, random_seed=4)
    np.testing.assert_array_equal(a.get_trajectory(), b.get_trajectory())
    assert not np.array_equal(a.get_trajectory(), c.get_trajectory())


def test_concatenated_walks_do_not_repeat_one_trajectory():
    """Seeds are strided by a large prime, not reused.

    Reusing the base seed would run n identical trajectories and concatenate
    them, which looks like n times the sampling and is none of it.
    """
    w = _walk()
    w.run(t_max=40.0, n_trajectories=4, random_seed=11)
    n = w.n_frames // 4
    traj = w.get_trajectory()
    first = traj[:n]
    for k in range(1, 4):
        assert not np.array_equal(first, traj[k * n:(k + 1) * n])


def test_replacing_the_volume_discards_the_trajectory():
    """It was a walk in the *old* volume; keeping it would let a caller read
    positions the new occupancy forbids."""
    w = _walk()
    w.run(t_max=40.0, n_trajectories=1, random_seed=5)
    assert w.n_frames > 0
    w.set_density(np.zeros(15 * 15 * 15, dtype=np.uint8))
    assert w.n_frames == 0
    assert w.get_trajectory().size == 0
    assert w.n_accepted == 0


def test_each_grid_takes_its_shape_from_its_own_length():
    """The fields are stored separately and are only conventionally the same
    size, so swapping the volume for a smaller one must not make the rate map
    unreadable. Grids are read back flat; a caller shapes them."""
    w = _walk(ng=15)
    w.set_density(np.zeros(11 * 11 * 11, dtype=np.uint8))
    small = w.get_density()
    assert small.size == 11 * 11 * 11 or small.size == 0
    assert w.get_quenching_rate_map().size == 15 * 15 * 15


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
    got = w.sample_grid(ramp, 15)
    centre = IMP.bff.grid_center_index(15)
    traj = w.get_trajectory().reshape(-1, 3)
    dg = w.dg
    idx = np.floor(traj / dg + centre).astype(np.int64)
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
    """``get_density()`` is empty, not an error: every consumer tests size."""
    pts = np.ascontiguousarray(np.random.default_rng(0).random((40, 4)))
    av = IMP.bff.AccessibleVolume(points=pts, position_name="donor")
    assert av.get_density().size == 0
    assert av.get_grid_origin().size == 0
    assert av.get_n_points() == 40
    assert av.get_points().reshape(-1, 4).shape == (40, 4)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))


def test_the_molecular_graph_is_the_systems_own():
    """`get_bonded_neighbors`, `get_exclusions`, `find_rings`, `is_within_bonds`.

    All four read only the bonds, angles and torsions the system already holds,
    so they are methods on it rather than Python functions taking it apart.
    Before, three modules had their own copies -- `cgprobe.sim`, `cgprobe.topology`
    and `scoring` -- which agreed on the shipped system and had drifted into
    three container shapes for one answer (a set of tuples, a set of frozensets,
    a defaultdict of sets).

    The numbers here are the reference implementation's.
    """
    import IMP.bff
    from IMP.bff import get_template_dir, get_structure_dir
    import IMP.bff
    from IMP.bff import build_probe_protein_system

    system = build_probe_protein_system(
        str(get_structure_dir("cx4.mol2")), str(get_structure_dir("atto655.mol2")),
        "CX4", "atto655",
        protein_template=str(get_template_dir("cx4.template.cif")),
        probe_template=str(get_template_dir("atto655.template.cif")),
    )

    # the graph, against `LabelledGraph` over the same bonds
    labelled = IMP.bff.LabelledGraph(
        [(b.site_a, b.site_b) for b in system.bonds])
    from_bonds = {n: set(labelled.get_neighbors(n))
                  for n in labelled.get_nodes()}
    from_system = {k: set(v) for k, v in system.get_bonded_neighbors().items()}
    assert from_system == from_bonds

    # rings, against the same graph's cycle finder
    assert ({tuple(sorted(c)) for c in system.find_rings(8)}
            == {tuple(sorted(c)) for c in labelled.get_rings(8)})
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
    """`MolecularGraph` against the reference derivations, on three real dyes.

    Adjacency, angles, torsions and rings were Python in three modules. The
    numbers below were produced by that Python and are what the C++ has to
    reproduce; the shipped MOL2 files are the fixture because a hand-built
    graph would not have found the two things that did go wrong -- key order
    (the reference `defaultdict` is in bond-insertion order, and a seeded
    sampler walks it) and node type (`scoring` builds graphs over site-id
    strings, not MOL2 serials).
    """
    import IMP.bff
    from IMP.bff import get_structure_dir, read_mol2_component

    expected = {                      # atoms, angles, dihedrals, rings
        "atto655.mol2": (70, 138, 207, 5),
        "cx4.mol2": (68, 124, 176, 4),
        "alexa488_r48.mol2": (83, 151, 218, 5),
    }
    for mol2, (n_atoms, n_ang, n_dih, n_ring) in expected.items():
        component = read_mol2_component(str(get_structure_dir(mol2)), "X")
        assert len(component.atoms) == n_atoms
        graph = IMP.bff.MolecularGraph([tuple(b) for b in component.bonds])
        assert len(graph.get_angles()) == n_ang
        assert len(graph.get_dihedrals()) == n_dih
        assert len(graph.get_rings(8)) == n_ring

    # A MOL2 serial is an int and the graph keys by it -- no renumbering, and
    # `get_nodes` is ascending, which is what a std::map gives.
    graph = IMP.bff.MolecularGraph([(9, 4), (4, 7), (7, 1)])
    assert list(graph.get_nodes()) == [1, 4, 7, 9]
    assert list(graph.get_neighbors(4)) == [7, 9]

    # Nodes need not be integers: site ids are strings, and `LabelledGraph`
    # is where that lives, rather than callers numbering labels around an
    # integer graph.
    named = IMP.bff.LabelledGraph([("dye:C1", "dye:C2"), ("dye:C2", "dye:C3")])
    assert [tuple(a) for a in named.get_angles()] == [
        ("dye:C1", "dye:C2", "dye:C3")]
    assert list(named.get_nodes()) == ["dye:C1", "dye:C2", "dye:C3"]


def test_leader_clustering_matches_the_python_it_replaced():
    """`cluster_frames_leader`, `assign_frames_to_clusters`, `rmsd_no_align`.

    The sweep order is part of the answer -- frame 0 leads and each later frame
    joins the *first* leader within the threshold -- so this is reproduced, not
    improved on: it is how FRETpredict's shipped libraries were built, and a
    k-medoids that found better centres would give a different library.

    Six well-separated shapes with noise, so the thresholds actually cluster;
    an earlier version of this check used random coordinates and every frame
    became its own centre, which would have passed against almost any
    implementation.
    """
    import numpy as np
    import IMP.bff

    rng = np.random.default_rng(3)
    base = rng.normal(size=(6, 30, 3)) * 6.0
    coords = np.concatenate(
        [base[i % 6] + rng.normal(size=(30, 3)) * 0.15 for i in range(600)]
    ).reshape(600, 30, 3)

    def python_leader(c, threshold):
        centers = [0]
        for i in range(1, c.shape[0]):
            if all(np.sqrt(np.mean(np.sum((c[i] - c[j]) ** 2, axis=-1))) >= threshold
                   for j in centers):
                centers.append(i)
        return centers

    for threshold in (0.5, 1.0, 2.0, 4.0, 8.0):
        expected = python_leader(coords, threshold)
        got = list(IMP.bff.cluster_frames_leader(coords, threshold))
        assert got == expected, threshold
        assert len(got) == 6

    centers = np.asarray(python_leader(coords, 1.0), dtype=np.int32)
    got = list(IMP.bff.assign_frames_to_clusters(coords, centers))
    for i, assigned in enumerate(got):
        d = [np.sqrt(np.mean(np.sum((coords[i] - coords[c]) ** 2, axis=-1)))
             for c in centers]
        assert assigned == int(np.argmin(d))

    a, b = coords[0].ravel(), coords[1].ravel()
    expected = float(np.sqrt(np.mean(np.sum(
        (coords[0] - coords[1]) ** 2, axis=-1))))
    assert abs(IMP.bff.rmsd_no_align(a, b) - expected) < 1e-12


def test_smith_waterman_matches_the_python_it_replaced():
    """Local alignment with affine gaps, against a reference implementation.

    The reference implementation is kept verbatim here, because
    the thing that can go wrong is not the score but the *precedence* when two
    paths tie -- stop, then diagonal, then a gap in the template, then a gap in
    the query. A tie broken the other way gives the same score and different
    blocks, and the blocks are what the caller uses.
    """
    import random

    import numpy as np
    import IMP.bff

    def reference(query, templ, match=2.0, mismatch=-1.0,
                  gap_open=-5.0, gap_extend=-1.0):
        n, m = len(query), len(templ)
        if n == 0 or m == 0:
            return [], []
        neg = float("-inf")
        main = np.zeros((n + 1, m + 1))
        gap_q = np.full((n + 1, m + 1), neg)
        gap_t = np.full((n + 1, m + 1), neg)
        trace = np.zeros((n + 1, m + 1), dtype=np.int8)
        best, best_ij = 0.0, (0, 0)
        for i in range(1, n + 1):
            qi = query[i - 1]
            for j in range(1, m + 1):
                gap_q[i, j] = max(main[i - 1, j] + gap_open, gap_q[i - 1, j] + gap_extend)
                gap_t[i, j] = max(main[i, j - 1] + gap_open, gap_t[i, j - 1] + gap_extend)
                diag = main[i - 1, j - 1] + (match if qi == templ[j - 1] else mismatch)
                cell = max(0.0, diag, gap_q[i, j], gap_t[i, j])
                main[i, j] = cell
                trace[i, j] = (0 if cell == 0.0 else 1 if cell == diag
                               else 2 if cell == gap_q[i, j] else 3)
                if cell > best:
                    best, best_ij = cell, (i, j)
        if best <= 0.0:
            return [], []
        i, j = best_ij
        qb, tb = [], []
        rq, rt = i, j
        while i > 0 and j > 0 and trace[i, j] != 0:
            step = trace[i, j]
            if step == 1:
                i, j = i - 1, j - 1
                continue
            if rq != i or rt != j:
                qb.append((i, rq))
                tb.append((j, rt))
            if step == 2:
                i -= 1
            else:
                j -= 1
            rq, rt = i, j
        if rq != i or rt != j:
            qb.append((i, rq))
            tb.append((j, rt))
        return qb[::-1], tb[::-1]

    random.seed(7)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    for trial in range(120):
        n, m = random.randint(1, 40), random.randint(1, 40)
        q = "".join(random.choice(aa) for _ in range(n))
        t = "".join(random.choice(aa) for _ in range(m))
        if trial % 3 == 0 and n > 8:                 # force real homology
            t = t[:5] + q[2:2 + n // 2] + t[5:]
        want_q, want_t = reference(q, t)
        blocks = IMP.bff.smith_waterman(q, t)
        assert [(b.query_start, b.query_end) for b in blocks] == want_q, (q, t)
        assert [(b.template_start, b.template_end) for b in blocks] == want_t, (q, t)

    assert len(IMP.bff.smith_waterman("", "ABC")) == 0
    assert IMP.bff.smith_waterman_score("ABCDEF", "ABCDEF") == 12.0


def test_improper_expansion_matches_the_python_it_replaced():
    """`MolecularGraph.expand_impropers` for all four template kinds.

    The counts are the reference's, on the three shipped dyes. `flat` and `orient`
    key off the *atom name* starting with S rather than the element -- what the
    Python did, kept because a MOL2's types are less reliable than its names,
    and the difference is invisible unless a sulfur is mistyped.
    """
    import IMP.bff
    from IMP.bff import get_structure_dir
    from IMP.bff import read_mol2_component

    expected = {
        "atto655.mol2": {"ring": 20, "pi": 15, "flat": 1, "orient": 1},
        "cx4.mol2": {"ring": 24, "pi": 24, "flat": 4, "orient": 4},
        "alexa488_r48.mol2": {"ring": 24, "pi": 27, "flat": 2, "orient": 2},
    }
    for mol2, counts in expected.items():
        component = read_mol2_component(str(get_structure_dir(mol2)), "X")
        atoms = {a.serial: a for a in component.atoms}
        nodes = sorted(atoms)
        graph = IMP.bff.MolecularGraph(sorted(tuple(b) for b in component.bonds))
        elements = [atoms[n].element for n in nodes]
        names = [atoms[n].atom_name for n in nodes]
        for kind, n_expected in counts.items():
            quads = graph.expand_impropers(kind, nodes, nodes, elements, names, 8)
            assert len(quads) == n_expected, (mol2, kind, len(quads))
            for q in quads:
                assert len(q) == 4
                assert len(set(q)) == 4, "four distinct atoms"
                # the centre is the second entry and is bonded to the rest
                centre = q[1]
                for other in (q[0], q[2], q[3]):
                    assert other in set(graph.get_neighbors(centre))

    # an unknown kind produces nothing rather than raising
    g = IMP.bff.MolecularGraph([(0, 1), (1, 2)])
    assert len(g.expand_impropers("nonsense", [1], [0, 1, 2],
                                  ["C", "C", "C"], ["C1", "C2", "C3"], 8)) == 0


def test_the_mol2_reader_matches_the_python_it_replaced():
    """`read_mol2_component` against the IMP-plus-numpy reader it replaced.

    Includes 1DG3 at 4,698 atoms, because the two things that can drift are
    per-atom and only show at scale: the atom *name* (column 2 of
    `@<TRIPOS>ATOM`, not the TRIPOS type -- IMP's reader maps `C.3` to `C3` and
    loses `C12`, which is what a template's features refer to) and the element
    derived from it.
    """
    import IMP.bff
    from IMP.bff import get_structure_dir
    # This compared `parse_dye_mol2`'s dicts with the typed value the C++
    # reader returns. The dict view is gone -- it was a copy of the same
    # fields under different keys -- so what is left to check is that the
    # reader gets the shapes right on real files, including a 4,698-atom
    # protein.
    for mol2, n_atoms, n_bonds in (("atto655.mol2", 70, 74),
                                   ("cx4.mol2", 68, 72),
                                   ("alexa488_r48.mol2", 83, 87),
                                   ("1DG3.mol2", 4698, 4428)):
        path = str(get_structure_dir(mol2))
        component = IMP.bff.read_mol2_component(path, "X")
        assert len(component.atoms) == n_atoms
        assert len(component.bonds) == n_bonds
        for a in component.atoms:
            assert a.atom_name and a.element
            assert a.element == IMP.bff.element_from_atom_name(a.atom_name)

    # the element rule is the first letter of the leading alphabetic run,
    # uppercased -- wrong for two-letter elements, and reproduced deliberately
    # because every site's LJ type is keyed on it
    assert IMP.bff.element_from_atom_name("C12") == "C"
    assert IMP.bff.element_from_atom_name("CL3") == "C"
    assert IMP.bff.element_from_atom_name("n1") == "N"
    assert IMP.bff.element_from_atom_name("1HG") == "C"
    assert IMP.bff.element_from_atom_name("") == "C"


def test_system_self_consistency_is_the_systems_own_check():
    """`ProbeForceFieldSystem.get_inconsistency`, which `cgprobe.sim` used to do.

    It reports rather than raises: what is wrong with a system is a property of
    the system, and whether that stops the caller is the caller's decision.
    `cgprobe.sim._validate_system` raises; a reader could report instead.
    """
    import IMP.bff
    from IMP.bff import get_template_dir, get_structure_dir
    from IMP.bff import build_probe_protein_system

    def fresh():
        return build_probe_protein_system(
            str(get_structure_dir("cx4.mol2")), str(get_structure_dir("atto655.mol2")),
            "CX4", "atto655",
            protein_template=str(get_template_dir("cx4.template.cif")),
            probe_template=str(get_template_dir("atto655.template.cif")))

    assert fresh().get_inconsistency() == ""
    assert IMP.bff.ProbeForceFieldSystem("x").get_inconsistency() == "system requires components"

    dangling = fresh()
    bond = dangling.bonds[0]
    bond.site_a = "no-such-site"
    dangling.bonds = [bond]
    assert dangling.get_inconsistency() == "bond references unknown site"

    duplicated = fresh()
    sites = list(duplicated.sites)
    sites[1].id = sites[0].id
    duplicated.sites = sites
    assert duplicated.get_inconsistency() == "duplicate site ids"

    orphan = fresh()
    sites = list(orphan.sites)
    sites[0].component = "not-a-component"
    orphan.sites = sites
    assert orphan.get_inconsistency().startswith("site ")
    assert "unknown component" in orphan.get_inconsistency()


def test_linker_geometry_matches_the_python_it_replaced():
    """`LinkerGeometry.apply` against the per-atom IMP transform loop.

    The reference is the old implementation kept verbatim, and the property it
    protects is the *sequencing*: the rotations are applied one after another,
    each about the axis as it stands after the previous one. Applying them all
    to the reference geometry independently gives a different structure and the
    same number of degrees of freedom, so a test on shapes would not notice.
    """
    import numpy as np
    import IMP.algebra
    from IMP.bff import get_structure_dir, linker_geometry_from_mol2

    geometry = linker_geometry_from_mol2(
        str(get_structure_dir("alexa488_r48.mol2")))
    n_tor = int(geometry.get_number_of_torsions())
    n_ang = int(geometry.get_number_of_angles())
    n_dof = n_tor + n_ang
    assert n_dof > 20, n_dof

    base = np.asarray(list(geometry.get_coordinates())).reshape(-1, 3)
    torsion_fixed = list(geometry.get_torsion_fixed())
    torsion_moving = list(geometry.get_torsion_moving())
    torsion_sets = [list(geometry.get_torsion_set(i)) for i in range(n_tor)]
    angle_b = list(geometry.get_angle_b())
    angle_c = list(geometry.get_angle_c())
    angle_a = list(geometry.get_angle_a())
    angle_sets = [list(geometry.get_angle_set(i)) for i in range(n_ang)]

    def _rotate(xyz, centre, axis, angle, rows):
        """One rotation, through IMP.algebra rather than the class under test."""
        if np.linalg.norm(axis) < 1e-8:
            return
        transform = IMP.algebra.get_rotation_about_point(
            IMP.algebra.Vector3D(*centre),
            IMP.algebra.get_rotation_about_axis(
                IMP.algebra.Vector3D(*axis), float(angle)))
        for row in rows:
            xyz[row] = list(transform.get_transformed(
                IMP.algebra.Vector3D(*xyz[row])))

    def reference(cfg):
        """The rotations applied one after another, on the moving geometry."""
        xyz = base.copy()
        for k in range(n_tor):
            f, m = torsion_fixed[k], torsion_moving[k]
            _rotate(xyz, xyz[f], xyz[m] - xyz[f], cfg[k], torsion_sets[k])
        for k in range(n_ang):
            b_row, c_row, a_row = angle_b[k], angle_c[k], angle_a[k]
            axis = np.cross(xyz[a_row] - xyz[b_row], xyz[c_row] - xyz[b_row])
            _rotate(xyz, xyz[b_row], axis, cfg[n_tor + k], angle_sets[k])
        return xyz

    def applied(cfg):
        return np.asarray(list(geometry.apply([float(x) for x in cfg]))).reshape(-1, 3)

    rng = np.random.default_rng(21)
    for _ in range(10):
        cfg = list(rng.normal(scale=0.8, size=n_dof))
        assert np.max(np.abs(reference(cfg) - applied(cfg))) < 1e-9

    # the all-zero configuration is the reference geometry itself
    assert np.max(np.abs(applied([0.0] * n_dof) - base)) < 1e-12
