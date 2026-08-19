"""The compiled AV distance kernels: registration, exactness, and what is only statistical.

``av/_kernels.py`` was numba until PRD-113 stage 3. Two properties of the port
are worth pinning rather than assuming:

* **``split_contact_volume`` is exact.** It is integer index arithmetic on a
  binary predicate, so the C++ and the numba it replaced agree voxel for voxel,
  and any future edit that moves a single voxel is a defect. The half-voxel and
  truncation defects this map once carried (fixed 2026-07-28, found first in
  QuEst) are what the registration test below exists to catch.
* **Everything built on ``random_distances`` is not.** It samples, and no C++
  generator reproduces numba's stream, so these are Monte-Carlo estimators that
  agree with anything else only to about :math:`1/\\sqrt{n}`. Testing them
  against recorded numbers would pin the generator, not the physics.

The coordinate kernels sit in between: exact in intent, but ``dg * i + r0``
contracts to an FMA under clang and does not under numba's LLVM, so they agree
to a few ULP and not bit for bit. That is a compiler artifact, ~7e-15 Å on a
grid, and the tolerance below says so rather than pretending otherwise.
"""

import numpy as np
import pytest

import IMP.bff
from IMP.bff.representation.av import _kernels


def test_split_is_exact_and_complementary():
    ng = 11
    rng = np.random.default_rng(3)
    density = np.where(rng.random((ng, ng, ng)) > 0.4, 1.0, 0.0)
    rs = np.array([[0.0, 0.0, 0.0], [2.0, -3.0, 1.0]])
    n_c, n_f, contact, free = _kernels.split_av_acv(
        density, 1.0, np.array([2.0, 1.5]), rs, np.zeros(3))

    occupied = int((density > 0).sum())
    assert n_c + n_f == occupied, "every occupied voxel is contact or free, and only one"
    assert not np.any(contact & free), "the two masks must not overlap"
    assert np.array_equal((contact | free).astype(bool), density > 0)
    assert contact.dtype == np.uint8 and free.dtype == np.uint8


def test_split_grid_registration_round_trips():
    """Every occupied voxel, taken to Angstrom and back, returns to itself.

    This is the check the header asks to repeat if the index map is ever
    touched. It fails for the float corner on even ``ng`` and for ``int()``
    truncation on the negative side -- the two defects the map used to have.
    """
    for ng in (8, 9, 20, 21):
        dg, r0 = 1.75, np.array([-3.0, 11.0, 0.5])
        half = (ng - 1) // 2
        idx = np.arange(ng)
        angstrom = r0[0] + (idx - half) * dg
        back = np.floor((angstrom - r0[0]) / dg).astype(int) + half
        assert np.array_equal(back, idx), f"index map is not an involution at ng={ng}"


def test_split_labels_are_the_enum():
    ng = 7
    density = np.ones((ng, ng, ng))
    density[0, 0, 0] = 0.0
    label = np.asarray(IMP.bff.split_contact_volume(
        density.ravel(), ng, 1.0, np.array([1.5]),
        np.zeros(3), np.zeros(3))).reshape(ng, ng, ng)
    assert set(np.unique(label)) <= {IMP.bff.AV_VOXEL_EMPTY,
                                     IMP.bff.AV_VOXEL_CONTACT,
                                     IMP.bff.AV_VOXEL_FREE}
    assert label[0, 0, 0] == IMP.bff.AV_VOXEL_EMPTY
    assert label[3, 3, 3] == IMP.bff.AV_VOXEL_CONTACT  # the anchor voxel itself


def test_density2points_places_voxels_where_the_grid_says():
    ng, dg = 5, 1.25
    r0 = np.array([2.0, -1.0, 7.5])
    density = np.zeros((ng, ng, ng))
    density[1, 2, 3] = 0.75
    n, pts = _kernels.density2points(ng, ng, ng, dg, density, r0)
    assert n == 1
    assert pts[0, 3] == 0.75
    np.testing.assert_allclose(pts[0, :3], r0 + dg * np.array([1, 2, 3]), atol=1e-12)


def test_weighted_mean_and_its_empty_case():
    pts = np.array([[0.0, 0.0, 0.0, 3.0], [4.0, 0.0, 0.0, 1.0]])
    np.testing.assert_allclose(_kernels.weighted_mean(pts, 2), [1.0, 0.0, 0.0], atol=1e-12)
    np.testing.assert_array_equal(_kernels.weighted_mean(np.empty((0, 4)), 0), np.zeros(3))


def test_random_distances_is_reproducible_but_seed_dependent():
    rng = np.random.default_rng(11)
    p1 = np.column_stack([rng.normal(0, 5, (300, 3)), np.ones(300)])
    p2 = np.column_stack([rng.normal(40, 5, (300, 3)), np.ones(300)])
    a = _kernels.random_distances(p1, p2, 5000, seed=4)
    assert a.shape == (5000, 2)
    np.testing.assert_array_equal(a, _kernels.random_distances(p1, p2, 5000, seed=4))
    assert not np.array_equal(a, _kernels.random_distances(p1, p2, 5000, seed=5))


def test_sampled_distances_converge_on_the_analytic_answer():
    """Two clouds a known distance apart, with the sampling error stated.

    A single point per cloud makes the estimator exact, which is the only way
    to pin a sampled kernel to a number. The spread case is checked against its
    own standard error instead.
    """
    a = np.array([[0.0, 0.0, 0.0, 1.0]])
    b = np.array([[30.0, 40.0, 0.0, 1.0]])
    assert _kernels.average_distance(a, 1, b, 1, n_samples=64) == pytest.approx(50.0)
    assert _kernels.mean_fret_distance(a, 1, b, 1, 50.0, n_samples=64) == pytest.approx(50.0)

    rng = np.random.default_rng(2)
    n = 2000
    c1 = np.column_stack([rng.normal(0, 3, (n, 3)), np.ones(n)])
    c2 = np.column_stack([rng.normal([50, 0, 0], 3, (n, 3)), np.ones(n)])
    exact = np.linalg.norm(
        c1[:, None, :3] - c2[rng.integers(0, n, 400), None, :3].reshape(-1, 1, 3).transpose(1, 0, 2),
        axis=-1).mean()
    got = _kernels.average_distance(c1, n, c2, n, n_samples=400000)
    assert abs(got - exact) < 0.1, f"{got} vs {exact}"


def test_fret_average_is_shorter_than_the_plain_average():
    """1/r^6 weights close pairs more heavily, so R_E <= <R_DA> always."""
    rng = np.random.default_rng(5)
    n = 1500
    c1 = np.column_stack([rng.normal(0, 8, (n, 3)), np.ones(n)])
    c2 = np.column_stack([rng.normal([52, 0, 0], 8, (n, 3)), np.ones(n)])
    r_da = _kernels.average_distance(c1, n, c2, n, n_samples=200000)
    r_e = _kernels.mean_fret_distance(c1, n, c2, n, 52.0, n_samples=200000)
    assert r_e < r_da


def test_fret_distance_limits():
    one = np.array([[0.0, 0.0, 0.0, 1.0]])
    assert _kernels.mean_fret_distance(one, 1, one, 1, 52.0, 100) == 0.0
    zero_weight = np.array([[0.0, 0.0, 0.0, 0.0]])
    assert np.isinf(_kernels.mean_fret_distance(one, 1, zero_weight, 1, 52.0, 100))


def test_kernels_carry_no_numba():
    """The code, not the prose -- the docstrings say "numba" precisely because
    they record what these used to be."""
    import ast
    import inspect
    import types

    tree = ast.parse(inspect.getsource(_kernels))
    imported = {
        n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.module
    } | {
        a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names
    }
    assert not any("numba" in m or m.endswith("_jit") for m in imported), imported
    assert not [
        n.name for n in ast.walk(tree)
        if isinstance(n, ast.FunctionDef) and n.decorator_list
    ], "a kernel still carries a decorator"
    for name in ("random_distances", "density2points", "weighted_mean",
                 "average_distance", "mean_fret_distance", "split_av_acv"):
        fn = getattr(_kernels, name)
        assert isinstance(fn, types.FunctionType), f"{name} is {type(fn)}, not a plain function"


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
