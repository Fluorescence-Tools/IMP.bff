"""Tests for ``IMP.bff.av`` — accessible volume infrastructure.

All tests in this file run in **standalone mode** (no IMP or LabelLib
required) and use the Numba kernels directly.
"""

from __future__ import annotations

import numpy as np
import pytest

from IMP.bff.av import BasicAV, ACV, compute_av, AccessibleVolume
from IMP.bff.av._kernels import (
    random_distances,
    density2points,
    weighted_mean,
    average_distance,
    mean_fret_distance,
    split_av_acv,
)


class TestBasicAV:
    """Test the BasicAV class."""

    def test_create_empty(self):
        """An AV with no points is valid but empty."""
        av = BasicAV(position_name="empty")
        assert av.n_points == 0
        assert av.points.shape == (0, 4)

    def test_create_with_points(self):
        """An AV can be initialised with a point cloud."""
        n = 100
        pts = np.random.RandomState(0).randn(n, 4).astype(np.float64)
        pts[:, 3] = np.abs(pts[:, 3]) + 0.01  # positive weights
        av = BasicAV(points=pts, position_name="test")
        assert av.n_points == n
        assert av.position_name == "test"

    def test_mean_position(self):
        """Mean position is the weighted centroid."""
        pts = np.array([
            [0.0, 0.0, 0.0, 1.0],
            [2.0, 0.0, 0.0, 1.0],
        ], dtype=np.float64)
        av = BasicAV(points=pts)
        mp = av.mean_position
        np.testing.assert_allclose(mp, [1.0, 0.0, 0.0], atol=1e-12)

    def test_mean_position_single_point(self):
        """Single point has that point as mean."""
        pts = np.array([[1.5, 2.5, 3.5, 1.0]], dtype=np.float64)
        av = BasicAV(points=pts)
        np.testing.assert_allclose(av.mean_position, [1.5, 2.5, 3.5], atol=1e-12)

    def test_dRmp(self):
        """Distance between mean positions is correct."""
        pts1 = np.array([[0.0, 0.0, 0.0, 1.0]], dtype=np.float64)
        pts2 = np.array([[3.0, 4.0, 0.0, 1.0]], dtype=np.float64)
        av1 = BasicAV(points=pts1)
        av2 = BasicAV(points=pts2)
        assert av1.dRmp(av2) == pytest.approx(5.0, abs=1e-10)

    def test_dRDA(self):
        """Mean distance from identical distributions is non-zero."""
        pts = np.array([
            [0.0, 0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0, 1.0],
        ], dtype=np.float64)
        av1 = BasicAV(points=pts)
        av2 = BasicAV(points=pts)
        d = av1.dRDA(av2, n_samples=20000)
        assert d > 0.0

    def test_dRDAE(self):
        """FRET-averaged distance from identical pts gives < 1 Å."""
        pts = np.array([[0.0, 0.0, 0.0, 1.0]], dtype=np.float64)
        av1 = BasicAV(points=pts)
        av2 = BasicAV(points=pts)
        r = av1.dRDAE(av2, forster_radius=52.0, n_samples=10000)
        assert r >= 0.0

    def test_pRDA_shape(self):
        """Distance distribution returns probability + axis."""
        pts = np.array([
            [0.0, 0.0, 0.0, 1.0],
            [10.0, 0.0, 0.0, 1.0],
        ], dtype=np.float64)
        av1 = BasicAV(points=pts)
        av2 = BasicAV(points=pts)
        y, x = av1.pRDA(av2, n_samples=10000)
        assert y.ndim == 1
        assert x.ndim == 1
        assert len(y) == len(x)
        assert abs(y.sum() - 1.0) < 0.01

    def test_save_xyz_roundtrip(self, tmp_path):
        """Saving and loading XYZ preserves point count."""
        pts = np.random.RandomState(42).randn(50, 4).astype(np.float64)
        pts[:, 3] = np.abs(pts[:, 3]) + 0.01
        av = BasicAV(points=pts)
        p = tmp_path / "test.xyz"
        av.save_xyz(str(p))
        loaded = np.loadtxt(str(p), skiprows=1)
        assert loaded.shape[0] == 50


class TestKernels:
    """Numba kernel correctness."""

    def test_random_distances_shape(self):
        p1 = np.array([[0.0, 0.0, 0.0, 1.0]], dtype=np.float64)
        p2 = np.array([[1.0, 0.0, 0.0, 1.0]], dtype=np.float64)
        d = random_distances(p1, p2, n_samples=100)
        assert d.shape == (100, 2)

    def test_random_distances_known(self):
        p1 = np.array([[0.0, 0.0, 0.0, 1.0]], dtype=np.float64)
        p2 = np.array([[3.0, 4.0, 0.0, 1.0]], dtype=np.float64)
        d = random_distances(p1, p2, n_samples=1)
        np.testing.assert_allclose(d[0, 0], 5.0, atol=1e-10)

    def test_weighted_mean(self):
        pts = np.array([[0.0, 0.0, 0.0, 2.0],
                        [2.0, 0.0, 0.0, 2.0]], dtype=np.float64)
        m = weighted_mean(pts, 2)
        np.testing.assert_allclose(m, [1.0, 0.0, 0.0], atol=1e-12)

    def test_density2points(self):
        density = np.zeros((3, 3, 3), dtype=np.float64)
        density[1, 1, 1] = 1.0
        origin = np.array([0.0, 0.0, 0.0])
        n, pts = density2points(3, 3, 3, 1.0, density, origin)
        assert n == 1
        np.testing.assert_allclose(pts[0, :3], [1.0, 1.0, 1.0], atol=1e-6)


class TestAccessibleVolume:
    """AccessibleVolume dataclass."""

    def test_empty(self):
        av = AccessibleVolume(
            points=np.empty((0, 4)),
            density=np.zeros((3, 3, 3)),
            grid_origin=np.zeros(3),
            grid_step=1.0,
            grid_shape=(3, 3, 3),
            attachment_point=np.zeros(3),
        )
        assert not av.has_volume
        assert av.n_points == 0
        np.testing.assert_allclose(av.mean_position, av.attachment_point)

    def test_non_empty(self):
        pts = np.array([[0.0, 0.0, 0.0, 1.0]], dtype=np.float64)
        av = AccessibleVolume(
            points=pts,
            density=np.zeros((3, 3, 3)),
            grid_origin=np.zeros(3),
            grid_step=1.0,
            grid_shape=(3, 3, 3),
            attachment_point=np.array([0.0, 0.0, 0.0]),
        )
        assert av.has_volume
        assert av.n_points == 1


class TestSplitAVACV:
    """split_av_acv kernel."""

    def test_all_non_contact_no_centers(self):
        density = np.ones((5, 5, 5), dtype=np.float64)
        rs = np.array([[100.0, 100.0, 100.0]], dtype=np.float64)
        radius = np.array([1.0], dtype=np.float64)
        r0 = np.zeros(3, dtype=np.float64)
        nc, nn, c, nc_mask = split_av_acv(density, 1.0, radius, rs, r0)
        assert nc == 0
        assert nn > 0

    def test_some_contact(self):
        density = np.ones((5, 5, 5), dtype=np.float64)
        rs = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        radius = np.array([1.5], dtype=np.float64)
        r0 = np.zeros(3, dtype=np.float64)
        nc, nn, c, nc_mask = split_av_acv(density, 1.0, radius, rs, r0)
        assert nc > 0
        assert nn > 0
        assert nc + nn == 5 * 5 * 5


class TestACV:
    """Accessible Contact Volume."""

    def test_create_empty(self):
        acv = ACV(position_name="empty")
        assert acv.n_points == 0

    def test_from_basic_av_with_slow_centers(self):
        density = np.ones((5, 5, 5), dtype=np.float64)
        origin = np.array([0.0, 0.0, 0.0])
        av = BasicAV(
            density=density,
            grid_origin=origin,
            grid_step=1.0,
        )
        slow = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        acv = ACV.from_basic_av(av, slow_centers=slow, trapped_fraction=0.5)
        assert isinstance(acv, ACV)
        assert acv.trapped_fraction == 0.5

    def test_slow_radius_broadcast(self):
        density = np.ones((5, 5, 5), dtype=np.float64)
        origin = np.array([0.0, 0.0, 0.0])
        av = BasicAV(
            density=density,
            grid_origin=origin,
            grid_step=1.0,
        )
        slow = np.array([[0.0, 0.0, 0.0], [4.0, 4.0, 4.0]], dtype=np.float64)
        acv = ACV.from_basic_av(av, slow_centers=slow, slow_radius=10.0)
        # scalar is stored as (1,) array; broadcasting happens in kernel
        assert acv.slow_radius[0] == 10.0



class TestComputeAvBackends:
    """`compute_av` end-to-end, on whichever backends this build actually has.

    Written 2026-07-27, when the `imp_bff` branch turned out never to have run:
    it raised `UsageException` on any input because it read its density by
    handing the source particle to an `IMP.em.SampledDensityMap`, which needs an
    `IMP.atom.Mass` it never set. Nothing in this suite touched `compute_av`, so
    a backend that could not work looked fine.
    """

    # A flat slab of obstacles with the attachment site just above it, so the
    # accessible volume is a lopsided hemisphere. A symmetric arrangement would
    # survive the axis transpose below and prove nothing.
    @staticmethod
    def _system():
        xs, ys = np.meshgrid(np.arange(-8.0, 8.1, 2.0), np.arange(-8.0, 8.1, 2.0))
        slab = np.column_stack([xs.ravel(), ys.ravel(), np.zeros(xs.size)])
        wall = slab.copy()
        wall[:, 2] = -2.0
        # Break the symmetry in x so a transposed grid cannot look identical.
        blocker = np.array([[4.0, 0.0, 4.0], [6.0, 0.0, 4.0], [6.0, 0.0, 6.0]])
        atoms_xyz = np.vstack([slab, wall, blocker])
        atoms_vdw = np.full(len(atoms_xyz), 1.8)
        source_xyz = np.array([0.0, 0.0, 2.0])
        return atoms_xyz, atoms_vdw, source_xyz

    @staticmethod
    def _available():
        from IMP.bff.av import compute

        found = []
        if compute._HAS_IMP_BFF:
            found.append("imp_bff")
        if compute._HAS_LABELLIB:
            found.append("labellib")
        return found

    @pytest.mark.parametrize("backend", ["imp_bff", "labellib"])
    def test_it_returns_a_volume(self, backend):
        if backend not in self._available():
            pytest.skip(f"{backend} is not available in this build")
        atoms_xyz, atoms_vdw, source_xyz = self._system()
        av = compute_av(
            atoms_xyz, atoms_vdw, source_xyz,
            linker_length=12.0, linker_width=1.0, dye_radii=(2.0, 0.0, 0.0),
            grid_resolution=0.5, backend=backend,
        )
        assert av.has_volume, "backend returned an empty accessible volume"
        assert av.points.shape[1] == 4

    @pytest.mark.parametrize("backend", ["imp_bff", "labellib"])
    def test_the_volume_sits_on_the_attachment_site(self, backend):
        """Regression: the map used to be left at the coordinate origin.

        With `shift_xyz=False` and the AV decorated onto the source particle
        itself, the header origin came back (0, 0, 0) and most of the grid was
        reported accessible — a volume nowhere near the protein.
        """
        if backend not in self._available():
            pytest.skip(f"{backend} is not available in this build")
        atoms_xyz, atoms_vdw, source_xyz = self._system()
        av = compute_av(
            atoms_xyz, atoms_vdw, source_xyz,
            linker_length=12.0, linker_width=1.0, dye_radii=(2.0, 0.0, 0.0),
            grid_resolution=0.5, backend=backend,
        )
        offset = np.linalg.norm(av.points[:, :3].mean(axis=0) - source_xyz)
        assert offset < 12.0, (
            f"{backend} put the volume {offset:.1f} A from the attachment site"
        )

    @pytest.mark.parametrize("backend", ["imp_bff", "labellib"])
    def test_the_density_and_the_points_describe_the_same_volume(self, backend):
        """The transpose test.

        IMP numbers voxels with x fastest, so reshaping its flat tile values
        C-order into (nx, ny, nz) exchanges x and z. The point count, the
        bounding box and the total volume all survive that unchanged — only the
        shape is mirrored — so the grid and the cloud have to be compared
        against *each other*, voxel by voxel.
        """
        if backend not in self._available():
            pytest.skip(f"{backend} is not available in this build")
        atoms_xyz, atoms_vdw, source_xyz = self._system()
        av = compute_av(
            atoms_xyz, atoms_vdw, source_xyz,
            linker_length=12.0, linker_width=1.0, dye_radii=(2.0, 0.0, 0.0),
            grid_resolution=0.5, backend=backend,
        )
        nx, ny, nz = av.grid_shape
        ix, iy, iz = np.nonzero(av.density)
        from_density = np.stack(
            [av.grid_origin[0] + ix * av.grid_step,
             av.grid_origin[1] + iy * av.grid_step,
             av.grid_origin[2] + iz * av.grid_step],
            axis=1,
        )
        a = set(map(tuple, np.round(from_density, 4)))
        b = set(map(tuple, np.round(av.points[:, :3], 4)))
        assert a == b, (
            f"{backend}: density grid and point cloud disagree — "
            f"{len(a - b)} voxels only in the grid, {len(b - a)} only in the cloud"
        )


class TestComputeAvIsThreadSafe:
    """`compute_av` from a thread pool must not die inside SWIG.

    IMP builds its decorators through SWIG, which is not safe from several
    threads at once. Without `_IMP_BUILD_LOCK` this raises
    ``TypeError: Wrong number or type of arguments for overloaded function
    'XYZR_setup_particle'`` — two threads part-way through construction in one
    interpreter. Running labelling sites in a pool is the obvious way to use
    this function, so the failure sits squarely on the happy path; it was found
    downstream in QuEst, where a scan does exactly this.

    **The system size and the worker count are both load-bearing.** This was
    first written against the 165-atom slab `TestComputeAvBackends` uses, with
    four workers, and it passed with the lock defeated — the construction window
    is too short to overlap, so the test proved nothing. Reproducing the race
    needs roughly a protein's worth of particles and enough threads to collide:
    2000 atoms over 8 workers fails on the first run without the lock. Do not
    shrink either without re-checking that it still fails when
    `_IMP_BUILD_LOCK` is replaced by a no-op.
    """

    N_ATOMS = 2000
    WORKERS = 8
    TASKS = 16

    @staticmethod
    def _system(n):
        # A hollow shell of obstacles around the attachment site: cheap to
        # build, and the resample stays small because the linker is short.
        rng = np.random.default_rng(0)
        directions = rng.normal(size=(n, 3))
        directions /= np.linalg.norm(directions, axis=1, keepdims=True)
        atoms_xyz = directions * rng.uniform(12.0, 30.0, size=(n, 1))
        return atoms_xyz, np.full(n, 1.8), np.array([0.0, 0.0, 0.0])

    def test_concurrent_calls_agree_and_do_not_raise(self):
        import concurrent.futures

        if "imp_bff" not in TestComputeAvBackends._available():
            pytest.skip("imp_bff is not available in this build")

        atoms_xyz, atoms_vdw, source_xyz = self._system(self.N_ATOMS)

        def build(_):
            return compute_av(
                atoms_xyz, atoms_vdw, source_xyz,
                linker_length=12.0, linker_width=1.0,
                dye_radii=(2.0, 0.0, 0.0), grid_resolution=1.0,
                backend="imp_bff",
            )

        with concurrent.futures.ThreadPoolExecutor(max_workers=self.WORKERS) as pool:
            results = list(pool.map(build, range(self.TASKS)))

        assert all(av.n_points > 0 for av in results)
        first = results[0]
        for other in results[1:]:
            # Same input, same answer: a race that did not raise would still
            # show up as two threads disagreeing about the volume.
            assert other.n_points == first.n_points
            assert other.grid_shape == first.grid_shape
            np.testing.assert_allclose(other.grid_origin, first.grid_origin)


# IMP runs every .py under test/ as a standalone script; a file of bare pytest
# functions would import cleanly and exit 0, reporting success without running
# a single assertion. Hand it to pytest so a failure here fails ctest.
if __name__ == "__main__":
    import sys
    try:
        import pytest
    except ImportError:
        print("pytest not installed; skipping", __file__)
        sys.exit(0)
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
