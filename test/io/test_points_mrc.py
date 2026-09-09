"""`write_points_mrc`: a weighted point cloud voxelised and written as MRC.

The voxelisation used to live in chisurf, which meant an application knew how
a cloud becomes a grid -- the rounding convention, the origin, the extent, the
suffix rule. All four are format questions, so they are here, and chisurf just
names the file.

No IMP: this is core, and it runs in the IMP-free lane.
"""

import numpy as np
import pytest

import IMP.bff as bff

mrcfile = pytest.importorskip("mrcfile")

#: Three weighted points on one plane, so the grid is 2 x 2 x 1 and every
#: value can be named. Deliberately not symmetric -- see the ordering test.
POINTS = np.asarray(
    [
        [0.0, 0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0, 2.0],
        [1.0, 1.0, 0.0, 3.0],
    ],
    dtype=np.float64,
)


def test_a_cloud_becomes_a_map(tmp_path):
    out = bff.write_points_mrc(str(tmp_path / "av.mrc"), POINTS, 1.0)
    with mrcfile.open(out) as m:
        assert m.data.size == 4
        assert float(m.data.max()) == pytest.approx(3.0)
        assert int(m.header.nversion) == 20140


def test_the_weights_land_on_the_voxels_they_belong_to(tmp_path):
    """MRC runs x fastest and a C-order array runs z fastest.

    A writer that skips the transposition still produces a valid file whose
    map is silently transposed, which for a symmetric cloud looks entirely
    correct. This cloud is not symmetric.
    """
    out = bff.write_points_mrc(str(tmp_path / "av.mrc"), POINTS, 1.0)
    with mrcfile.open(out) as m:
        # mrcfile presents the map as (nz, ny, nx); the cloud is one plane at
        # z = 0, so index [0, y, x].
        data = np.asarray(m.data)
        assert data[0, 0, 0] == pytest.approx(1.0)
        assert data[0, 0, 1] == pytest.approx(2.0)
        assert data[0, 1, 1] == pytest.approx(3.0)
        assert data[0, 1, 0] == pytest.approx(0.0)


def test_points_without_weights_count_once(tmp_path):
    xyz = np.ascontiguousarray(POINTS[:, :3])
    out = bff.write_points_mrc(str(tmp_path / "plain.mrc"), xyz, 1.0)
    with mrcfile.open(out) as m:
        assert float(m.data.sum()) == pytest.approx(3.0)


def test_coincident_points_add_up(tmp_path):
    """Two points in one voxel are one voxel with the sum, not the last one."""
    pts = np.asarray([[0.0, 0.0, 0.0, 1.0], [0.1, 0.0, 0.0, 2.0]])
    out = bff.write_points_mrc(str(tmp_path / "same.mrc"), pts, 1.0)
    with mrcfile.open(out) as m:
        assert float(m.data.sum()) == pytest.approx(3.0)


def test_the_origin_is_the_clouds_own_corner(tmp_path):
    shifted = POINTS.copy()
    shifted[:, :3] += 100.0
    out = bff.write_points_mrc(str(tmp_path / "far.mrc"), shifted, 1.0)
    with mrcfile.open(out) as m:
        assert float(m.header.origin.x) == pytest.approx(100.0)
        # Shifting the cloud may not change its shape.
        assert m.data.shape == (1, 2, 2)


@pytest.mark.parametrize(
    "given, expected",
    [("av", "av.mrc"), ("av.mrc", "av.mrc"), ("av.map", "av.map"),
     ("av.ccp4", "av.ccp4"), ("av.txt", "av.mrc"), ("av.MRC", "av.MRC")],
)
def test_the_suffix_rule(tmp_path, given, expected):
    out = bff.write_points_mrc(str(tmp_path / given), POINTS, 1.0)
    assert out.endswith(expected), out


def test_a_dotted_directory_does_not_lose_the_name(tmp_path):
    """`find_last_of('.')` finds the one in `v1.2/av`, not a suffix."""
    d = tmp_path / "v1.2"
    d.mkdir()
    out = bff.write_points_mrc(str(d / "av"), POINTS, 1.0)
    assert out.endswith("av.mrc"), out


@pytest.mark.parametrize("bad", [np.zeros((0, 4)), np.zeros((3, 2))])
def test_a_cloud_this_cannot_read_is_refused(tmp_path, bad):
    with pytest.raises(Exception):
        bff.write_points_mrc(str(tmp_path / "x.mrc"), bad, 1.0)


@pytest.mark.parametrize("step", [0.0, -1.0, float("nan")])
def test_a_grid_step_that_is_not_a_length_is_refused(tmp_path, step):
    with pytest.raises(Exception):
        bff.write_points_mrc(str(tmp_path / "x.mrc"), POINTS, step)
