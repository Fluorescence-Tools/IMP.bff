"""The density and the point cloud must describe the same volume.

``IMP.bff.compute_av`` returns both a voxel ``density`` and a ``points`` cloud.
They are two readings of one accessible volume and they have to agree — but a
mirrored grid keeps the right voxel count, the right bounding box and the right
total volume, so nothing but a voxel-by-voxel comparison catches it.

Nothing compared them here until PRD-113 stage 3, and they did not agree. IMP
numbers voxels with *x* fastest (``i = x + nx*y + nx*ny*z``), so the C-order
reshape into ``(nx, ny, nz)`` — which makes the *last* axis fastest — returned
the volume transposed. Of 2180 cloud points on T4L A132 at 1.5 Å, 1548 (71 %)
landed on a voxel the density called occupied. Downstream, on the same site at
2.0 Å, the mean donor lifetime came out 3.3974 ns instead of 3.5626 ns.

The array path in ``IMP.bff.compute_av`` had found and fixed this already, and
its test is the model for this one: it reads the point cloud from IMP rather
than deriving it from the grid, precisely so the two are independent readings.
Two readings can disagree; one cannot.
"""

import numpy as np
import pytest

import IMP.bff
from IMP.bff import compute_av_from_structure

_SITE = dict(chain_identifier="A", residue_seq_number=132, atom_name="CB",
             simulation_type="AV1", linker_length=20.0, linker_width=0.5,
             radius1=3.5, allowed_sphere_radius=2.1)


@pytest.fixture(scope="module")
def pdb_path():
    return IMP.bff.get_example_path("structure/T4L/3GUN.pdb")


def _build(pdb_path, resolution):
    return compute_av_from_structure(pdb_path, dict(_SITE), disc_step=resolution)


def _voxel_indices(av):
    points = np.asarray(av.get_points(), dtype=float).reshape(-1, 4)[:, :3]
    origin = np.asarray(av.get_grid_origin(), dtype=float)
    return np.rint((points - origin) / float(av.get_grid_step())).astype(int)


def _density(av):
    ng = av.get_ng()
    return np.asarray(av.get_density(), dtype=float).reshape(ng, ng, ng)


@pytest.mark.parametrize("resolution", [1.5, 2.0, 2.5])
def test_every_cloud_point_indexes_to_an_occupied_voxel(pdb_path, resolution):
    av = _build(pdb_path, resolution)
    density = _density(av)
    idx = _voxel_indices(av)
    inside = np.all((idx >= 0) & (idx < np.array(density.shape)), axis=1)
    assert inside.all(), f"{(~inside).sum()} points index outside the grid"
    occupied = density[idx[:, 0], idx[:, 1], idx[:, 2]] > 0
    assert occupied.all(), (
        f"{(~occupied).sum()} of {len(idx)} cloud points land on an empty "
        "voxel: the density and the point cloud are misregistered")


def test_a_transposed_density_would_be_caught(pdb_path):
    """The guard on the guard.

    A mirrored volume overlaps the truth substantially — 71 % when this was
    live — so the test above only means something if a transpose actually
    fails it. Asserting that keeps it from silently becoming vacuous.
    """
    av = _build(pdb_path, 1.5)
    density = _density(av)
    idx = _voxel_indices(av)
    mirrored = np.ascontiguousarray(density.transpose(2, 1, 0))
    occupied = mirrored[idx[:, 0], idx[:, 1], idx[:, 2]] > 0
    assert not occupied.all(), (
        "a transposed density satisfies the registration test, so the test "
        "cannot detect the defect it exists for")


def test_the_occupied_voxel_count_matches_the_cloud(pdb_path):
    av = _build(pdb_path, 2.0)
    density = _density(av)
    assert int((density > 0).sum()) == av.get_n_points()


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
