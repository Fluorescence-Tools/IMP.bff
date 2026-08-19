"""The C++ kernels behind ``asa`` and ``pet``, and the coverage they lacked.

Two of these functions had **no test at all**, which is how a missing
``import IMP`` in ``pet.py`` survived a green 139-test run: nothing called
``quenching_rate_per_frame``. The same gap let ``worm_like_chain_linker`` raise
for a whole commit. A kernel with no caller in the suite is a kernel the suite
cannot vouch for, whatever the total says.
"""

import numpy as np
import pytest

import IMP.bff
from IMP.bff.quenching import _asa, sphere_points
from IMP.bff.quenching import quenching_rate_per_frame


class TestSpherePoints:

    @pytest.mark.parametrize("n", [1, 7, 200, 961])
    def test_the_points_are_on_the_unit_sphere(self, n):
        pts = sphere_points(n)
        assert pts.shape == (n, 3)
        np.testing.assert_allclose(np.linalg.norm(pts, axis=1), 1.0, atol=1e-12)

    def test_it_is_float64_now(self):
        """The Python it replaced built the spiral in float32.

        That cost up to 2.4e-6 against an exact float64 construction; the C++
        agrees with one to 1e-15. The port is the more accurate of the two, so
        this is not a tolerance to relax later.
        """
        n = 961
        pts = sphere_points(n)
        assert pts.dtype == np.float64
        inc = np.pi * (3.0 - np.sqrt(5.0))
        off = 2.0 / float(n)
        k = np.arange(n, dtype=np.float64)
        y = k * off - 1.0 + off / 2.0
        rd = np.sqrt(np.maximum(0.0, 1.0 - y * y))
        ref = np.stack([np.cos(k * inc) * rd, y, np.sin(k * inc) * rd], axis=1)
        np.testing.assert_allclose(pts, ref, atol=1e-14)

    def test_they_spread_over_the_whole_sphere(self):
        pts = sphere_points(400)
        assert np.linalg.norm(pts.mean(axis=0)) < 0.05   # no hemisphere bias


class TestSolventAccessibleSurface:

    def test_a_lone_atom_is_fully_exposed(self):
        """With nothing to occlude it, the area is the full sphere at `radius`."""
        xyz = np.zeros((1, 3))
        vdw = np.array([2.0])
        pts = sphere_points(2000)
        area = _asa(xyz, vdw, np.array([0], dtype=np.uint32), pts, 1.4, 2.5)
        assert area[0] == pytest.approx(4.0 * np.pi * vdw[0] ** 2, rel=1e-12)

    def test_a_buried_atom_has_no_area(self):
        """Shells of neighbours at contact leave no accessible sample."""
        shell = sphere_points(300) * 3.0
        xyz = np.vstack([np.zeros(3), shell])
        vdw = np.full(len(xyz), 2.0)
        area = _asa(xyz, vdw, np.array([0], dtype=np.uint32),
                    sphere_points(500), 1.4, 2.5)
        assert area[0] == 0.0

    def test_the_neighbour_cutoff_reaches_far_enough(self):
        """The bug this kernel was ported with: a 2.45 A cutoff, not 6.0 A.

        An occluder at 5 A is inside ``2*radius + probe = 6.4`` and must be
        seen. A cutoff of ``sqrt(2*(vdw + probe)) = 2.45`` would miss it and
        report the atom fully exposed.
        """
        xyz = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]])
        vdw = np.array([2.0, 2.0])
        area = _asa(xyz, vdw, np.array([0], dtype=np.uint32),
                    sphere_points(2000), 1.4, 2.5)
        full = 4.0 * np.pi * vdw[0] ** 2
        assert area[0] < full, "the occluder at 5 A was not seen"


class TestQuenchingRatePerFrame:
    """The function whose missing ``import IMP`` a green suite did not catch."""

    def test_rates_add_over_the_quenchers_in_contact(self):
        collided = np.array([[1, 0, 1], [0, 0, 0], [1, 1, 1]], dtype=bool)
        k = np.array([1.0, 2.0, 3.0])
        np.testing.assert_allclose(quenching_rate_per_frame(collided, k),
                                   [4.0, 0.0, 6.0])

    def test_no_contact_is_zero_not_the_radiative_floor(self):
        collided = np.zeros((5, 4), dtype=bool)
        out = quenching_rate_per_frame(collided, np.arange(4.0))
        np.testing.assert_allclose(out, np.zeros(5))

    def test_it_returns_one_value_per_frame(self):
        collided = np.ones((17, 3), dtype=bool)
        assert quenching_rate_per_frame(collided, np.ones(3)).shape == (17,)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
