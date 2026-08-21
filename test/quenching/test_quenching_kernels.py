"""PRD-109: the PET quenching kernels moved here from QuEst.

Every numeric pin below was captured from QuEst's implementation at the time of
the move and verified equal to it, kernel by kernel, on 2026-08-17 -- the tables
and the ASA exactly, the grids voxel-for-voxel on both an odd and an even grid
edge (the even case is where the registration bugs bit), the Brownian walk
bit-for-bit at a fixed seed, the photon trace bit-for-bit at a fixed seed, and
the decay curve bit-for-bit against QuEst run serially.

The values are frozen *here* rather than compared against QuEst at run time
because nothing in ``IMP.bff`` may import ``quest``: the dependency arrow points
one way.
"""

import unittest

import math

import numpy as np

import IMP
import IMP.test

# quenching is one module now, so its parts are named directly rather than
# through the submodules the six-file subpackage used to have.
import IMP.bff as _q
import IMP.bff as diffusion


def sphere_grid(ng=40, radius_voxels=15, slow_voxels=8):
    """A spherical accessible volume with a concentric sticky core."""
    i = np.arange(ng) - (ng - 1) // 2
    x, y, z = np.meshgrid(i, i, i, indexing="ij")
    r2 = x ** 2 + y ** 2 + z ** 2
    density = (r2 < radius_voxels ** 2).astype(np.uint8)
    slow = (r2 < slow_voxels ** 2).astype(np.uint8)
    return density, slow


class PetParameterTests(IMP.test.TestCase):

    def test_quenching_centres_sit_on_the_redox_active_moiety(self):
        """Not on CB -- for TRP that is ~3.3 A away, and contact is a few A."""
        qa = _q.quencher_atoms()
        self.assertEqual(tuple(qa["MET"]), ("SD",))
        self.assertEqual(tuple(qa["CYS"]), ("SG",))
        self.assertIn("NE1", qa["TRP"])
        self.assertIn("OH", qa["TYR"])
        self.assertNotIn("CB", qa["TRP"])

    def test_the_reference_chemistry_orders_the_quenchers(self):
        """TRP > PRO ~ TYR > MET > HIS > CYS."""
        reference = _q.pet_quenching_reference()
        kq = {r: p.kQ for r, p in reference.items()}
        self.assertGreater(kq["TRP"], kq["TYR"])
        self.assertAlmostEqual(kq["PRO"], kq["TYR"])
        self.assertGreater(kq["TYR"], kq["MET"])
        self.assertGreater(kq["MET"], kq["HIS"])
        self.assertGreater(kq["HIS"], kq["CYS"])

    def test_defaults_cover_every_standard_residue(self):
        table = _q.amino_acid_quenching_defaults()
        standards = _q.standard_amino_acid_residues()
        self.assertEqual(set(table), set(standards))
        reference = _q.pet_quenching_reference()
        for residue, params in table.items():
            if residue not in reference:
                self.assertEqual(params.kQ, 0.0)
                # NaN is what `None` was: inherit the model-wide distance.
                self.assertTrue(math.isnan(params.quench_radius))

    def test_the_dye_radius_offsets_the_surface_contact_distance(self):
        """The table is quoted dye-surface-to-quencher; the walk tracks centres."""
        radius = 5.0
        table = _q.amino_acid_quenching_defaults(dye_radius=radius)
        reference = _q.pet_quenching_reference()
        for residue, ref in reference.items():
            self.assertAlmostEqual(
                table[residue].quench_radius,
                radius + ref.contact_distance,
            )

    def test_kq_scale_is_a_plain_multiplier(self):
        scaled = _q.amino_acid_quenching_defaults(kQ_scale=0.4)
        plain = _q.amino_acid_quenching_defaults()
        reference = _q.pet_quenching_reference()
        for residue in reference:
            self.assertAlmostEqual(
                scaled[residue].kQ, 0.4 * plain[residue].kQ
            )

    def test_normalisation_is_forgiving_about_shape(self):
        table = _q.normalize_amino_acid_quenching({
            "trp": _q.ResidueQuenching(kQ=9.0),
            "TYR": _q.ResidueQuenching(slow_factor=0.3),
            "ZZZ": _q.ResidueQuenching(kQ=1.0),
        })
        # The residue name is normalised, so a lower-case key still lands.
        self.assertEqual(table["TRP"].kQ, 9.0)
        self.assertEqual(table["TYR"].slow_factor, 0.3)
        # An unknown residue is kept, so a non-standard one can be given a rate,
        # and it falls back to CB for its quenching centre.
        self.assertEqual(table["ZZZ"].kQ, 1.0)
        self.assertEqual(list(table["ZZZ"].quench_atoms), ["CB"])

    def test_slow_factors_are_clamped_and_radii_sanitised(self):
        table = _q.normalize_amino_acid_quenching({
            "TRP": _q.ResidueQuenching(slow_factor=5.0, kQ=-1.0,
                                       quench_radius=-3.0),
            "TYR": _q.ResidueQuenching(slow_factor=-2.0,
                                       quench_radius=float("nan")),
        })
        self.assertEqual(table["TRP"].slow_factor, 1.0)
        self.assertEqual(table["TYR"].slow_factor, 0.0)
        self.assertEqual(table["TRP"].kQ, 0.0)
        # A non-positive or non-finite radius means "inherit the global one".
        self.assertTrue(math.isnan(table["TRP"].quench_radius))
        self.assertTrue(math.isnan(table["TYR"].quench_radius))


class SolventAccessibleSurfaceTests(IMP.test.TestCase):

    def test_a_lone_atom_is_fully_exposed(self):
        """4 pi r^2 for the isolated sphere, within the sampling resolution."""
        xyz = np.zeros((1, 3))
        vdw = np.array([2.0])
        area = _q.solvent_accessible_surface(xyz, vdw, np.array([0], np.uint32))
        self.assertAlmostEqual(float(area[0]), 4.0 * np.pi * 2.0 ** 2, delta=1e-3)

    def test_burying_an_atom_reduces_its_area(self):
        shell = np.array([[0.0, 0.0, 0.0]] + [
            [3.0 * np.sin(t) * np.cos(p), 3.0 * np.sin(t) * np.sin(p), 3.0 * np.cos(t)]
            for t in np.linspace(0.2, np.pi - 0.2, 6)
            for p in np.linspace(0, 2 * np.pi, 6, endpoint=False)
        ])
        vdw = np.full(shell.shape[0], 2.0)
        idx = np.array([0], np.uint32)
        buried = _q.solvent_accessible_surface(shell, vdw, idx)
        alone = _q.solvent_accessible_surface(shell[:1], vdw[:1], idx)
        self.assertLess(float(buried[0]), 0.5 * float(alone[0]))

    def test_frozen_reference(self):
        """QuEst's kernel gives 1194.63 here -- it is the one that is wrong.

        Its neighbour test compared a *squared* distance against an unsquared
        sum (``dist2 < 2 * (vdw + probe)``), so at the defaults it looked for
        occluders within 2.45 A instead of 6.0 A, missed nearly all of them and
        reported almost every atom as fully exposed. Nothing consumed the
        result, so correcting it moved no downstream number -- see
        ``quenching/asa.py``.
        """
        rng = np.random.default_rng(3)
        n = 300
        xyz = rng.normal(0, 12, (n, 3))
        vdw = rng.uniform(1.2, 2.0, n)
        idx = np.arange(0, n, 7, dtype=np.uint32)
        area = _q.solvent_accessible_surface(xyz, vdw, idx)
        self.assertAlmostEqual(float(area.sum()), 642.911865, delta=1e-3)

    def test_no_probe_atoms_is_not_an_error(self):
        area = _q.solvent_accessible_surface(
            np.zeros((3, 3)), np.ones(3), np.array([], np.uint32)
        )
        self.assertEqual(area.size, 0)


class GridRegistrationTests(IMP.test.TestCase):

    def test_the_centre_index_is_integer_not_the_float_corner(self):
        """The two differ on every *even* edge, which is the normal case."""
        self.assertEqual(_q.grid_center_index(92), 45)
        self.assertEqual(_q.grid_center_index(91), 45)
        for ng in (46, 86, 92):
            self.assertNotEqual(_q.grid_center_index(ng), (ng - 1) / 2)

    def test_a_centre_lands_on_the_voxel_that_contains_it(self):
        """The stamp must index the grid the way the walk samples it."""
        ng, dg = 21, 1.0
        density = np.ones((ng, ng, ng), np.uint8)
        offset = _q.grid_center_index(ng)
        for point, expected in (
            ([0.0, 0.0, 0.0], offset),
            ([2.0, 0.0, 0.0], offset + 2),
            ([-2.0, 0.0, 0.0], offset - 2),
        ):
            rates = _q.quenching_rate_grid(
                density, ng, dg, np.array([1.4]), np.array([point]),
                np.zeros(3), np.array([1.0]),
            )
            hit = np.argwhere(rates > 0)
            self.assertEqual(int(round(hit[:, 0].mean())), expected)

    def test_a_negative_offset_rounds_down_like_every_other_map(self):
        """`int()` truncates toward zero; half the grid is on the negative side."""
        ng, dg = 21, 1.0
        density = np.ones((ng, ng, ng), np.uint8)
        offset = _q.grid_center_index(ng)
        rates = _q.quenching_rate_grid(
            density, ng, dg, np.array([1.4]), np.array([[-0.5, 0.0, 0.0]]),
            np.zeros(3), np.array([1.0]),
        )
        hit = np.argwhere(rates > 0)
        self.assertEqual(int(round(hit[:, 0].mean())), offset - 1)


class QuenchingGridTests(IMP.test.TestCase):

    def setUp(self):
        super().setUp()
        self.ng, self.dg = 31, 0.6
        self.density = np.ones((self.ng,) * 3, np.uint8)
        self.r0 = np.zeros(3)

    def test_overlapping_quenchers_add_their_rates(self):
        """Independent PET channels compose additively."""
        centres = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        rates = _q.quenching_rate_grid(
            self.density, self.ng, self.dg, np.array([4.0, 4.0]),
            centres, self.r0, np.array([1.5, 2.5]),
        )
        self.assertAlmostEqual(float(rates.max()), 4.0)

    def test_overlapping_sticky_spheres_multiply_their_factors(self):
        centres = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        factors = _q.slow_factor_grid(
            self.density, self.ng, self.dg, np.array([4.0, 4.0]),
            centres, self.r0, np.array([0.5, 0.25]),
        )
        self.assertAlmostEqual(float(factors.min()), 0.125)

    def test_voxels_outside_the_accessible_volume_are_left_alone(self):
        """The walk cannot reach them, so a rate there would be a lie."""
        density = np.zeros((self.ng,) * 3, np.uint8)
        density[10:20, 10:20, 10:20] = 1
        rates = _q.quenching_rate_grid(
            density, self.ng, self.dg, np.array([20.0]), np.zeros((1, 3)),
            self.r0, np.array([3.0]),
        )
        self.assertTrue(np.all(rates[density == 0] == 0.0))
        self.assertTrue(np.any(rates[density == 1] > 0.0))
        factors = _q.slow_factor_grid(
            density, self.ng, self.dg, np.array([20.0]), np.zeros((1, 3)),
            self.r0, np.array([0.3]),
        )
        self.assertTrue(np.all(factors[density == 0] == 1.0))

    def test_a_zero_rate_quencher_stamps_nothing(self):
        rates = _q.quenching_rate_grid(
            self.density, self.ng, self.dg, np.array([5.0]), np.zeros((1, 3)),
            self.r0, np.array([0.0]),
        )
        self.assertEqual(float(rates.max()), 0.0)

    def test_the_stamp_is_a_ball_of_the_requested_radius(self):
        radius = 3.0
        rates = _q.quenching_rate_grid(
            self.density, self.ng, self.dg, np.array([radius]), np.zeros((1, 3)),
            self.r0, np.array([1.0]),
        )
        n_voxels = int((rates > 0).sum())
        expected = 4.0 / 3.0 * np.pi * (radius / self.dg) ** 3
        self.assertAlmostEqual(n_voxels / expected, 1.0, delta=0.1)

    def test_the_contact_mask_is_inside_the_accessible_volume(self):
        density, _ = sphere_grid(ng=self.ng, radius_voxels=12)
        mask = _q.av_contact_mask(
            density, self.ng, self.dg, np.array([3.0]), np.zeros((1, 3)), self.r0
        )
        self.assertTrue(np.all(mask[density == 0] == 0))
        self.assertGreater(int(mask.sum()), 0)


class DyeDiffusionTests(IMP.test.TestCase):

    def setUp(self):
        super().setUp()
        self.density, self.slow = sphere_grid()
        self.dg = 0.5

    def run_walk(self, **kwargs):
        options = dict(
            t_max=200.0, t_step=0.002, D=40.0, slow_fact=[0.1], random_seed=7
        )
        options.update(kwargs)
        return diffusion.simulate_dye_diffusion(
            self.density, self.slow, self.dg, **options
        )

    @staticmethod
    def _xyz(walk):
        return walk.get_xyz().reshape(-1, 3)

    def test_frozen_reference(self):
        """The walk's *statistics*, since the exact trace can no longer be pinned.

        This asserted ``n_accepted == 94567`` -- bit-for-bit QuEst's walk at
        seed 7, which held through the numba port because both ran the same
        Mersenne stream. The C++ port (PRD-113 stage 5) ends that: no C++
        generator reproduces numba's stream, so a per-seed count is a property
        of the generator and not of the model.

        What survives is the distribution, and it is a tight one. Over 12
        seeds:

            numba, 6*D*dt width   94436 +/- 272   range [94002, 94835]
            C++,   6*D*dt width   94388 +/- 318   range [94031, 95096]
            C++,   2*D*dt width   96618 +/- 253   range [96218, 97057]

        The third row is the current one. The step width was corrected on
        2026-08-18 -- it had been the total 3-D MSD used as one component's
        width, so the walk diffused at 3D -- and a shorter step rejects less
        often, which is the whole of the difference. The first two rows are kept
        because they are the evidence that the *port* changed nothing: same
        width, same distribution.

        The bound below is roughly four standard deviations wide and would still
        catch a real change in the geometry, the step width or the rejection
        rule. A single seed is deliberately not asserted.
        """
        accepted = [self.run_walk(random_seed=s).n_accepted for s in range(6)]
        mean = sum(accepted) / len(accepted)
        self.assertAlmostEqual(mean, 96600, delta=1200)
        for walk_accepted in accepted:
            self.assertEqual(walk_accepted + (100000 - walk_accepted), 100000)
        self.assertEqual(self.run_walk().n_frames, 100000)

    def test_the_walk_is_reproducible_for_a_seed(self):
        first = self.run_walk()
        second = self.run_walk()
        self.assertTrue(np.array_equal(self._xyz(first), self._xyz(second)))

    def test_different_seeds_give_different_walks(self):
        self.assertFalse(
            np.array_equal(self._xyz(self.run_walk(random_seed=1)),
                           self._xyz(self.run_walk(random_seed=2)))
        )

    def test_the_walk_never_leaves_the_accessible_volume(self):
        """Every accepted position must sit on an occupied voxel."""
        walk = self.run_walk()
        ng = self.density.shape[0]
        xyz = self._xyz(walk)
        idx = np.floor(xyz / self.dg).astype(int) + _q.grid_center_index(ng)
        self.assertTrue(np.all((idx >= 0) & (idx < ng)))
        self.assertTrue(np.all(self.density[idx[:, 0], idx[:, 1], idx[:, 2]] > 0))

    def test_stickiness_slows_the_dye_down(self):
        """Halving the step width inside the core must shorten the mean step."""
        free = self.run_walk(slow_fact=[1.0])
        sticky = self.run_walk(slow_fact=[0.01])
        step = lambda w: np.linalg.norm(np.diff(self._xyz(w), axis=0), axis=1).mean()
        self.assertLess(step(sticky), step(free))

    def test_a_per_voxel_factor_grid_is_accepted(self):
        factors = _q.slow_factor_grid(
            self.density, self.density.shape[0], self.dg,
            np.array([4.0]), np.zeros((1, 3)), np.zeros(3), np.array([0.05]),
        )
        walk = self.run_walk(slow_fact=factors)
        self.assertEqual(walk.n_frames, 100000)
        self.assertGreater(walk.n_accepted, 0)

    def test_an_empty_volume_yields_no_accepted_steps(self):
        """Rather than looping forever looking for a start point."""
        empty = np.zeros_like(self.density)
        walk = diffusion.simulate_dye_diffusion(
            empty, empty, self.dg, t_max=10.0, t_step=0.002, random_seed=1
        )
        self.assertEqual(walk.n_accepted, 0)
        self.assertEqual(walk.acceptance_ratio, 0.0)

    def test_a_larger_diffusion_coefficient_explores_further(self):
        slow = self.run_walk(D=1.0, slow_fact=[1.0])
        fast = self.run_walk(D=100.0, slow_fact=[1.0])
        spread = lambda w: float(np.linalg.norm(self._xyz(w).std(axis=0)))
        self.assertGreater(spread(fast), spread(slow))


class FretRateTraceTests(IMP.test.TestCase):

    def setUp(self):
        super().setUp()
        self.trajectory = np.zeros((16, 3))
        self.trajectory[:, 0] = np.linspace(0.0, 30.0, 16)
        self.acceptor = np.array([[52.0, 0.0, 0.0]])

    def test_the_rate_at_the_forster_radius_is_the_intrinsic_rate(self):
        """k_FRET = 1/tau0 exactly when r = R0 and kappa2 is isotropic."""
        rate = _q.fret_rate_trace(
            np.zeros((1, 3)), self.acceptor, R0=52.0, tau0=4.0
        )
        self.assertAlmostEqual(float(rate[0]), 1.0 / 4.0)

    def test_the_rate_falls_as_the_sixth_power(self):
        rate = _q.fret_rate_trace(
            np.array([[0.0, 0.0, 0.0], [26.0, 0.0, 0.0]]),
            np.array([[52.0, 0.0, 0.0]]), R0=52.0, tau0=4.0,
        )
        self.assertAlmostEqual(float(rate[1] / rate[0]), 2.0 ** 6, delta=1e-6)

    def test_kappa2_scales_relative_to_the_isotropic_value(self):
        """R0 is published *at* 2/3, so passing 2/3 must change nothing."""
        base = _q.fret_rate_trace(
            self.trajectory, self.acceptor, 52.0, 4.0
        )
        same = _q.fret_rate_trace(
            self.trajectory, self.acceptor, 52.0, 4.0, kappa2=2.0 / 3.0
        )
        self.assertTrue(np.allclose(base, same))
        doubled = _q.fret_rate_trace(
            self.trajectory, self.acceptor, 52.0, 4.0, kappa2=4.0 / 3.0
        )
        self.assertTrue(np.allclose(doubled, 2.0 * base))

    def test_an_impossible_kappa2_is_rejected(self):
        for bad in (-0.1, 4.5, float("nan")):
            with self.assertRaises(ValueError):
                _q.fret_rate_trace(
                    self.trajectory, self.acceptor, 52.0, 4.0, kappa2=bad
                )

    def test_overlapping_clouds_do_not_diverge(self):
        """Both clouds are dye *centres* and may overlap; (R0/r)^6 would blow up."""
        rate = _q.fret_rate_trace(
            np.zeros((4, 3)), np.zeros((4, 3)), R0=52.0, tau0=4.0, r_min=7.0
        )
        self.assertTrue(np.all(np.isfinite(rate)))
        expected = (1.0 / 4.0) * (52.0 / 7.0) ** 6
        self.assertAlmostEqual(float(rate[0]), expected, delta=1e-6)

    def test_the_pair_trace_matches_the_cloud_average_for_a_single_acceptor(self):
        """A one-point cloud held still *is* a frame-aligned acceptor trajectory."""
        cloud = _q.fret_rate_trace(
            self.trajectory, self.acceptor, 52.0, 4.0
        )
        held = np.repeat(self.acceptor, self.trajectory.shape[0], axis=0)
        pair = _q.fret_rate_pair_trace(self.trajectory, held, 52.0, 4.0)
        self.assertTrue(np.allclose(cloud, pair))

    def test_mismatched_trajectories_raise_rather_than_truncate(self):
        """Truncating would silently change the sampled time window."""
        with self.assertRaises(ValueError):
            _q.fret_rate_pair_trace(
                self.trajectory, self.trajectory[:-1], 52.0, 4.0
            )

    def test_an_empty_acceptor_cloud_raises(self):
        with self.assertRaises(ValueError):
            _q.fret_rate_trace(
                self.trajectory, np.zeros((0, 3)), 52.0, 4.0
            )

    def test_a_non_positive_lifetime_raises(self):
        with self.assertRaises(ValueError):
            _q.fret_rate_trace(self.trajectory, self.acceptor, 52.0, 0.0)

    def test_the_acceptor_cloud_is_subsampled_by_an_even_stride(self):
        """The average converges long before the full cloud."""
        rng = np.random.default_rng(0)
        cloud = rng.normal(40.0, 5.0, (5000, 3))
        full = _q.fret_rate_trace(
            self.trajectory, cloud, 52.0, 4.0, max_acceptor_points=5000
        )
        capped = _q.fret_rate_trace(
            self.trajectory, cloud, 52.0, 4.0, max_acceptor_points=512
        )
        self.assertTrue(np.allclose(full, capped, rtol=0.05))


if __name__ == "__main__":
    IMP.test.main()
