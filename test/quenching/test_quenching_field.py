"""PRD-109: the field formulation, folded in from ChiSurf.

Rate and mobility maps on the AV grid plus an explicit solver for
``dp/dt = div(D grad p) - k p``, the deterministic counterpart of the Brownian
walk in ``test_quenching_kernels.py``.

Pinned against analytic diffusion rather than against ChiSurf's numbers, because
three of ChiSurf's kernels were wrong (see :class:`PortedDefectTests`).
"""

import numpy as np

import IMP
import IMP.test

from IMP.bff.quenching import maps
from IMP.bff.sampling.smoluchowski import (
    GridDiffusionSolver,
    diffusion_stability_limit,
    equilibrium_occupancy,
)


def open_box(ng=41):
    """A domain filling the grid except the outer shell the stencil cannot reach."""
    bounds = np.zeros((ng,) * 3)
    bounds[1:-1, 1:-1, 1:-1] = 1.0
    return bounds


def point_source(ng=41):
    p = np.zeros((ng,) * 3)
    c = (ng - 1) // 2
    p[c, c, c] = 1.0
    return p


class GridAxisTests(IMP.test.TestCase):

    def test_the_axis_uses_the_same_centre_as_every_other_map(self):
        from IMP.bff.quenching.grids import grid_center_index
        for ng in (31, 46, 92):
            axis = maps.grid_axis(ng, 0.5)
            self.assertEqual(axis.size, ng)
            self.assertAlmostEqual(axis[grid_center_index(ng)], 0.0)

    def test_the_spacing_is_the_voxel_edge(self):
        axis = maps.grid_axis(21, 0.75)
        self.assertAlmostEqual(float(np.diff(axis).max()), 0.75)


class AtomicQuenchingParameterTests(IMP.test.TestCase):

    def setUp(self):
        super().setUp()
        self.atoms = np.zeros(
            4, dtype=[("res_name", "U4"), ("atom_name", "U4"), ("coord", "f8", 3)]
        )
        self.atoms["res_name"] = ["TRP", "TRP", "ALA", "MET"]
        self.atoms["atom_name"] = ["NE1", "CB", "CB", "SD"]

    def test_only_named_atoms_quench(self):
        kQ, rC = maps.atomic_quenching_parameters(
            self.atoms, {"TRP": {"NE1": (3.5, 1.0)}, "MET": {"SD": (1.67, 0.8)}}
        )
        self.assertTrue(np.allclose(kQ, [3.5, 0.0, 0.0, 1.67]))
        self.assertTrue(np.allclose(rC, [1.0, 0.0, 0.0, 0.8]))

    def test_an_empty_table_quenches_nothing(self):
        kQ, rC = maps.atomic_quenching_parameters(self.atoms, {})
        self.assertEqual(float(kQ.sum()), 0.0)
        self.assertEqual(float(rC.sum()), 0.0)


class RateMapTests(IMP.test.TestCase):

    def setUp(self):
        super().setUp()
        self.ng, self.dg = 21, 1.0
        self.density = np.ones((self.ng,) * 3)
        self.r0 = np.zeros(3)

    def test_the_intrinsic_rate_is_the_floor_without_quenchers(self):
        rates = maps.quenching_rate_map(
            self.density, self.r0, self.dg, np.zeros((0, 3)),
            np.zeros(0), np.zeros(0), tau0=4.0, dye_radius=3.5,
        )
        self.assertTrue(np.allclose(rates, 0.25))

    def test_quenching_falls_off_exponentially(self):
        """kQ * exp(-(d - r_dye) / rC), measured from the dye surface."""
        kQ, rC, r_dye, tau0 = 3.5, 2.0, 3.5, 4.0
        rates = maps.quenching_rate_map(
            self.density, self.r0, self.dg, np.zeros((1, 3)),
            np.array([kQ]), np.array([rC]), tau0=tau0, dye_radius=r_dye,
        )
        centre = (self.ng - 1) // 2
        for offset in (2, 4, 6):
            got = rates[centre + offset, centre, centre]
            expected = 1.0 / tau0 + kQ * np.exp(-(offset - r_dye) / rC)
            self.assertAlmostEqual(float(got), expected, places=9)

    def test_a_zero_rate_or_length_atom_is_skipped(self):
        rates = maps.quenching_rate_map(
            self.density, self.r0, self.dg, np.zeros((2, 3)),
            np.array([0.0, 3.5]), np.array([2.0, 0.0]), tau0=4.0, dye_radius=3.5,
        )
        self.assertTrue(np.allclose(rates, 0.25))

    def test_voxels_outside_the_volume_carry_no_rate(self):
        density = np.zeros((self.ng,) * 3)
        density[5:15, 5:15, 5:15] = 1.0
        rates = maps.quenching_rate_map(
            density, self.r0, self.dg, np.zeros((1, 3)),
            np.array([3.5]), np.array([2.0]), tau0=4.0, dye_radius=3.5,
        )
        self.assertTrue(np.all(rates[density == 0] == 0.0))
        self.assertTrue(np.all(rates[density > 0] > 0.0))

    def test_the_whole_grid_is_written(self):
        """ChiSurf's `range(-npm, npm)` never reached the outer slab."""
        rates = maps.quenching_rate_map(
            self.density, self.r0, self.dg, np.zeros((1, 3)),
            np.array([3.5]), np.array([2.0]), tau0=4.0, dye_radius=3.5,
        )
        self.assertTrue(np.all(rates > 0.0))
        for face in (rates[0], rates[-1], rates[:, 0], rates[:, -1],
                     rates[..., 0], rates[..., -1]):
            self.assertTrue(np.all(face > 0.0))


class MobilityMapTests(IMP.test.TestCase):

    def setUp(self):
        super().setUp()
        self.ng, self.dg = 21, 1.0
        self.density = np.ones((self.ng,) * 3)

    def test_far_from_every_atom_the_dye_is_free(self):
        d_map = maps.diffusion_coefficient_map(
            self.density, np.zeros(3), self.dg, np.array([[100.0, 100.0, 100.0]]),
            free_diffusion=8.0, min_distance=5.0, slow_factor=0.5,
        )
        self.assertTrue(np.allclose(d_map, 8.0))

    def test_each_contacting_atom_slows_the_dye_again(self):
        """The factor compounds -- two atoms in contact means slow_factor^2."""
        atoms = np.zeros((2, 3))
        d_map = maps.diffusion_coefficient_map(
            self.density, np.zeros(3), self.dg, atoms,
            free_diffusion=8.0, min_distance=3.0, slow_factor=0.5,
        )
        centre = (self.ng - 1) // 2
        self.assertAlmostEqual(float(d_map[centre, centre, centre]), 8.0 * 0.25)

    def test_the_map_is_zero_where_the_dye_cannot_go(self):
        """A non-zero coefficient there would let the solver leak into the protein."""
        density = np.zeros((self.ng,) * 3)
        density[5:15, 5:15, 5:15] = 1.0
        d_map = maps.diffusion_coefficient_map(
            density, np.zeros(3), self.dg, np.array([[100.0, 0.0, 0.0]]),
            free_diffusion=8.0, min_distance=3.0, slow_factor=0.5,
        )
        self.assertTrue(np.all(d_map[density == 0] == 0.0))

    def test_a_radial_profile_replaces_the_constant_base(self):
        d_map = maps.radial_diffusion_map(
            self.density, self.dg, lambda r: 1.0 + r
        )
        centre = (self.ng - 1) // 2
        self.assertAlmostEqual(float(d_map[centre, centre, centre]), 1.0)
        self.assertAlmostEqual(float(d_map[centre + 3, centre, centre]), 4.0)

    def test_the_whole_grid_is_written(self):
        """ChiSurf's variant 1 left the last slab as uninitialised memory."""
        d_map = maps.diffusion_coefficient_map(
            self.density, np.zeros(3), self.dg, np.array([[100.0, 100.0, 100.0]]),
            free_diffusion=8.0, min_distance=3.0, slow_factor=0.5,
        )
        self.assertTrue(np.all(d_map == 8.0))


class FretMapTests(IMP.test.TestCase):

    def single_voxel_grids(self, separation, ng=11, dg=1.0):
        donor = np.zeros((ng,) * 3)
        acceptor = np.zeros((ng,) * 3)
        c = (ng - 1) // 2
        donor[c, c, c] = 1.0
        acceptor[c, c, c] = 1.0
        return donor, acceptor, np.zeros(3), np.array([separation, 0.0, 0.0]), dg

    def test_one_donor_and_one_acceptor_give_the_forster_rate(self):
        d, a, r0d, r0a, dg = self.single_voxel_grids(52.0)
        rates = maps.fret_rate_map(d, a, r0d, r0a, dg, dg, 52.0, 0.25,
                                   acceptor_step=1)
        c = (d.shape[0] - 1) // 2
        # r = R0 -> k = kf
        self.assertAlmostEqual(float(rates[c, c, c]), 0.25, places=9)

    def test_the_rate_falls_as_the_sixth_power(self):
        c = None
        rates = []
        for separation in (52.0, 104.0):
            d, a, r0d, r0a, dg = self.single_voxel_grids(separation)
            m = maps.fret_rate_map(d, a, r0d, r0a, dg, dg, 52.0, 0.25,
                                   acceptor_step=1)
            c = (d.shape[0] - 1) // 2
            rates.append(float(m[c, c, c]))
        self.assertAlmostEqual(rates[0] / rates[1], 2.0 ** 6, delta=1e-6)

    def test_it_is_the_harmonic_mean_of_the_rates(self):
        """1/<1/k>, not <k> -- the mean transfer *time* is what is averaged."""
        ng, dg = 11, 1.0
        c = (ng - 1) // 2
        donor = np.zeros((ng,) * 3)
        donor[c, c, c] = 1.0
        acceptor = np.zeros((ng,) * 3)
        acceptor[c - 2, c, c] = 1.0
        acceptor[c + 2, c, c] = 1.0
        r0a = np.array([60.0, 0.0, 0.0])
        kf, R0 = 0.25, 52.0
        got = maps.fret_rate_map(
            donor, acceptor, np.zeros(3), r0a, dg, dg, R0, kf, acceptor_step=1
        )[c, c, c]
        distances = np.array([60.0 - 2.0, 60.0 + 2.0])
        k = kf * (R0 / distances) ** 6
        harmonic = 1.0 / np.mean(1.0 / k)
        self.assertAlmostEqual(float(got), float(harmonic), places=9)
        self.assertNotAlmostEqual(float(got), float(np.mean(k)), places=6)

    def test_a_zero_radiative_rate_raises(self):
        d, a, r0d, r0a, dg = self.single_voxel_grids(52.0)
        with self.assertRaises(ValueError):
            maps.fret_rate_map(d, a, r0d, r0a, dg, dg, 52.0, 0.0)


class GridDiffusionSolverTests(IMP.test.TestCase):

    def setUp(self):
        super().setUp()
        self.ng, self.dg, self.D = 41, 1.0, 1.0
        self.bounds = open_box(self.ng)
        self.d_map = np.full((self.ng,) * 3, self.D)
        self.t_step = 0.5 * diffusion_stability_limit(self.D, self.dg)

    def solver(self, density, rate_map=None):
        return GridDiffusionSolver(
            self.d_map, self.bounds, density, rate_map,
            t_step=self.t_step, dg=self.dg,
        )

    def test_free_diffusion_spreads_as_two_D_t(self):
        """The pin that caught ChiSurf's half-rate ping-pong."""
        n_steps = 200
        result = self.solver(point_source(self.ng)).run(n_steps, n_out=n_steps)
        axis = maps.grid_axis(self.ng, self.dg)
        density = result.density
        variance = float((density * axis[:, None, None] ** 2).sum() / density.sum())
        expected = 2.0 * self.D * n_steps * self.t_step
        self.assertAlmostEqual(variance / expected, 1.0, delta=0.01)

    def test_mass_is_conserved_without_decay(self):
        result = self.solver(point_source(self.ng)).run(200, n_out=50)
        self.assertAlmostEqual(float(result.density.sum()), 1.0, places=9)
        self.assertTrue(np.allclose(result.fluorescence, 1.0, atol=1e-9))

    def test_a_uniform_rate_gives_a_single_exponential(self):
        rate = 0.25
        result = self.solver(
            self.bounds.copy(), np.full((self.ng,) * 3, rate)
        ).run(400, n_out=20)
        slope = np.polyfit(result.time, np.log(result.fluorescence), 1)[0]
        # Explicit Euler is first order in dt, hence the 2% tolerance.
        self.assertAlmostEqual(-slope, rate, delta=0.02 * rate)

    def test_nothing_leaks_outside_the_domain(self):
        result = self.solver(point_source(self.ng)).run(300, n_out=300)
        self.assertTrue(np.all(result.density[self.bounds == 0] == 0.0))

    def test_a_domain_touching_the_grid_edge_raises(self):
        """The stencil cannot be evaluated there, so population would vanish."""
        solver = GridDiffusionSolver(
            self.d_map, np.ones((self.ng,) * 3), np.ones((self.ng,) * 3),
            t_step=self.t_step, dg=self.dg,
        )
        with self.assertRaises(ValueError):
            solver.run(1)

    def test_an_unstable_time_step_raises_rather_than_diverging(self):
        limit = diffusion_stability_limit(self.D, self.dg)
        solver = GridDiffusionSolver(
            self.d_map, self.bounds, point_source(self.ng),
            t_step=2.0 * limit, dg=self.dg,
        )
        with self.assertRaises(ValueError):
            solver.run(1)

    def test_the_stability_limit_is_the_explicit_three_d_one(self):
        self.assertAlmostEqual(diffusion_stability_limit(2.0, 1.0), 1.0 / 12.0)
        self.assertEqual(diffusion_stability_limit(0.0, 1.0), float("inf"))

    def test_equilibrium_is_flat_where_the_mobility_is_flat(self):
        solver = self.solver(point_source(self.ng))
        equilibrium = solver.equilibrium(n_steps=6000, tolerance=1e-10)
        interior = equilibrium[self.bounds > 0]
        self.assertAlmostEqual(float(equilibrium.sum()), 1.0, places=9)
        self.assertLess(float(interior.std() / interior.mean()), 0.05)

    def test_ito_equilibrium_favours_the_slow_region(self):
        """The inherited behaviour, kept as the ``"ito"`` branch and pinned there.

        This asserted the *default* until 2026-08-18, when Smoluchowski became
        the default: a mobility field must not decide where the dye sits (see
        :class:`FluxFormTests`). Occupancy differing from the AV density is a
        real effect, but it needs an attractive potential to produce it, not
        friction.
        """
        d_map = np.full((self.ng,) * 3, self.D)
        half = self.ng // 2
        d_map[:half] = self.D / 20.0
        solver = GridDiffusionSolver(
            d_map, self.bounds, self.bounds.copy(),
            t_step=0.5 * diffusion_stability_limit(float(d_map.max()), self.dg),
            dg=self.dg, flux_form="ito",
        )
        equilibrium = solver.equilibrium(n_steps=20000, tolerance=1e-12)
        slow = float(equilibrium[:half][self.bounds[:half] > 0].mean())
        fast = float(equilibrium[half:][self.bounds[half:] > 0].mean())
        self.assertGreater(slow, fast)


class EquilibriumOccupancyTests(IMP.test.TestCase):
    """The stationary state has a closed form, and it is `p ∝ 1/D`.

    The flux is discretised as `d[i]*p[i] - d[j]*p[j]`, i.e. `∂p/∂t = ∇²(Dp)`
    (Itô), whose stationary state is `D p = const` — **not** `∇·(D∇p)`, whose
    stationary state with no flux is uniform. Getting this wrong is not
    cosmetic: it is the difference between a dye that accumulates where it moves
    slowly and one that does not, and accumulating is the claim the mobility
    field exists to make.
    """

    def setUp(self):
        super().setUp()
        self.ng, self.dg = 25, 1.0
        self.bounds = open_box(self.ng)
        axis = np.arange(self.ng)[:, None, None]
        self.d_map = np.where(self.bounds > 0, 2.0 + 1.5 * np.cos(axis * 0.3), 0.0)

    def test_the_closed_form_is_normalised_and_masked(self):
        occupancy = equilibrium_occupancy(self.d_map, self.bounds)
        self.assertAlmostEqual(float(occupancy.sum()), 1.0, places=12)
        self.assertTrue(np.all(occupancy[self.bounds == 0] == 0.0))

    def test_d_times_p_is_constant(self):
        """The defining property of the ``"ito"`` branch, and only of it."""
        occupancy = equilibrium_occupancy(self.d_map, self.bounds, "ito")
        inside = self.bounds > 0
        product = self.d_map[inside] * occupancy[inside]
        self.assertLess(float(product.std() / product.mean()), 1e-12)

    def test_it_agrees_with_propagating_to_equilibrium(self):
        """The iterative solver is the independent check on the algebra."""
        solver = GridDiffusionSolver(
            self.d_map, self.bounds, self.bounds.copy(),
            t_step=0.5 * diffusion_stability_limit(self.d_map.max(), self.dg),
            dg=self.dg, flux_form="ito",
        )
        iterated = solver.equilibrium(n_steps=200000, tolerance=1e-14, n_check=1000)
        closed = equilibrium_occupancy(self.d_map, self.bounds, "ito")
        inside = self.bounds > 0
        deviation = np.max(np.abs(iterated[inside] - closed[inside]) / closed[inside])
        self.assertLess(float(deviation), 1e-10)

    def test_a_uniform_mobility_gives_a_uniform_occupancy(self):
        uniform = np.where(self.bounds > 0, 3.0, 0.0)
        occupancy = equilibrium_occupancy(uniform, self.bounds)
        inside = occupancy[self.bounds > 0]
        self.assertLess(float(inside.std() / inside.mean()), 1e-12)

    def test_slower_regions_hold_more_of_the_dye_under_ito(self):
        """Under ``"ito"`` only. Under the default this ratio is exactly 1."""
        d_map = np.where(self.bounds > 0, 4.0, 0.0)
        half = self.ng // 2
        d_map[:half] = np.where(self.bounds[:half] > 0, 0.25, 0.0)
        uniform = equilibrium_occupancy(d_map, self.bounds)
        self.assertAlmostEqual(
            float(uniform[self.bounds > 0].std()), 0.0, delta=1e-18)
        occupancy = equilibrium_occupancy(d_map, self.bounds, "ito")
        slow = occupancy[:half][self.bounds[:half] > 0].mean()
        fast = occupancy[half:][self.bounds[half:] > 0].mean()
        self.assertAlmostEqual(float(slow / fast), 16.0, delta=1e-9)

    def test_an_empty_domain_is_not_an_error(self):
        empty = np.zeros((self.ng,) * 3)
        self.assertEqual(float(equilibrium_occupancy(empty, empty).sum()), 0.0)


class PortedDefectTests(IMP.test.TestCase):
    """The three ChiSurf defects this port fixes, each stated as a property.

    Recorded as tests rather than only in the commit message, because each one
    was invisible in ChiSurf: the field still looked like a plausible field.
    """

    def test_the_map_builders_write_every_voxel(self):
        """`range(-npm, npm)` is one short; on an even `ng` it is two short."""
        for ng in (21, 22):
            density = np.ones((ng,) * 3)
            d_map = maps.diffusion_coefficient_map(
                density, np.zeros(3), 1.0, np.array([[500.0, 0.0, 0.0]]),
                free_diffusion=3.0, min_distance=1.0, slow_factor=0.5,
            )
            self.assertTrue(
                np.all(d_map == 3.0),
                f"ng={ng}: {int((d_map != 3.0).sum())} voxels unwritten",
            )

    def test_the_solver_advances_at_the_requested_rate(self):
        """ChiSurf swapped buffers only on odd steps and lost half the evolution."""
        ng, dg, D = 41, 1.0, 1.0
        t_step = 0.5 * diffusion_stability_limit(D, dg)
        n_steps = 200
        result = GridDiffusionSolver(
            np.full((ng,) * 3, D), open_box(ng), point_source(ng),
            t_step=t_step, dg=dg,
        ).run(n_steps, n_out=n_steps)
        axis = maps.grid_axis(ng, dg)
        variance = float(
            (result.density * axis[:, None, None] ** 2).sum() / result.density.sum()
        )
        ratio = variance / (2.0 * D * n_steps * t_step)
        self.assertAlmostEqual(ratio, 1.0, delta=0.01)
        # The defect gave 0.503 -- far outside that band, so this pins it.

    def test_the_outer_shell_holds_no_stale_population(self):
        """The stencil skips it; with ping-pong buffers it kept two-steps-ago data."""
        ng = 31
        bounds = open_box(ng)
        solver = GridDiffusionSolver(
            np.full((ng,) * 3, 1.0), bounds, point_source(ng),
            t_step=0.5 * diffusion_stability_limit(1.0, 1.0), dg=1.0,
        )
        result = solver.run(101, n_out=101)  # odd, so the buffers end swapped
        shell = result.density.copy()
        shell[1:-1, 1:-1, 1:-1] = 0.0
        self.assertEqual(float(shell.sum()), 0.0)


class FluxFormTests(IMP.test.TestCase):
    """Where the dye sits at equilibrium, and why it must not depend on ``D``.

    The inherited ChiSurf kernel discretises the flux as ``d[i]p[i] - d[j]p[j]``,
    i.e. ``∂p/∂t = ∇²(Dp)``, whose stationary state is ``p ∝ 1/D``: the dye
    accumulates wherever it moves slowly. That makes a *friction* field act as an
    *attractive potential*, which is not a convention -- it is wrong.
    Equilibrium is thermodynamics and mobility is kinetics; a dye slowed near the
    surface with no attraction is still found uniformly across its accessible
    volume, it merely takes longer to get around.

    The field's own canonical treatment agrees: the **Haas-Steinberg** equation
    for diffusion-modulated FRET is written as
    ``D d/dr [ p(r) d/dr ( N/p(r) ) ]`` precisely so its stationary state is the
    *given* ``p(r)`` for any ``D``, with ``p(r)`` coming from chain statistics
    and ``D`` a separate kinetic parameter fitted against it.

    ``flux_form="smoluchowski"`` is therefore the default; ``"ito"`` is kept so
    the inherited behaviour can be reproduced and compared.
    """

    def slab(self, ng=21, fast=8.0, slow=0.5):
        """A domain with a strong mobility contrast down one side."""
        bounds = open_box(ng)
        d_map = np.full((ng,) * 3, fast)
        d_map[:, :, 1:6] = slow
        return bounds, d_map, fast / slow

    def test_smoluchowski_equilibrium_is_uniform_whatever_the_mobility(self):
        bounds, d_map, _ = self.slab()
        occupancy = equilibrium_occupancy(d_map, bounds, "smoluchowski")
        inside = occupancy[bounds > 0]
        self.assertAlmostEqual(float(inside.max() / inside.min()), 1.0, delta=1e-12)
        self.assertAlmostEqual(float(occupancy.sum()), 1.0, delta=1e-12)
        self.assertEqual(float(occupancy[bounds == 0].sum()), 0.0)

    def test_ito_equilibrium_is_one_over_D(self):
        bounds, d_map, contrast = self.slab()
        occupancy = equilibrium_occupancy(d_map, bounds, "ito")
        inside = occupancy[bounds > 0]
        self.assertAlmostEqual(float(inside.max() / inside.min()), contrast, delta=1e-9)

    def test_the_two_forms_agree_when_the_mobility_is_uniform(self):
        ng = 15
        bounds = open_box(ng)
        d_map = np.full((ng,) * 3, 4.0)
        a = equilibrium_occupancy(d_map, bounds, "smoluchowski")
        b = equilibrium_occupancy(d_map, bounds, "ito")
        self.assertLess(float(np.abs(a - b).max()), 1e-15)

    def test_the_iterated_solver_reproduces_each_closed_form(self):
        """The closed forms are checked against the kernel, not asserted."""
        for form, tolerance in (("smoluchowski", 1e-12), ("ito", 1e-9)):
            bounds, d_map, _ = self.slab()
            start = bounds / bounds.sum()
            solver = GridDiffusionSolver(
                d_map, bounds, start, dg=1.0, flux_form=form,
                t_step=0.5 * diffusion_stability_limit(float(d_map.max()), 1.0),
            )
            iterated = solver.equilibrium(n_steps=60000, tolerance=1e-13)
            closed = equilibrium_occupancy(d_map, bounds, form)
            inside = bounds > 0
            deviation = float(
                np.abs(iterated[inside] - closed[inside]).max() / closed[inside].max())
            self.assertLess(deviation, tolerance, f"{form}: {deviation:.3g}")

    def test_smoluchowski_conserves_mass_without_decay(self):
        ng = 21
        bounds, d_map, _ = self.slab(ng)
        solver = GridDiffusionSolver(
            d_map, bounds, point_source(ng), dg=1.0, flux_form="smoluchowski",
            t_step=0.5 * diffusion_stability_limit(float(d_map.max()), 1.0),
        )
        result = solver.run(400, n_out=40)
        self.assertAlmostEqual(float(result.density.sum()), 1.0, delta=1e-10)

    def test_smoluchowski_free_diffusion_still_spreads_as_two_D_t(self):
        """A uniform mobility must give textbook diffusion under either form."""
        ng = 61
        d, dg = 1.0, 1.0
        t_step = 0.5 * diffusion_stability_limit(d, dg)
        n_steps = 600
        solver = GridDiffusionSolver(
            np.full((ng,) * 3, d), open_box(ng), point_source(ng),
            t_step=t_step, dg=dg, flux_form="smoluchowski",
        )
        density = solver.run(n_steps, n_out=n_steps).density
        axis = (np.arange(ng) - (ng - 1) // 2) * dg
        profile = density.sum(axis=(1, 2))
        variance = float((profile * axis ** 2).sum() / profile.sum())
        self.assertAlmostEqual(variance / (2.0 * d * n_steps * t_step), 1.0, delta=0.02)

    def test_an_unknown_flux_form_raises(self):
        bounds, d_map, _ = self.slab()
        self.assertRaises(
            ValueError, equilibrium_occupancy, d_map, bounds, "stratonovich")
        self.assertRaises(
            ValueError, GridDiffusionSolver, d_map, bounds, bounds,
            None, 1e-3, 1.0, True, "stratonovich")


class RateStabilityTests(IMP.test.TestCase):
    """The rate term constrains the step too, and it used to be ignored.

    ``diffusion_stability_limit`` validated only ``dt <= dg^2/(6D)``. On a real
    site the quenching term is far larger: T4L site 19 at 2.5 A, with a 25 A^2/ns
    ceiling on ``D`` setting the step, contributes 2.01 to the update coefficient
    against diffusion's 0.16. Their sum of 2.17 is past the divergence threshold
    of 2 and the decay reached 7e36, while the solver called the step safe. Site
    124, at 2.10, returned a smooth, finite, entirely plausible decay that was
    2.6 % wrong -- which is why this is pinned rather than left to inspection.

    Note the two thresholds: the coefficient goes negative at a sum of 1 and the
    scheme diverges at 2. PRD-111 stages 0-1 ran at sums of 0.92-0.98, inside
    both, and re-measuring after this fix reproduced them to about 1 %.

    The fix is to integrate the rate exactly, ``exp(-k dt)`` rather than
    ``1 - k dt``, so it contributes no stability constraint at all.
    """

    def setUp(self):
        super().setUp()
        self.ng = 21
        self.bounds = open_box(self.ng)
        self.dg = 1.0

    def test_the_limit_accounts_for_the_rate(self):
        without = diffusion_stability_limit(8.0, 1.0)
        with_rate = diffusion_stability_limit(8.0, 1.0, k_max=96.0)
        self.assertAlmostEqual(without, 1.0 / 48.0, delta=1e-12)
        self.assertAlmostEqual(with_rate, 1.0 / (48.0 + 96.0), delta=1e-12)
        self.assertLess(with_rate, without)

    def test_no_rate_reproduces_the_plain_limit(self):
        self.assertAlmostEqual(
            diffusion_stability_limit(4.0, 2.0),
            2.0 ** 2 / (6.0 * 4.0), delta=1e-15)

    def test_a_fast_rate_no_longer_diverges(self):
        """``k dt`` far above 1: fatal when subtracted, harmless as a factor."""
        d_map = np.full((self.ng,) * 3, 8.0)
        rate = np.full((self.ng,) * 3, 96.0)
        t_step = 0.5 * diffusion_stability_limit(8.0, self.dg)   # k*dt = 1.0
        solver = GridDiffusionSolver(
            d_map, self.bounds, self.bounds / self.bounds.sum(), rate,
            t_step=t_step, dg=self.dg)
        result = solver.run(400, n_out=40)
        self.assertTrue(np.all(np.isfinite(result.fluorescence)))
        self.assertTrue(np.all(result.fluorescence >= 0.0))
        self.assertTrue(np.all(np.diff(result.fluorescence) <= 1e-12))
        self.assertGreaterEqual(float(result.density.min()), 0.0)

    def test_a_uniform_rate_is_exact_whatever_the_step(self):
        """With ``k`` uniform the answer is ``exp(-k t)`` -- to round-off."""
        d_map = np.full((self.ng,) * 3, 8.0)
        k = 96.0
        rate = np.full((self.ng,) * 3, k)
        t_step = 0.5 * diffusion_stability_limit(8.0, self.dg)
        n_steps = 200
        solver = GridDiffusionSolver(
            d_map, self.bounds, self.bounds / self.bounds.sum(), rate,
            t_step=t_step, dg=self.dg)
        result = solver.run(n_steps, n_out=n_steps)
        expected = np.exp(-k * n_steps * t_step)
        self.assertAlmostEqual(
            float(result.fluorescence[-1]) / expected, 1.0, delta=1e-10)


if __name__ == "__main__":
    IMP.test.main()
