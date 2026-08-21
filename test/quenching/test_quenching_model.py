"""PRD-109 stage 3: the particle model's objects, moved here from QuEst.

``DyeDiffusionSimulation`` and ``QuenchedDonorDecay`` -- what QuEst called
``SimulateDiffusion`` and ``DonorDecay``, less the application layer.
"""

import numpy as np

import IMP
import IMP.test

import IMP.bff as sites
from IMP.bff import grid_center_index, quenching_rate_grid
from IMP.bff import (
    AccessibleVolume,
    DyeDiffusionSimulation,
    ObstacleAtoms,
    QuenchedDonorDecay,
    resolve_trajectory_count,
    trajectory_seeds,
)
from IMP.bff import ResidueQuenching
from IMP.bff import amino_acid_quenching_defaults


ATOM_DTYPE = [
    ("chain", "U4"), ("res_id", "i8"), ("res_name", "U4"),
    ("atom_name", "U4"), ("coord", "f8", 3),
]


def make_atoms(rows):
    """The obstacles as an :class:`ObstacleAtoms` value (parallel arrays)."""
    atoms = ObstacleAtoms()
    atoms.chains = [r[0] for r in rows]
    atoms.res_ids = [int(r[1]) for r in rows]
    atoms.res_names = [r[2] for r in rows]
    atoms.atom_names = [r[3] for r in rows]
    atoms.coords = [float(v) for r in rows for v in r[4]]
    return atoms


def sphere_density(ng=41, radius_voxels=15):
    i = np.arange(ng) - (ng - 1) // 2
    x, y, z = np.meshgrid(i, i, i, indexing="ij")
    return ((x ** 2 + y ** 2 + z ** 2) < radius_voxels ** 2).astype(np.uint8)


def _sites(atoms, table=None):
    """`residue_sites` over an :class:`ObstacleAtoms` value."""
    kwargs = {}
    if table is not None:
        kwargs["table"] = table
    return sites.residue_sites(
        atoms.chains, atoms.res_ids, atoms.res_names, atoms.atom_names,
        np.ascontiguousarray(atoms.coords, dtype=np.float64).reshape(-1, 3),
        **kwargs)


class ResidueSiteTests(IMP.test.TestCase):

    def test_the_slow_centre_is_cb_and_the_quench_centre_the_moiety(self):
        atoms = make_atoms([
            ("A", 1, "TRP", "CB", [0.0, 0.0, 0.0]),
            ("A", 1, "TRP", "NE1", [4.0, 0.0, 0.0]),
            ("A", 1, "TRP", "CD2", [6.0, 0.0, 0.0]),
        ])
        found = _sites(atoms).get_slow_centers().reshape(-1, 3)
        slow = found
        quench = _sites(atoms).get_quench_centers().reshape(-1, 3)
        self.assertEqual(list(_sites(atoms).get_residue_names()), ["TRP"])
        self.assertTrue(np.allclose(slow[0], [0.0, 0.0, 0.0]))
        # The centroid of the PET-active atoms present, not CB.
        self.assertTrue(np.allclose(quench[0], [5.0, 0.0, 0.0]))

    def test_ca_stands_in_when_there_is_no_cb(self):
        atoms = make_atoms([("A", 1, "GLY", "CA", [1.0, 2.0, 3.0])])
        sites_value = _sites(atoms)
        self.assertTrue(
            np.allclose(sites_value.get_slow_centers()[0:3], [1.0, 2.0, 3.0]))

    def test_a_residue_with_no_active_atom_falls_back_to_its_slow_centre(self):
        """So a user-supplied rate still has a well-defined centre."""
        atoms = make_atoms([("A", 1, "TRP", "CB", [1.0, 1.0, 1.0])])
        found = _sites(atoms)
        self.assertTrue(
            np.allclose(found.get_quench_centers(),
                        found.get_slow_centers())
        )

    def test_residues_are_keyed_by_chain_as_well_as_number(self):
        """Numbers restart per chain: in a homodimer every one occurs twice.

        QuEst keyed on ``res_id`` alone and folded both chains' side chains into
        one centre on every multi-chain structure.
        """
        atoms = make_atoms([
            ("A", 1, "TRP", "CB", [0.0, 0.0, 0.0]),
            ("B", 1, "TRP", "CB", [50.0, 0.0, 0.0]),
        ])
        found = _sites(atoms)
        self.assertEqual(found.size(), 2)
        slow = found.get_slow_centers().reshape(-1, 3)
        self.assertAlmostEqual(
            float(np.linalg.norm(slow[0] - slow[1])),
            50.0,
        )

    def test_a_custom_quench_atom_selection_is_honoured(self):
        atoms = make_atoms([
            ("A", 1, "TRP", "CB", [0.0, 0.0, 0.0]),
            ("A", 1, "TRP", "NE1", [4.0, 0.0, 0.0]),
        ])
        table = {"TRP": ResidueQuenching(quench_atoms=["CB"], kQ=1.0)}
        found = _sites(atoms, table)
        self.assertTrue(
            np.allclose(found.get_quench_centers()[0:3], [0.0, 0.0, 0.0]))

    def test_the_per_residue_lookups_follow_the_name_order(self):
        table = amino_acid_quenching_defaults()
        names = ["TRP", "ALA", "MET"]
        rates = sites.quenching_rates_for_residues(names, table)
        self.assertGreater(rates[0], 0.0)
        self.assertEqual(rates[1], 0.0)
        self.assertGreater(rates[2], 0.0)
        self.assertGreater(rates[0], rates[2])

    def test_an_unset_radius_inherits_the_critical_distance(self):
        table = {"TRP": ResidueQuenching(kQ=1.0),          # radius: inherit
                 "TYR": ResidueQuenching(kQ=1.0, quench_radius=9.0)}
        radii = sites.quench_radii_for_residues(["TRP", "TYR"], table, 6.5)
        self.assertAlmostEqual(radii[0], 6.5)
        self.assertAlmostEqual(radii[1], 9.0)


class TrajectorySeedTests(IMP.test.TestCase):

    def test_a_base_seed_gives_distinct_reproducible_seeds(self):
        first = list(trajectory_seeds(11, 4))
        self.assertEqual(first, list(trajectory_seeds(11, 4)))
        self.assertEqual(len(set(first)), 4)

    def test_no_seed_gives_fresh_ones_each_time(self):
        # A negative seed is "draw freely" -- the C++ spelling of ``None``.
        self.assertNotEqual(list(trajectory_seeds(-1, 4)),
                            list(trajectory_seeds(-1, 4)))

    def test_no_seed_and_one_trajectory_stays_unseeded(self):
        self.assertEqual(list(trajectory_seeds(-1, 1)), [-1])

    def test_the_trajectory_count_is_capped(self):
        from IMP.bff import MAX_PARALLEL_TRAJECTORIES
        self.assertEqual(resolve_trajectory_count(100),
                         MAX_PARALLEL_TRAJECTORIES)
        self.assertEqual(resolve_trajectory_count(2), 2)
        self.assertGreaterEqual(resolve_trajectory_count(-1), 1)
        self.assertEqual(resolve_trajectory_count(0), 1)


class DyeDiffusionSimulationTests(IMP.test.TestCase):

    def setUp(self):
        super().setUp()
        self.density = sphere_density()
        self.dg = 1.0
        self.x0 = np.array([10.0, -5.0, 2.0])

    def simulation(self, rate_map=None):
        return DyeDiffusionSimulation(
            self.density, self.dg, self.x0, quenching_rate_map=rate_map
        )

    def test_the_trajectory_comes_back_in_the_structure_frame(self):
        """The walk runs on the grid; the anchor puts it back on the protein."""
        simulation = self.simulation()
        simulation.run(t_max=100.0, t_step=0.002, n_trajectories=1, random_seed=3)
        self.assertAlmostEqual(
            float(np.linalg.norm(simulation.get_mean_position() - self.x0)), 0.0,
            delta=3.0)

    def test_trajectories_are_concatenated_not_averaged(self):
        one = self.simulation()
        one.run(t_max=100.0, t_step=0.002, n_trajectories=1, random_seed=3)
        four = self.simulation()
        four.run(t_max=100.0, t_step=0.002, n_trajectories=4, random_seed=3)
        self.assertEqual(four.n_frames, 4 * one.n_frames)

    def test_a_seed_pins_the_walk(self):
        first = self.simulation()
        first.run(t_max=100.0, t_step=0.002, n_trajectories=2, random_seed=5)
        second = self.simulation()
        second.run(t_max=100.0, t_step=0.002, n_trajectories=2, random_seed=5)
        self.assertTrue(
            np.array_equal(first.get_trajectory(), second.get_trajectory()))

    def test_an_empty_volume_yields_no_trajectory(self):
        simulation = DyeDiffusionSimulation(
            np.zeros_like(self.density), self.dg, self.x0
        )
        self.assertEqual(simulation.run(t_max=10.0, random_seed=1), 0)
        self.assertEqual(simulation.n_frames, 0)

    def test_no_rate_map_means_no_quenching(self):
        simulation = self.simulation()
        simulation.run(t_max=100.0, t_step=0.002, n_trajectories=1, random_seed=3)
        self.assertEqual(float(simulation.get_k_quench().sum()), 0.0)
        self.assertEqual(simulation.collision_fraction, 0.0)

    def test_the_rate_is_read_from_where_the_quencher_was_stamped(self):
        """The sampling map must invert the stamping map, or the walk reads
        rates from beside the quenchers. Two conventions carry that, and both
        were wrong in QuEst once: the integer centre offset, and floor-not-trunc.
        """
        ng = self.density.shape[0]
        rate_map = quenching_rate_grid(
            self.density, ng, self.dg, np.array([30.0]),
            np.array([self.x0]), self.x0, np.array([3.0]),
        )
        simulation = self.simulation(rate_map)
        simulation.run(t_max=200.0, t_step=0.002, n_trajectories=1, random_seed=3)
        # The sphere of influence covers the whole volume, so every frame sees it.
        self.assertAlmostEqual(float(simulation.get_k_quench().min()), 3.0, places=5)
        self.assertTrue(np.all(simulation.get_k_quench() > 0.0))
        self.assertAlmostEqual(simulation.collision_fraction, 1.0)

    def test_sampling_matches_a_hand_computed_voxel(self):
        ng = self.density.shape[0]
        field = np.arange(ng ** 3, dtype=np.float64).reshape((ng,) * 3)
        simulation = self.simulation()
        simulation.run(t_max=50.0, t_step=0.002, n_trajectories=1, random_seed=3)
        sampled = simulation.sample_grid(field, ng)
        centre = grid_center_index(ng)
        expected = []
        for position in simulation.get_trajectory().reshape(-1, 3):
            idx = np.floor((position - self.x0) / self.dg + centre).astype(int)
            expected.append(field[tuple(idx)])
        self.assertTrue(np.allclose(sampled, expected))

    def test_reading_results_before_running_raises(self):
        simulation = self.simulation()
        for call in (lambda: simulation.get_mean_position(),
                     lambda: simulation.get_k_quench(),
                     lambda: simulation.sample_grid(
                         np.zeros(self.density.size, dtype=np.float64),
                         self.density.shape[0])):
            with self.assertRaises(ValueError):
                call()


class DonorModelFixture:
    """Shared setup. Not a TestCase -- subclassing one re-runs all its tests."""

    def setUp(self):
        super().setUp()
        self.x0 = np.zeros(3)
        self.av = AccessibleVolume(
            density=sphere_density(41, 15).astype(np.float64), grid_step=1.0,
            attachment_point=self.x0)
        # One tryptophan sitting inside the volume, and an inert alanine.
        self.atoms = make_atoms([
            ("A", 1, "TRP", "CB", [5.0, 0.0, 0.0]),
            ("A", 1, "TRP", "NE1", [7.0, 0.0, 0.0]),
            ("A", 1, "TRP", "CD2", [7.0, 1.0, 0.0]),
            ("A", 2, "ALA", "CB", [-6.0, 0.0, 0.0]),
        ])

    def model(self, **kwargs):
        options = dict(
            tau0=4.0,
            quenching_table=amino_acid_quenching_defaults(),
            critical_distance=7.0, slow_radius=10.0,
            t_max=400.0, t_step=0.004, n_photons=40000,
            n_trajectories=1, random_seed=7,
        )
        options.update(kwargs)
        return QuenchedDonorDecay(self.av, self.atoms, **options)


class QuenchedDonorDecayTests(DonorModelFixture, IMP.test.TestCase):
    """Driven on a synthetic AV so the test needs no structure file."""

    def test_only_the_quenching_residue_stamps_a_rate(self):
        model = self.model()
        model.update_grids()
        rate_map = model.get_quenching_rate_map()
        self.assertGreater(float(rate_map.max()), 0.0)
        # TRP's rate, not ALA's zero.
        table = amino_acid_quenching_defaults()
        self.assertAlmostEqual(float(rate_map.max()), table["TRP"].kQ, places=6)

    def test_quenching_lowers_the_quantum_yield_and_the_lifetime(self):
        quenched = self.model()
        free = self.model(quenching_table={})
        self.assertLess(quenched.quantum_yield, free.quantum_yield)
        self.assertLess(quenched.fluorescence_lifetime,
                        free.fluorescence_lifetime)

    def test_without_quenchers_the_lifetime_is_tau0(self):
        free = self.model(quenching_table={}, n_photons=200000)
        self.assertAlmostEqual(free.quantum_yield, 1.0, places=6)
        self.assertAlmostEqual(free.fluorescence_lifetime, 4.0, delta=0.05)

    def test_a_stronger_quencher_quenches_more(self):
        weak = self.model(quenching_table={"TRP": ResidueQuenching(kQ=0.1)})
        strong = self.model(quenching_table={"TRP": ResidueQuenching(kQ=20.0)})
        self.assertGreater(weak.quantum_yield, strong.quantum_yield)

    def test_the_histogram_counts_only_emitted_photons(self):
        """A quenched excitation returns dt = 0 and would pile into bin 0."""
        model = self.model()
        emitted = int(model.get_emitted().sum())
        self.assertLess(emitted, model.n_photons)
        flat = model.decay_histogram(512, 0.0, 400.0)
        edges, counts = flat[:513], flat[513:].astype(int)
        # The curve holds QY * N photons, not N. Histogramming the whole trace
        # would give N, because a quenched excitation comes back as dt = 0 and
        # piles into the first bin -- a spike of photons that never existed.
        self.assertEqual(int(counts.sum()), emitted)
        self.assertLess(int(counts.sum()), model.n_photons)

    def test_the_model_is_reproducible(self):
        self.assertAlmostEqual(
            self.model().quantum_yield, self.model().quantum_yield, places=12
        )

    def test_zero_photons_is_an_empty_simulation_not_a_crash(self):
        self.assertEqual(self.model(n_photons=0).quantum_yield, 0.0)

    def test_the_photon_seed_differs_from_the_walk_seed(self):
        """Or the photon draws replay the trajectory they are scored against."""
        model = self.model(random_seed=11)
        self.assertNotEqual(model.get_photon_seed(), 11)
        # -1 is the C++ spelling of "draw freely": an unseeded model must not
        # derive a seed from the one it does not have.
        self.assertEqual(self.model(random_seed=-1).get_photon_seed(), -1)

    def test_sites_are_cached_across_uses(self):
        model = self.model()
        a = model.get_sites()
        b = model.get_sites()
        self.assertEqual(a.size(), b.size())


class QuenchedDonorFretTests(DonorModelFixture, IMP.test.TestCase):

    def test_a_far_acceptor_transfers_almost_nothing(self):
        import math
        model = self.model()
        far = np.array([[400.0, 0.0, 0.0, 1.0]])
        rates = model.fret_rate_trace_cloud(
            np.ascontiguousarray(far[:, :3]), 52.0, 2.0 / 3.0, 7.0)
        # (1/tau0) * (R0/r)^6 = 0.25 * (52/400)^6 ~ 1.9e-6 /ns, i.e. nothing
        # against the 0.25 /ns radiative rate.
        self.assertTrue(np.all(rates < 1e-5))
        self.assertLess(float(rates.max()) * model.tau0, 1e-4)

    def test_efficiency_rises_as_the_acceptor_approaches(self):
        model = self.model()
        near = np.array([[40.0, 0.0, 0.0, 1.0]])
        far = np.array([[90.0, 0.0, 0.0, 1.0]])
        r_near = model.fret_rate_trace_cloud(
            np.ascontiguousarray(near[:, :3]), 52.0, 2.0 / 3.0, 7.0)
        r_far = model.fret_rate_trace_cloud(
            np.ascontiguousarray(far[:, :3]), 52.0, 2.0 / 3.0, 7.0)
        self.assertGreater(model.fret_efficiency(r_near),
                           model.fret_efficiency(r_far))

    def test_efficiency_is_bounded(self):
        """E = 1 - QY_DA/QY_D, and both terms are sampled.

        The two traces share a seed, so they are strongly correlated and the
        difference is far quieter than either yield -- but it is still a
        difference of two Monte-Carlo estimates, so a genuinely-zero efficiency
        scatters about zero rather than sitting on it. The tolerance says that,
        instead of pretending the estimator is exact.
        """
        model = self.model(n_photons=200000)
        for cloud in (np.array([[10.0, 0.0, 0.0, 1.0]]),
                      np.array([[300.0, 0.0, 0.0, 1.0]])):
            rates = model.fret_rate_trace_cloud(
                np.ascontiguousarray(cloud[:, :3]), 52.0, 2.0 / 3.0, 7.0)
            efficiency = model.fret_efficiency(rates)
            self.assertGreaterEqual(efficiency, -0.01)
            self.assertLessEqual(efficiency, 1.0 + 1e-9)

    def test_a_very_close_acceptor_transfers_almost_everything(self):
        model = self.model(n_photons=200000)
        close = np.array([[8.0, 0.0, 0.0, 1.0]])
        rates = model.fret_rate_trace_cloud(
            np.ascontiguousarray(close[:, :3]), 52.0, 2.0 / 3.0, 7.0)
        self.assertGreater(model.fret_efficiency(rates), 0.95)

    def test_a_paired_acceptor_model_is_frame_paired(self):
        donor = self.model()
        acceptor = self.model()
        rates = donor.fret_rate_trace_paired(acceptor, 52.0, 2.0 / 3.0, 7.0)
        self.assertEqual(rates.shape[0], donor.n_frames)


if __name__ == "__main__":
    IMP.test.main()
