"""`Minimizer::compute_objective_batch`: many candidates, one crossing.

A chi-square surface, a support-plane interval, a population sampler and a
random restart all ask one question of a great many points, and asked one at
a time each point pays a boundary crossing for arithmetic that is often
shorter than the crossing itself. The batch walks them in C++.

The numbers must be the optimiser's own -- not a second implementation that
agrees to five digits -- so the tests here compare against driving the graph
by hand, and require exact equality rather than a tolerance.
"""

import unittest

import numpy as np

from IMP import bff


def decay_fit(n=128):
    x = np.linspace(0.1, 15.0, n)
    y = 2.0 * np.exp(-x / 3.0)

    curve = bff.GraphExpression("model")
    curve.set_expression("a*exp(-x/t)")
    ports = {}
    for v in curve.get_variable_names():
        if v == "x":
            continue
        p = bff.GraphPort(1.0)
        curve.add_input_port(v, p)
        ports[v] = p
    axis = bff.GraphPort([0.0])
    axis.set_values_array(np.ascontiguousarray(x))
    curve.add_input_port("x", axis)
    out = bff.GraphPort([0.0], False, True)
    curve.add_output_port("model", out)

    chi2 = bff.ChiSquared("chi2")
    chi2.set_data_arrays(np.ascontiguousarray(y), np.ones(n))
    model_in = bff.GraphPort([0.0])
    model_in.link = out
    chi2.add_input_port("model", model_in)
    chi2.add_output_port("chi2", bff.GraphPort(0.0, False, True))
    chi2.add_output_port("residuals", bff.GraphPort([0.0], False, True))
    chi2._graph = (curve, out, model_in, axis)

    free = [ports["a"], ports["t"]]
    for port, v in zip(free, (2.0, 3.0)):
        port.value = v
    m = bff.Minimizer()
    m.set_parameter_ports(free)
    m.set_objective(chi2, "residuals")
    m._graph = (chi2, curve, ports)
    return m, chi2, free


def one_at_a_time(chi2, free, candidates):
    out = []
    for row in candidates:
        for port, v in zip(free, row):
            port.value = float(v)
        chi2.update()
        out.append(chi2.get_chi2())
    return np.asarray(out)


class ObjectiveBatchTests(unittest.TestCase):
    def setUp(self):
        self.m, self.chi2, self.free = decay_fit()
        self.candidates = np.array([
            [2.0, 3.0],      # the truth: zero
            [1.8, 2.7],
            [2.2, 3.4],
            [2.0, 3.0],      # again, to show the walk is stateless
            [0.5, 12.0],
        ])

    def test_the_numbers_are_the_graphs_own(self):
        batch = self.m.compute_objective_batch(self.candidates)
        by_hand = one_at_a_time(self.chi2, self.free, self.candidates)
        np.testing.assert_array_equal(batch, by_hand)

    def test_it_returns_one_value_per_candidate(self):
        batch = self.m.compute_objective_batch(self.candidates)
        self.assertEqual(batch.shape, (len(self.candidates),))
        self.assertEqual(batch[0], 0.0)          # scored at the truth
        self.assertEqual(batch[0], batch[3])     # the same point twice

    def test_the_ports_come_back_where_they_were(self):
        """A surface scan must not quietly move a fit somebody is holding."""
        for port, v in zip(self.free, (1.234, 5.678)):
            port.value = v
        self.chi2.update()
        held = self.chi2.get_chi2()
        self.m.compute_objective_batch(self.candidates)
        self.assertEqual([p.value for p in self.free], [1.234, 5.678])
        # and the graph is consistent with them, not with the last candidate
        self.assertEqual(self.chi2.get_chi2(), held)

    def test_a_wrong_number_of_columns_is_refused(self):
        with self.assertRaises(ValueError):
            self.m.compute_objective_batch(np.zeros((3, 5)))

    def test_no_candidates_is_an_empty_answer_not_an_error(self):
        batch = self.m.compute_objective_batch(np.zeros((0, 2)))
        self.assertEqual(batch.shape, (0,))

    def test_a_sanitising_parameter_port_never_delivers_the_nan(self):
        """Worth pinning before the next test, because it is why a NaN
        candidate is so rarely seen: a `GraphPort` sanitises by default, and a NaN
        written into one arrives as the smallest positive double. The
        objective is then finite and merely very bad, which is the behaviour
        a fit wants and not something this batch changes."""
        batch = self.m.compute_objective_batch(np.array([[np.nan, 3.0]]))
        self.assertTrue(np.isfinite(batch[0]))

    def test_a_nan_that_does_reach_the_objective_scores_infinite(self):
        """With sanitising off the NaN travels, and then it must come back as
        `+inf` -- `ChiSquared`'s convention, so a sampler rejects the
        candidate instead of carrying a NaN into the posterior."""
        for port in self.free:
            port.set_sanitize(False)
        batch = self.m.compute_objective_batch(
            np.array([[np.nan, 3.0], [2.0, 3.0]]))
        self.assertTrue(np.isinf(batch[0]))
        self.assertFalse(np.isnan(batch[0]))
        self.assertEqual(batch[1], 0.0)

    def test_it_agrees_with_run_at_the_solution(self):
        """The batch and the optimiser are the same objective, so the batch
        scored at the fitted point is the fit's own chi-square."""
        for port, v in zip(self.free, (1.0, 1.0)):
            port.value = v
        info = self.m.run()
        self.assertIn(info, (1, 2, 3, 4))
        batch = self.m.compute_objective_batch(np.asarray([self.m.x]))
        self.assertAlmostEqual(batch[0], self.m.get_chi2(), places=12)

    def test_a_surface_scan(self):
        """What it is for, at the size it is for: a 60 x 60 grid in one call,
        with the minimum where the fit says it is."""
        a = np.linspace(1.5, 2.5, 60)
        t = np.linspace(2.5, 3.5, 60)
        grid = np.stack(np.meshgrid(a, t, indexing="ij"), axis=-1).reshape(-1, 2)
        surface = self.m.compute_objective_batch(
            np.ascontiguousarray(grid)).reshape(60, 60)
        i, j = np.unravel_index(np.argmin(surface), surface.shape)
        self.assertAlmostEqual(a[i], 2.0, delta=0.02)
        self.assertAlmostEqual(t[j], 3.0, delta=0.02)
        # the grid does not land exactly on (2, 3), so the floor is small
        # rather than zero; what matters is that it is where the fit is
        self.assertLess(surface.min(), 1e-3)

    def test_a_batch_without_an_objective_is_refused(self):
        m = bff.Minimizer()
        with self.assertRaises(ValueError):
            m.compute_objective_batch(np.zeros((2, 2)))


if __name__ == "__main__":
    unittest.main()
