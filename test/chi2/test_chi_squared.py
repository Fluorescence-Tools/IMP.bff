"""bff.FitChiSquared against ChiSurf's own residual arithmetic.

The objective is the one number a sampler trusts, so it is checked against
the Python it was ported from rather than against hand-written expectations.
ChiSurf is imported when it is available; the parity tests skip without it,
while the arithmetic tests always run.
"""

import unittest

import numpy as np

from IMP import bff

try:
    from chisurf.core.fitting import (
        calculate_weighted_residuals as cs_wres,
        deviance_residuals as cs_deviance,
    )
    HAVE_CHISURF = True
except Exception:  # pragma: no cover - depends on the environment
    HAVE_CHISURF = False


def make_chi2(y, ey, xmin=0, xmax=-1, noise="default", mask=None):
    node = bff.FitChiSquared("chi2")
    node.set_data(list(map(float, y)), list(map(float, ey)))
    node.set_fit_range(xmin, xmax)
    node.set_noise_model_name(noise)
    if mask is not None:
        node.set_mask(list(map(float, mask)))
    return node


class ChiSquaredArithmeticTests(unittest.TestCase):
    def test_neyman_residuals(self):
        c = make_chi2([10.0, 12.0, 9.0], [1.0, 2.0, 1.0])
        self.assertEqual(list(c.compute_weighted_residuals([9.0, 10.0, 9.0])),
                         [1.0, 1.0, 0.0])
        self.assertEqual(c.compute_chi2([9.0, 10.0, 9.0]), 2.0)

    def test_fit_range_restricts_the_residuals(self):
        c = make_chi2([1.0, 2.0, 3.0, 4.0], [1.0] * 4, xmin=1, xmax=3)
        wres = c.compute_weighted_residuals([0.0, 0.0, 0.0, 0.0])
        self.assertEqual(list(wres), [2.0, 3.0])

    def test_open_ended_range_runs_to_the_end(self):
        c = make_chi2([1.0, 2.0, 3.0], [1.0] * 3, xmin=1, xmax=-1)
        self.assertEqual(len(c.compute_weighted_residuals([0.0] * 3)), 2)

    def test_short_model_truncates_rather_than_failing(self):
        """ChiSurf takes the shorter of data and model, so a model that
        stops early is a shorter residual vector, not an error."""
        c = make_chi2([1.0, 2.0, 3.0, 4.0], [1.0] * 4)
        self.assertEqual(len(c.compute_weighted_residuals([0.0, 0.0])), 2)

    def test_nan_becomes_an_infinite_misfit(self):
        c = make_chi2([1.0, 2.0], [1.0, 1.0])
        self.assertEqual(c.compute_chi2([float("nan"), 0.0]), float("inf"))

    def test_mask_multiplies_the_residuals(self):
        c = make_chi2([1.0, 2.0, 3.0], [1.0] * 3, mask=[1.0, 0.0, 1.0])
        self.assertEqual(c.compute_chi2([0.0, 0.0, 0.0]), 1.0 + 9.0)

    def test_mask_of_the_wrong_length_is_ignored(self):
        """ChiSurf's guard: a mask that no longer spans the window is
        skipped outright rather than applied to the wrong points."""
        c = make_chi2([1.0, 2.0, 3.0], [1.0] * 3, mask=[0.0, 0.0])
        self.assertEqual(c.compute_chi2([0.0, 0.0, 0.0]), 1.0 + 4.0 + 9.0)

    def test_reduced_chi2_uses_the_free_parameter_count(self):
        c = make_chi2([1.0] * 10, [1.0] * 10)
        node_out = bff.GraphPort(0.0, name="chi2")
        c.add_input_port("model", bff.GraphPort([0.0] * 10))
        c.add_output_port("chi2", node_out)
        c.evaluate()
        self.assertEqual(c.get_chi2(), 10.0)
        self.assertAlmostEqual(c.get_chi2r(2), 10.0 / (10 - 2 - 1))

    def test_unknown_noise_model_is_refused(self):
        c = bff.FitChiSquared("chi2")
        with self.assertRaises(ValueError):
            c.set_noise_model_name("student-t")

    def test_mismatched_errors_are_refused(self):
        c = bff.FitChiSquared("chi2")
        with self.assertRaises(ValueError):
            c.set_data([1.0, 2.0], [1.0])


class ChiSquaredAsANodeTests(unittest.TestCase):
    def _graph(self, y, ey, model):
        node = make_chi2(y, ey)
        node.add_input_port("model", bff.GraphPort(list(model)))
        node.add_output_port("chi2", bff.GraphPort(0.0, name="chi2"))
        return node

    def test_evaluating_the_node_writes_chi2_to_its_output(self):
        node = self._graph([10.0, 12.0], [1.0, 1.0], [9.0, 11.0])
        node.evaluate()
        self.assertEqual(node.get_output_port("chi2").value, 2.0)

    def test_updating_the_model_port_changes_chi2(self):
        node = self._graph([10.0, 12.0], [1.0, 1.0], [9.0, 11.0])
        node.evaluate()
        node.get_input_port("model").set_value_vector([10.0, 12.0])
        node.evaluate()
        self.assertEqual(node.get_output_port("chi2").value, 0.0)

    def test_a_missing_output_port_is_refused(self):
        node = make_chi2([1.0], [1.0])
        node.add_input_port("model", bff.GraphPort([0.0]))
        with self.assertRaises(ValueError):
            node.evaluate()

    def test_a_missing_model_port_is_refused(self):
        node = make_chi2([1.0], [1.0])
        node.add_output_port("chi2", bff.GraphPort(0.0, name="chi2"))
        with self.assertRaises(ValueError):
            node.evaluate()


@unittest.skipUnless(HAVE_CHISURF, "chisurf not importable")
class ChiSurfParityTests(unittest.TestCase):
    """The port must agree with the Python it came from, point for point."""

    def _curves(self, seed, n=64):
        rng = np.random.default_rng(seed)
        y = rng.uniform(5.0, 500.0, n)
        ey = np.sqrt(y)
        model = y + rng.normal(0.0, 3.0, n)
        return y, ey, model

    def _cs_wres(self, y, ey, model, xmin, xmax, noise):
        import chisurf.core.data
        import chisurf.core.curve
        data = chisurf.core.data.DataCurve(
            x=np.arange(y.size, dtype=np.float64), y=y, ey=ey)
        curve = chisurf.core.curve.Curve(
            x=np.arange(model.size, dtype=np.float64), y=model)
        return np.asarray(
            cs_wres(data, curve, xmin, xmax, noise_model=noise))

    def test_neyman_matches_chisurf(self):
        for seed in range(6):
            y, ey, model = self._curves(seed)
            ref = self._cs_wres(y, ey, model, 0, y.size, "default")
            got = np.asarray(
                make_chi2(y, ey).compute_weighted_residuals(model.tolist()))
            np.testing.assert_allclose(got, ref, rtol=1e-12, atol=1e-12)

    def test_poisson_matches_chisurf(self):
        for seed in range(6):
            y, ey, model = self._curves(seed)
            ref = self._cs_wres(y, ey, model, 0, y.size, "poisson")
            got = np.asarray(
                make_chi2(y, ey, noise="poisson")
                .compute_weighted_residuals(model.tolist()))
            np.testing.assert_allclose(got, ref, rtol=1e-12, atol=1e-12)

    def test_poisson_matches_chisurf_at_zero_counts(self):
        """The y == 0 branch, where the y log(y/mu) term is defined to vanish."""
        y = np.array([0.0, 0.0, 3.0, 10.0])
        mu = np.array([1e-30, 2.0, 3.0, 4.0])
        ref = np.asarray(cs_deviance(y, mu))
        got = np.asarray(
            make_chi2(y, np.ones_like(y), noise="poisson")
            .compute_weighted_residuals(mu.tolist()))
        np.testing.assert_allclose(got, ref, rtol=1e-12, atol=1e-12)

    def test_fit_range_matches_chisurf(self):
        y, ey, model = self._curves(7, n=48)
        for xmin, xmax in ((0, 48), (5, 40), (10, 11), (0, 1)):
            ref = self._cs_wres(y, ey, model, xmin, xmax, "default")
            got = np.asarray(
                make_chi2(y, ey, xmin=xmin, xmax=xmax)
                .compute_weighted_residuals(model.tolist()))
            np.testing.assert_allclose(got, ref, rtol=1e-12, atol=1e-12,
                                       err_msg=f"window ({xmin}, {xmax})")

    def test_chi2_matches_chisurf(self):
        for seed in range(6):
            y, ey, model = self._curves(seed)
            ref = float((self._cs_wres(y, ey, model, 0, y.size,
                                       "default") ** 2).sum())
            got = make_chi2(y, ey).compute_chi2(model.tolist())
            self.assertAlmostEqual(got, ref, places=9)


if __name__ == "__main__":
    unittest.main()


class WholeFitInCppTests(unittest.TestCase):
    """A polynomial fit whose model *and* objective are bff nodes.

    This is the end of the chain the port exists for: Sampler drives
    parameters, the graph computes the curve, FitChiSquared scores it, and no
    step of a move enters Python.
    """

    def setUp(self):
        import sys
        sys.path.insert(0, __file__.rsplit("/", 1)[0])
        from bench_fit_objective import build_graph
        rng = np.random.default_rng(0)
        self.x = np.linspace(1.0, 2.0, 96)
        self.truth = (2.0, 0.5, 1.0)  # a, b, c
        self.y = (1.0 + 2.0 * self.x + 0.5 * self.x ** 2
                  + rng.normal(0.0, 0.1, self.x.size))
        self.ey = np.ones_like(self.y) * 0.1
        self.params, self.objective = build_graph(self.x, self.y, self.ey)

    def _chi2_at(self, a, b, c):
        self.params[0].value, self.params[1].value, self.params[2].value = a, b, c
        self.objective.update()
        return self.objective.get_output_port("chi2").value

    def test_graph_chi2_matches_numpy(self):
        """The graph's chi-square is the textbook one, not an approximation."""
        for a, b, c in ((2.0, 0.5, 1.0), (1.0, 0.0, 0.0), (2.3, 0.4, 0.8)):
            model = c + a * self.x + b * self.x ** 2
            expected = float((((self.y - model) / self.ey) ** 2).sum())
            self.assertAlmostEqual(self._chi2_at(a, b, c), expected, places=6,
                                   msg=f"at ({a}, {b}, {c})")

    def test_truth_fits_better_than_a_wrong_parameter_vector(self):
        self.assertLess(self._chi2_at(*self.truth), self._chi2_at(0.0, 0.0, 0.0))

    def test_sampling_the_graph_recovers_the_curve(self):
        """The three terms are collinear by construction, so the marginals
        are wide while the *prediction* is tight. The curve is what the fit
        actually determines, so that is what is asserted."""
        sampler = bff.Sampler("stretch", 7)
        sampler.set_parameter_ports(self.params)
        sampler.set_objective(self.objective, "chi2")
        sampler.set_number_of_walkers(16)
        sampler.run(3000)
        chain = np.asarray(sampler.chain).reshape(-1, 3)
        a, b, c = chain.mean(axis=0)
        predicted = c + a * self.x + b * self.x ** 2
        truth_curve = (self.truth[2] + self.truth[0] * self.x
                       + self.truth[1] * self.x ** 2)
        rms = float(np.sqrt(np.mean((predicted - truth_curve) ** 2)))
        self.assertLess(rms, 0.1, f"predicted curve is off by {rms:.4f}")
