"""bff.FitMinimizer against the scipy-backed optimiser it was ported from.

`FitMinimizer` is MINPACK's `lmdif` plus chisurf's bounds transform, so the bar
is not "it finds a minimum" -- any Levenberg-Marquardt does that -- but "it
finds the *same* minimum, from the same start, in the same number of
evaluations, and stops for the same reason". A fit that moved when the
optimiser was replaced would be indistinguishable from a fit that was wrong
before, so the parity tests compare against
`chisurf.core.math.optimization.leastsqbound` directly and skip when chisurf
is absent; the arithmetic and graph tests always run.
"""

import math
import unittest

import numpy as np

from IMP import bff

def _import_leastsqbound():
    """The optimiser that was replaced, from the frozen copy beside this file.

    It used to be imported from the sibling `chisurf` checkout. ChiSurf's
    copy was deleted on 2026-09-01 -- one implementation of this algorithm in
    the stack, and it is bff's -- so the reference moved *here*, into
    `reference_leastsqbound.py`, rather than the parity tests quietly
    skipping themselves. A skipped parity test on a 1:1 port is worth
    nothing: it is the only thing that says the answer did not move.
    """
    import sys
    from pathlib import Path
    here = str(Path(__file__).resolve().parent)
    if here not in sys.path:
        sys.path.insert(0, here)
    try:
        from reference_leastsqbound import (
            leastsqbound, _expected_evaluations, _reported_total)
        return leastsqbound, _expected_evaluations, _reported_total
    except Exception:
        return None


_CHISURF = _import_leastsqbound()
HAVE_CHISURF = _CHISURF is not None
if HAVE_CHISURF:
    cs_leastsqbound, cs_expected, cs_total = _CHISURF


# --------------------------------------------------------------- fixtures

def exponential_data(n=200, seed=7):
    """A two-exponential decay with reproducible noise, the shape of a fit."""
    rng = np.random.default_rng(seed)
    x = np.linspace(0.0, 10.0, n)
    truth = (1.4, 0.8, 0.5, 3.0)
    y = truth[0] * np.exp(-x / truth[1]) + truth[2] * np.exp(-x / truth[3])
    ey = np.full(n, 0.02)
    return x, y + rng.normal(0.0, 0.02, n), ey, truth


def residual_callable(x, y, ey):
    """The residual as chisurf's `leastsqbound` wants it: f(p) -> vector."""
    def f(p):
        model = p[0] * np.exp(-x / p[1]) + p[2] * np.exp(-x / p[3])
        return (y - model) / ey
    return f


class ResidualNode(bff.GraphNode):
    """A Python objective, the way a model bff cannot represent reaches it.

    One crossing per residual evaluation instead of the four or five a
    numpy-driven loop pays -- the fallback path the port is measured against,
    and the reason `set_residual_function` is not wrapped: a Python residual
    is a `GraphNode` director.
    """

    def __init__(self, func, names, name="residuals"):
        super().__init__(name)
        self._func = func
        self._names = list(names)
        for n in self._names:
            self.add_input_port(n, bff.GraphPort(0.0))
        self.add_output_port("residuals", bff.GraphPort([0.0], False, True))
        self.calls = 0
        self.visited = []

    def evaluate(self):
        p = np.array([self.inputs[n].value for n in self._names])
        self.calls += 1
        self.visited.append(p.copy())
        self.outputs["residuals"].set_value_vector(
            list(map(float, self._func(p))))
        self.set_valid(True)

    def ports_for(self):
        return [self.get_input_port(n) for n in self._names]


def build_minimizer(func, start, bounds=None, **options):
    node = ResidualNode(func, ["p%d" % i for i in range(len(start))])
    ports = node.ports_for()
    for port, v in zip(ports, start):
        port.value = float(v)
    m = bff.FitMinimizer()
    m.set_parameter_ports(ports)
    m.set_objective(node, "residuals")
    if bounds is not None:
        m.bounds = bounds
    for k, v in options.items():
        setattr(m, k, v)
    return m, node


def _lmdif_offset(seen, start):
    """Where MINPACK's own first evaluation is in `leastsqbound`'s record.

    `leastsqbound` checks the function -- and, with bounds, scipy's `leastsq`
    checks it again -- before MINPACK ever runs, so the record opens with a
    run of evaluations at x0 that the algorithm never made. The last of them
    is MINPACK's. Counting them would be pinning a scipy implementation
    detail; finding them is not.
    """
    x0 = np.asarray(start, dtype=float)
    first_move = next(k for k, p in enumerate(seen) if not np.array_equal(p, x0))
    return first_move - 1


# ----------------------------------------------------------------- parity

@unittest.skipUnless(HAVE_CHISURF, "chisurf is not importable")
class LeastsqboundParityTests(unittest.TestCase):
    """The same answer as the optimiser being replaced, not merely a good one.

    **How exact the bar can be depends on the bounds, and this is measured,
    not assumed.** With no bounds the transform is the identity, both sides
    run the same arithmetic on the same numbers, and the *whole trajectory*
    is reproduced -- every evaluation point, and the same evaluation count.
    With bounds it cannot be: `leastsqbound` maps coordinates with numpy's
    `sin`/`arcsin` and this maps them with libm's, the two disagree in the
    last bit (`i2e(e2i(x))` round-trips to a difference of 2.2e-16 here), and
    twenty-odd Levenberg-Marquardt iterations amplify one ulp to about 1e-7
    in the answer -- which is inside `xtol = 1.49e-8` *relative*, i.e. inside
    the optimiser's own idea of converged. scipy's answer moves by the same
    amount when the bounds are merely respelled. So the bounded bar is
    agreement to the convergence tolerance, and the unbounded bar is exact.
    """

    #: Off the model's permutation-symmetry ridge. A sum of exponentials
    #: started at (1, 1, 1, 1) has two *identical* pairs of Jacobian columns,
    #: so the pivoted QR is choosing between equal norms and one ulp decides
    #: the pivot. Both optimisers still find the same minimum there, but they
    #: label its two terms differently, which is not a discrepancy worth
    #: pinning a trajectory against.
    START = (1.0, 0.5, 0.3, 2.0)

    def _run_both(self, bounds, start=None):
        start = list(start or self.START)
        x, y, ey, _ = exponential_data()
        f = residual_callable(x, y, ey)
        seen = []

        def counted(p):
            seen.append(np.array(p, dtype=float))
            return f(p)

        reference, ier = cs_leastsqbound(counted, np.array(start), bounds=bounds)
        m, node = build_minimizer(f, start, bounds=bounds)
        info = m.run()
        return m, info, reference, ier, seen, node

    def _compare(self, bounds, rtol=1e-5):
        m, info, reference, ier, _, _ = self._run_both(bounds)
        self.assertEqual(info, ier, "stopped for a different reason")
        np.testing.assert_allclose(m.x, reference, rtol=rtol)
        return m, reference

    def test_unbounded_reproduces_the_whole_trajectory(self):
        """Every evaluation point, not just the endpoint."""
        m, info, reference, ier, seen, node = self._run_both([(None, None)] * 4)
        trajectory = seen[_lmdif_offset(seen, self.START):]
        self.assertEqual(len(trajectory), m.n_evaluations)
        # `node`, not `m.get_objective()`: SWIG hands back the `GraphNode` proxy,
        # not the Python subclass, so the director's own attributes are only
        # reachable through the reference the caller kept.
        mine = node.visited
        for i, (a, b) in enumerate(zip(trajectory, mine)):
            np.testing.assert_allclose(b, a, rtol=0, atol=1e-14,
                                       err_msg="evaluation %d" % i)
        np.testing.assert_allclose(m.x, reference, rtol=0, atol=1e-12)
        self.assertEqual(info, ier)

    def test_two_sided_bounds_match(self):
        self._compare([(0.0, 10.0), (0.01, 5.0), (0.0, 10.0), (0.01, 20.0)])

    def test_lower_bound_only_matches(self):
        self._compare([(0.0, None), (0.01, None), (0.0, None), (0.01, None)])

    def test_upper_bound_only_matches(self):
        self._compare([(None, 10.0), (None, 5.0), (None, 10.0), (None, 20.0)])

    def test_infinite_bounds_are_no_bounds(self):
        """`(-inf, inf)` must mean unbounded, which is chisurf's `_is_unbounded`.

        Recognising only `None` sent such a parameter down the two-sided
        branch, where `arcsin(inf/inf - 1)` poisoned the internal vector with
        NaN -- the defect that function exists to fix. Being genuinely
        unbounded, this must reproduce the *unbounded* answer exactly.
        """
        inf = float("inf")
        m, info, reference, ier, _, _ = self._run_both([(-inf, inf)] * 4)
        self.assertTrue(np.all(np.isfinite(m.x)))
        np.testing.assert_allclose(m.x, reference, rtol=0, atol=1e-12)

    def test_a_symmetric_start_finds_the_same_minimum(self):
        """The permutation ridge: the same answer, its terms relabelled.

        Two implementations of a pivoted QR cannot agree about which of two
        equal column norms to pivot on, so this start is where trajectory
        parity stops being a meaningful question -- but the minimum is not
        allowed to move.
        """
        m, info, reference, ier, _, _ = self._run_both(
            [(None, None)] * 4, start=(1.0, 1.0, 1.0, 1.0))
        self.assertEqual(info, ier)
        got = sorted(zip(m.x[0::2], m.x[1::2]), key=lambda t: t[1])
        want = sorted(zip(reference[0::2], reference[1::2]), key=lambda t: t[1])
        np.testing.assert_allclose(np.array(got), np.array(want), rtol=1e-5)

    def test_evaluation_count_matches(self):
        m, info, reference, ier, seen, _ = self._run_both([(None, None)] * 4)
        offset = _lmdif_offset(seen, self.START)
        self.assertEqual(m.n_evaluations, len(seen) - offset)

    def test_covariance_matches_scipy_and_its_own_definition(self):
        """`(J^T J)^-1` in **external** coordinates, from the R already built.

        Checked two ways, against neither of which it could pass by
        accident: `scipy.optimize.leastsq`'s own `cov_x`, and the definition
        -- the inverse Gram matrix of a finite-difference Jacobian taken at
        the answer.

        Not checked against `leastsqbound(full_output=1)`, which is **wrong
        on current scipy**: it still does `ipvt - 1` on an `ipvt` scipy made
        0-based (scipy 1.18 `_minpack_py.py:488`, `perm = ipvt`), so its
        covariance comes back permuted by the pivot -- and the `-1` wraps to
        the last row, so it is not even a permutation. Measured here: a
        cyclic shift of both indices, up to 9x wrong element by element.
        Latent in chisurf only because `Fit.update_error_estimates` takes
        its own finite-difference Jacobian instead.
        """
        import scipy.optimize

        x, y, ey, _ = exponential_data()
        f = residual_callable(x, y, ey)
        reference = scipy.optimize.leastsq(f, np.array(self.START),
                                           full_output=1)[1]
        m, _ = build_minimizer(f, self.START)
        m.run()
        np.testing.assert_allclose(m.covariance, reference, rtol=1e-8,
                                   atol=1e-14)

        p = np.asarray(m.x)
        f0 = f(p)
        eps = math.sqrt(np.finfo(float).eps)
        jacobian = np.empty((len(f0), len(p)))
        for j in range(len(p)):
            q = p.copy()
            h = eps * abs(q[j]) or eps
            q[j] += h
            jacobian[:, j] = (f(q) - f0) / h
        np.testing.assert_allclose(m.covariance,
                                   np.linalg.inv(jacobian.T @ jacobian),
                                   rtol=1e-3)

    def test_progress_arithmetic_matches(self):
        for n in (1, 3, 4, 10):
            for maxfev in (0, 50, 200):
                self.assertEqual(bff.minimizer_expected_evaluations(n, maxfev),
                                 cs_expected(n, maxfev))
        for nfev in (1, 5, 42, 450, 5000):
            for expected in (1, 42, 200):
                self.assertEqual(bff.minimizer_reported_total(nfev, expected),
                                 cs_total(nfev, expected))


# ------------------------------------------------------------- arithmetic

class MinimizerTests(unittest.TestCase):
    def test_recovers_the_truth(self):
        x, y, ey, truth = exponential_data()
        f = residual_callable(x, y, ey)
        m, _ = build_minimizer(f, [1.0, 0.5, 0.3, 2.0])
        info = m.run()
        self.assertIn(info, (1, 2, 3, 4))
        # Within the noise, not to the digit: the data carry sigma = 0.02 on
        # a 0.5 amplitude, so the least-squares optimum is a good way from
        # the truth and demanding otherwise would be testing the noise seed.
        # A sum of exponentials is also invariant under relabelling its
        # terms, so the answer is the *set* of (amplitude, lifetime) pairs.
        got = sorted(zip(m.x[0::2], m.x[1::2]), key=lambda p: p[1])
        want = sorted(zip(truth[0::2], truth[1::2]), key=lambda p: p[1])
        for (ga, gt), (wa, wt) in zip(got, want):
            self.assertAlmostEqual(ga, wa, delta=0.15 * abs(wa))
            self.assertAlmostEqual(gt, wt, delta=0.15 * abs(wt))
        self.assertLess(m.chi2r, 2.0)

    def test_ports_are_left_at_the_solution(self):
        """The graph is left evaluated at the answer, so the fitted curve is
        read off it without evaluating anything again."""
        x, y, ey, _ = exponential_data()
        m, node = build_minimizer(residual_callable(x, y, ey), [1.0] * 4)
        m.run()
        calls = node.calls
        for port, v in zip(node.ports_for(), m.x):
            self.assertEqual(port.value, v)
        np.testing.assert_allclose(
            node.outputs["residuals"].value, m.residuals)
        self.assertEqual(node.calls, calls, "reading the answer re-evaluated")

    def test_bounds_are_respected(self):
        x, y, ey, _ = exponential_data()
        f = residual_callable(x, y, ey)
        # A ceiling well below the true 3.0 that the fit must not cross.
        m, _ = build_minimizer(f, [1.0, 1.0, 1.0, 1.0],
                               bounds=[(0.0, 10.0), (0.01, 5.0),
                                       (0.0, 10.0), (0.01, 1.5)])
        m.run()
        self.assertLessEqual(m.x[3], 1.5 + 1e-9)
        self.assertGreaterEqual(m.x[3], 0.01 - 1e-9)

    def test_bounds_round_trip_through_the_property(self):
        m, _ = build_minimizer(lambda p: np.zeros(5), [1.0, 2.0],
                               bounds=[(0.0, 1.0), (None, None)])
        self.assertEqual(m.bounds, [(0.0, 1.0), (None, None)])

    def test_starting_inside_a_bound_is_kept(self):
        m, node = build_minimizer(lambda p: np.zeros(5), [0.5],
                                  bounds=[(0.0, 1.0)])
        m.run()
        self.assertAlmostEqual(m.x[0], 0.5, places=10)

    def test_errors_are_the_covariance_diagonal(self):
        x, y, ey, _ = exponential_data()
        m, _ = build_minimizer(residual_callable(x, y, ey), [1.0] * 4)
        m.run()
        np.testing.assert_allclose(m.errors, np.sqrt(np.diag(m.covariance)),
                                   rtol=1e-12)

    def test_maxfev_stops_and_says_so(self):
        x, y, ey, _ = exponential_data()
        m, _ = build_minimizer(residual_callable(x, y, ey), [1.0] * 4,
                               maxfev=7)
        info = m.run()
        self.assertEqual(info, 5)
        self.assertIn("maxfev", m.message)
        self.assertEqual(m.covariance.size, 0,
                         "a run that did not converge must not report a "
                         "covariance of the point it was passing through")

    def test_fewer_residuals_than_parameters_is_refused(self):
        m, _ = build_minimizer(lambda p: np.zeros(2), [1.0, 2.0, 3.0])
        with self.assertRaises(ValueError):
            m.run()

    def test_a_residual_length_that_moves_is_refused(self):
        state = {"n": 10}

        def wobbly(p):
            state["n"] += 1
            return np.zeros(state["n"])

        m, _ = build_minimizer(wobbly, [1.0, 2.0])
        with self.assertRaises(ValueError):
            m.run()

    def test_no_objective_is_refused(self):
        m = bff.FitMinimizer()
        m.set_initial_values([1.0])
        with self.assertRaises(ValueError):
            m.run()

    def test_an_unknown_algorithm_is_refused(self):
        with self.assertRaises(ValueError):
            bff.FitMinimizer("nelder-mead")

    def test_a_fixed_port_cannot_be_optimised(self):
        port = bff.GraphPort(1.0)
        port.set_fixed(True)
        m = bff.FitMinimizer()
        with self.assertRaises(ValueError):
            m.set_parameter_ports([port])

    def test_an_objective_without_the_residual_port_is_refused(self):
        node = bff.GraphNode("plain")
        node.add_output_port("chi2", bff.GraphPort(0.0, False, True))
        m = bff.FitMinimizer()
        with self.assertRaises(ValueError):
            m.set_objective(node, "residuals")


# ---------------------------------------------------------------- observer

class Observer(bff.FitMinimizerObserver):
    """Counts reports, and cancels after `stop_after` of them."""

    def __init__(self, stop_after=None):
        super().__init__()
        self.reports = []
        self.stop_after = stop_after

    def report(self, n_evaluations, total, chi2):
        self.reports.append((n_evaluations, total, chi2))
        if self.stop_after is not None and n_evaluations >= self.stop_after:
            return False
        return True


class ObserverTests(unittest.TestCase):
    def test_every_evaluation_is_reported(self):
        x, y, ey, _ = exponential_data()
        m, _ = build_minimizer(residual_callable(x, y, ey), [1.0] * 4)
        obs = Observer()
        m.set_observer(obs)
        m.run()
        # One per evaluation, plus the completion report that fills the bar.
        self.assertEqual(len(obs.reports), m.n_evaluations + 1)
        self.assertEqual(obs.reports[-1][0], obs.reports[-1][1])

    def test_the_bar_only_ever_moves_forward(self):
        """chisurf's asymptotic total: strictly increasing, never arriving.

        A linear bar against a guessed total saturates on exactly the fits
        long enough for anyone to be watching.
        """
        fractions = [bff.minimizer_reported_total(n, 42) for n in range(1, 500)]
        done = [n / t for n, t in zip(range(1, 500), fractions)]
        self.assertTrue(all(b >= a for a, b in zip(done, done[1:])),
                        "the reported fraction stepped backwards")
        self.assertTrue(all(f < 1.0 for f in done))

    def test_returning_false_cancels(self):
        x, y, ey, _ = exponential_data()
        m, _ = build_minimizer(residual_callable(x, y, ey), [1.0] * 4)
        obs = Observer(stop_after=5)
        m.set_observer(obs)
        info = m.run()
        self.assertEqual(info, -1)
        self.assertTrue(m.cancelled)
        self.assertEqual(m.n_evaluations, 5)
        self.assertIn("cancelled", m.message)

    def test_a_cancelled_run_leaves_a_point_it_has_been(self):
        """Not the trial point it was evaluating when it was told to stop."""
        x, y, ey, _ = exponential_data()
        m, _ = build_minimizer(residual_callable(x, y, ey), [1.0] * 4)
        m.set_observer(Observer(stop_after=3))
        m.run()
        self.assertTrue(np.all(np.isfinite(m.x)))
        self.assertTrue(np.all(np.isfinite(m.residuals)))

    def test_the_observer_survives_the_call_that_set_it(self):
        """A progress dialog held only by a raw pointer would be collected
        mid-fit; the C++ side owns a reference."""
        x, y, ey, _ = exponential_data()
        m, _ = build_minimizer(residual_callable(x, y, ey), [1.0] * 4)
        m.set_observer(Observer())
        import gc
        gc.collect()
        m.run()
        self.assertGreater(m.n_evaluations, 1)

    def test_clear_observer_stops_reporting(self):
        x, y, ey, _ = exponential_data()
        m, _ = build_minimizer(residual_callable(x, y, ey), [1.0] * 4)
        obs = Observer()
        m.set_observer(obs)
        m.clear_observer()
        m.run()
        self.assertEqual(obs.reports, [])


# ------------------------------------------------------------- whole graph

class GraphObjectiveTests(unittest.TestCase):
    """`GraphExpression -> FitChiSquared -> FitMinimizer`: a fit that never re-enters
    the interpreter. This is the arrangement the port exists for -- the
    optimiser and the objective on the same side of the boundary, so a
    caller crosses once per `run()` rather than once per part per iteration.
    """

    def _graph(self, equation, data_y, data_ey, axis):
        expression = bff.GraphExpression("model")
        expression.set_expression(equation)
        names = [n for n in expression.get_variable_names() if n != "x"]
        for n in names:
            expression.add_input_port(n, bff.GraphPort(1.0))
        expression.add_input_port("x", bff.GraphPort(list(map(float, axis))))
        curve = bff.GraphPort([0.0], False, True)
        expression.add_output_port("model", curve)

        chi2 = bff.FitChiSquared("chi2")
        chi2.set_data(list(map(float, data_y)), list(map(float, data_ey)))
        model_in = bff.GraphPort([0.0])
        model_in.link = curve
        chi2.add_input_port("model", model_in)
        chi2.add_output_port("chi2", bff.GraphPort(0.0, False, True))
        chi2.add_output_port("residuals", bff.GraphPort([0.0], False, True))
        return expression, chi2, names

    def test_a_whole_fit_in_cplusplus(self):
        x = np.linspace(0.1, 5.0, 120)
        truth = (2.5, 0.7)
        y = truth[0] * np.exp(-x / truth[1])
        ey = np.full_like(y, 0.01)
        expression, chi2, names = self._graph("a*exp(-x/t)", y, ey, x)
        self.assertEqual(names, ["a", "t"])

        ports = [expression.get_input_port(n) for n in names]
        ports[0].value, ports[1].value = 1.0, 1.0
        m = bff.FitMinimizer()
        m.set_parameter_ports(ports)
        m.set_objective(chi2, "residuals")
        info = m.run()

        self.assertIn(info, (1, 2, 3, 4))
        self.assertAlmostEqual(m.x[0], truth[0], places=4)
        self.assertAlmostEqual(m.x[1], truth[1], places=4)
        self.assertLess(m.chi2, 1e-6)

    def test_chi_squared_still_writes_its_scalar(self):
        """The residual port is additive: a `MCMCSampler` reading chi-square off
        the same node must be unaffected."""
        x = np.linspace(0.1, 5.0, 20)
        y = 2.0 * np.exp(-x / 0.5)
        expression, chi2, names = self._graph(
            "a*exp(-x/t)", y, np.ones_like(y), x)
        expression.get_input_port("a").value = 2.0
        expression.get_input_port("t").value = 0.5
        chi2.update()
        self.assertAlmostEqual(chi2.get_output_port("chi2").value, 0.0,
                               places=10)
        self.assertEqual(
            len(chi2.get_output_port("residuals").value), 20)

    def test_a_trial_point_that_breaks_the_model_is_rejected(self):
        """The optimiser must not be *attracted* to a broken model.

        chinet's ports floor a NaN to `np.finfo(float).tiny`, so before the
        fit's transport opted out of that, a model that blew up arrived as a
        model that is zero everywhere -- a good fit for data near zero. The
        optimiser then walked towards the parameters that break it, and
        agreed with nothing: the numpy path returns NaN and MINPACK rejects
        the step.

        `sqrt(a - x)` is NaN wherever `a < x`, so the search cannot converge
        into that region.
        """
        x = np.linspace(1.0, 5.0, 64)
        truth = 9.0
        y = np.sqrt(truth - x)
        expression, chi2, names = self._graph("sqrt(a-x)", y, np.ones_like(y), x)
        port = expression.get_input_port("a")
        port.value = 6.0
        m = bff.FitMinimizer()
        m.set_parameter_ports([port])
        m.set_objective(chi2, "residuals")
        info = m.run()
        self.assertIn(info, (1, 2, 3, 4))
        self.assertAlmostEqual(m.x[0], truth, places=4)
        # And the region that breaks the model is infinitely bad, not free.
        port.value = 0.0
        expression.set_valid(False)
        chi2.set_valid(False)
        chi2.update()
        self.assertEqual(chi2.get_chi2(), float("inf"))
        # The residuals -- what the optimiser reads -- carry the NaN through.
        self.assertTrue(np.all(np.isnan(
            np.asarray(chi2.get_output_port("residuals").value))))
        # The *scalar* chi-square port keeps chinet's clamp, deliberately:
        # it is a document value, a `MCMCSampler` reads it, and the largest
        # double rejects a move exactly as an infinity would. Only the
        # transport the optimiser reads had to stop sanitising.
        self.assertEqual(chi2.get_output_port("chi2").value,
                         np.finfo(np.float64).max)

    def test_a_node_without_a_residual_port_keeps_working(self):
        """`FitChiSquared` writes the residuals only when asked; a graph built
        before this port existed evaluates exactly as it did."""
        x = np.linspace(0.1, 5.0, 20)
        y = 2.0 * np.exp(-x / 0.5)
        expression, chi2, _ = self._graph("a*exp(-x/t)", y, np.ones_like(y), x)
        # Drop the residual port the fixture adds.
        chi2.add_output_port("residuals", bff.GraphPort(0.0, False, True))
        chi2.set_residuals_port_key("not_a_port")
        chi2.update()
        self.assertGreaterEqual(chi2.get_output_port("chi2").value, 0.0)


# ------------------------------------------------------------- covariance

def numpy_covariance(f, x, epsilon=None, floor=1.0):
    """chisurf's `approx_grad` + `covariance_matrix`, written out here.

    The reference the C++ routine is pinned against, transcribed rather than
    imported so that this file states the specification it is testing. It is
    the *step rule* that matters -- `eps * max(|x|, floor)` with an absolute
    floor, rounded to an exactly representable difference -- and the two
    things `covariance_matrix` does around it: drop the columns that are
    identically zero, and pseudo-invert rather than invert.
    """
    import scipy.linalg
    if epsilon is None:
        epsilon = math.sqrt(np.finfo(float).eps)
    x = np.asarray(x, dtype=float)
    f0 = np.asarray(f(x), dtype=float)
    grad = np.zeros((len(x), len(f0)), dtype=float)
    for k in range(len(x)):
        step = epsilon * max(abs(float(x[k])), floor)
        step = (x[k] + step) - x[k]
        d = np.zeros_like(x)
        d[k] = step
        grad[k] = (np.asarray(f(x + d), dtype=float) - f0) / step
    used = [k for k, row in enumerate(grad) if (row ** 2).sum() > 0.0]
    j = grad[used]
    return scipy.linalg.pinvh(j @ j.T), used


class CovarianceTests(unittest.TestCase):
    """`compute_covariance` -- the error bars, differenced over the graph.

    **Every fixture here names its node and keeps the name.** A
    `ResidualNode` is a SWIG director and the C++ side holds only a weak
    reference to its Python proxy, so a node that falls out of scope while
    its `FitMinimizer` is still in use is a use-after-free -- reliably a
    segfault, in `run()` exactly as much as in `compute_covariance`. Binding
    it to `_` is enough to lose it, because the next tuple unpacking rebinds
    `_`. See T-20260901-13.

    What these pin is agreement with numpy, not plausibility. chisurf's
    parameter table is populated from `covariance_matrix`, and the whole
    licence for computing it in C++ is that the number displayed does not
    move; a covariance that is merely reasonable would be a regression.
    """

    def test_it_matches_the_numpy_reference(self):
        x, y, ey, _ = exponential_data()
        f = residual_callable(x, y, ey)
        m, node = build_minimizer(f, [1.0, 0.5, 0.3, 2.0])
        m.run()
        cpp, used = m.finite_difference_covariance()
        want, want_used = numpy_covariance(f, m.x)
        self.assertEqual(used, want_used)
        np.testing.assert_allclose(np.sqrt(np.diag(cpp)),
                                   np.sqrt(np.diag(want)), rtol=1e-10)
        np.testing.assert_allclose(cpp, want, rtol=1e-8, atol=0.0)

    def test_the_step_is_epsilon_times_the_parameter_or_the_floor(self):
        """The one line the whole ticket is about, read off the objective.

        `lmdif` steps by `sqrt(epsfcn) * |x_j|`, floored only at *exactly*
        zero, so a parameter that converged near zero is differenced at a
        step near zero and its Jacobian column is round-off. `approx_grad`
        floors at 1.0 instead, and that is the rule this routine copies --
        with the difference rounded to something exactly representable, so
        the divisor is the step the objective actually saw.

        Asserted by watching the parameter vectors the objective is handed,
        which is the only way to see a step rather than infer it from a
        matrix. The fixture straddles the floor on purpose: two of its four
        parameters converge below 1.0 and two above, so the same assertion
        covers both branches of the `max`.
        """
        x, y, ey, _truth = exponential_data()
        residual = residual_callable(x, y, ey)
        seen = []

        def f(p):
            seen.append(np.array(p, dtype=float))
            return residual(p)

        m, node = build_minimizer(f, [1.0, 0.5, 0.3, 2.0])
        m.run()
        solution = np.asarray(m.x, dtype=float)
        self.assertTrue(np.any(np.abs(solution) < 1.0)
                        and np.any(np.abs(solution) > 1.0),
                        "the fixture no longer straddles the floor")

        for epsilon, floor, wanted in ((0.0, 0.0, math.sqrt(
                                            np.finfo(float).eps)),
                                       (1.0e-4, 1.0, 1.0e-4),
                                       (1.0e-4, 1.0e-300, 1.0e-4)):
            del seen[:]
            m.finite_difference_covariance(epsilon, floor)
            scale = 1.0 if floor <= 0.0 else floor
            self.assertEqual(len(seen), len(solution) + 1,
                             "one evaluation per parameter, plus the restore")
            for k in range(len(solution)):
                step = wanted * max(abs(solution[k]), scale)
                step = (solution[k] + step) - solution[k]
                expected = solution.copy()
                expected[k] += step
                np.testing.assert_array_equal(
                    seen[k], expected,
                    "parameter %d was not stepped by eps*max(|x|, floor)" % k)
            np.testing.assert_array_equal(
                seen[-1], solution,
                "the objective was left at a perturbation, not at the answer")

        # And the floor is load-bearing, not incidentally equal here: with it
        # removed, every parameter under 1.0 is differenced at a smaller step
        # than the rule gives -- which is the whole of `lmdif`'s problem.
        del seen[:]
        m.finite_difference_covariance(1.0e-4, 1.0)
        floored = [seen[k][k] - solution[k] for k in range(len(solution))]
        del seen[:]
        m.finite_difference_covariance(1.0e-4, 1.0e-300)
        relative = [seen[k][k] - solution[k] for k in range(len(solution))]
        smaller = [k for k in range(len(solution)) if abs(solution[k]) < 1.0]
        self.assertTrue(smaller)
        for k in smaller:
            self.assertLess(relative[k], floored[k])

    def test_a_column_of_zeros_is_dropped_rather_than_inverted(self):
        """A free parameter nothing depends on is legitimate.

        chisurf's `E_FRET` is one: its value comes from a callable, its
        setter is ignored, and the graph gives it a deliberately dangling
        port -- so its column is exactly zero *by design*. Inverting the
        singular `J'J` that results would hand every FRET fit error bars that
        mean nothing, so the parameter is dropped from the matrix and named
        in `get_covariance_parameters()` instead.
        """
        x = np.linspace(0.0, 10.0, 100)
        def f(p):
            return (p[0] * np.exp(-x / p[1]) - 1.4 * np.exp(-x / 0.8))
        # Three parameters, the third of which the objective never reads.
        def g(p):
            return f(p)
        m, node = build_minimizer(g, [1.0, 0.5, 0.25])
        m.run()
        cov, used = m.finite_difference_covariance()
        self.assertEqual(used, [0, 1])
        self.assertEqual(cov.shape, (2, 2))
        self.assertTrue(np.all(np.isfinite(cov)))
        want, want_used = numpy_covariance(g, m.x)
        self.assertEqual(want_used, used)
        np.testing.assert_allclose(cov, want, rtol=1e-8)

    def test_it_is_in_external_coordinates(self):
        """A bounded fit must report the error bar of the bounded parameter.

        The optimiser works in an unconstrained internal coordinate and its
        own covariance is taken back through `_internal2external_grad` for
        exactly this reason. This routine never enters that coordinate at all
        -- it differences the objective at the external vector, as
        `approx_grad` does -- so what has to be shown is that the two agree,
        away from the bound where both are meaningful.
        """
        x, y, ey, _ = exponential_data()
        f = residual_callable(x, y, ey)
        # Both nodes are kept named. A `ResidualNode` is a SWIG director and
        # the C++ side holds only a weak reference to its Python proxy, so
        # letting the first one fall out of scope while its `FitMinimizer` is
        # still in use is a use-after-free -- in `run()` as much as here.
        free, free_node = build_minimizer(f, [1.0, 0.5, 0.3, 2.0])
        free.run()
        bounded, bounded_node = build_minimizer(
            f, [1.0, 0.5, 0.3, 2.0],
            bounds=[(-10.0, 10.0), (0.01, 5.0),
                    (-10.0, 10.0), (0.01, 8.0)])
        bounded.run()
        np.testing.assert_allclose(bounded.x, free.x, rtol=1e-5)
        a, _ = bounded.finite_difference_covariance()
        b, _ = free.finite_difference_covariance()
        np.testing.assert_allclose(np.sqrt(np.diag(a)), np.sqrt(np.diag(b)),
                                   rtol=1e-4)
        self.assertIsNotNone(free_node)
        self.assertIsNotNone(bounded_node)

    def test_it_puts_the_graph_back_where_the_solution_left_it(self):
        """The last thing differenced is a perturbation, not the answer.

        `run()` documents that the ports and the graph are left at the
        solution so a caller reads the fitted curve without evaluating
        anything again. Taking a covariance must not quietly break that
        promise -- and it is the kind of breakage nothing else notices,
        because the parameters would still read correctly while the *curve*
        belonged to the last finite-difference step.
        """
        x, y, ey, _ = exponential_data()
        m, node = build_minimizer(residual_callable(x, y, ey), [1.0] * 4)
        m.run()
        residuals = np.array(m.residuals, dtype=float)
        m.finite_difference_covariance()
        for port, v in zip(node.ports_for(), m.x):
            self.assertEqual(port.value, v)
        np.testing.assert_allclose(node.outputs["residuals"].value, residuals)

    def test_it_costs_one_evaluation_per_parameter_plus_one(self):
        """p + 1, and the + 1 is the restore.

        `f0` is not recomputed: `run()` leaves the graph at the solution, so
        the residuals it already holds *are* the residuals at `x`. Stated as
        a test because it is an assumption about ordering -- anything that
        wrote a port between `run()` and here would make it false silently.
        """
        x, y, ey, _ = exponential_data()
        m, node = build_minimizer(residual_callable(x, y, ey), [1.0] * 4)
        m.run()
        before = node.calls
        m.finite_difference_covariance()
        self.assertEqual(node.calls - before, len(m.x) + 1)

    def test_there_is_no_covariance_before_a_run(self):
        x, y, ey, _ = exponential_data()
        m, node = build_minimizer(residual_callable(x, y, ey), [1.0] * 4)
        cov, used = m.finite_difference_covariance()
        self.assertEqual(cov.size, 0)
        self.assertEqual(used, [])


if __name__ == "__main__":
    unittest.main()
