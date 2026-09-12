"""The grouping: one misfit over several datasets, in C++.

A joint fit is two things. **Sharing a parameter** is `GraphPort::set_link`, which
bff has had since the GraphPort runtime; **one objective over every dataset** is
`JointChiSquared`. Together they make a group a graph a `Minimizer` can
optimise without the caller crossing the SWIG boundary per iteration.

The tests that matter are the ones separating a joint fit from N separate
fits: a shared parameter has to come out at the *compromise* the datasets
force, not at either dataset's own answer, and it has to be reached in one
optimisation rather than several.
"""

import math
import unittest

import numpy as np

from IMP import bff


def make_member(equation, x, y, ey, name):
    """One dataset: `GraphExpression -> ChiSquared`, residuals on a port."""
    curve = bff.GraphExpression(name + "_model")
    curve.set_expression(equation)
    parameters = {}
    for v in curve.get_variable_names():
        if v == "x":
            continue
        port = bff.GraphPort(1.0)
        curve.add_input_port(v, port)
        parameters[v] = port
    axis = bff.GraphPort([0.0])
    axis.set_values_array(np.ascontiguousarray(x, dtype=float))
    curve.add_input_port("x", axis)
    out = bff.GraphPort([0.0], False, True)
    curve.add_output_port(name + "_model", out)

    chi2 = bff.ChiSquared(name)
    chi2.set_data_arrays(np.ascontiguousarray(y, dtype=float),
                         np.ascontiguousarray(ey, dtype=float))
    model_in = bff.GraphPort([0.0])
    model_in.link = out
    chi2.add_input_port("model", model_in)
    chi2.add_output_port(name, bff.GraphPort(0.0, False, True))
    chi2.add_output_port("residuals", bff.GraphPort([0.0], False, True))
    # The Python proxies have to outlive this function or the C++ graph is
    # left holding nodes whose directors have been collected.
    chi2._graph = (curve, out, model_in, axis)
    return chi2, parameters


def two_datasets(t_shared=3.0, a1=2.0, a2=5.0, n=256, sigma=0.01, seed=11):
    """Two decays with a **shared lifetime** and their own amplitudes."""
    rng = np.random.default_rng(seed)
    x = np.linspace(0.1, 15.0, n)
    y1 = a1 * np.exp(-x / t_shared) + rng.normal(0, sigma, n)
    y2 = a2 * np.exp(-x / t_shared) + rng.normal(0, sigma, n)
    ey = np.full(n, sigma)
    return x, y1, y2, ey


def build_group(x, y1, y2, ey):
    """Two members, one shared `t`, one joint objective."""
    m1, p1 = make_member("a*exp(-x/t)", x, y1, ey, "d1")
    m2, p2 = make_member("a*exp(-x/t)", x, y2, ey, "d2")
    # The sharing: dataset 2's `t` follows dataset 1's. This is the whole
    # coupling, and it is a port link -- already on the C++ side.
    p2["t"].link = p1["t"]

    joint = bff.JointChiSquared("joint")
    joint.add_output_port("joint", bff.GraphPort(0.0, False, True))
    joint.add_output_port("residuals", bff.GraphPort([0.0], False, True))
    joint.add_member(m1, "residuals")
    joint.add_member(m2, "residuals")

    free = [p1["a"], p2["a"], p1["t"]]
    for port, v in zip(free, (1.0, 1.0, 1.0)):
        port.value = v
    m = bff.Minimizer()
    m.set_parameter_ports(free)
    m.set_objective(joint, "residuals")
    m._graph = (m1, m2, joint, p1, p2)
    return m, joint, free, (m1, m2)


class GroupingTests(unittest.TestCase):
    def test_a_group_fits_a_shared_parameter(self):
        x, y1, y2, ey = two_datasets()
        m, joint, free, _ = build_group(x, y1, y2, ey)
        info = m.run()
        self.assertIn(info, (1, 2, 3, 4))
        a1, a2, t = m.x
        self.assertAlmostEqual(a1, 2.0, delta=0.02)
        self.assertAlmostEqual(a2, 5.0, delta=0.02)
        self.assertAlmostEqual(t, 3.0, delta=0.02)

    def test_the_shared_parameter_is_genuinely_one_number(self):
        """Not two that happen to agree: the second dataset's port is a
        follower, so it holds whatever the first holds."""
        x, y1, y2, ey = two_datasets()
        m, joint, free, (m1, m2) = build_group(x, y1, y2, ey)
        m.run()
        d1_t = m1._graph[0].get_input_port("t")
        d2_t = m2._graph[0].get_input_port("t")
        self.assertEqual(d1_t.value, d2_t.value)
        self.assertEqual(d2_t.value, m.x[2])

    def test_the_joint_residual_is_the_blocks_end_to_end(self):
        x, y1, y2, ey = two_datasets()
        m, joint, free, (m1, m2) = build_group(x, y1, y2, ey)
        m.run()
        sizes = joint.get_block_sizes()
        self.assertEqual(list(sizes), [len(x), len(x)])
        self.assertEqual(list(joint.get_block_offsets()), [0, len(x)])
        residuals = np.asarray(m.residuals)
        self.assertEqual(residuals.size, 2 * len(x))
        np.testing.assert_allclose(residuals[:len(x)],
                                   m1.get_weighted_residuals(), rtol=0,
                                   atol=0)
        np.testing.assert_allclose(residuals[len(x):],
                                   m2.get_weighted_residuals(), rtol=0,
                                   atol=0)

    def test_chi2_is_the_sum_of_the_members(self):
        x, y1, y2, ey = two_datasets()
        m, joint, free, (m1, m2) = build_group(x, y1, y2, ey)
        m.run()
        self.assertAlmostEqual(joint.get_chi2(),
                               m1.get_chi2() + m2.get_chi2(), places=9)
        self.assertAlmostEqual(m.chi2, joint.get_chi2(), places=9)

    def test_the_group_is_not_two_separate_fits(self):
        """The point of a joint fit, stated as a test.

        Give the two datasets *different* true lifetimes. Fitted separately
        each recovers its own; fitted jointly the shared parameter must land
        between them -- at the inverse-variance compromise, not at either.
        A group that quietly optimised its members one at a time would give
        one of the two, and would pass every other test here.
        """
        n, sigma = 256, 0.01
        rng = np.random.default_rng(5)
        x = np.linspace(0.1, 15.0, n)
        ey = np.full(n, sigma)
        y1 = 2.0 * np.exp(-x / 2.6) + rng.normal(0, sigma, n)
        y2 = 2.0 * np.exp(-x / 3.4) + rng.normal(0, sigma, n)

        def alone(y):
            member, p = make_member("a*exp(-x/t)", x, y, ey, "solo")
            member.add_output_port("residuals_solo", bff.GraphPort([0.0], False, True))
            for port in p.values():
                port.value = 1.0
            one = bff.Minimizer()
            one.set_parameter_ports([p["a"], p["t"]])
            one.set_objective(member, "residuals")
            one._graph = member
            one.run()
            return one.x[1]

        t1, t2 = alone(y1), alone(y2)
        self.assertAlmostEqual(t1, 2.6, delta=0.05)
        self.assertAlmostEqual(t2, 3.4, delta=0.05)

        m, joint, free, _ = build_group(x, y1, y2, ey)
        m.run()
        t_joint = m.x[2]
        self.assertGreater(t_joint, min(t1, t2) + 1e-3)
        self.assertLess(t_joint, max(t1, t2) - 1e-3)

    def test_members_keep_their_own_window(self):
        """A group whose members share no fit range is the normal case."""
        x, y1, y2, ey = two_datasets()
        m, joint, free, (m1, m2) = build_group(x, y1, y2, ey)
        m1.set_fit_range(0, 100)
        m2.set_fit_range(50, 200)
        m.run()
        self.assertEqual(list(joint.get_block_sizes()), [100, 150])
        self.assertEqual(np.asarray(m.residuals).size, 250)

    def test_one_member_is_a_plain_fit(self):
        x, y1, y2, ey = two_datasets()
        m1, p1 = make_member("a*exp(-x/t)", x, y1, ey, "only")
        joint = bff.JointChiSquared("joint")
        joint.add_output_port("joint", bff.GraphPort(0.0, False, True))
        joint.add_output_port("residuals", bff.GraphPort([0.0], False, True))
        joint.add_member(m1, "residuals")
        for port in p1.values():
            port.value = 1.0
        m = bff.Minimizer()
        m.set_parameter_ports([p1["a"], p1["t"]])
        m.set_objective(joint, "residuals")
        m._graph = (m1, joint, p1)
        self.assertIn(m.run(), (1, 2, 3, 4))
        self.assertAlmostEqual(m.x[1], 3.0, delta=0.02)

    def test_a_member_without_residuals_is_refused(self):
        joint = bff.JointChiSquared("joint")
        plain = bff.GraphNode("plain")
        plain.add_output_port("chi2", bff.GraphPort(0.0, False, True))
        with self.assertRaises(ValueError):
            joint.add_member(plain, "residuals")

    def test_a_node_without_its_own_output_port_is_refused(self):
        x, y1, y2, ey = two_datasets()
        m1, _ = make_member("a*exp(-x/t)", x, y1, ey, "d1")
        joint = bff.JointChiSquared("joint")
        joint.add_member(m1, "residuals")
        with self.assertRaises(ValueError):
            joint.update()

    def test_members_are_reported_in_the_order_they_were_added(self):
        x, y1, y2, ey = two_datasets()
        m, joint, free, (m1, m2) = build_group(x, y1, y2, ey)
        self.assertEqual(joint.get_number_of_members(), 2)
        self.assertEqual(list(joint.get_member_names()), ["d1", "d2"])
        self.assertEqual(joint.get_member(0).get_name(), "d1")
        with self.assertRaises(ValueError):
            joint.get_member(2)
        self.assertIn("d1", joint.describe())

    def test_a_nan_in_one_member_makes_the_group_infinite(self):
        """A model that blows up must be rejected, not rewarded.

        chinet's ports floor a NaN to ``np.finfo(float).tiny`` so a stored
        value stays serialisable, and for a fit that is exactly backwards: a
        NaN residual floored to ~0 reads as a *perfect* block, so the
        optimiser walks towards the parameters that break the model. The
        numpy path returns NaN and MINPACK rejects the step, so the graph
        has to as well -- the fit's own transport does not sanitise.
        """
        n = 64
        x = np.linspace(1.0, 5.0, n)
        ey = np.ones(n)
        m1, p1 = make_member("sqrt(a-x)", x, np.ones(n), ey, "d1")
        m2, p2 = make_member("a*exp(-x/t)", x, np.ones(n), ey, "d2")
        joint = bff.JointChiSquared("joint")
        joint.add_output_port("joint", bff.GraphPort(0.0, False, True))
        joint.add_output_port("residuals", bff.GraphPort([0.0], False, True))
        joint.add_member(m1, "residuals")
        joint.add_member(m2, "residuals")

        p1["a"].value = 10.0          # finite everywhere
        joint.update()
        healthy = joint.get_chi2()
        self.assertTrue(math.isfinite(healthy))

        p1["a"].value = 0.0           # sqrt of a negative at every point
        m1._graph[0].set_valid(False)
        m1.set_valid(False)
        joint.set_valid(False)
        joint.update()
        self.assertEqual(joint.get_chi2(), float("inf"))
        self.assertTrue(np.all(np.isnan(
            np.asarray(m1._graph[1].value))),
            "the model curve was floored before the misfit could see it")

    def test_an_ordinary_port_still_sanitises(self):
        """Only the fit's transport opts out; a document value is unchanged."""
        p = bff.GraphPort([1.0])
        p.set_value_vector([float("nan")])
        self.assertGreater(float(np.asarray(p.value)[0]), 0.0)
        self.assertTrue(p.get_sanitize())


if __name__ == "__main__":
    unittest.main()
