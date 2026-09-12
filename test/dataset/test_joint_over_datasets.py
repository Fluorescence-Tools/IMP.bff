"""PRD-140: a group whose members disagree about noise, through the graph.

`test_objective_over_datasets.py` proves the arithmetic -- a Poisson member
and a Gaussian member each weight their own way and sum to the total -- but
it drives `compute_chi2` directly, so it says nothing about whether the
*node* composes them. That is the part a joint fit actually runs: nobody
calls `compute_chi2` per member during a minimisation, `GraphNode::update()` walks
the group and every crossing into Python is one that a sampler pays for
thousands of times.

So these tests hold the objective at arm's length. They evaluate the joint
node and read the ports, which is the path `Minimizer` takes.
"""

import unittest

import numpy as np

from IMP import bff


def poisson_dataset(y):
    d = bff.Dataset()
    d.set_values_array(np.ascontiguousarray(y, dtype=float))
    d.set_noise_family(bff.NOISE_FAMILY_POISSON)
    return d


def gaussian_dataset(y, variance):
    d = bff.Dataset()
    d.set_values_array(np.ascontiguousarray(y, dtype=float))
    d.set_noise_family(bff.NOISE_FAMILY_STORED)
    d.set_stored_variance_array(np.ascontiguousarray(variance, dtype=float))
    return d


def make_member(equation, x, dataset, name):
    """One dataset as a node: `GraphExpression -> ChiSquared`, residuals on a port.

    The same shape as `test/minimizer/test_joint.py`'s member, except the
    misfit is *given the dataset* rather than told a noise model name.
    """
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
    chi2.set_dataset(dataset)
    model_in = bff.GraphPort([0.0])
    model_in.link = out
    chi2.add_input_port("model", model_in)
    chi2.add_output_port(name, bff.GraphPort(0.0, False, True))
    chi2.add_output_port("residuals", bff.GraphPort([0.0], False, True))
    # The Python proxies have to outlive this function or the C++ graph is
    # left holding nodes whose directors have been collected.
    chi2._graph = (curve, out, model_in, axis)
    return chi2, parameters


class HeterogeneousGroupTests(unittest.TestCase):
    """A decay scored as counts beside a correlation curve scored against a
    measured variance -- the case PRD-140 exists for."""

    def setUp(self):
        self.x = np.linspace(0.5, 8.0, 24)
        # counts: a decay, Poisson
        self.y_counts = np.round(400.0 * np.exp(-self.x / 3.0)) + 1.0
        # a correlation curve: small numbers, a measured variance, and
        # emphatically not Poisson -- its values are near 1, so counting
        # statistics would weight it by nothing.
        self.y_corr = 1.0 + 0.2 * np.exp(-self.x / 3.0)
        self.var_corr = np.full(self.x.size, 1e-4)

        self.d_counts = poisson_dataset(self.y_counts)
        self.d_corr = gaussian_dataset(self.y_corr, self.var_corr)

    def build(self):
        m1, p1 = make_member("a*exp(-x/t)", self.x, self.d_counts, "decay")
        m2, p2 = make_member("1.0+a*exp(-x/t)", self.x, self.d_corr, "corr")
        # the coupling: one lifetime, two datasets, two noise families
        p2["t"].link = p1["t"]

        joint = bff.JointChiSquared("joint")
        joint.add_output_port("joint", bff.GraphPort(0.0, False, True))
        joint.add_output_port("residuals", bff.GraphPort([0.0], False, True))
        joint.add_member(m1, "residuals")
        joint.add_member(m2, "residuals")
        joint._graph = (m1, m2, p1, p2)
        return joint, (m1, m2), (p1, p2)

    def test_the_node_composes_members_that_disagree_about_noise(self):
        """The claim the arithmetic test could not make: evaluating the
        *graph* gives each member the weighting its own dataset carries."""
        joint, (m1, m2), (p1, p2) = self.build()
        p1["a"].value = 380.0
        p1["t"].value = 2.8
        p2["a"].value = 0.19
        joint.update()

        # what each member's dataset says, computed outside the graph
        mu_counts = 380.0 * np.exp(-self.x / 2.8)
        mu_corr = 1.0 + 0.19 * np.exp(-self.x / 2.8)
        expected = (self.d_counts.objective(list(mu_counts))
                    + self.d_corr.objective(list(mu_corr)))

        self.assertAlmostEqual(
            joint.get_output_port("joint").value / expected, 1.0, places=10)
        self.assertAlmostEqual(joint.get_chi2() / expected, 1.0, places=10)

    def test_each_block_is_its_own_family(self):
        """Not merely that the total matches: the blocks are separable, and
        each equals what its dataset alone would report."""
        joint, (m1, m2), (p1, p2) = self.build()
        p1["a"].value = 380.0
        p1["t"].value = 2.8
        p2["a"].value = 0.19
        joint.update()

        offsets = joint.get_block_offsets()
        sizes = joint.get_block_sizes()
        self.assertEqual(list(sizes), [self.x.size, self.x.size])

        res = np.asarray(joint.get_output_port("residuals").get_value_view())
        mu_counts = 380.0 * np.exp(-self.x / 2.8)
        mu_corr = 1.0 + 0.19 * np.exp(-self.x / 2.8)
        # the counts block carries the deviance (the family's own residual),
        # the correlation block Pearson -- two definitions in one vector
        np.testing.assert_allclose(
            res[offsets[0]:offsets[0] + sizes[0]],
            np.asarray(self.d_counts.residuals(list(mu_counts),
                                               bff.RESIDUAL_DEVIANCE)),
            rtol=1e-12)
        np.testing.assert_allclose(
            res[offsets[1]:offsets[1] + sizes[1]],
            np.asarray(self.d_corr.residuals(list(mu_corr),
                                             bff.RESIDUAL_PEARSON)),
            rtol=1e-12)

    def test_the_weighting_is_not_incidental(self):
        """Scoring the correlation curve as counts gives a different number,
        so the family reaching the graph is doing work rather than agreeing
        with whatever it replaced."""
        joint, (m1, m2), (p1, p2) = self.build()
        p1["a"].value = 380.0
        p1["t"].value = 2.8
        p2["a"].value = 0.19
        joint.update()
        proper = joint.get_chi2()

        m2.set_dataset(poisson_dataset(self.y_corr))
        m2.set_valid(False)
        joint.set_valid(False)
        joint.update()
        self.assertNotAlmostEqual(joint.get_chi2() / proper, 1.0, places=3)

    def test_a_minimizer_moves_the_shared_parameter_over_both_families(self):
        """The point of the exercise: one step uses the curvature of a
        Poisson dataset and a Gaussian one at the same time."""
        joint, (m1, m2), (p1, p2) = self.build()
        free = [p1["a"], p1["t"], p2["a"]]
        for port, v in zip(free, (300.0, 2.0, 0.1)):
            port.value = v
        m = bff.Minimizer()
        m.set_parameter_ports(free)
        m.set_objective(joint, "residuals")
        m._graph = (joint, m1, m2)
        info = m.run()
        self.assertIn(info, (1, 2, 3, 4))
        a_counts, t, a_corr = m.x
        # the data were generated at t = 3 with amplitudes 400 and 0.2; the
        # counts were rounded and floored by 1, so t comes out near but not
        # exactly at 3
        self.assertAlmostEqual(t, 3.0, delta=0.05)
        self.assertAlmostEqual(a_counts, 400.0, delta=8.0)
        self.assertAlmostEqual(a_corr, 0.2, delta=0.01)
        # and the shared lifetime really is one number
        self.assertEqual(p1["t"].value, p2["t"].value)


if __name__ == "__main__":
    unittest.main()
