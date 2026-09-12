"""A Python GraphNode held only by its FitMinimizer/MCMCSampler must survive (T-20260901-13).

The C++ side keeps a shared_ptr to the objective node, but a SWIG director
keeps only a weak pointer back to the Python proxy: before the fix, a
ResidualNode bound to `_` was collected when the next tuple unpacking
rebound `_`, and the next Node_update dispatched into a dead object — a
segfault on the fallback path every graph-less chisurf fit takes.

The fix is a %pythonappend on set_objective stashing the proxy on the
wrapper that holds the C++ reference (IMP_SWIG_DIRECTOR does not work here:
its registry silently refuses anything without IMP's get_ref_count, and
GraphNode is a plain shared_ptr class).
"""
import gc
import unittest

import numpy as np

import IMP.bff as bff


class Residual(bff.GraphNode):
    def __init__(self):
        super().__init__("resid")
        self._p = [bff.GraphPort(1.0), bff.GraphPort(2.0)]
        for i, p in enumerate(self._p):
            self.add_input_port("x%d" % i, p)
        self.add_output_port("residuals", bff.GraphPort([0.0], False, True))

    def evaluate(self):
        x = [self.get_input_port("x%d" % i).value for i in range(2)]
        r = [x[0] - 3.0, x[1] - 5.0, (x[0] - 3.0) * 0.5]
        self.get_output_port("residuals").set_values_array(np.asarray(r))


class TestNodeLifetime(unittest.TestCase):

    def test_a_dropped_python_node_survives_a_minimizer_run(self):
        node = Residual()
        m = bff.FitMinimizer()
        m.set_parameter_ports(node._p)
        m.set_objective(node, "residuals")
        # The exact pattern that lost the proxy: bind to _, rebind _.
        _ = node
        node = None
        _ = (1, 2)
        gc.collect()
        m.run()
        np.testing.assert_allclose(list(m.get_x()), [3.0, 5.0], atol=1e-6)

    def test_a_dropped_python_node_survives_a_sampler_run(self):
        node = Residual()
        s = bff.MCMCSampler("metropolis", 7)
        s.set_parameter_ports(node._p)
        s.set_objective(node, "residuals")
        node = None
        gc.collect()
        s.run(50, 1)
        self.assertEqual(int(s.iteration), 50)

    def test_replacing_the_objective_releases_the_old_proxy(self):
        """The keepalive mirrors the C++ reference — it must not leak: a
        node replaced as the objective loses its stash and may die."""
        import weakref
        node = Residual()
        m = bff.FitMinimizer()
        m.set_parameter_ports(node._p)
        m.set_objective(node, "residuals")
        ref = weakref.ref(node)
        replacement = Residual()
        m.set_parameter_ports(replacement._p)
        m.set_objective(replacement, "residuals")
        node = None
        gc.collect()
        self.assertIsNone(ref(), "the replaced proxy must be collectable")


if __name__ == '__main__':
    unittest.main()
