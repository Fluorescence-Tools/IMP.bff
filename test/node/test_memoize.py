"""`GraphNode::set_memoize`: skip the work when no input actually moved.

The invalidation contract is about *reachability* -- a port write invalidates
every node the change can reach -- and that is the right default, because it
is the only one that never returns a stale number. But reachable is not the
same as changed, and there is one shape where the difference is most of the
work: a value written back over itself.

That is not a contrived case. A sampler that rejects a move restores the
parameters it came from; a line search that overshoots puts the old point
back; a UI re-evaluates after an event that touched nothing; and
`FitMinimizer::compute_objective_batch` restores the ports it borrowed. Each of
those writes the values already there, invalidates the whole graph, and pays
for a full evaluation to arrive at the number it already had.

So a node may be told it is a *function* of its inputs, and then it compares
them before evaluating. Off by default, because that claim is not true of
every node and the class cannot check it.
"""

import unittest

import numpy as np

from IMP import bff


def decay_graph(n=512):
    """`GraphExpression -> FitChiSquared`, the shape a fit actually runs."""
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

    chi2 = bff.FitChiSquared("chi2")
    chi2.set_data_arrays(np.ascontiguousarray(y), np.ones(n))
    model_in = bff.GraphPort([0.0])
    model_in.link = out
    chi2.add_input_port("model", model_in)
    chi2.add_output_port("chi2", bff.GraphPort(0.0, False, True))
    chi2.add_output_port("residuals", bff.GraphPort([0.0], False, True))
    chi2._graph = (curve, out, model_in, axis)
    return chi2, curve, ports


class MemoizeTests(unittest.TestCase):
    def setUp(self):
        self.chi2, self.curve, self.ports = decay_graph()
        self.ports["a"].value = 2.0
        self.ports["t"].value = 3.0
        self.chi2.update()

    def set_same_values(self):
        self.ports["a"].value = 2.0
        self.ports["t"].value = 3.0

    def test_off_by_default(self):
        self.assertFalse(self.curve.get_memoize())
        self.assertEqual(self.curve.get_memo_hit_count(), 0)

    def test_a_value_written_over_itself_costs_a_full_evaluation_by_default(self):
        """The behaviour being improved on, pinned so the improvement is
        measured against something rather than asserted."""
        before = self.curve.get_evaluation_count()
        for _ in range(10):
            self.set_same_values()
            self.chi2.update()
        self.assertEqual(self.curve.get_evaluation_count() - before, 10)

    def test_memoized_it_costs_one(self):
        for node in (self.curve, self.chi2):
            node.set_memoize(True)
        self.set_same_values()
        self.chi2.update()          # the evaluation that records the inputs
        before = self.curve.get_evaluation_count()
        hits = self.curve.get_memo_hit_count()
        for _ in range(10):
            self.set_same_values()
            self.chi2.update()
        self.assertEqual(self.curve.get_evaluation_count() - before, 0)
        self.assertEqual(self.curve.get_memo_hit_count() - hits, 10)

    def test_the_saving_travels_downstream_as_silence(self):
        """A node that does not evaluate does not write its output, so the
        node below is never invalidated and does not even reach its own memo.
        The second saving is larger than the first and costs nothing."""
        for node in (self.curve, self.chi2):
            node.set_memoize(True)
        self.set_same_values()
        self.chi2.update()
        evals = self.chi2.get_evaluation_count()
        hits = self.chi2.get_memo_hit_count()
        for _ in range(10):
            self.set_same_values()
            self.chi2.update()
        self.assertEqual(self.chi2.get_evaluation_count() - evals, 0)
        self.assertEqual(self.chi2.get_memo_hit_count() - hits, 0)

    def test_the_answer_is_the_same_answer(self):
        """The whole point: identical numbers, fewer evaluations."""
        self.ports["a"].value = 1.7
        self.ports["t"].value = 2.6
        self.chi2.update()
        plain = self.chi2.get_chi2()
        plain_residuals = np.asarray(
            self.chi2.get_output_port("residuals").get_value_view())

        for node in (self.curve, self.chi2):
            node.set_memoize(True)
        for _ in range(5):
            self.ports["a"].value = 1.7
            self.ports["t"].value = 2.6
            self.chi2.update()
        self.assertEqual(self.chi2.get_chi2(), plain)
        np.testing.assert_array_equal(
            np.asarray(self.chi2.get_output_port("residuals").get_value_view()),
            plain_residuals)

    def test_a_real_change_still_gets_through(self):
        """The failure that would make this unusable: a memo that holds when
        an input genuinely moved."""
        for node in (self.curve, self.chi2):
            node.set_memoize(True)
        self.chi2.update()
        first = self.chi2.get_chi2()
        self.ports["t"].value = 2.0
        self.chi2.update()
        self.assertNotEqual(self.chi2.get_chi2(), first)
        # and back again, to the number it had
        self.ports["t"].value = 3.0
        self.chi2.update()
        self.assertEqual(self.chi2.get_chi2(), first)

    def test_a_changed_vector_input_is_a_miss(self):
        """Not only the scalars: the axis is an input too, and a memo that
        compared only what is cheap to compare would be wrong here."""
        for node in (self.curve, self.chi2):
            node.set_memoize(True)
        self.set_same_values()
        self.chi2.update()
        before = self.curve.get_evaluation_count()
        axis = self.chi2._graph[3]
        axis.set_values_array(np.linspace(0.2, 16.0, 512))
        self.chi2.update()
        self.assertEqual(self.curve.get_evaluation_count() - before, 1)

    def test_an_explicit_invalidation_defeats_the_memo(self):
        """`set_valid(False)` means "something changed that you cannot see in
        your inputs" -- a swapped dataset, an internal coefficient, a caller
        who simply wants it run again. A memo that survived it would make
        that call a no-op, which is the one way this feature could silently
        return a stale number."""
        for node in (self.curve, self.chi2):
            node.set_memoize(True)
        self.set_same_values()
        self.chi2.update()
        before = self.curve.get_evaluation_count()

        self.set_same_values()
        self.chi2.update()
        self.assertEqual(self.curve.get_evaluation_count() - before, 0)

        self.curve.set_valid(False)
        self.chi2.update()
        self.assertEqual(self.curve.get_evaluation_count() - before, 1)

    def test_turning_it_off_and_on_does_not_leak_a_stale_record(self):
        for node in (self.curve, self.chi2):
            node.set_memoize(True)
        self.set_same_values()
        self.chi2.update()

        # off: a phase where the node is not a function of its inputs
        self.curve.set_memoize(False)
        self.ports["t"].value = 2.5
        self.chi2.update()
        expected = self.chi2.get_chi2()

        # on again, at the same point -- the record must describe *now*
        self.curve.set_memoize(True)
        self.ports["t"].value = 2.5
        self.chi2.update()
        self.assertEqual(self.chi2.get_chi2(), expected)
        self.ports["t"].value = 3.0
        self.chi2.update()
        self.assertNotEqual(self.chi2.get_chi2(), expected)

    def test_an_off_phase_cannot_leave_a_record_that_hits_later(self):
        """The sequence that would return a genuinely stale number, and the
        reason `set_memoize(False)` drops the record rather than parking it:
        record at t = 3, evaluate elsewhere with the memo off so the record
        is not updated, then come back to t = 3. A parked record would match
        and keep the outputs computed at the *other* point."""
        for node in (self.curve, self.chi2):
            node.set_memoize(True)
        self.ports["t"].value = 3.0
        self.chi2.update()
        at_three = self.chi2.get_chi2()

        self.curve.set_memoize(False)
        self.chi2.set_memoize(False)
        self.ports["t"].value = 2.5
        self.chi2.update()
        self.assertNotEqual(self.chi2.get_chi2(), at_three)

        for node in (self.curve, self.chi2):
            node.set_memoize(True)
        self.ports["t"].value = 3.0
        self.chi2.update()
        self.assertEqual(self.chi2.get_chi2(), at_three)


if __name__ == "__main__":
    unittest.main()
