"""Assemble without evaluating, run once, and recompute only what moved.

PRD-139's evaluation-graph requirement, and the owner's wording for it
(2026-09-08): "the user can define an evaluation graph that will produce the
posterior densities... should have a run option, so that it is not eval
automatically".

Three of the four properties were already true for nodes whose callback is one
of the built-in C++ operators. **None of the last two were true for a node
whose `evaluate()` is overridden in Python**, which is the only kind a
consumer supplying its own physics can use, and that is what these tests
pin:

* `GraphNode::evaluate()` sets the validity flag as its last statement, so an
  override -- which replaces the method wholesale -- left the node permanently
  invalid and `update()` re-ran it on every pass;
* `update()` copies each linked input from its source before evaluating, and
  that copy went through the ordinary write path, which invalidates the
  node. So even a node that had just run was invalid again by the end of the
  next `update()`.

A counter on the callback is the only honest way to assert "did not
recompute", which is why these read the way they do.
"""

import numpy as np
import pytest

import IMP.bff


class Doubler(IMP.bff.GraphNode):
    """y = 2x, and a tally of how often it actually ran."""

    def __init__(self, name):
        IMP.bff.GraphNode.__init__(self, name)
        self.calls = 0

    def evaluate(self):
        self.calls += 1
        x = np.asarray(self.get_input_port("x").get_value_view(), dtype=float)
        self.get_output_port("y").set_value_vector((x * 2.0).tolist())


def _chain(n):
    """`n` doublers, each fed by the one before."""
    nodes = []
    for i in range(n):
        node = Doubler("n%d" % i)
        node.add_input_port("x", IMP.bff.GraphPort([0.0]))
        node.add_output_port("y", IMP.bff.GraphPort([0.0]))
        if nodes:
            node.get_input_port("x").set_link(nodes[-1].get_output_port("y"))
        nodes.append(node)
    return nodes


def test_assembling_evaluates_nothing():
    nodes = _chain(3)
    assert [n.calls for n in nodes] == [0, 0, 0]
    assert not nodes[-1].is_valid()
    # ...and neither does writing an input: a non-reactive port invalidates,
    # it does not compute.
    nodes[0].get_input_port("x").set_value_vector([1.0, 2.0, 3.0])
    assert [n.calls for n in nodes] == [0, 0, 0]


def test_one_run_produces_the_answer():
    nodes = _chain(3)
    nodes[0].get_input_port("x").set_value_vector([1.0, 2.0, 3.0])
    nodes[-1].update()
    assert [n.calls for n in nodes] == [1, 1, 1]
    got = np.asarray(nodes[-1].get_output_port("y").get_value_view())
    np.testing.assert_allclose(got, np.array([1.0, 2.0, 3.0]) * 8.0)


def test_running_again_with_nothing_changed_recomputes_nothing():
    nodes = _chain(3)
    nodes[0].get_input_port("x").set_value_vector([1.0, 2.0, 3.0])
    nodes[-1].update()
    before = [n.calls for n in nodes]
    nodes[-1].update()
    nodes[-1].update()
    assert [n.calls for n in nodes] == before


def test_only_what_is_downstream_of_the_change_recomputes():
    nodes = _chain(3)
    nodes[0].get_input_port("x").set_value_vector([1.0, 2.0, 3.0])
    nodes[-1].update()
    base = [n.calls for n in nodes]

    # the head: everything below it
    nodes[0].get_input_port("x").set_value_vector([5.0, 6.0, 7.0])
    nodes[-1].update()
    assert [n.calls for n in nodes] == [b + 1 for b in base]
    np.testing.assert_allclose(
        np.asarray(nodes[-1].get_output_port("y").get_value_view()),
        np.array([5.0, 6.0, 7.0]) * 8.0)

    # the middle: itself and the tail, not the head
    base = [n.calls for n in nodes]
    nodes[1].get_input_port("x").set_value_vector([9.0])
    nodes[-1].update()
    assert [n.calls for n in nodes] == [base[0], base[1] + 1, base[2] + 1]


def test_a_node_with_nothing_to_do_stays_invalid():
    """chinet's rule, and the reason `update()` asks whether anything came out
    rather than whether a callback exists: a bare node writes nothing."""
    bare = IMP.bff.GraphNode("bare")
    bare.update()
    assert not bare.get_node_valid()


def test_a_vector_port_needs_get_value_view():
    """`get_value()` is the scalar accessor and returns the first element of a
    vector port -- quietly, which is how a graph comes to compute one number
    where it should compute an array. Pinned because the consumer will meet
    it."""
    p = IMP.bff.GraphPort([0.0])
    p.set_value_vector([1.0, 2.0, 3.0])
    assert p.current_size() == 3
    assert p.get_value() == 1.0
    np.testing.assert_allclose(np.asarray(p.get_value_view()), [1.0, 2.0, 3.0])
