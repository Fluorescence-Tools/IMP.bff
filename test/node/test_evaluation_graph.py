"""PRD-139: a graph with named outputs, run on demand.

The consumer's requirements, in their words and pinned here:

* a **label** naming a (node, port), not a node name and not a port name, so a
  settings file or a notebook can say `p_R` and a node can still be renamed or
  re-plumbed;
* `run()` for everything registered, `run(outputs=[...])` for exactly what
  those labels need **and no more**;
* a label may name **any** port, not only a terminal one -- pinning an
  intermediate is the normal mode of use, not a corner case: an instrument
  basis several histograms share and a weight scan should not touch, or the
  expected counts, which one wants to compare against the data for two
  milliseconds rather than compute a posterior for twenty seconds;
* two labels sharing an upstream subgraph evaluate it **once**;
* the run returns a report, not values -- values stay on the ports.

The graph below is that shape in miniature: one shared `basis`, two consumers
of it, and a `total` that joins them.
"""

import numpy as np
import pytest

import IMP.bff


class Add(IMP.bff.Node):
    """out = sum of every input port, elementwise, with a tally."""

    def __init__(self, name, inputs=("x",)):
        IMP.bff.Node.__init__(self, name)
        self.calls = 0
        for i in inputs:
            self.add_input_port(i, IMP.bff.Port([0.0]))
        self.add_output_port("out", IMP.bff.Port([0.0]))
        self._inputs = inputs

    def evaluate(self):
        self.calls += 1
        acc = None
        for i in self._inputs:
            v = np.asarray(self.get_input_port(i).get_value_view(), dtype=float)
            acc = v.copy() if acc is None else acc + v
        self.get_output_port("out").set_value_vector((acc + 1.0).tolist())


@pytest.fixture
def graph():
    """source -> basis -> {left, right} -> total"""
    source = Add("source")
    basis = Add("basis")
    left = Add("left")
    right = Add("right")
    total = Add("total", ("a", "b"))
    basis.get_input_port("x").set_link(source.get_output_port("out"))
    left.get_input_port("x").set_link(basis.get_output_port("out"))
    right.get_input_port("x").set_link(basis.get_output_port("out"))
    total.get_input_port("a").set_link(left.get_output_port("out"))
    total.get_input_port("b").set_link(right.get_output_port("out"))
    source.get_input_port("x").set_value_vector([1.0, 2.0])

    g = IMP.bff.EvaluationGraph()
    g.add_output("basis", basis, "out")
    g.add_output("left", left, "out")
    g.add_output("right", right, "out")
    g.add_output("total", total, "out")
    nodes = {n.get_name(): n for n in (source, basis, left, right, total)}
    return g, nodes


def test_registering_outputs_evaluates_nothing(graph):
    g, nodes = graph
    assert g.get_number_of_outputs() == 4
    assert list(g.get_output_labels()) == ["basis", "left", "right", "total"]
    assert all(n.calls == 0 for n in nodes.values())


def test_a_run_of_one_label_stops_at_that_label(graph):
    """`basis` is an intermediate. Asking for it must not compute anything
    that consumes it -- that is the two-millisecond diagnostic instead of the
    twenty-second posterior."""
    g, nodes = graph
    rep = g.run(["basis"])
    assert nodes["source"].calls == 1 and nodes["basis"].calls == 1
    assert nodes["left"].calls == 0 and nodes["right"].calls == 0
    assert nodes["total"].calls == 0
    assert rep.nodes_visited == 2 and rep.nodes_evaluated == 2


def test_a_shared_subgraph_is_evaluated_once(graph):
    """`left` and `right` both read `basis`. Asking for both must run the
    basis once -- the case this graph hits constantly."""
    g, nodes = graph
    rep = g.run(["left", "right"])
    assert nodes["basis"].calls == 1
    assert nodes["source"].calls == 1
    assert nodes["left"].calls == 1 and nodes["right"].calls == 1
    assert rep.nodes_visited == 4 and rep.nodes_evaluated == 4


def test_running_everything_produces_the_answer(graph):
    g, nodes = graph
    g.run()
    # source: x+1 = [2,3]; basis: +1 = [3,4]; left and right: +1 = [4,5];
    # total: a+b+1 = [9,11]
    np.testing.assert_allclose(
        np.asarray(g.get_output_port("total").get_value_view()), [9.0, 11.0])
    np.testing.assert_allclose(
        np.asarray(g.get_output_port("basis").get_value_view()), [3.0, 4.0])


def test_a_second_run_with_nothing_changed_evaluates_nothing(graph):
    g, nodes = graph
    g.run()
    rep = g.run()
    assert rep.nodes_evaluated == 0
    assert rep.nodes_visited == 5


def test_only_what_the_change_reached_is_re_evaluated(graph):
    g, nodes = graph
    g.run()
    before = {k: n.calls for k, n in nodes.items()}
    # a change at the head of the shared subgraph reaches everything
    nodes["source"].get_input_port("x").set_value_vector([5.0, 6.0])
    rep = g.run()
    assert rep.nodes_evaluated == 5
    # a change in one branch reaches that branch and the join, not the other
    before = {k: n.calls for k, n in nodes.items()}
    nodes["left"].get_input_port("x").set_value_vector([0.0, 0.0])
    rep = g.run()
    assert nodes["basis"].calls == before["basis"]
    assert nodes["right"].calls == before["right"]
    assert nodes["left"].calls == before["left"] + 1
    assert nodes["total"].calls == before["total"] + 1


def test_the_dependencies_of_a_label_are_its_upstream(graph):
    g, _ = graph
    assert list(g.get_dependencies("basis")) == ["source", "basis"]
    deps = list(g.get_dependencies("total"))
    assert deps[-1] == "total"
    assert set(deps) == {"source", "basis", "left", "right", "total"}


def test_an_unregistered_label_is_refused_not_ignored(graph):
    """A misspelt label that produced nothing would look exactly like a graph
    with nothing to do."""
    g, _ = graph
    with pytest.raises(ValueError):
        g.run(["p_R"])
    with pytest.raises(ValueError):
        g.get_dependencies("p_R")


def test_a_label_is_not_silently_rebound(graph):
    g, nodes = graph
    with pytest.raises(ValueError):
        g.add_output("basis", nodes["left"], "out")
    g.remove_output("basis")
    g.add_output("basis", nodes["left"], "out")
    assert g.get_output_node("basis").get_name() == "left"


def test_a_port_that_does_not_exist_is_refused(graph):
    g, nodes = graph
    with pytest.raises(ValueError):
        g.add_output("nope", nodes["basis"], "no_such_port")


def test_the_labels_round_trip_onto_rebuilt_nodes(graph):
    """The nodes are not written -- a node's work is its callback, which is
    the caller's code. What round-trips is the naming."""
    g, nodes = graph
    text = g.to_json()
    assert '"label": "total"' in text and '"port": "out"' in text
    back = IMP.bff.EvaluationGraph()
    back.from_json(text, nodes)
    assert list(back.get_output_labels()) == list(g.get_output_labels())
    for label in g.get_output_labels():
        assert back.get_output_node(label).get_name() == \
               g.get_output_node(label).get_name()


def test_loading_against_nodes_that_are_missing_says_so(graph):
    g, nodes = graph
    text = g.to_json()
    back = IMP.bff.EvaluationGraph()
    with pytest.raises(ValueError):
        back.from_json(text, {"source": nodes["source"]})
