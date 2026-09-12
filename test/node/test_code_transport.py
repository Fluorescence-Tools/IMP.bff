"""Carrying a Python node's code in the document, so a graph can travel.

Python callbacks are the prototyping path, not the intended one -- what a
model needs for speed belongs in C++ nodes (`okf/reference-cxx-nodes.md`).
But a prototype still has to be shareable, and the naive attempt at that
fails in a way worth pinning: the class source alone rebuilds into a
NameError, because a callback closes over its imports and its module
constants and **the source is not the closure**.

So the document carries the imports and the referenced module-level values
with it, and then it does work. It is gated on `trusted=True` because loading
such a document executes it.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

import IMP.bff

sys.path.insert(0, str(Path(__file__).resolve().parent / "codefixture"))
import mymodel  # noqa: E402


def _run(cls, values):
    n = cls()
    n.add_input_port("x", IMP.bff.GraphPort([0.0]))
    n.add_output_port("y", IMP.bff.GraphPort([0.0]))
    n.get_input_port("x").set_value_vector(list(values))
    n.update()
    return np.asarray(n.get_output_port("y").get_value_view())


def test_the_capture_carries_the_imports_and_the_constants():
    code = IMP.bff.node_code(mymodel.Scaled)
    assert code["class"] == "Scaled"
    assert "import numpy as np" in code["imports"]
    assert any("IMP.bff as bff" in i for i in code["imports"]), \
        "the base class alias must travel; it is evaluated at class creation"
    assert code["globals"]["SCALE"] == 3.0
    assert code["globals"]["OFFSET"] == {"__ndarray__": [0.5]}
    assert code["unresolved"] == [], "nothing was left behind"


def test_a_node_rebuilt_from_the_document_computes_the_same_thing():
    code = json.loads(json.dumps(IMP.bff.node_code(mymodel.Scaled)))
    rebuilt = IMP.bff.rebuild_node_class(code, trusted=True)
    assert rebuilt is not mymodel.Scaled
    np.testing.assert_allclose(_run(rebuilt, [1.0, 2.0]),
                               _run(mymodel.Scaled, [1.0, 2.0]))
    np.testing.assert_allclose(_run(rebuilt, [1.0, 2.0]), [3.5, 6.5])


def test_the_source_alone_is_not_enough():
    """The naive version, pinned so nobody re-derives it: strip what the
    callback closes over and the rebuild dies on the first name it wants."""
    code = IMP.bff.node_code(mymodel.Scaled)
    stripped = dict(code, imports=[], globals={})
    with pytest.raises(NameError):
        IMP.bff.rebuild_node_class(stripped, trusted=True)


def test_rebuilding_will_not_happen_by_accident():
    """A document that carries code executes when it is loaded, and whether a
    file is trustworthy is the caller's knowledge."""
    code = IMP.bff.node_code(mymodel.Scaled)
    with pytest.raises(ValueError):
        IMP.bff.rebuild_node_class(code)


def test_a_class_with_no_readable_source_can_supply_its_own():
    """A node written in a notebook cell has no source to introspect;
    passing it explicitly is the way through."""
    ns = {"bff": IMP.bff}
    src = ("class Cell(bff.GraphNode):\n"
           "    def evaluate(self):\n"
           "        self.get_output_port('y').set_value_vector([7.0])\n")
    exec(src, ns)
    with pytest.raises((OSError, TypeError)):
        IMP.bff.node_code(ns["Cell"])
    # ...and neither is there a module to take the imports from, so those are
    # given too. Introspection gets what it can; the caller supplements.
    code = IMP.bff.node_code(ns["Cell"], source=src,
                             imports=["import IMP.bff as bff"])
    rebuilt = IMP.bff.rebuild_node_class(code, trusted=True)
    n = rebuilt()
    n.add_output_port("y", IMP.bff.GraphPort([0.0]))
    n.update()
    np.testing.assert_allclose(np.asarray(n.get_output_port("y").get_value_view()), [7.0])


def test_a_whole_graph_document_carries_its_nodes():
    a = mymodel.Scaled()
    a.add_input_port("x", IMP.bff.GraphPort([0.0]))
    a.add_output_port("y", IMP.bff.GraphPort([0.0]))
    a.set_name("a")
    g = IMP.bff.GraphEvaluation()
    g.add_output("out", a, "y")
    doc = IMP.bff.graph_code(g, {"a": a})
    assert "Scaled" in doc
    classes = IMP.bff.rebuild_node_classes(doc, trusted=True)
    assert set(classes) == {"a"}
    np.testing.assert_allclose(_run(classes["a"], [1.0]), [3.5])
