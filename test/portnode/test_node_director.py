"""Python callbacks on Node, through the SWIG director.

chinet's nodes accept Python callables: set_python_callback_function wraps a
function into a node's ports and evaluate() calls it with the input-port
values. bff's Node is a director, so a Python subclass overrides evaluate()
(or update()) and the C++ machinery drives it: a reactive port write, a
linked follower's update(), a session-registered node -- every path that
calls Node::evaluate() from C++ crosses into the Python override.

The set_python_callback_function ergonomics themselves (inspect the
signature, auto-create input ports from the parameters and output ports from
a dict return) are a chisurf-side helper (chisurf/core/nodes.py); these
tests cover the director mechanics it is built on.
"""
import pytest

from IMP.bff import Node, Port, get_session


class CallbackNode(Node):
    """A director Node running ``func`` over its input-port values.

    The evaluate() body mirrors chinet's: kwargs from the inputs, a dict
    return written port by port, a tuple/list return to out_00.., anything
    else to out_00, a TypeError falling back to the ports-maps calling
    convention, other nodes sharing an output port invalidated, and the
    node valid afterwards.
    """

    def __init__(self, func, name=""):
        super().__init__(name)
        self._func = func

    def evaluate(self):
        args = {k: p.value for k, p in self.inputs.items()}
        outs = self.outputs
        try:
            res = self._func(**args)
        except TypeError:
            self._func(self.inputs, self.outputs)
        else:
            if isinstance(res, dict):
                for k, v in res.items():
                    if k in outs:
                        outs[k].value = v
            elif isinstance(res, (tuple, list)) and len(res) > 1:
                for i, v in enumerate(res):
                    oname = "out_%02d" % i
                    if oname in outs:
                        outs[oname].value = v
            elif "out_00" in outs:
                outs["out_00"].value = res
        for p in outs.values():
            n = p.get_node()
            if n is not None and n.get_uid() != self.get_uid():
                n.set_valid(False)
        self.set_valid(True)


def two_outputs(a=0.0, b=0.0):
    return {"sum": a + b, "prod": a * b}


def test_director_subclass_reads_inputs_and_writes_outputs():
    node = CallbackNode(two_outputs, "fn")
    node.add_input_port("a", Port(3.0))
    node.add_input_port("b", Port(4.0))
    node.add_output_port("sum", Port(0.0, False, True))
    node.add_output_port("prod", Port(0.0, False, True))
    assert not node.is_valid()
    node.evaluate()
    assert node.outputs["sum"].value == pytest.approx(7.0)
    assert node.outputs["prod"].value == pytest.approx(12.0)
    assert node.is_valid()


def test_cpp_update_dispatches_into_python_evaluate():
    """A C++ update() pull crosses the director boundary.

    The follower node is plain C++ (operator callback); its linked input
    pulls from the Python node's output, evaluating the director upstream.
    """
    py = CallbackNode(lambda x=1.0: {"y": 2.0 * x}, "doubler")
    py.add_input_port("x", Port(1.0))
    py.add_output_port("y", Port(0.0, False, True))

    follower = Node("follower")
    follower.add_input_port("x", Port(0.0))
    follower.inputs["x"].link = py.outputs["y"]
    follower.add_input_port("one", Port(1.0))
    follower.add_output_port("follower", Port(0.0, False, True))
    follower.set_callback("multiply_double", "C")

    py.inputs["x"].value = 5.0
    follower.update()
    assert py.outputs["y"].value == pytest.approx(10.0)
    assert follower.outputs["follower"].value == pytest.approx(10.0)


def test_reactive_port_write_evaluates_python_node_from_cpp():
    node = CallbackNode(lambda x=1.0: {"y": 3.0 * x}, "reactor")
    node.add_input_port("x", Port(2.0, fixed=False, is_reactive=True))
    node.add_output_port("y", Port(0.0, False, True, True))
    node.inputs["x"].value = 7.0
    assert node.outputs["y"].value == pytest.approx(21.0)
    assert node.is_valid()


def test_director_registered_with_session_stays_driven():
    """A session-registered director keeps evaluating from C++ calls."""
    node = CallbackNode(two_outputs, "adder")
    node.add_input_port("a", Port(1.0))
    node.add_input_port("b", Port(1.0))
    node.add_output_port("sum", Port(0.0, False, True))
    node.add_output_port("prod", Port(0.0, False, True))
    get_session().add_node("director_adder", node)
    node.inputs["a"].value = 10.0
    node.update()
    assert node.outputs["sum"].value == pytest.approx(11.0)
    assert node.outputs["prod"].value == pytest.approx(10.0)


def test_update_override_is_dispatched():
    """update() is virtual too; a Python override can intercept the pull."""
    calls = []

    class UpdateSpy(CallbackNode):
        def update(self):
            calls.append("update")
            super().update()

    node = UpdateSpy(lambda x=0.0: {"y": x + 1.0}, "spy")
    node.add_input_port("x", Port(4.0))
    node.add_output_port("y", Port(0.0, False, True))
    node.update()
    assert calls == ["update"]
    assert node.outputs["y"].value == pytest.approx(5.0)


def test_type_error_falls_back_to_ports_maps_convention():
    """chinet's second calling convention: func(inputs, outputs)."""
    seen = {}

    def maps_cb(inputs, outputs):
        seen["keys"] = sorted(k for k in inputs)
        outputs["out_00"].value = 42.0

    node = CallbackNode(maps_cb, "maps")
    node.add_input_port("a", Port(1.0))
    node.add_output_port("out_00", Port(0.0, False, True))
    node.evaluate()
    assert seen["keys"] == ["a"]
    assert node.outputs["out_00"].value == pytest.approx(42.0)


def test_director_node_sharing_an_output_port_is_invalidated():
    """chinet: writing an output invalidates every other node reading it.

    A second node adopts the first's output port (add_port of a shared
    port, chinet's set_ports); evaluating the director must mark that node
    invalid while the director itself becomes valid.
    """
    producer = CallbackNode(lambda x=0.0: {"out": x}, "producer")
    producer.add_input_port("x", Port(1.0))
    producer.add_output_port("out", Port(0.0, False, True))

    reader = Node("reader")
    # An input, so is_valid() consults node_valid_ (an input-less node is
    # trivially valid in both chinet and bff and the flag is unobservable).
    reader.add_input_port("seed", Port(0.0))
    reader.add_output_port("out", producer.outputs["out"])

    producer.evaluate()
    assert producer.is_valid()
    assert not reader.is_valid()
    assert reader.outputs["out"].value == pytest.approx(1.0)
