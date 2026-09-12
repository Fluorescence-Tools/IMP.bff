"""The reactive GraphPort/GraphNode runtime, in bff (ported from chinet).

chinet is the parameter/node layer chisurf's models are built on: a GraphPort is
a value cell (scalar or array) with a fixed flag, hard bounds, an optional
prior, and a link to a port it follows; a GraphNode owns named ports and computes
its outputs from its inputs, lazily, driven by invalidation. The link graph
must stay acyclic -- linking is the only operation that adds dependency
edges, and set_link refuses a cycle with Kahn's algorithm.

These tests port the scenarios of chinet's test/test_port.py and
test/test_node.py onto the bff C++ implementation. Differences from the
chinet suite, all deliberate (see include/GraphPort.h):

- chinet's GraphLinkCycleError subclasses ValueError; here it crosses the
  wrapper as IMP.ValueException, also a ValueError -- so the tests catch
  ValueError, the contract chisurf's callers rely on.
- values are stored as double; the int/float distinction lives in the
  value-type code (0/1 scalar, 2/3 vector) and integral ports round.
- chinet's Python-callable node callbacks (callback_class,
  callback_function) run through the SWIG director and are covered by
  test_node_director.py; the node tests here drive the string operator
  callbacks ("multiply_double", "C" ...) instead,
  which compute over the first two input ports and write the output port
  keyed by the node's name. A pass-through is a multiply by a port of 1.0.
- the DB-backed tests (write_to_db/read_from_db, JSON round-trip) are
  phase 2 and are not ported; prior round-trip is covered directly.
"""
import numpy as np
import pytest

from IMP.bff import GraphExpression, GraphNode, GraphPort


def passthrough_node(name, value):
    """A node that copies its 'inA' input to its output, chinet-style.

    Stands in for chinet's CallbackNodePassOn: multiply by a constant 1.0
    through the C operator protocol (output port keyed by the node name).
    """
    node = GraphNode(name)
    node.add_input_port("inA", GraphPort(value))
    node.add_input_port("one", GraphPort(1.0))
    node.add_output_port(name, GraphPort(0.0, False, True))
    node.set_callback("multiply_double", "C")
    return node


# ---------------------------------------------------------------------------
# GraphPort: construction, value, fixed (chinet test_port_init_singelton)
# ---------------------------------------------------------------------------


def test_port_init_singleton():
    v1 = 23.0
    v2 = 29.0
    p1 = GraphPort(v1)
    p2 = GraphPort()
    p2.value = v1
    assert p1.value == pytest.approx(p2.value)

    p3 = GraphPort(v1, True)  # fixed, positionally
    p4 = GraphPort(v1, False)
    assert p3.fixed is True
    assert p4.fixed is False

    p5 = GraphPort(v2)
    p5.link = p4
    assert p5.value == pytest.approx(p4.value)


def test_port_init_singleton_bounded():
    """chinet's second singleton test: bounds clip the initial value."""
    p6 = GraphPort(0.0, False, False, False, True, 2, 5)
    assert p6.value <= 5
    assert p6.value >= 2
    # the lower bound is not part of the [lb, ub) half-open interval
    assert p6.value == pytest.approx(2)


def test_port_bounds():
    v1 = np.array(
        [1, 2, 3, 6, 5.5, -3, -2, -6.1, -10000, 10000], dtype=np.double
    )
    p1 = GraphPort()
    p1.value = v1
    assert np.allclose(p1.value, v1)

    p1.bounds = 0, 1
    p1.bounded = True

    assert (p1.value <= 1).all()
    assert (p1.value >= 0).all()


def test_port_get_set_value():
    v1 = [1, 2, 3]
    p1 = GraphPort()
    p1.value = v1

    p2 = GraphPort()
    p2.value = v1
    assert (p1.value == p2.value).all()


def test_port_init_vector():
    v1 = [1, 2, 3, 5, 8]
    v2 = [1, 2, 4, 8, 16]

    p1 = GraphPort()
    p1.value = v1
    p2 = GraphPort(v1)
    assert list(p2.value) == list(p1.value)

    p3 = GraphPort(v1, True)
    p4 = GraphPort(v1, False)
    assert p3.fixed is True
    assert p4.fixed is False

    p5 = GraphPort(v2, False)
    p5.link = p4
    assert np.allclose(p5.value, p4.value)


def test_port_init_array():
    array = np.array([1, 2, 3, 5, 8, 13], dtype=np.double)
    p1 = GraphPort()
    p1.value = array
    p2 = GraphPort(array)
    assert list(p1.value) == list(p2.value)


def test_set_get_value_1():
    value = 23.0
    port = GraphPort(value)
    assert port.value == value


def test_set_get_value_2():
    """A one-element tuple is a vector write; it reads back as one."""
    value = (1,)
    port = GraphPort()
    port.value = value
    assert np.allclose(port.value, value)


# ---------------------------------------------------------------------------
# GraphPort: links
# ---------------------------------------------------------------------------


def test_port_link_value():
    value1 = np.array([12], dtype=np.double)
    value2 = np.array([6], dtype=np.double)
    p1 = GraphPort(value1)
    p2 = GraphPort(value2)

    assert np.allclose(p1.value, value1)
    assert np.allclose(p2.value, value2)
    assert not np.allclose(p1.value, p2.value)

    p2.link = p1
    assert np.allclose(p1.value, p2.value)
    p2.unlink()
    assert np.allclose(p2.value, value2)


def test_port_link_dag_enforced():
    # The link graph must remain acyclic; set_link runs Kahn's algorithm to
    # reject any link that would create a cycle.
    a = GraphPort(1.0)
    b = GraphPort(2.0)
    c = GraphPort(3.0)

    # Build the chain c -> b -> a (each follows the previous).
    b.link = a
    c.link = b

    # would_create_cycle is a side-effect-free predicate.
    assert a.would_create_cycle(c)  # a -> c closes the loop
    e = GraphPort(5.0)
    assert not a.would_create_cycle(e)  # independent port is fine

    # Closing the loop must raise and leave the graph untouched.
    with pytest.raises(ValueError):
        a.link = c
    assert not a.is_linked()

    # Linking a port to itself is a degenerate cycle.
    with pytest.raises(ValueError):
        a.link = a

    # A non-cyclic re-link is still allowed.
    d = GraphPort(4.0)
    c.link = d
    assert c.is_linked()


def test_node_link_dag_enforced():
    # Cycles across nodes (n1 depends on n2 depends on n1) are rejected too.
    n1 = GraphNode("n1")
    n1.add_input_port("in", GraphPort(0.0))
    n1.add_output_port("out", GraphPort(0.0, False, True))
    n2 = GraphNode("n2")
    n2.add_input_port("in", GraphPort(0.0))
    n2.add_output_port("out", GraphPort(0.0, False, True))

    # n2 depends on n1.
    n2.get_input_port("in").link = n1.get_output_port("out")
    # n1 depending back on n2 would create a node-level cycle.
    with pytest.raises(ValueError):
        n1.get_input_port("in").link = n2.get_output_port("out")


def test_port_link_same_node_allowed():
    """Two ports of the same node may be linked (a node reads its output)."""
    node = GraphNode("nd")
    a = GraphPort(1.0)
    b = GraphPort(2.0)
    node.add_input_port("in", a)
    node.add_output_port("out", b)
    a.link = b  # same vertex, not a cycle
    assert a.is_linked()


def test_update_does_not_recurse_on_a_self_link():
    """A node whose input follows another of its *own* ports must not recurse.

    Two variables of one expression that are one number -- which is how a
    parameter shared inside a single model is spelt -- link input to input on
    the same node. `GraphNode::update()` used to ask the source port's node to
    update, find itself, and take the stack with it: a **segfault**, not an
    exception. `inputs_valid()` had always skipped the self-link; `update()`
    had not.
    """
    x = np.linspace(0.1, 10.0, 8)
    curve = GraphExpression("model")
    curve.set_expression("a*exp(-x/t)+b")
    a, t, b = GraphPort(1.0), GraphPort(1.0), GraphPort(1.0)
    curve.add_input_port("a", a)
    curve.add_input_port("t", t)
    curve.add_input_port("b", b)
    axis = GraphPort([0.0])
    axis.set_values_array(np.ascontiguousarray(x))
    curve.add_input_port("x", axis)
    out = GraphPort([0.0], False, True)
    curve.add_output_port("model", out)

    b.link = a                      # the offset *is* the amplitude
    curve.update()
    np.testing.assert_allclose(np.asarray(out.value),
                               1.0 * np.exp(-x / 1.0) + 1.0)

    # And the follower follows: writing the master moves both terms.
    a.value = 2.0
    curve.update()
    np.testing.assert_allclose(np.asarray(out.value),
                               2.0 * np.exp(-x / 1.0) + 2.0)


def test_port_write_propagates_to_follower():
    p4 = GraphPort(23.0)
    p5 = GraphPort(29.0)
    p5.link = p4
    p4.value = 31.0
    # chinet's quirk, reproduced: the write pushes the source's array into
    # every follower, so the follower now reads as a vector.
    assert np.allclose(p5.value, [31.0])
    assert p5.get_is_vector()

    # The follower keeps the last pushed value once unlinked.
    p5.unlink()
    assert np.allclose(p5.value, [31.0])


# ---------------------------------------------------------------------------
# GraphPort: flags, type codes, sanitisation, prior, metadata
# ---------------------------------------------------------------------------


def test_port_fixed():
    p1 = GraphPort(12)
    p1.fixed = True
    assert p1.fixed is True
    p1.value = 55
    assert p1.value == pytest.approx(12)  # a fixed port ignores writes

    p1.fixed = False
    assert p1.fixed is False
    p1.value = 55
    assert p1.value == pytest.approx(55)


def test_port_reactive():
    p1 = GraphPort(12)
    p1.reactive = True
    assert p1.reactive is True

    p1.reactive = False
    assert p1.reactive is False


def test_port_value_type():
    # the constructor argument declares the type; nothing later changes it.
    assert GraphPort(12).get_value_type() == 0
    assert GraphPort(12.0).get_value_type() == 1
    # a vector write sets the vector code (the constructor leaves 1).
    p0 = GraphPort()
    p0.value = [1.0, 2.0]
    assert p0.get_value_type() == 3
    assert p0.get_is_vector()

    # Writing a float into an int port does NOT convert the port: since
    # 2026-09-02 the element type is declared at construction and a write of
    # another kind is accepted and coerced. chinet inferred the dtype from
    # every write and so let a port become something else mid-session.
    p = GraphPort(12)
    p.value = 1.5
    assert p.get_value_type() == 0
    assert p.value == 1

    # Writing an int into a float port keeps it float.
    p2 = GraphPort(1.5)
    p2.value = 2
    assert p2.get_value_type() == 1
    assert p2.value == pytest.approx(2.0)

    # set_value_type converts the stored data (astype semantics: truncate).
    p3 = GraphPort(5.5)
    p3.set_value_type(0)
    assert p3.get_value_type() == 0
    assert p3.value == pytest.approx(5)


def test_port_value_sanitised():
    """chinet sanitises float writes: optimisers emit NaN and inf."""
    p = GraphPort(1.0)
    p.value = float("nan")
    assert p.value == pytest.approx(2.2250738585072014e-308)  # np tiny
    p.value = float("inf")
    assert p.value == pytest.approx(1.7976931348623157e308)
    p.value = float("-inf")
    assert p.value == pytest.approx(-1.7976931348623157e308)

    pv = GraphPort([1.0, 1.0])
    pv.value = [1.0, float("nan")]  # sanitisation is a write-path behaviour
    assert pv.value[1] == pytest.approx(2.2250738585072014e-308)


def test_port_prior():
    p = GraphPort(1.0)
    assert p.prior is None

    spec = {"kind": "gaussian", "mu": 1.0, "sigma": 2.0}
    p.prior = spec
    assert p.prior == spec

    p.prior = None
    assert p.prior is None


def test_port_bounds_accessors():
    p = GraphPort(1.0)
    p.set_bounds(0.0, 2.0)
    assert p.get_lower_bound() == pytest.approx(0.0)
    assert p.get_upper_bound() == pytest.approx(2.0)

    # Unenforced bounds report nan, standing in for chinet's (None, None).
    assert np.isnan(p.bounds[0]) and np.isnan(p.bounds[1])
    p.bounded = True
    assert p.bounds == (0.0, 2.0)

    p.value = 10.0
    assert p.value == pytest.approx(2.0)  # clipped on write

    # Tightening the bounds clips the stored data.
    p.set_bounds(0.0, 1.0)
    assert p.value == pytest.approx(1.0)


def test_port_metadata():
    p = GraphPort(1.0)
    p.name = "gamma"
    assert p.get_name() == "gamma"
    assert p.name == "gamma"

    q = GraphPort(1.0)
    assert p.uid != q.uid  # every object gets its own id
    q.set_uid("patched")
    assert q.uid == "patched"


def test_port_size_and_validity():
    p = GraphPort([1.0, 2.0, 3.0])
    assert p.current_size() == 3
    assert p.is_valid()  # ports are always self-consistent
    assert not p.is_linked()


# ---------------------------------------------------------------------------
# GraphNode: construction and port management (chinet test_node_init)
# ---------------------------------------------------------------------------


def test_node_init():
    node_with_ports = GraphNode.make_graph_node(
        "NodeName",
        {
            "inA": GraphPort(7.0),
            "inB": GraphPort(13.0),
            "outC": GraphPort(0.0, False, True),
        },
    )
    assert list(node_with_ports.get_input_ports().keys()) == ["inA", "inB"]
    values = [p.value for p in node_with_ports.get_ports().values()]
    assert np.allclose(sorted(values), [0.0, 7.0, 13.0])


def test_node_ports():
    node = GraphNode()
    node.add_input_port("portA", GraphPort(17))
    node.add_output_port("portB", GraphPort(23))

    assert node.get_input_port("portA").value == pytest.approx(17)
    assert node.get_output_port("portB").value == pytest.approx(23)
    assert node.get_port("portA") is not None
    assert node.get_port("missing") is None
    # add_port sets the output flag; set_ports names ports by their key.
    assert node.get_port("portA").is_output is False
    assert node.get_port("portB").is_output is True


def test_node_set_ports_names_ports():
    node = GraphNode.make_graph_node(
        "nd", {"x": GraphPort(1.0), "y": GraphPort(2.0, False, True)}
    )
    assert node.get_port("x").name == "x"
    assert node.get_port("y").name == "y"
    assert node.name == "nd"


# ---------------------------------------------------------------------------
# GraphNode: callbacks
# ---------------------------------------------------------------------------


def test_node_callback_types():
    node = GraphNode()
    node.set_callback("multiply_double", "C")
    assert node.get_callback() == "multiply_double"
    assert node.get_callback_type() == 0
    node.set_callback("multiply_double", "class")
    assert node.get_callback_type() == 1
    node.set_callback("multiply_double", "bogus")
    assert node.get_callback_type() == -1


def test_node_c_operator_callback():
    """The four chinet operator callbacks, over vector and scalar ports."""
    node = GraphNode("outA")
    node.add_input_port("inA", GraphPort([2.0, 3.0, 4.0]))
    node.add_input_port("inB", GraphPort([2.0, 3.0, 4.0]))
    node.add_output_port("outA", GraphPort(0.0, False, True))
    node.set_callback("multiply_double", "C")
    node.evaluate()
    assert np.allclose(
        node.get_output_port("outA").value,
        node.get_input_port("inA").value * node.get_input_port("inB").value,
    )


def test_node_c_operator_callback_addition():
    node = GraphNode("res")
    node.add_input_port("inA", GraphPort(2.5))
    node.add_input_port("inB", GraphPort(3.5))
    node.add_output_port("res", GraphPort(0.0, False, True))
    node.set_callback("addition_double", "C")
    node.evaluate()
    assert node.get_output_port("res").value == pytest.approx(6.0)


def test_node_c_operator_callback_int_truncates():
    node = GraphNode("res")
    node.add_input_port("inA", GraphPort(2.5))
    node.add_input_port("inB", GraphPort(2.5))
    node.add_output_port("res", GraphPort(0.0, False, True))
    node.set_callback("addition_int", "C")
    node.evaluate()
    assert node.get_output_port("res").value == pytest.approx(5)


def test_node_c_operator_scalar_broadcast():
    node = GraphNode("res")
    node.add_input_port("inA", GraphPort(2.0))
    node.add_input_port("inB", GraphPort([1.0, 2.0, 3.0]))
    node.add_output_port("res", GraphPort(0.0, False, True))
    node.set_callback("multiply_double", "C")
    node.evaluate()
    assert np.allclose(node.get_output_port("res").value, [2.0, 4.0, 6.0])


# ---------------------------------------------------------------------------
# GraphNode: validity and reactivity (chinet test_node_valid*)
# ---------------------------------------------------------------------------


def test_node_valid():
    """(inA)-(node)-(outA), non-reactive: invalid until evaluated."""
    node_1 = passthrough_node("node_1", 3.0)

    out_node_1 = node_1.get_output_port("node_1")
    in_node_1 = node_1.get_input_port("inA")
    assert out_node_1.value == pytest.approx(0.0)
    assert in_node_1.value == pytest.approx(3.0)
    assert node_1.is_valid() is False
    node_1.evaluate()
    assert node_1.is_valid() is True
    assert out_node_1.value == pytest.approx(3.0)


def test_node_valid_reactive_port():
    """A reactive input evaluates the node when its value changes."""
    node_1 = passthrough_node("node_1", 3.0)
    node_1.get_input_port("inA").reactive = True
    assert node_1.is_valid() is False

    node_1.get_input_port("inA").value = 12
    assert node_1.is_valid() is True
    assert node_1.get_output_port("node_1").value == pytest.approx(12)


def test_node_valid_connected_nodes():
    """Two chained nodes, non-reactive: writes invalidate, evaluate restores."""
    node_1 = passthrough_node("node_1", 3.0)
    in_node_1 = node_1.get_input_port("inA")
    out_node_1 = node_1.get_output_port("node_1")
    node_1.evaluate()
    assert out_node_1.value == pytest.approx(3.0)

    node_2 = passthrough_node("node_2", 13.0)
    in_node_2 = node_2.get_input_port("inA")
    out_node_2 = node_2.get_output_port("node_2")
    in_node_2.link = out_node_1

    assert node_2.is_valid() is False
    node_2.evaluate()
    assert node_2.is_valid() is True
    assert out_node_2.value == pytest.approx(3.0)  # follows the link

    # A change upstream invalidates both nodes.
    in_node_1.value = 13
    assert node_1.is_valid() is False
    assert node_2.is_valid() is False  # its input's source node is invalid

    node_1.evaluate()
    assert out_node_1.value == pytest.approx(13.0)

    node_2.evaluate()
    assert out_node_2.value == pytest.approx(13.0)


def test_node_valid_connected_nodes_reactive_ports():
    """A reactive chain: one write cascades evaluation through both nodes."""
    node_1 = passthrough_node("node_1", 3.0)
    node_2 = passthrough_node("node_2", 1.0)
    in_node_1 = node_1.get_input_port("inA")
    out_node_1 = node_1.get_output_port("node_1")
    out_node_2 = node_2.get_output_port("node_2")
    node_1.get_input_port("inA").reactive = True
    node_2.get_input_port("inA").reactive = True
    node_2.get_input_port("inA").link = out_node_1

    in_node_1.value = 13

    assert node_1.is_valid() is True
    assert node_2.is_valid() is True
    assert out_node_2.value == pytest.approx(13.0)


def test_node_update():
    """update() pulls linked inputs to their sources, evaluating upstream."""
    node_1 = passthrough_node("node_1", 3.0)
    node_2 = passthrough_node("node_2", 1.0)
    out_node_1 = node_1.get_output_port("node_1")
    out_node_2 = node_2.get_output_port("node_2")
    node_2.get_input_port("inA").link = out_node_1

    node_1.get_input_port("inA").value = 5  # non-reactive: both now stale
    node_2.update()  # recursively updates node_1, then itself

    assert node_1.is_valid() is True
    assert node_2.is_valid() is True
    assert out_node_1.value == pytest.approx(5.0)
    assert out_node_2.value == pytest.approx(5.0)


def test_node_no_callback_stays_invalid():
    """chinet: a node with neither callback object nor operator does nothing."""
    node = GraphNode("nd")
    node.add_input_port("inA", GraphPort(1.0))
    node.add_output_port("nd", GraphPort(0.0, False, True))
    node.evaluate()
    assert node.is_valid() is False


def test_node_without_inputs_is_valid():
    """chinet quirk: a node with no input ports reports valid."""
    node = GraphNode("empty")
    node.add_output_port("empty", GraphPort(0.0, False, True))
    assert node.is_valid() is True


def test_node_output_write_invalidates_sharer():
    """A port shared by two nodes: evaluating one invalidates the other."""
    n1 = passthrough_node("n1", 2.0)
    n2 = passthrough_node("n2", 3.0)
    shared = n1.get_output_port("n1")
    # n2 adopts n1's output port under a second key (the port's node
    # back-pointer now names n2, as in chinet's set_node).
    n2.add_output_port("n1", shared)
    n2.evaluate()
    assert n2.is_valid() is True
    n1.evaluate()  # writes the shared port
    assert n1.is_valid() is True
    assert n2.is_valid() is False  # its shared output moved under it


def test_node_metadata_and_describe():
    node = GraphNode("comp")
    node.add_input_port("inA", GraphPort(1.0))
    assert node.name == "comp"
    assert node.uid != GraphNode("other").uid
    assert "comp" in node.describe()
    assert list(node.inputs.keys()) == ["inA"]
    assert len(node.outputs) == 0
