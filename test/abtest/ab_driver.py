"""Differential A/B driver: frozen Python chinet vs. IMP.bff's C++ re-implementation.

chinet was originally C++, crashed constantly and was rewritten in pure Python.
That Python is therefore the semantic truth, and bff's C++ port has to match it.
This module drives both runtimes through identical operation sequences and
compares every observable after every step.

Two documented deviations are mapped here rather than silently skipped:

* **Unbounded bounds.** chinet reports ``(None, None)`` when a port is not
  bounded; bff stores the bounds regardless and reports ``(nan, nan)``. Both
  are normalised to ``(None, None)`` when ``is_bounded`` is false.
* **Cycle errors.** chinet raises ``LinkCycleError``, bff raises the SWIG
  ``ValueError`` subclass. Both are recorded as ``"cycle-rejected"``.

Anything else that differs is a bug in bff and must be fixed there.
"""

import math

import numpy as np

import reference_chinet as rc
from IMP import bff

CYCLE_REJECTED = "cycle-rejected"


def _scalars(value):
    """Normalise a port value to a plain list of floats.

    SWIG hands back tuples, chinet hands back numpy arrays or bare scalars;
    neither difference is a semantic one.
    """
    arr = np.atleast_1d(np.asarray(value, dtype=np.float64))
    return [float(v) for v in arr.ravel()]


def _read(port, name):
    """Read an attribute the way whichever runtime exposes it.

    bff has ``get_<name>()`` accessors, chinet plain properties. Attribute
    lookup goes through the *type* because chinet's ``Port.__getattr__``
    forwards unknown names to the underlying numpy array, which would turn a
    missing accessor into a confusing numpy error instead of a clean miss.
    """
    if hasattr(type(port), "get_" + name):
        return getattr(port, "get_" + name)()
    value = getattr(port, name)
    return value() if callable(value) else value


def _write(port, name, value):
    """Write an attribute the way whichever runtime exposes it."""
    if hasattr(type(port), "set_" + name):
        getattr(port, "set_" + name)(value)
    else:
        setattr(port, name, value)


def _bound(value):
    if value is None:
        return None
    value = float(value)
    return None if math.isnan(value) else value


class PortAdapter:
    """One port, in one runtime, behind a runtime-neutral surface."""

    def __init__(self, port, label):
        self.port = port
        self.label = label

    # -- writes ---------------------------------------------------------
    def set_value(self, value):
        self.port.value = value

    def set_fixed(self, flag):
        _write(self.port, "fixed", bool(flag))

    def set_reactive(self, flag):
        _write(self.port, "is_reactive", bool(flag))

    def set_bounds(self, lb, ub):
        """bff takes two arguments, chinet a single pair."""
        try:
            self.port.set_bounds(lb, ub)
        except TypeError:
            self.port.set_bounds((lb, ub))

    def set_is_bounded(self, flag):
        _write(self.port, "is_bounded", bool(flag))

    def link(self, other):
        """Link to another port, reporting a refused cycle uniformly."""
        try:
            self.port.link = other.port
        except (ValueError, rc.LinkCycleError):
            return CYCLE_REJECTED
        return None

    def unlink(self):
        self.port.unlink()

    # -- reads ----------------------------------------------------------
    def _is_bounded(self):
        return bool(_read(self.port, "is_bounded"))

    def _is_linked(self):
        return bool(_read(self.port, "is_linked"))

    def _link_label(self, labels):
        """Identify the upstream by name.

        Identity is useless here: SWIG hands back a fresh proxy object for the
        same underlying C++ port on every call, so ``id()`` never matches.
        Names are assigned per label by the driver and compare across runtimes.
        """
        target = self.port.get_link()
        if target is None:
            return None
        return target.name

    def state(self, labels):
        bounded = self._is_bounded()
        lb, ub = (None, None)
        if bounded:
            lb, ub = _bound(self.port.bounds[0]), _bound(self.port.bounds[1])
        return {
            "name": self.port.name,
            "value": _scalars(self.port.value),
            "fixed": bool(self.port.fixed),
            "is_bounded": bounded,
            "bounds": (lb, ub),
            "is_output": bool(_read(self.port, "is_output")),
            "is_reactive": bool(_read(self.port, "is_reactive")),
            "is_linked": self._is_linked(),
            "link": self._link_label(labels),
        }


class NodeAdapter:
    def __init__(self, node, label):
        self.node = node
        self.label = label

    def add_input_port(self, key, port):
        self.node.add_input_port(key, port.port)

    def add_output_port(self, key, port):
        self.node.add_output_port(key, port.port)

    def state(self):
        return {"name": self.node.name, "valid": bool(self.node.is_valid())}


class Runtime:
    """A whole runtime (reference or bff) behind one neutral surface."""

    def __init__(self, kind):
        assert kind in ("ref", "bff")
        self.kind = kind
        self.ports = []
        self.nodes = []
        self._labels = {}

    def _module(self):
        return rc if self.kind == "ref" else bff

    def make_port(self, label, value=0.0, name=None, **kwargs):
        port = self._module().Port(value, name=name or label, **kwargs)
        adapter = PortAdapter(port, label)
        self.ports.append(adapter)
        self._labels[id(port)] = label
        return adapter

    def make_node(self, label, name=None):
        node = self._module().Node()
        node.name = name or label
        adapter = NodeAdapter(node, label)
        self.nodes.append(adapter)
        return adapter

    def port(self, label):
        return next(p for p in self.ports if p.label == label)

    def node(self, label):
        return next(n for n in self.nodes if n.label == label)

    def state(self):
        """Full observable state: every port, every node, keyed by label."""
        return {
            "ports": {p.label: p.state(self._labels) for p in self.ports},
            "nodes": {n.label: n.state() for n in self.nodes},
        }


class Pair:
    """The two runtimes, driven in lockstep."""

    def __init__(self):
        self.ref = Runtime("ref")
        self.bff = Runtime("bff")

    def both(self, method, *args, **kwargs):
        """Apply the same call to both runtimes, resolving port/node labels."""

        def resolve(runtime, value):
            if isinstance(value, str) and value.startswith("@"):
                label = value[1:]
                try:
                    return runtime.port(label)
                except StopIteration:
                    return runtime.node(label)
            return value

        results = []
        for runtime in (self.ref, self.bff):
            target = runtime
            if "." in method:
                holder, method_name = method.split(".", 1)
                target = resolve(runtime, holder)
            else:
                method_name = method
            bound = getattr(target, method_name)
            call_args = [resolve(runtime, a) for a in args]
            results.append(bound(*call_args, **kwargs))
        return results

    def compare(self, context=""):
        """Assert the two runtimes are observably identical."""
        ref_state, bff_state = self.ref.state(), self.bff.state()
        if ref_state != bff_state:
            raise AssertionError(_diff_message(ref_state, bff_state, context))


def _diff_message(ref_state, bff_state, context):
    lines = ["A/B mismatch" + (f" after {context}" if context else "")]
    for kind in ("ports", "nodes"):
        for label in sorted(set(ref_state[kind]) | set(bff_state[kind])):
            a = ref_state[kind].get(label)
            b = bff_state[kind].get(label)
            if a == b:
                continue
            lines.append(f"  {kind[:-1]} {label!r}:")
            keys = sorted(set(a or {}) | set(b or {}))
            for key in keys:
                if (a or {}).get(key) != (b or {}).get(key):
                    lines.append(
                        f"    {key}: ref={(a or {}).get(key)!r} bff={(b or {}).get(key)!r}"
                    )
    return "\n".join(lines)
