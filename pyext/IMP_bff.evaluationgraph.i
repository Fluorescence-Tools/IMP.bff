/*
 * The evaluation graph, and how a saved one is checked against the code that
 * rebuilds it.
 *
 * `EvaluationGraph` stores a provenance string per output and never looks
 * inside it. What a useful one contains is a Python question -- a node's
 * module, its class, and a hash of its source -- so it is answered here.
 */

%include "IMP/bff/EvaluationGraph.h"

%pythoncode %{
import hashlib as _eg_hashlib
import inspect as _eg_inspect


def node_provenance(node):
    """Identify the code behind a node, as well as it can be identified.

    Returns ``"<module>:<qualname>:<sha1 of the source, 12 hex>"``, with the
    hash left empty when the source cannot be read. It is a *fingerprint*, not
    a recipe: a saved graph records what it was computed with, so that loading
    it against changed code can say so instead of quietly computing something
    else.

    **Why not save the source and rebuild from it.** Measured, and it does not
    work. `inspect.getsource` raises "is a built-in class" for a class defined
    in a notebook cell or by ``exec``, which is where models of this kind are
    usually written. Where it does succeed, executing the text in a fresh
    namespace fails on the first name it does not carry -- the source is not
    the closure, and a callback closes over its imports, its module constants
    and everything it calls. Capturing all of that is what ``cloudpickle`` and
    ``dill`` do, and they break across library versions.

    And a document that carries code executes when it is loaded. These are
    shared between colleagues and copied off servers, so a graph file is
    deliberately inert.
    """
    cls = type(node)
    try:
        source = _eg_inspect.getsource(cls)
        digest = _eg_hashlib.sha1(source.encode("utf-8")).hexdigest()[:12]
    except (OSError, TypeError):
        # A class from a notebook cell, a REPL, an exec, or a byte-compiled
        # module with no source beside it. The name still identifies it.
        digest = ""
    return "%s:%s:%s" % (getattr(cls, "__module__", "?"),
                         getattr(cls, "__qualname__", cls.__name__), digest)


def add_output_with_provenance(graph, label, node, port_name):
    """`add_output`, recording what the node's code was at the time."""
    graph.add_output(label, node, port_name, node_provenance(node))


def check_provenance(graph, nodes):
    """Labels whose node's code differs from what the graph recorded.

    `nodes` maps label to **your** node object. Returns a dict of
    ``label -> (recorded, now)`` for the ones that disagree, empty when they
    all match -- so ``assert not check_provenance(g, nodes)`` is the whole
    check. A label with no recorded provenance, or one absent from `nodes`, is
    skipped rather than reported: nothing was promised about it.

    **Your objects, not the graph's.** `graph.get_output_node(label)` hands
    back a `Node` -- the same C++ node, same uid, and the Python override
    still runs when it is evaluated -- but *not* the Python object you
    registered, because a director's Python identity is not carried back out
    through the binding. So its class is `Node` and its provenance would be
    meaningless. This function therefore asks for the nodes you rebuilt, which
    you have, since rebuilding them is what made the graph loadable at all.
    """
    out = {}
    for label in graph.get_output_labels():
        recorded = graph.get_output_provenance(label)
        if not recorded:
            continue
        node = nodes.get(label) if hasattr(nodes, "get") else None
        if node is None:
            continue
        now = node_provenance(node)
        if now != recorded:
            out[label] = (recorded, now)
    return out
%}
