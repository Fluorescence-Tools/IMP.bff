/*
 * The evaluation graph, and how a saved one is checked against the code that
 * rebuilds it.
 *
 * `GraphEvaluation` stores a provenance string per output and never looks
 * inside it. What a useful one contains is a Python question -- a node's
 * module, its class, and a hash of its source -- so it is answered here.
 */

%include "IMP/bff/GraphEvaluation.h"

%pythoncode %{
import hashlib as _eg_hashlib
import inspect as _eg_inspect
import sys as _eg_sys
import types as _eg_types


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
    back a `GraphNode` -- the same C++ node, same uid, and the Python override
    still runs when it is evaluated -- but *not* the Python object you
    registered, because a director's Python identity is not carried back out
    through the binding. So its class is `GraphNode` and its provenance would be
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

%pythoncode %{
# ---------------------------------------------------------------------------
# Carrying the code, for a graph that has to travel
# ---------------------------------------------------------------------------
#
# `node_provenance` above records *what the code was*. This records *the code*,
# so a graph of Python nodes can be rebuilt somewhere that does not have the
# module it came from.
#
# It is not the default and it is gated, because a document that carries code
# executes when it is loaded, and these files are shared and copied. The gate
# is one argument, not an obstacle: the caller says it trusts the document,
# which is a thing only the caller can know.
#
# The naive version does not work and it is worth saying why, because it is
# what one tries first: the class source alone rebuilds into a NameError on
# the first thing it uses. A callback closes over its imports and its module
# constants, and the source is not the closure. So the imports travel with it.


def _eg_names_used(cls):
    """Global names the class's methods reference."""
    out = set()
    for _, fn in vars(cls).items():
        code = getattr(fn, "__code__", None)
        if code is None:
            continue
        out.update(code.co_names)
        for const in code.co_consts:
            if isinstance(const, _eg_types.CodeType):
                out.update(const.co_names)
    return out


def node_code(node, source=None, imports=(), values=None):
    """What it takes to rebuild this node's class somewhere else.

    Returns a dict with `class`, `source`, `imports`, `globals` and
    `unresolved`.

    Introspection gets what it can and the caller supplements the rest. A
    class written in a notebook cell or by ``exec`` has neither readable
    source nor a module to harvest imports from, so `source`, `imports` and
    `values` are there to be given: they are added to whatever was found, and
    they are the difference between this working in a script and working
    everywhere.

    `imports` is every module the defining module imported, under the alias it
    used. Wholesale rather than only what the methods name, because a base
    class is evaluated when the class is created: ``class N(bff.GraphNode)`` never
    puts `bff` in any method's names, and capturing only those rebuilds into a
    NameError on the class statement itself.

    `globals` carries the module-level values the methods reference, as long
    as they are JSON-safe; a numpy array travels as a list and comes back as
    an array. **`unresolved` names what could not travel** -- a module-level
    object, an open file, a loaded 24 MB table -- and an empty list is the
    only claim that this is complete. Check it.
    """
    cls = type(node) if not isinstance(node, type) else node
    src = source if source is not None else _eg_inspect.getsource(cls)
    module = _eg_sys.modules.get(getattr(cls, "__module__", ""), None)
    found_imports, found_values, unresolved = list(imports), dict(values or {}), []
    used = _eg_names_used(cls)
    for name, value in sorted(vars(module).items()) if module else []:
        if name.startswith("__"):
            continue
        if isinstance(value, _eg_types.ModuleType):
            found_imports.append("import %s as %s" % (value.__name__, name))
            continue
        if name not in used:
            continue
        try:
            import numpy as _np
        except ImportError:
            _np = None
        if name in found_values:
            continue                       # the caller has already said
        if _np is not None and isinstance(value, _np.ndarray):
            found_values[name] = {"__ndarray__": value.tolist()}
        elif isinstance(value, (int, float, str, bool, list, dict, type(None))):
            found_values[name] = value
        elif not callable(value):
            unresolved.append(name)
    return {"class": cls.__name__, "source": src,
            "imports": sorted(set(found_imports)), "globals": found_values,
            "unresolved": unresolved}


def rebuild_node_class(code, trusted=False):
    """Recreate a class from `node_code`'s output.

    **This executes the code in the document.** `trusted` must be given
    explicitly, because whether a file is trustworthy is the caller's
    knowledge and nothing here can infer it: a graph read from your own run
    directory and one downloaded from a colleague look identical.
    """
    if not trusted:
        raise ValueError(
            "rebuild_node_class executes the source carried in the document, "
            "so it will not run unless you pass trusted=True. Rebuild the "
            "nodes yourself and use GraphEvaluation.from_json if you would "
            "rather the document stayed inert.")
    namespace = {}
    for statement in code["imports"]:
        exec(statement, namespace)
    for name, value in code.get("globals", {}).items():
        if isinstance(value, dict) and "__ndarray__" in value:
            import numpy as _np
            namespace[name] = _np.asarray(value["__ndarray__"])
        else:
            namespace[name] = value
    exec(code["source"], namespace)
    return namespace[code["class"]]


def graph_code(graph, nodes, sources=None):
    """The graph's naming, plus what it takes to rebuild each node's class.

    `nodes` maps node name to your node object -- your objects, because one
    fetched back from the graph has lost its Python class. `sources` may map a
    node name to its source for classes introspection cannot read.
    """
    sources = sources or {}
    import json as _json
    document = _json.loads(graph.to_json())
    document["code"] = {
        name: node_code(node, sources.get(name))
        for name, node in nodes.items()
    }
    return _json.dumps(document, indent=2)


def rebuild_node_classes(document, trusted=False):
    """GraphNode name -> class, from a `graph_code` document."""
    import json as _json
    doc = _json.loads(document) if isinstance(document, str) else document
    return {name: rebuild_node_class(code, trusted=trusted)
            for name, code in doc.get("code", {}).items()}
%}
