"""The subpackage import graph: acyclic, standalone, and no reaching into privates.

This is the test that should have existed already. `import IMP.bff.representation.av` as a
process's first import raised ImportError at every commit before 2026-08-18 --
``representation/__init__`` imports ``distribution``, which imports ``av``,
which imports ``representation`` -- and nothing in 676 tests caught it, because
every test imports ``IMP.bff`` first and that happens to resolve the packages in
an order the cycle survives.

Three properties, and the first two are different:

* **Acyclic at module scope.** A cycle can be broken by deferring an import into
  a function, which is a legitimate fix when the layering is genuinely one-way
  and only the package ``__init__`` closes the loop. What must not exist is a
  cycle among *module-level* imports.
* **Standalone importable.** Even an acyclic graph breaks if a package
  ``__init__`` pulls in a module that imports back through the package. Only
  starting a fresh interpreter proves this, so that is what the test does.
* **No reaching into another domain's privates.** ``from IMP.bff.quenching.pet
  import _foo`` couples two domains through something neither promised to keep.
"""

import ast
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pytest

SOURCE = Path(__file__).resolve().parent.parent / "pyext" / "src"

#: The domains. Everything else is a loose module at the package root.
DOMAINS = sorted(p.name for p in SOURCE.iterdir()
                 if p.is_dir() and (p / "__init__.py").exists())


def _is_type_checking(node) -> bool:
    """``if TYPE_CHECKING:`` -- an annotation-only import, not a runtime edge.

    It is still worth *writing* such imports, which is why the graph ignores
    them rather than the code avoiding them: a name used only in a string
    annotation costs nothing at import time and everything to a type checker.
    """
    if not isinstance(node, ast.If):
        return False
    test = node.test
    return ((isinstance(test, ast.Name) and test.id == "TYPE_CHECKING")
            or (isinstance(test, ast.Attribute) and test.attr == "TYPE_CHECKING"))


def _module_level_imports(path: Path):
    """``IMP.bff.*`` targets imported at module scope, with the alias names."""
    tree = ast.parse(path.read_text())
    out = []
    for node in tree.body:                      # body, not walk: module scope only
        if isinstance(node, ast.ImportFrom) and node.module:
            if node.level:                      # relative: resolve against the file
                parts = path.relative_to(SOURCE).parent.parts
                base = parts[:len(parts) - (node.level - 1)] if node.level > 1 else parts
                target = ".".join(("IMP", "bff") + tuple(base) + (node.module,))
            else:
                target = node.module
            out.append((target, [a.name for a in node.names]))
        elif isinstance(node, ast.Import):
            for a in node.names:
                out.append((a.name, []))
        elif isinstance(node, (ast.Try, ast.If)):
            if _is_type_checking(node):
                continue                        # provably never executed
            for sub in ast.walk(node):          # guarded imports are still module scope
                if isinstance(sub, ast.ImportFrom) and sub.module and not sub.level:
                    out.append((sub.module, [a.name for a in sub.names]))
                elif isinstance(sub, ast.Import):
                    out.extend((a.name, []) for a in sub.names)
    return out


def _relative_imports(path: Path):
    """Relative ``from .. import`` targets, which the resolver above approximates."""
    tree = ast.parse(path.read_text())
    return [n for n in tree.body if isinstance(n, ast.ImportFrom) and n.level]


def _domain_of(dotted: str):
    parts = dotted.split(".")
    if len(parts) >= 3 and parts[0] == "IMP" and parts[1] == "bff":
        return parts[2] if parts[2] in DOMAINS else None
    return None


def _graph():
    edges = defaultdict(set)
    for path in sorted(SOURCE.rglob("*.py")):
        rel = path.relative_to(SOURCE)
        if not rel.parts or rel.parts[0] not in DOMAINS:
            continue
        here = rel.parts[0]
        for target, _names in _module_level_imports(path):
            there = _domain_of(target)
            if there and there != here:
                edges[here].add(there)
    return edges


def test_the_domains_are_what_we_think_they_are():
    """A guard on the test itself: if this list empties, the rest passes vacuously.

    The four stages -- representation, scoring, sampling, analysis -- plus the
    things that cut across every model: the species, where it is attached, the
    photophysics, the output contract, the formats, and scoring against data.
    """
    assert len(DOMAINS) >= 6, DOMAINS
    for expected in ("representation", "scoring", "sampling", "analysis",
                     "dye", "label", "photophysics", "observables",
                     "io", "restraints"):
        assert expected in DOMAINS, (expected, DOMAINS)


def test_the_subpackage_graph_is_acyclic():
    edges = _graph()
    colour = {}
    stack = []

    def visit(node):
        colour[node] = "grey"
        stack.append(node)
        for nxt in sorted(edges.get(node, ())):
            if colour.get(nxt) == "grey":
                cycle = stack[stack.index(nxt):] + [nxt]
                raise AssertionError("import cycle: " + " -> ".join(cycle))
            if colour.get(nxt) is None:
                visit(nxt)
        colour[node] = "black"
        stack.pop()

    for node in sorted(edges):
        if colour.get(node) is None:
            visit(node)


@pytest.mark.parametrize("domain", DOMAINS)
def test_each_subpackage_imports_on_its_own(domain):
    """In a fresh interpreter, as the *first* IMP.bff import.

    Importing ``IMP.bff`` first hides this, which is exactly how the
    ``av``/``representation`` cycle survived every test in the suite.
    """
    result = subprocess.run(
        [sys.executable, "-c", f"import IMP.bff.{domain}; print('ok')"],
        capture_output=True, text=True, timeout=300)
    assert result.returncode == 0, result.stderr[-2000:]
    assert "ok" in result.stdout


def test_no_domain_imports_another_domains_private_names():
    offenders = []
    for path in sorted(SOURCE.rglob("*.py")):
        rel = path.relative_to(SOURCE)
        if not rel.parts or rel.parts[0] not in DOMAINS:
            continue
        here = rel.parts[0]
        for target, names in _module_level_imports(path):
            there = _domain_of(target)
            if there is None or there == here:
                continue
            private = [n for n in names if n.startswith("_")]
            # a private *module* is as much of a promise as a private name
            if target.split(".")[-1].startswith("_"):
                private.append(target)
            if private:
                offenders.append(f"{rel} imports {private} from {target}")
    assert not offenders, offenders


def test_relative_imports_do_not_climb_out_of_their_domain():
    """``from ... import x`` reaches the package root and hides the coupling.

    One level (``from .foo``) is a sibling; two (``from ..foo``) is the domain's
    own parent, which for a nested module is still its domain. Beyond that the
    import is cross-domain but does not look it.
    """
    offenders = []
    for path in sorted(SOURCE.rglob("*.py")):
        rel = path.relative_to(SOURCE)
        if not rel.parts or rel.parts[0] not in DOMAINS:
            continue
        depth = len(rel.parts) - 1          # directories below pyext/src
        for node in _relative_imports(path):
            if node.level > depth:
                offenders.append(
                    f"{rel}: 'from {'.' * node.level}{node.module or ''} import ...' "
                    f"climbs past its domain")
    assert not offenders, offenders


def test_every_module_imports():
    """Import every file in the package, in one fresh interpreter.

    This is the test that was missing. Four modules carried stale
    sibling-relative imports from before ``fret/`` was decomposed --
    ``observables/pair_distribution.py`` broke at module scope and could not be
    imported at all, and ``restraints/docking.py`` had three that would have
    failed at call time. **808 tests passed over them**, because nothing in the
    suite imports those two modules.

    A module that no test exercises still has to import. That is the cheapest
    possible check and it costs one subprocess.
    """
    modules = []
    for path in sorted(SOURCE.rglob("*.py")):
        if "__pycache__" in path.parts:
            continue
        rel = path.relative_to(SOURCE).with_suffix("")
        parts = [p for p in rel.parts if p != "__init__"]
        modules.append("IMP.bff" + ("." + ".".join(parts) if parts else ""))
    assert len(modules) > 50, "the source scan found almost nothing"

    code = (
        "import importlib, sys\n"
        f"broken = []\n"
        f"for name in {modules!r}:\n"
        "    try:\n"
        "        importlib.import_module(name)\n"
        "    except Exception as exc:\n"
        "        broken.append(f'{name}: {type(exc).__name__}: {exc}')\n"
        "sys.stdout.write(chr(10).join(broken))\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True,
                            text=True, timeout=600)
    assert result.returncode == 0, result.stderr[-2000:]
    assert not result.stdout.strip(), result.stdout


def test_test_module_basenames_are_unique():
    """pytest imports test modules by basename, so a collision silently drops one.

    ``test/io/test_io.py`` collided with ``test/cgdye/test_io.py`` and was
    reported as a collection ERROR at the very bottom of a 780-test run -- easy
    to read past, and the file's assertions simply did not run.
    """
    from collections import defaultdict

    root = Path(__file__).resolve().parent
    by_name = defaultdict(list)
    for path in sorted(root.rglob("test_*.py")):
        if "__pycache__" in path.parts:
            continue
        by_name[path.name].append(str(path.relative_to(root)))
    clashes = {name: paths for name, paths in by_name.items() if len(paths) > 1}
    assert not clashes, clashes


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
