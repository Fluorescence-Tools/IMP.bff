"""Every shipped Python source compiles: the example scripts and the programs.

A script that does not parse is invisible to the test suite unless something
compiles every file -- ``examples/structure/rotamer_hgbp1_site481.py`` shipped
with a SyntaxError for a while. The manifest check this file used to carry died
with ``pyext/src``: the package's Python is the `.i` files' ``%pythoncode``
now, which the build compiles into ``__init__.py`` and imports with every test.
"""

import py_compile
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]


def _python_sources():
    out = []
    for base in ("examples", "bin"):
        root = REPO / base
        if not root.is_dir():
            continue
        for path in sorted(root.rglob("*")):
            if "__pycache__" in path.parts or path.is_dir():
                continue
            if path.suffix == ".py" or (base == "bin" and path.is_file()
                                        and not path.name.startswith(".")
                                        and "CMakeLists" not in path.name
                                        and "Files.cmake" not in path.name):
                out.append(path)
    assert out, "no shipped Python sources found"
    return out


@pytest.mark.parametrize("path", _python_sources(), ids=str)
def test_source_compiles(path):
    py_compile.compile(str(path), doraise=True)


def test_the_programs_import_the_flat_surface():
    """The programs reach every name flat -- no submodule survived pyext/src."""
    import ast
    offenders = []
    for path in _python_sources():
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            target = None
            if isinstance(node, ast.ImportFrom) and node.module:
                target = node.module
            elif isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name.startswith("IMP.bff."):
                        offenders.append(f"{path.name}: import {alias.name}")
            if target and target.startswith("IMP.bff."):
                suffix = target[len("IMP.bff."):]
                if "." in suffix or suffix != "":
                    offenders.append(f"{path.name}: from {target} import ...")
    assert not offenders, offenders


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q", "-p no:cacheprovider"]))
