"""Shared fixtures.

``imp_bff_program`` loads ``bin/imp_bff`` as a module. The command tree lives
there rather than in the package because a click command is a decorated
function, and a library module carrying one cannot be imported without click.
Tests that exercise commands therefore load the program, and need an explicit
loader -- IMP's installed programs have no file extension, and
``spec_from_file_location`` cannot infer a loader without one.
"""

import importlib.machinery
import importlib.util
import re
import sys
from pathlib import Path

import pytest

_PROGRAM = Path(__file__).resolve().parent.parent / "bin" / "imp_bff"

# --- the standalone lane -----------------------------------------------------
# The core builds without IMP (PRD-137). Under that build `import IMP.bff`
# works and nothing else of IMP exists, so a test file that names another IMP
# module -- IMP.atom, IMP.Model, IMP.test -- cannot even be imported. Those
# files are the connection layer's tests; they are deselected here, not
# forked, so the one suite serves both builds. Two rules, both cheap: any
# `IMP.<something other than bff>` in the text, or a `from IMP.bff import`
# of a name the standalone module does not have (a layer function).
_STANDALONE = importlib.util.find_spec("IMP.atom") is None
_NEEDS_IMP = re.compile(r"\bIMP\.(?!bff\b)[A-Za-z_]")
_FROM_BFF = re.compile(r"^\s*from IMP\.bff import \(?([^)\n]*(?:\n[^)]*)?)\)?", re.M)


_BFF_ALIAS = re.compile(r"^\s*(?:import IMP\.bff as (\w+)|from IMP import bff(?: as (\w+))?)", re.M)


def _uses_a_layer_name(text):
    """True when the file names something IMP.bff only has under IMP: a
    `from IMP.bff import x`, or an attribute `IMP.bff.x` / `bff.x` through
    any alias the file declares, where x is not in the standalone module."""
    import IMP.bff
    for m in _FROM_BFF.finditer(text):
        for name in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", m.group(1).split("#")[0]):
            if name != "as" and not hasattr(IMP.bff, name):
                return True
    aliases = {"IMP.bff"}
    for m in _BFF_ALIAS.finditer(text):
        aliases.add(m.group(1) or m.group(2) or "bff")
    for alias in aliases:
        for name in set(re.findall(r"\b" + re.escape(alias) + r"\.([A-Za-z_]\w*)", text)):
            if not hasattr(IMP.bff, name):
                return True
    return False


def pytest_ignore_collect(collection_path, config):
    if not _STANDALONE or collection_path.suffix != ".py":
        return None
    try:
        text = collection_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    if _NEEDS_IMP.search(text) or _uses_a_layer_name(text):
        return True
    return None


def pytest_report_header(config):
    if _STANDALONE:
        return "IMP.bff lane: standalone (no IMP; files naming other IMP modules are not collected)"
    return "IMP.bff lane: IMP build"


@pytest.fixture(scope="session")
def imp_bff_program():
    if "imp_bff_program" not in sys.modules:
        loader = importlib.machinery.SourceFileLoader("imp_bff_program", str(_PROGRAM))
        spec = importlib.util.spec_from_file_location(
            "imp_bff_program", _PROGRAM, loader=loader)
        module = importlib.util.module_from_spec(spec)
        sys.modules["imp_bff_program"] = module
        try:
            spec.loader.exec_module(module)
        except ModuleNotFoundError as e:
            del sys.modules["imp_bff_program"]
            if _STANDALONE and e.name and e.name.startswith("IMP"):
                pytest.skip("bin/imp_bff imports %s; the standalone lane has no IMP" % e.name)
            raise
    return sys.modules["imp_bff_program"]
