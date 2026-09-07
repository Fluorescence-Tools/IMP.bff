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
import sys
from pathlib import Path

import pytest

_PROGRAM = Path(__file__).resolve().parent.parent / "bin" / "imp_bff"


@pytest.fixture(scope="session")
def imp_bff_program():
    if "imp_bff_program" not in sys.modules:
        loader = importlib.machinery.SourceFileLoader("imp_bff_program", str(_PROGRAM))
        spec = importlib.util.spec_from_file_location(
            "imp_bff_program", _PROGRAM, loader=loader)
        module = importlib.util.module_from_spec(spec)
        sys.modules["imp_bff_program"] = module
        spec.loader.exec_module(module)
    return sys.modules["imp_bff_program"]
