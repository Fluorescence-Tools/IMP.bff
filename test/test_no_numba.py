"""numba is gone: it must not be imported, and blocking it must change nothing.

PRD-113's rule is that what runs in CI is what runs for users -- no optional
accelerator that silently swaps in a slower path, and no numerics outside C++.
``IMP.bff._jit`` was the shim that made numba optional; the shim is deleted
because there is nothing left to make optional.

The check that matters is not "numba is absent from the environment" -- it is
installed here, and on many users' machines. It is that ``IMP.bff`` never
reaches for it.
"""

import ast
import importlib
import subprocess
import sys
from pathlib import Path

import pytest

import IMP.bff


def _sources():
    """The checked-in sources, not the build tree.

    The build tree is a directory of symlinks, and a deleted source leaves a
    dangling one behind -- ``spectroscopy/`` did, months after stage 0 removed
    it. Scanning the source is both the right question and the robust one.
    """
    root = Path(__file__).resolve().parent.parent / "pyext" / "src"
    assert root.is_dir(), root
    return sorted(root.rglob("*.py"))


def test_the_jit_shim_is_gone():
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("IMP.bff._jit")


def test_no_module_imports_numba_or_the_shim():
    offenders = []
    for path in _sources():
        tree = ast.parse(path.read_text())
        names = {n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.module}
        names |= {a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
        if any("numba" in n or n.endswith("_jit") for n in names):
            offenders.append(path.name)
    assert not offenders, offenders


def test_no_jit_decorator_survives():
    offenders = []
    for path in _sources():
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            for dec in node.decorator_list:
                text = ast.dump(dec)
                if "njit" in text or "numba" in text or "'jit'" in text:
                    offenders.append(f"{path.name}::{node.name}")
    assert not offenders, offenders


def test_the_package_works_with_numba_blocked():
    """A fresh interpreter that cannot import numba at all."""
    code = (
        "import sys\n"
        "sys.modules['numba'] = None\n"
        "import numpy as np\n"
        "import IMP.bff\n"
        "from IMP.bff.av import _kernels\n"
        "from IMP.bff.photophysics import orientation\n"
        "from IMP.bff.quenching import diffusion, photon\n"
        "from IMP.bff.representation import distance\n"
        "assert 'numba' not in [m for m in sys.modules if sys.modules[m] is not None "
        "and m == 'numba']\n"
        # exercise one kernel from each ported domain
        "pts = np.array([[0., 0., 0., 1.], [2., 0., 0., 1.]])\n"
        "assert _kernels.weighted_mean(pts, 2)[0] == 1.0\n"
        "assert abs(orientation.kappasq(0.0, 0.0, 0.0, 0.0, 0.0) - 2/3) < 1e-12\n"
        "d, _ = photon.simulate_photon_trace(100, np.zeros(10), 0.01, 4.0, random_seed=1)\n"
        "assert (d >= 0).all()\n"
        "print('ok')\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True,
                            text=True, timeout=300)
    assert result.returncode == 0, result.stderr
    assert "ok" in result.stdout


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
