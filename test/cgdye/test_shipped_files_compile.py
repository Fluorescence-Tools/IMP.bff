"""Every shipped cgdye source compiles, is listed in the manifest, and imports.

A script that does not parse (scripts/rotamer_hgbp1_site481.py shipped with a
SyntaxError for a while) is invisible to the test suite unless something
compiles every file. The manifest check keeps pyext/src/Files.cmake honest
in both directions.
"""

import importlib
import os
import py_compile
import re
import subprocess
import sys
from pathlib import Path

import pytest

import IMP.bff.cgdye

PKG_DIR = Path(IMP.bff.cgdye.__file__).resolve().parent


def _manifest_path():
    for parent in Path(__file__).resolve().parents:
        cand = parent / "pyext" / "src" / "Files.cmake"
        if cand.exists():
            return cand
    return None


def _manifest_cgdye_files():
    manifest = _manifest_path()
    if manifest is None:
        pytest.skip("pyext/src/Files.cmake not found (installed tree)")
    text = manifest.read_text()
    m = re.search(r'set\(pyfiles "(.*?)"\)', text, re.S)
    assert m, "Files.cmake has no set(pyfiles ...)"
    return sorted(f for f in m.group(1).split(";") if f.startswith("cgdye/"))


def _disk_cgdye_files():
    out = []
    for path in PKG_DIR.rglob("*.py"):
        if "__pycache__" in path.parts:
            continue
        if path.name == "__init__.py":
            # IMP's build regenerates Files.cmake from disk and deliberately
            # leaves package __init__.py files out of it.
            continue
        out.append("cgdye/" + path.relative_to(PKG_DIR).as_posix())
    return sorted(out)


def test_every_cgdye_source_compiles():
    for path in PKG_DIR.rglob("*.py"):
        if "__pycache__" not in path.parts:
            py_compile.compile(str(path), doraise=True)


def test_manifest_matches_disk():
    manifest = set(_manifest_cgdye_files())
    disk = set(_disk_cgdye_files())
    assert manifest == disk, (
        f"missing from Files.cmake: {sorted(disk - manifest)}; "
        f"listed but not on disk: {sorted(manifest - disk)}")


def _importable_modules():
    mods = []
    for rel in _disk_cgdye_files():
        parts = rel[:-3].split("/")
        if "scripts" in parts:
            continue  # scripts have side effects; compile-only
        mods.append("IMP.bff." + ".".join(parts))
        mods.append("IMP.bff." + ".".join(parts[:-1]))  # the package itself
    return sorted(set(mods))


@pytest.mark.parametrize("module", _importable_modules())
def test_module_imports(module):
    importlib.import_module(module)


def test_library_imports_without_click():
    """click is a CLI dependency, not a library one (declared in conda-recipe)."""
    code = (
        "import sys; sys.modules['click'] = None\n"
        "import IMP.bff.cgdye\n"
        "import IMP.bff.cgdye.topology, IMP.bff.cgdye.sim\n"
        "print('ok')\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    assert "ok" in result.stdout


# IMP runs every .py under test/ as a standalone script, and a file of bare
# pytest functions would import cleanly and exit 0 -- reporting success without
# running a single assertion. Hand the file to pytest explicitly so a failure
# here is a failure in ctest.
if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
