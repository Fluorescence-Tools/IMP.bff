"""tttrlib is a *soft* dependency of IMP.bff: most of IMP.bff works without it.

IMP.bff sits above tttrlib in the layering and may reach down into it — tttrlib
is the lower-level library, and delegating photon-level work there is better
than duplicating it. But IMP.bff ships through conda-forge as part of IMP, and
its job is structure, dye simulation, photophysics and scoring, none of which
needs photons. So tttrlib is optional, and the core API has to work with it
absent.

Softness of this kind hardens quietly. One module-level ``import tttrlib`` on a
path the core API touches is enough, it passes every test on a developer machine
where tttrlib is installed, and it surfaces only for a user who does not have
it. Hence a test rather than a convention.

The rule is written down in ../chisurf/okf/references/imp-ecosystem.md.
"""

import ast
import os
import subprocess
import sys
import unittest

import IMP
import IMP.test


def _module_roots():
    """Every directory holding IMP.bff's own Python sources."""
    here = os.path.dirname(os.path.abspath(__file__))
    return [os.path.join(os.path.dirname(here), "pyext", "src")]


def _unguarded_imports(path):
    """Imports of tttrlib that are *not* inside a try/except or a function.

    A guarded import (inside ``try:`` or inside a function body) is fine: it is
    reached only when the caller asked for the feature. A module-level bare
    import is not, because merely importing IMP.bff would then require tttrlib.
    """
    with open(path, "r", errors="replace") as fh:
        try:
            tree = ast.parse(fh.read())
        except SyntaxError:
            return []

    guarded = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Try, ast.FunctionDef, ast.AsyncFunctionDef)):
            for child in ast.walk(node):
                guarded.add(id(child))

    bad = []
    for node in ast.walk(tree):
        if id(node) in guarded:
            continue
        names = []
        if isinstance(node, ast.Import):
            names = [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            names = [node.module]
        if any(n.split(".")[0] == "tttrlib" for n in names):
            bad.append(getattr(node, "lineno", 0))
    return bad


class Tests(IMP.test.TestCase):

    def test_no_unguarded_tttrlib_import(self):
        """IMP.bff sources must not import tttrlib at module level"""
        offenders = {}
        for root in _module_roots():
            for dirpath, dirnames, filenames in os.walk(root):
                dirnames[:] = [d for d in dirnames if d != "__pycache__"]
                for fn in filenames:
                    if not fn.endswith(".py"):
                        continue
                    p = os.path.join(dirpath, fn)
                    lines = _unguarded_imports(p)
                    if lines:
                        offenders[os.path.relpath(p, root)] = lines
        self.assertEqual(
            offenders, {},
            "tttrlib is an optional dependency of IMP.bff: most of IMP.bff has "
            "to work without it. Put the import inside the function that needs "
            "it, or behind try/except. Offenders (file: lines): %s" % offenders)

    def test_imports_with_tttrlib_hidden(self):
        """IMP.bff imports and scores with tttrlib made unavailable"""
        # A fresh interpreter with tttrlib blocked at the finder level, so this
        # holds even on a machine that has it installed.
        script = (
            "import sys\n"
            "class Block:\n"
            "    def find_module(self, name, path=None):\n"
            "        return self if name.split('.')[0] == 'tttrlib' else None\n"
            "    def find_spec(self, name, path=None, target=None):\n"
            "        if name.split('.')[0] == 'tttrlib':\n"
            "            raise ImportError('tttrlib blocked for this test')\n"
            "        return None\n"
            "sys.meta_path.insert(0, Block())\n"
            "import IMP, IMP.bff\n"
            "m = IMP.Model()\n"
            "assert IMP.bff.ProbeAccessibleVolumeDecorator is not None\n"
            "assert IMP.bff.PathMap is not None\n"
            "print('ok')\n"
        )
        out = subprocess.run([sys.executable, "-c", script],
                             capture_output=True, text=True)
        self.assertEqual(
            out.returncode, 0,
            "IMP.bff must import and expose its core API without tttrlib.\n"
            "stdout: %s\nstderr: %s" % (out.stdout, out.stderr))
        self.assertIn("ok", out.stdout)


if __name__ == '__main__':
    IMP.test.main()
