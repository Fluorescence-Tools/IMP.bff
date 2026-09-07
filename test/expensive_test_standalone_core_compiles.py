"""The core compiles with no IMP on the include path.

PRD-137's independent core: every source outside `src/imp/` (and the layer's
bridge, `src/ImpLayer.cpp`) is compiled as one translation unit against
`standalone/include/` -- the module's own config header and the `IMP::Object`,
`IMP::Pointer`, `IMP::algebra::VectorD` and `IMP::PI` shims -- plus the
externals the core does use (Eigen, cereal, Boost.Random headers, RMF's
headers, python-ihm's C parser). IMP's own headers are nowhere on the path,
so a core source that reaches for one fails here first.

Expensive (a full syntax check of the core, a minute or two) and
environment-bound: skipped where a compiler or an external is not found.
"""
import glob
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _find_dir(candidates):
    for c in candidates:
        if c and os.path.isdir(c):
            return c
    return None


def _externals():
    prefix = sys.prefix
    eigen = _find_dir([os.path.join(prefix, "include", "eigen3"), "/usr/include/eigen3", "/opt/homebrew/include/eigen3"])
    inc = _find_dir([os.path.join(prefix, "include")])
    ok = inc and all(os.path.isdir(os.path.join(inc, d)) for d in ("cereal", "boost", "RMF"))
    ihm = _find_dir([os.path.join(_ROOT, "..", "imp", "modules", "core", "dependency", "python-ihm", "src"),
                     os.path.join(_ROOT, "src", "standalone", "ihm")])
    return (eigen, inc, ihm) if (eigen and ok and ihm) else None


class TestStandaloneCoreCompiles(unittest.TestCase):

    def test_the_core_needs_nothing_of_imp(self):
        cxx = shutil.which("c++") or shutil.which("clang++") or shutil.which("g++")
        ext = _externals()
        if not cxx or not ext:
            self.skipTest("no compiler, or Eigen/cereal/Boost/RMF/ihm not found in this environment")
        eigen, inc, ihm = ext
        tmp = tempfile.mkdtemp()
        try:
            shim = os.path.join(tmp, "inc", "IMP")
            os.makedirs(shim)
            os.symlink(os.path.join(_ROOT, "include"), os.path.join(shim, "bff"))
            tu = os.path.join(tmp, "core.cpp")
            with open(tu, "w") as fh:
                for src in sorted(glob.glob(os.path.join(_ROOT, "src", "*.cpp")) +
                                  glob.glob(os.path.join(_ROOT, "src", "standalone", "*.cpp"))):
                    if os.path.basename(src) == "ImpLayer.cpp":
                        continue
                    fh.write("#include <%s>\n" % src)
            cmd = [cxx, "-std=c++17", "-fsyntax-only", "-DIMPBFF_STANDALONE", "-DIMPBFF_COMPILATION",
                   "-I", os.path.join(_ROOT, "standalone", "include"), "-I", os.path.join(tmp, "inc"),
                   "-isystem", inc, "-isystem", eigen, "-I", ihm, "-Xclang", "-fopenmp", tu]
            r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            self.assertEqual(r.returncode, 0, "\n".join(l for l in r.stderr.splitlines() if "error" in l)[-4000:])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
