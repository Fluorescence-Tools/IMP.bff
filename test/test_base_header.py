"""`IMPCompatibility.h`: the module's vocabulary compiles with no IMP on the include path.

Nearly every header here says `IMP_THROW`, `IMP_SHOWABLE_INLINE` or
`IMP_VALUES`. `include/IMPCompatibility.h` is the one door those come through: with IMP
present it forwards to IMP's definitions, and with `IMPBFF_STANDALONE`
defined it supplies equivalents of its own. That second branch is the whole
basis of the independent core, and nothing in the ordinary build ever
compiles it -- so this test does, with the system compiler, against an
include path that contains this module's headers and **nothing of IMP's**.

Three things are pinned:

* the standalone branch compiles and runs -- the plural typedef, the
  `show()` method and `IMP_THROW` all work, and the exceptions are catchable
  as `std::exception`, which is what a caller who has never heard of IMP
  expects;
* it is the flag that switches, not luck -- the same TU **without**
  `IMPBFF_STANDALONE` and without IMP's headers must fail to compile;
* the header carries no IMP include in its standalone branch, checked at the
  source level so it holds even where no compiler is available.

Skipped, not failed, where there is no C++ compiler: this is a property of
the header, and a machine that cannot compile C++ cannot refute it.
"""

import os
import re
import shutil
import subprocess
import tempfile
import unittest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_BASE_H = os.path.join(_ROOT, "include", "IMPCompatibility.h")
_SHIMS = os.path.join(_ROOT, "standalone", "include")  # the standalone build's IMP/ tree

_CONFIG_STUB = """\
#ifndef IMPBFF_CONFIG_H
#define IMPBFF_CONFIG_H
#define IMPBFFEXPORT
#define IMPBFF_BEGIN_NAMESPACE namespace IMP { namespace bff {
#define IMPBFF_END_NAMESPACE } }
#endif
"""

_TU = r"""
#include <IMP/bff/IMPCompatibility.h>
#include <iostream>
#include <vector>

IMPBFF_BEGIN_NAMESPACE
struct Thing {
    int a;
    Thing(int a = 0) : a(a) {}
    IMP_SHOWABLE_INLINE(Thing, out << "Thing(" << a << ")");
};
IMP_VALUES(Thing, Things);

inline void boom(int x) {
    if (x < 0) IMP_THROW("x is " << x << ", which is negative", ValueException);
}
inline void readfail() { IMP_THROW("no such file", IOException); }
IMPBFF_END_NAMESPACE

int main() {
    IMP::bff::Things ts;
    ts.push_back(IMP::bff::Thing(7));
    std::cout << ts[0] << " n=" << ts.size() << "\n";
    try { IMP::bff::boom(-3); }
    catch (const IMP::ValueException &e) { std::cout << "VE:" << e.what(); }
    try { IMP::bff::readfail(); }
    catch (const std::exception &e) { std::cout << "SE:" << e.what(); }
    return 0;
}
"""


def _compiler():
    for c in ("c++", "clang++", "g++"):
        p = shutil.which(c)
        if p:
            return p
    return None


class TestBaseHeaderStandalone(unittest.TestCase):

    def setUp(self):
        self.cxx = _compiler()
        self.tmp = tempfile.mkdtemp()
        inc = os.path.join(self.tmp, "IMP", "bff")
        os.makedirs(inc)
        shutil.copy(_BASE_H, os.path.join(inc, "IMPCompatibility.h"))
        # the shim tree the standalone build puts on the path: IMP/Object.h,
        # IMP/Pointer.h, IMP/constants.h, IMP/algebra/ -- and nothing of IMP's
        for name in os.listdir(os.path.join(_SHIMS, "IMP")):
            if name == "bff":
                continue  # the stub config above stands in for the shim tree's
            src = os.path.join(_SHIMS, "IMP", name)
            dst = os.path.join(self.tmp, "IMP", name)
            if os.path.isdir(src):
                shutil.copytree(src, dst)
            elif name.endswith(".h"):
                shutil.copy(src, dst)
        with open(os.path.join(inc, "bff_config.h"), "w") as fh:
            fh.write(_CONFIG_STUB)
        self.src = os.path.join(self.tmp, "tu.cpp")
        with open(self.src, "w") as fh:
            fh.write(_TU)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _compile(self, *extra):
        exe = os.path.join(self.tmp, "tu")
        cmd = [self.cxx, "-std=c++14", "-I", self.tmp, "-o", exe, self.src]
        cmd.extend(extra)
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           text=True)
        return r.returncode, r.stderr, exe

    def test_standalone_branch_compiles_and_runs(self):
        if not self.cxx:
            self.skipTest("no C++ compiler on PATH")
        rc, err, exe = self._compile("-DIMPBFF_STANDALONE")
        self.assertEqual(rc, 0, err[-2000:])
        out = subprocess.run([exe], stdout=subprocess.PIPE, text=True).stdout
        self.assertIn("Thing(7) n=1", out)          # IMP_VALUES + show()
        self.assertIn("VE:x is -3, which is negative", out)  # IMP_THROW, typed
        self.assertIn("SE:no such file", out)       # catchable as std::exception

    def test_without_the_flag_and_without_imp_it_does_not_compile(self):
        # Proves the flag is the switch. If this ever starts passing, IMPCompatibility.h
        # has grown an IMP-free path that is on by default -- which would be a
        # different design, and one that should be chosen, not stumbled into.
        if not self.cxx:
            self.skipTest("no C++ compiler on PATH")
        rc, err, _ = self._compile()
        self.assertNotEqual(rc, 0)
        self.assertIn("IMP/exception.h", err)

    def test_standalone_branch_includes_only_the_shims(self):
        # Source-level, so it holds on a machine with no compiler too. The
        # branch is everything from `#else` to the closing `#endif`. Every
        # `<IMP/...>` it includes must be a file of the standalone shim tree
        # (`standalone/include/`), never one of IMP's own headers.
        src = open(_BASE_H).read()
        m = re.search(r"#ifndef IMPBFF_STANDALONE(.*?)#else(.*?)#endif  // IMPBFF_STANDALONE",
                      src, re.S)
        self.assertIsNotNone(m, "IMPCompatibility.h no longer has the two-branch shape")
        imp_branch, standalone_branch = m.group(1), m.group(2)
        self.assertRegex(imp_branch, r"#include <IMP/exception\.h>")
        for inc in re.findall(r"#include <(IMP/[^>]+)>", standalone_branch):
            self.assertTrue(os.path.isfile(os.path.join(_SHIMS, inc)),
                            "%s is not a standalone shim" % inc)


if __name__ == "__main__":
    unittest.main()
