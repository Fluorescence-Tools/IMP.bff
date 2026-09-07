"""The vendored expression engine must be byte-identical to tttrlib's.

`include/internal/ExpressionEngine.h` is a copy of tttrlib's
`modules/core/include/ExpressionEngine.h`. Copied rather than depended on:
the engine is header-only and pure arithmetic, so there is nothing to link,
and taking it from an installed tttrlib would still have made this repository
need tttrlib present to build a file that has no dependencies of its own.

The obvious hazard of a copy is that it drifts, and a copy that has silently
drifted is worse than either a dependency or a fork -- it looks like one
implementation while behaving as two. So this is checked rather than trusted.

**tttrlib owns the original.** A change goes there and is copied here:

    cp ../tttrlib/modules/core/include/ExpressionEngine.h \\
       include/internal/ExpressionEngine.h

Skipped, not failed, when the sibling checkout is absent -- a release tarball
or a CI job that builds only imp.bff cannot run this, and should not be
reported as broken for it.
"""

import hashlib
import pathlib
import unittest

HERE = pathlib.Path(__file__).resolve()
VENDORED = HERE.parents[2] / "include" / "internal" / "ExpressionEngine.h"
UPSTREAM = (HERE.parents[3] / "tttrlib" / "modules" / "core" / "include"
            / "ExpressionEngine.h")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class EngineCopyTests(unittest.TestCase):

    def test_the_vendored_copy_exists(self):
        self.assertTrue(VENDORED.is_file(), f"missing: {VENDORED}")

    @unittest.skipUnless(UPSTREAM.is_file(),
                         "sibling tttrlib checkout not present")
    def test_the_copy_has_not_drifted_from_tttrlib(self):
        ours, theirs = digest(VENDORED), digest(UPSTREAM)
        self.assertEqual(
            ours, theirs,
            "include/internal/ExpressionEngine.h has drifted from tttrlib's "
            f"copy.\n  here:     {VENDORED}\n  upstream: {UPSTREAM}\n"
            "tttrlib owns the original; re-copy it rather than editing here:\n"
            f"  cp {UPSTREAM} {VENDORED}")


if __name__ == "__main__":
    unittest.main()
