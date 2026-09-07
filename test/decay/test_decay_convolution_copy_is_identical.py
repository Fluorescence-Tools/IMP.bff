"""The vendored decay kernels must be byte-identical to tttrlib's.

`include/internal/DecayConvolution.h` is a copy of tttrlib's
`modules/spectroscopy/decay/include/DecayConvolution.h`, and `TcspcDecay`
evaluates a fit's model curve with two of the header-only pieces in it:
``fconv_per_cs_ad<double>`` for the periodic reconvolution and
``shift_lamp_ad<double>`` for the timeshift.

Copied rather than depended on, for the same reason the expression engine is
(see `test/expression/test_engine_copy_is_identical.py`): the pieces used
here are header-only, so there is nothing to link, and taking them from an
installed tttrlib would make this repository need tttrlib present to build.

The hazard of a copy is that it drifts, and a drifted copy is worse than
either a dependency or a fork -- it looks like one implementation while
behaving as two, and here the two would be two different opinions about what
a lifetime is. So it is checked rather than trusted.

**tttrlib owns the original.** A change goes there and is copied here:

    cp ../tttrlib/modules/spectroscopy/decay/include/DecayConvolution.h \\
       include/internal/DecayConvolution.h

Skipped, not failed, when the sibling checkout is absent -- a release tarball
or a CI job that builds only imp.bff cannot run this, and should not be
reported as broken for it.
"""

import hashlib
import pathlib
import unittest

HERE = pathlib.Path(__file__).resolve()
VENDORED = HERE.parents[2] / "include" / "internal" / "DecayConvolution.h"
UPSTREAM = (HERE.parents[3] / "tttrlib" / "modules" / "spectroscopy" / "decay"
            / "include" / "DecayConvolution.h")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class DecayConvolutionCopyTests(unittest.TestCase):

    def test_the_vendored_copy_exists(self):
        self.assertTrue(VENDORED.is_file(), f"missing: {VENDORED}")

    def test_it_carries_the_two_header_only_kernels(self):
        """The copy is only worth having for what can be used from a header.

        Both are templated cores that tttrlib's own ``fconv_per_cs`` and
        ``shift_lamp`` call, so there is one implementation and not two.
        """
        text = VENDORED.read_text()
        self.assertIn("fconv_per_cs_ad", text)
        self.assertIn("shift_lamp_ad", text)

    @unittest.skipUnless(UPSTREAM.is_file(),
                         "sibling tttrlib checkout not present")
    def test_the_copy_has_not_drifted_from_tttrlib(self):
        ours, theirs = digest(VENDORED), digest(UPSTREAM)
        self.assertEqual(
            ours, theirs,
            "include/internal/DecayConvolution.h has drifted from tttrlib's "
            f"copy.\n  here:     {VENDORED}\n  upstream: {UPSTREAM}\n"
            "tttrlib owns the original; re-copy it rather than editing here:\n"
            f"  cp {UPSTREAM} {VENDORED}")


if __name__ == "__main__":
    unittest.main()
