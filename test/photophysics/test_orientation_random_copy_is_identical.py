"""The vendored centralized RNG must be byte-identical to tttrlib's.

`include/internal/Random.h` is a copy of tttrlib's
`modules/math/include/Random.h`, and the two Monte-Carlo samplers in
`FRETOrientationFactor.cpp` (`wobbling_kappa2_distribution`,
`sample_kappa2_diffusion_with_traps`) draw through it instead of a local
`std::mt19937_64`, so their randomness follows the rest of the app's
`TTTR_RNG_SEED` / `TTTR_RNG_ENGINE` discipline rather than a third,
independent generator.

Copied rather than depended on, for the same reason the decay kernels are
(see `test/decay/test_decay_convolution_copy_is_identical.py`): the pieces
used here are header-only, so there is nothing to link, and taking them from
an installed tttrlib would make this repository need tttrlib present to
build.

`include/internal/info.h` is NOT pinned the same way: it is a deliberately
reduced subset of tttrlib's `modules/util/include/info.h` (just the three
environment-variable helpers `Random.h` needs), not a full copy -- the rest
of that file drives AVX/NEON dispatch through macros tttrlib's own CMake
defines, which do not exist in this build. Its three functions are still
transcribed verbatim; there is no automated pin for that, because there is
no single upstream file to hash against.

**tttrlib owns the original.** A change goes there and is copied here:

    cp ../tttrlib/modules/math/include/Random.h \\
       include/internal/Random.h

Skipped, not failed, when the sibling checkout is absent -- a release
tarball or a CI job that builds only imp.bff cannot run this, and should not
be reported as broken for it.
"""

import hashlib
import pathlib
import unittest

HERE = pathlib.Path(__file__).resolve()
VENDORED = HERE.parents[2] / "include" / "internal" / "Random.h"
UPSTREAM = (HERE.parents[3] / "tttrlib" / "modules" / "math" / "include"
            / "Random.h")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class OrientationRandomCopyTests(unittest.TestCase):

    def test_the_vendored_copy_exists(self):
        self.assertTrue(VENDORED.is_file(), f"missing: {VENDORED}")

    def test_it_carries_the_centralized_rng(self):
        text = VENDORED.read_text()
        self.assertIn("class Random", text)
        self.assertIn("double normal()", text)

    @unittest.skipUnless(UPSTREAM.is_file(),
                         "sibling tttrlib checkout not present")
    def test_the_copy_has_not_drifted_from_tttrlib(self):
        ours, theirs = digest(VENDORED), digest(UPSTREAM)
        self.assertEqual(
            ours, theirs,
            "include/internal/Random.h has drifted from tttrlib's copy.\n"
            f"  here:     {VENDORED}\n  upstream: {UPSTREAM}\n"
            "tttrlib owns the original; re-copy it rather than editing here:\n"
            f"  cp {UPSTREAM} {VENDORED}")


if __name__ == "__main__":
    unittest.main()
