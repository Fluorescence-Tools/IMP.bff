"""``include/internal/MlpCore.h`` is a verbatim copy of tttrlib's, and stays one.

The differentiable MLP kernels (forward, backward, the Taylor-augmented passes
that give ``dy/dx`` and ``d2y/dx2`` and their adjoint) are written once, in
tttrlib ``modules/math/include/MlpCore.h`` — header-only and std-only for
exactly this reason. IMP.bff cannot link tttrlib (it is a soft, Python-level
dependency; see ``test_tttrlib_is_optional.py``), so it carries a copy the way
it carries pcg and nlohmann/json, and a network trained in tttrlib evaluates
and differentiates here bit for bit.

A copy diverges silently — someone fixes a derivative on one side, the other
keeps training to the wrong minimum. Hence a test, not a convention: when the
sibling checkout is present the two files must be identical. To refresh::

    cp ../tttrlib/modules/math/include/MlpCore.h include/internal/MlpCore.h

The direction is one-way; tttrlib is the source. Never edit the copy.
"""

import hashlib
import os
import unittest

import IMP
import IMP.test


def _sha(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


class Tests(IMP.test.TestCase):
    def _paths(self):
        here = os.path.dirname(os.path.abspath(__file__))
        repo = os.path.dirname(here)
        ours = os.path.join(repo, "include", "internal", "MlpCore.h")
        theirs = os.path.join(os.path.dirname(repo), "tttrlib", "modules", "math",
                              "include", "MlpCore.h")
        return ours, theirs

    def test_copy_exists_and_is_std_only(self):
        ours, _ = self._paths()
        self.assertTrue(os.path.exists(ours), ours)
        with open(ours) as fh:
            text = fh.read()
        # The contract that makes the copy possible: nothing but the standard
        # library, and no dependency on the rest of either repository.
        for forbidden in ("Mat.h", "nlohmann", "Registry.h", "SimPcgRandom", "Eigen", "IMP/"):
            self.assertNotIn('#include "' + forbidden, text)
            self.assertNotIn("#include <" + forbidden, text)
        self.assertIn("TTTRLIB_MLPCORE_H", text)

    def test_copy_matches_tttrlib_when_the_checkout_is_present(self):
        ours, theirs = self._paths()
        if not os.path.exists(theirs):
            self.skipTest("../tttrlib checkout not present; cannot compare")
        self.assertEqual(
            _sha(ours), _sha(theirs),
            "include/internal/MlpCore.h differs from ../tttrlib/modules/math/include/MlpCore.h; "
            "tttrlib is the source -- refresh with "
            "`cp ../tttrlib/modules/math/include/MlpCore.h include/internal/MlpCore.h`")


if __name__ == "__main__":
    IMP.test.main()
