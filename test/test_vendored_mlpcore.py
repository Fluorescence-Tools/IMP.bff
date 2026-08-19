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

    def test_bff_evaluates_and_differentiates_a_tttrlib_trained_model(self):
        """A network trained by tttrlib runs, and differentiates, here.

        Trains a small net with scalers in tttrlib, writes the JSON, compiles
        ``test/cpp_snippets/mlpcore_eval.cpp`` against bff's vendored headers
        only (``internal/json.h`` + ``internal/MlpCore.h``, under the
        ``IMP::bff::internal`` namespace) and checks that the predictions
        agree to 1e-12 and that the C++ side's ``dL/dparams`` and ``dL/dx``
        match central differences. This is the contract PRD-115 builds on.
        """
        import shutil
        import subprocess
        import tempfile

        import numpy as np
        try:
            import tttrlib
        except ImportError:
            self.skipTest("tttrlib not installed")
        cxx = shutil.which("c++") or shutil.which("clang++") or shutil.which("g++")
        if cxx is None:
            self.skipTest("no C++ compiler on PATH")
        here = os.path.dirname(os.path.abspath(__file__))
        repo = os.path.dirname(here)
        src = os.path.join(here, "cpp_snippets", "mlpcore_eval.cpp")

        rng = np.random.default_rng(3)
        X = rng.uniform(-1, 1, size=(300, 2))
        Y = np.column_stack([np.sin(2 * X[:, 0]) + X[:, 1] ** 2, X[:, 0] * X[:, 1]])
        opt = tttrlib.TrainOptions()
        opt.hidden_layer_sizes = tttrlib.VectorInt32([8, 8])
        opt.activation = tttrlib.activation_from_string("tanh")
        opt.max_iter = 30
        opt.seed = 1
        net = tttrlib.NeuralNet.train_np(X, Y, opt)
        Xt = rng.uniform(-1, 1, size=(5, 2))
        expect = net.predict_batch_np(Xt)

        with tempfile.TemporaryDirectory() as tmp:
            # the header must compile with IMP/bff/internal/... include paths
            inc = os.path.join(tmp, "IMP", "bff")
            os.makedirs(inc)
            os.symlink(os.path.join(repo, "include", "internal"), os.path.join(inc, "internal"))
            exe = os.path.join(tmp, "mlpcore_eval")
            subprocess.check_call([cxx, "-std=c++17", "-O2", "-I", tmp, src, "-o", exe])
            model_json = os.path.join(tmp, "model.json")
            net.to_json_file(model_json)
            xfile = os.path.join(tmp, "X.txt")
            with open(xfile, "w") as fh:
                fh.write("%d %d\n" % Xt.shape)
                for row in Xt:
                    fh.write(" ".join("%.17g" % v for v in row) + "\n")
            out = subprocess.check_output([exe, model_json, xfile], text=True).strip().splitlines()
        got = np.array([[float(v) for v in line.split()] for line in out[: Xt.shape[0]]])
        self.assertLess(np.abs(got - expect).max(), 1e-12)
        checks = dict(line.split() for line in out[Xt.shape[0]:])
        self.assertLess(float(checks["fd_check"]), 1e-6)
        self.assertLess(float(checks["dx_check"]), 1e-6)


if __name__ == "__main__":
    IMP.test.main()
