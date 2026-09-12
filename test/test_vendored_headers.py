"""The headers vendored from tttrlib are verbatim copies, and stay so.

Two kernels are written once, in tttrlib, header-only and std-only for
exactly this reason, and carried here the way pcg and nlohmann/json are:

* ``include/internal/MlpCore.h`` — the differentiable MLP (forward, backward,
  the Taylor-augmented passes for ``dy/dx`` and ``d²y/dx²`` and their adjoint,
  the whole-model struct with scalers and the JSON format);
* ``include/internal/LatticeDiffusion.h`` — the masked-lattice diffusion
  solver behind ``GridDiffusionSolver`` and its adjoint.

IMP.bff cannot link tttrlib (it is a soft, Python-level dependency; see
``test_tttrlib_is_optional.py``), so it carries copies, and a network trained
in tttrlib evaluates and differentiates here bit for bit; a gradient of the
field solver computed here is the one tttrlib's test validated.

A copy diverges silently — someone fixes a derivative on one side, the other
keeps training to the wrong minimum. Hence a test, not a convention: when the
sibling checkout is present each pair must be identical. To refresh::

    cp ../tttrlib/modules/math/include/MlpCore.h include/internal/MlpCore.h
    cp ../tttrlib/modules/math/include/LatticeDiffusion.h include/internal/LatticeDiffusion.h

The direction is one-way; tttrlib is the source. Never edit the copies.
"""

import hashlib
import json
import os
import unittest

import IMP
import IMP.test


def _sha(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


class Tests(IMP.test.TestCase):
    # name -> the tttrlib module it was taken from. A map rather than one
    # hardcoded directory because the copies do not all come from `math`:
    # DecayConvolution.h is the spectroscopy module's, and while it said in
    # its own header that it was "kept in step", nothing checked that it was.
    VENDORED = {
        "Dual.h": ("math",),
        "MlpCore.h": ("math",),
        "LatticeDiffusion.h": ("math",),
        "DecayConvolution.h": ("spectroscopy", "decay"),
    }

    def _paths(self, name="MlpCore.h"):
        here = os.path.dirname(os.path.abspath(__file__))
        repo = os.path.dirname(here)
        ours = os.path.join(repo, "include", "internal", name)
        theirs = os.path.join(os.path.dirname(repo), "tttrlib", "modules",
                              *self.VENDORED[name], "include", name)
        return ours, theirs

    def test_copies_exist_and_are_std_only(self):
        for name in self.VENDORED:
            ours, _ = self._paths(name)
            self.assertTrue(os.path.exists(ours), ours)
            with open(ours) as fh:
                text = fh.read()
            # The contract that makes the copy possible: nothing but the standard
            # library, and no dependency on the rest of either repository.
            for forbidden in ("Mat.h", "nlohmann", "Registry.h", "SimPcgRandom", "Eigen", "IMP/"):
                self.assertNotIn('#include "' + forbidden, text, name)
                self.assertNotIn("#include <" + forbidden, text, name)
            guard = ("TTTRLIB_FSCONV_H" if name == "DecayConvolution.h"
                     else "TTTRLIB_" + name[:-2].upper() + "_H")
            self.assertIn(guard, text)

    def test_copies_match_tttrlib_when_the_checkout_is_present(self):
        for name in self.VENDORED:
            ours, theirs = self._paths(name)
            if not os.path.exists(theirs):
                self.skipTest("../tttrlib checkout not present; cannot compare")
            self.assertEqual(
                _sha(ours), _sha(theirs),
                f"include/internal/{name} differs from ../tttrlib/modules/math/include/{name}; "
                "tttrlib is the source -- refresh with "
                f"`cp ../tttrlib/modules/math/include/{name} include/internal/{name}`")

    def _ptolib_paths(self):
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        return (os.path.join(repo, "thirdparty", "ptolib"),
                os.environ.get("PTOLIB_CHECKOUT", os.path.join(os.path.dirname(repo), "ptolib")))

    def test_ptolib_manifest_matches_vendored_sources(self):
        """Validate every vendored source even without a sibling checkout."""
        vendor = self._ptolib_paths()[0]
        with open(os.path.join(vendor, "VENDORING.json")) as fh:
            manifest = json.load(fh)
        self.assertIn("CMakeLists.txt", manifest)
        self.assertIn("include/ptolib/ptolib.h", manifest)
        self.assertIn("src/ptolib.cpp", manifest)
        for rel, expected in manifest.items():
            self.assertFalse(os.path.isabs(rel), rel)
            self.assertNotIn("..", rel.split("/"), rel)
            path = os.path.join(vendor, rel)
            self.assertTrue(os.path.isfile(path), rel)
            self.assertFalse(os.path.islink(path), rel)
            self.assertEqual(_sha(path), expected, rel + ": refresh ptolib")

    def test_ptolib_sources_match_checkout_when_present(self):
        vendor, checkout = self._ptolib_paths()
        if not os.path.isdir(os.path.join(checkout, "include", "ptolib")):
            self.skipTest("ptolib checkout not present; cannot compare")
        with open(os.path.join(vendor, "VENDORING.json")) as fh:
            manifest = json.load(fh)
        for rel in manifest:
            self.assertEqual(_sha(os.path.join(vendor, rel)),
                             _sha(os.path.join(checkout, rel)),
                             rel + ": refresh the ptolib source package")

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

    def test_bff_loads_an_onnx_model_without_tttrlib(self):
        """The same bff-only program reads an ONNX file written by PyTorch --
        tttrlib's committed fixture, with PyTorch's own outputs -- so a network
        trained anywhere runs inside bff with nothing but the vendored header."""
        import shutil
        import subprocess
        import tempfile

        import numpy as np
        cxx = shutil.which("c++") or shutil.which("clang++") or shutil.which("g++")
        if cxx is None:
            self.skipTest("no C++ compiler on PATH")
        here = os.path.dirname(os.path.abspath(__file__))
        repo = os.path.dirname(here)
        fixtures = os.path.join(os.path.dirname(repo), "tttrlib", "test", "python", "misc", "fixtures", "nn")
        if not os.path.exists(os.path.join(fixtures, "expected.json")):
            self.skipTest("../tttrlib fixtures not present")
        import json
        with open(os.path.join(fixtures, "expected.json")) as fh:
            e = json.load(fh)
        X = np.asarray(e["X"])
        src = os.path.join(here, "cpp_snippets", "mlpcore_eval.cpp")
        with tempfile.TemporaryDirectory() as tmp:
            inc = os.path.join(tmp, "IMP", "bff")
            os.makedirs(inc)
            os.symlink(os.path.join(repo, "include", "internal"), os.path.join(inc, "internal"))
            exe = os.path.join(tmp, "mlpcore_eval")
            subprocess.check_call([cxx, "-std=c++17", "-O2", "-I", tmp, src, "-o", exe])
            xfile = os.path.join(tmp, "X.txt")
            with open(xfile, "w") as fh:
                fh.write("%d %d\n" % X.shape)
                for row in X:
                    fh.write(" ".join("%.17g" % v for v in row) + "\n")
            for name, key, tol in (("mlp_torch_legacy.onnx", "Y_float32", 1e-6),
                                   ("mlp_matmul_add.onnx", "Y_double", 1e-12)):
                out = subprocess.check_output([exe, os.path.join(fixtures, name), xfile], text=True).strip().splitlines()
                got = np.array([[float(v) for v in line.split()] for line in out[: X.shape[0]]])
                self.assertLess(np.abs(got - np.asarray(e[key])).max(), tol, name)


if __name__ == "__main__":
    IMP.test.main()
