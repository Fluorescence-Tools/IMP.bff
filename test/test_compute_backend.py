"""Where a kernel runs, and what the library says about it.

The accelerator is a plugin: looked for at import, used if it is there, and
otherwise not mentioned again beyond one warning. That makes two things worth
pinning. The **honesty** of the report -- `get_compute_backend_name()` must
say `cpu` when the kernels are on the CPU, because a silent fallback is how
somebody comes to believe they measured a GPU. And the **refusal** path: a
plugin that is missing, broken or built against another ABI must leave the
library working on the CPU rather than take it down with it.

None of this needs a GPU, and none of it is skipped where there is none.
"""

import os
import subprocess
import sys
import unittest

import IMP.bff


class TestTheReport(unittest.TestCase):

    def test_it_names_a_backend(self):
        name = IMP.bff.get_compute_backend_name()
        self.assertIsInstance(name, str)
        self.assertTrue(name)

    def test_without_a_plugin_it_is_the_cpu(self):
        # `cpu`, or a real backend if one was found and loaded -- but never
        # empty and never a lie about which one answered.
        name = IMP.bff.get_compute_backend_name()
        if IMP.bff.get_compute_backend() is None:
            self.assertEqual(name, "cpu")
        else:
            self.assertNotEqual(name, "cpu")


class TestRefusal(unittest.TestCase):

    def setUp(self):
        self.before = IMP.bff.get_compute_backend_name()

    def test_a_library_that_is_not_there_is_declined_with_a_reason(self):
        self.assertFalse(IMP.bff.load_compute_backend("/no/such/plugin.dylib"))
        self.assertIn("/no/such/plugin", IMP.bff.get_compute_backend_error())
        self.assertEqual(IMP.bff.get_compute_backend_name(), self.before)

    def test_an_empty_path_is_declined(self):
        self.assertFalse(IMP.bff.load_compute_backend(""))
        self.assertTrue(IMP.bff.get_compute_backend_error())
        self.assertEqual(IMP.bff.get_compute_backend_name(), self.before)

    def test_a_library_without_the_entry_point_is_declined(self):
        # any real shared library that is not a backend will do
        import ctypes.util
        path = ctypes.util.find_library("m") or ctypes.util.find_library("c")
        if not path or not os.path.isabs(path):
            self.skipTest("no system library to point at")
        self.assertFalse(IMP.bff.load_compute_backend(path))
        self.assertIn("imp_bff_compute_backend", IMP.bff.get_compute_backend_error())
        self.assertEqual(IMP.bff.get_compute_backend_name(), self.before)


class TestDiscovery(unittest.TestCase):

    def test_the_wgpu_library_is_a_file_or_nothing(self):
        path = IMP.bff._wgpu_library()
        if path is not None:
            self.assertTrue(os.path.isfile(path), path)
            self.assertIn("wgpu", os.path.basename(path))

    def test_the_plugin_is_a_file_or_nothing(self):
        path = IMP.bff._gpu_plugin()
        if path is not None:
            self.assertTrue(os.path.isfile(path), path)

    def test_off_is_honoured(self):
        out = subprocess.run(
            [sys.executable, "-c",
             "import IMP.bff; print(IMP.bff.get_compute_backend_name())"],
            env=dict(os.environ, IMP_BFF_GPU="off"), capture_output=True, text=True)
        self.assertEqual(out.returncode, 0, out.stderr[-2000:])
        self.assertEqual(out.stdout.strip().splitlines()[-1], "cpu")


class TestTheKernelStillWorks(unittest.TestCase):
    """Whatever answered, the numbers have to be the numbers."""

    def test_a_short_propagation_conserves_what_it_should(self):
        import numpy as np
        ng = 11
        i = np.arange(ng) - ng // 2
        X, Y, Z = np.meshgrid(i, i, i, indexing="ij")
        inside = (X ** 2 + Y ** 2 + Z ** 2) <= (ng // 2 - 1) ** 2
        cur = (inside / inside.sum()).ravel()
        d = (np.full(ng ** 3, 0.05) * inside.ravel())
        decay = np.ones(ng ** 3)          # nothing decays
        bounds = inside.ravel().astype(float)
        fl, dens = IMP.bff.diffusion_propagate(
            cur.tolist(), d.tolist(), decay.tolist(), bounds.tolist(),
            ng, 0, 200, 50)
        fl = np.asarray(fl)
        self.assertTrue(np.all(np.isfinite(fl)))
        # with no decay and a closed boundary the population is preserved
        np.testing.assert_allclose(fl, fl[0], rtol=1e-9)


if __name__ == "__main__":
    unittest.main()
