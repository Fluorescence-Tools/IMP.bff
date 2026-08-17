"""The hGBP1 site-481 Langevin/Brownian script runs end to end (md and bd)."""

import os
import subprocess
import sys
import tempfile
import unittest

import RMF


class TestLangevinHGBP1Script(unittest.TestCase):
    def _run(self, integrator):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_rmf = os.path.join(tmpdir, "traj.rmf3")
            cmd = [sys.executable, "-m", "IMP.bff.cgdye.scripts.langevin_hgbp1_site481",
                   "--integrator", integrator, "--n-steps", "1000", "--write-every", "100",
                   "--output-rmf", out_rmf]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            self.assertIn(f"{integrator} dynamics finished", result.stdout)
            fh = RMF.open_rmf_file_read_only(out_rmf)
            self.assertGreaterEqual(fh.get_number_of_frames(), 10)
            fh.close()

    def test_md(self):
        self._run("md")

    def test_bd(self):
        self._run("bd")


if __name__ == "__main__":
    unittest.main()
