#!/usr/bin/env python3

import os
import sys
import subprocess
import tempfile
import unittest
from pathlib import Path

import RMF

from IMP.bff.cgdye.utils import get_structure_dir

def _cgdye_file(*parts):
    """A path inside the installed IMP.bff.cgdye package.

    These tests used to walk up from __file__ into imp-tricks' src/IMP/bff
    layout, which stopped existing when cgdye moved into imp.bff. Asking the
    package where it is works wherever it is installed from.
    """
    import IMP.bff.cgdye
    from pathlib import Path as _P
    return _P(IMP.bff.cgdye.__file__).parent.joinpath(*parts)



class TestHGBP1RRTScript(unittest.TestCase):
    def _run_case(self, extra_args):
        script = _cgdye_file("scripts", "rrt_hgbp1_site481.py")
        if not os.path.exists(script):
            self.skipTest(f"{script} not found")

        required_inputs = [
            str(get_structure_dir("1DG3.pdb")),
            str(get_structure_dir("alexa488_r48.pdb")),
            str(get_structure_dir("alexa488_r48.mol2")),
        ]
        missing = [p for p in required_inputs if not os.path.exists(p)]
        if missing:
            self.skipTest(f"Missing required input files: {missing}")

        with tempfile.TemporaryDirectory() as tmpdir:
            out_rmf = os.path.join(tmpdir, "trajectory.rmf3")
            cmd = [
                sys.executable,
                script,
                "--n-iter",
                "4",
                "--output-rmf",
                out_rmf,
            ] + list(extra_args)
            env = os.environ.copy()
            # Prepend, do not replace: overwriting PYTHONPATH with "." threw away
            # everything the parent had, so the child could not import IMP.bff.
            env["PYTHONPATH"] = os.pathsep.join(
                [p for p in (".", env.get("PYTHONPATH", "")) if p])
            result = subprocess.run(
                cmd,
                env=env,
                capture_output=True,
                text=True,
                timeout=180,
            )

            self.assertEqual(result.returncode, 0, msg=result.stderr)
            self.assertIn("RRT(linker torsions only) finished", result.stdout)
            self.assertTrue(os.path.exists(out_rmf))

            fh = RMF.open_rmf_file_read_only(out_rmf)
            fh.close()

    def test_basic_rrt(self):
        self._run_case([])

    def test_rrt_with_seed(self):
        self._run_case(["--seed", "42"])

    def test_rrt_with_goal_bias(self):
        self._run_case(["--goal-bias", "0.5"])


if __name__ == "__main__":
    unittest.main()
