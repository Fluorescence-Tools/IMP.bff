import unittest
import subprocess
import os
import sys
import tempfile
import shutil
from pathlib import Path

from IMP.bff.cgdye.utils import get_output_dir, get_structure_dir

def _cgdye_file(*parts):
    """A path inside the installed IMP.bff.cgdye package.

    These tests used to walk up from __file__ into imp-tricks' src/IMP/bff
    layout, which stopped existing when cgdye moved into imp.bff. Asking the
    package where it is works wherever it is installed from.
    """
    import IMP.bff.cgdye
    from pathlib import Path as _P
    return _P(IMP.bff.cgdye.__file__).parent.joinpath(*parts)


class TestIntegration(unittest.TestCase):
    def setUp(self):
        self.system_cif = str(get_output_dir("systems", "cx4_atto655.system.cif"))
        self.nmr_cif = str(Path(__import__("IMP.bff", fromlist=["x"]).get_data_path("cgdye")) / "inputs" / "restraints" / "nmr.restraints.cif")
        self.script = _cgdye_file("sim", "runner.py")
        
        # Ensure dependencies exist
        if not os.path.exists(self.system_cif):
            self.skipTest(f"{self.system_cif} not found. Run pixi run build-system first.")

    def run_sim(self, args):
        env = os.environ.copy()
        # Prepend, do not replace: overwriting PYTHONPATH with "." threw away
        # everything the parent had, so the child could not import IMP.bff.
        env["PYTHONPATH"] = os.pathsep.join(
            [p for p in (".", env.get("PYTHONPATH", "")) if p])
        cmd = [sys.executable, self.script] + args
        # Use a reasonable timeout for integration tests
        return subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=60)

    def test_hybrid_md_mc_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = [
                "--system-cif", self.system_cif,
                "--output-root", tmpdir,
                "--sampling-mode", "hybrid_md_mc",
                "--md-steps", "10",
                "--write-every", "5",
                "--nmr-cif", self.nmr_cif,
                "--com-pull-k", "1.0"
            ]
            result = self.run_sim(args)
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            
            run_dir = os.path.join(tmpdir, "CX4_atto655_imp")
            self.assertTrue(os.path.exists(run_dir))
            # Hybrid mode writes to rmfs/0.rmf3
            self.assertTrue(os.path.exists(os.path.join(run_dir, "rmfs", "0.rmf3")))

    def test_multi_restart_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = [
                "--system-cif", self.system_cif,
                "--output-root", tmpdir,
                "--sampling-mode", "multi_restart",
                "--n-restarts", "2",
                "--md-steps", "10",
                "--write-every", "5",
                "--nmr-cif", self.nmr_cif,
                "--convergence-threshold", "999999.0" 
            ]
            result = self.run_sim(args)
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            
            run_dir = os.path.join(tmpdir, "CX4_atto655_imp")
            self.assertTrue(os.path.exists(run_dir))
            self.assertTrue(os.path.exists(os.path.join(run_dir, "rmfs", "0.rmf3")))
            self.assertTrue(os.path.exists(os.path.join(run_dir, "stat.0.out")))

    def test_simple_md_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = [
                "--system-cif", self.system_cif,
                "--output-root", tmpdir,
                "--sampling-mode", "simple_md",
                "--md-steps", "10",
                "--write-every", "5"
            ]
            result = self.run_sim(args)
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            
            run_dir = os.path.join(tmpdir, "CX4_atto655_imp")
            self.assertTrue(os.path.exists(run_dir))
            self.assertTrue(os.path.exists(os.path.join(run_dir, "rmfs", "0.rmf3")))

if __name__ == "__main__":
    unittest.main()
