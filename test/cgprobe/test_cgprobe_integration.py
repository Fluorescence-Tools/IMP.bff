import unittest
import subprocess
import os
import sys
import tempfile
import shutil
from pathlib import Path

from IMP.bff import get_output_dir, get_structure_dir

class TestIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Build the CX4+atto655 force-field system into a temporary directory
        # instead of expecting a pre-built output/systems/*.system.cif: the
        # test used to skip everywhere the build-system task had not run.
        from IMP.bff import write_probe_forcefield_cif
        from IMP.bff import build_probe_protein_system
        from IMP.bff import get_template_dir
        cls._tmp = tempfile.TemporaryDirectory()
        system = build_probe_protein_system(
            str(get_structure_dir("cx4.mol2")),
            str(get_structure_dir("atto655.mol2")),
            "CX4",
            "atto655",
            protein_template=str(get_template_dir("cx4.template.cif")),
            probe_template=str(get_template_dir("atto655.template.cif")),
        )
        cls.system_cif = os.path.join(cls._tmp.name, "cx4_atto655.system.cif")
        write_probe_forcefield_cif(cls.system_cif, system)

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def run_sim(self, args):
        env = os.environ.copy()
        # Prepend, do not replace: overwriting PYTHONPATH with "." threw away
        # everything the parent had, so the child could not import IMP.bff.
        env["PYTHONPATH"] = os.pathsep.join(
            [p for p in (".", env.get("PYTHONPATH", "")) if p])
        # The command is `imp_bff simulate` now. `IMP.bff.cgprobe.sim` kept the
        # function and lost the click wrapper, so `-m` on it exits 0 having
        # done nothing -- which is how this test failed: a returncode of 0 and
        # no output directory.
        program = Path(__file__).resolve().parents[2] / "bin" / "imp_bff"
        cmd = [sys.executable, str(program), "simulate"] + args
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
