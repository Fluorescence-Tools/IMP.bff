#!/usr/bin/env python3

import os
import sys
import subprocess
import tempfile
import unittest
from pathlib import Path

import IMP
import IMP.algebra
import IMP.atom
import IMP.core
import IMP.rmf
import RMF

def _cgdye_file(*parts):
    """A path inside the installed IMP.bff.cgdye package.

    These tests used to walk up from __file__ into imp-tricks' src/IMP/bff
    layout, which stopped existing when cgdye moved into imp.bff. Asking the
    package where it is works wherever it is installed from.
    """
    import IMP.bff.cgdye
    from pathlib import Path as _P
    return _P(IMP.bff.cgdye.__file__).parent.joinpath(*parts)



class TestLangevinHGBP1Script(unittest.TestCase):
    def test_langevin_run_writes_dynamic_rmf(self):
        script = _cgdye_file("scripts", "langevin_hgbp1_site481.py")
        if not script.exists():
            self.skipTest(f"{script} not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            out_rmf = os.path.join(tmpdir, "traj.rmf3")
            cmd = [
                sys.executable,
                script,
                "--n-steps",
                "100",
                "--write-every",
                "20",
                "--output-rmf",
                out_rmf,
                "--step-size",
                "0.5",  # large step to ensure displacement
            ]
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
                timeout=300,
            )
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            self.assertIn("Langevin (Internal-DOF stochastic walk) finished:", result.stdout)

            fh = RMF.open_rmf_file_read_only(out_rmf)
            # 100/20 = 5 frames + init = 6 frames
            self.assertGreaterEqual(fh.get_number_of_frames(), 2)

            m = IMP.Model()
            hs = IMP.rmf.create_hierarchies(fh, m)
            IMP.rmf.load_frame(fh, RMF.FrameID(0))

            # Find dye atoms using Selection
            sel = IMP.atom.Selection(hs[0], residue_index=2)
            dye_atoms = sel.get_selected_particles()
            # Sort dye atoms by index to be sure they match the Mol2 order
            dye_atoms = sorted(dye_atoms, key=lambda p: p.get_index())

            # 1-based index in dump means 0-based index 0 and 11
            p1 = dye_atoms[0]
            p10 = dye_atoms[11]
            targets = {"C1": p1, "C10": p10}

            def get_dye_atom_coords(frame_idx):
                IMP.rmf.load_frame(fh, RMF.FrameID(frame_idx))
                return {name: IMP.core.XYZ(p).get_coordinates() for name, p in targets.items()}

            n_frames = fh.get_number_of_frames()

            distances = []
            for f in range(n_frames):
                c = get_dye_atom_coords(f)
                dist = IMP.algebra.get_distance(c["C1"], c["C10"])
                distances.append(dist)

            d0 = distances[0]
            d1 = distances[-1]

            # Internal distance between core atoms MUST be preserved exactly in internal-DOF sampling.
            self.assertAlmostEqual(d0, d1, places=4)

            # Check for overall displacement of the dye core center
            def get_center(c):
                return [(c["C1"][i] + c["C10"][i]) / 2 for i in range(3)]

            c0 = get_dye_atom_coords(0)
            c1 = get_dye_atom_coords(n_frames - 1)

            ctr0 = get_center(c0)
            ctr1 = get_center(c1)

            disp = IMP.algebra.get_distance(ctr0, ctr1)
            # Since step-size is 0.5 rad and site 481 is exposed, it should move.
            self.assertGreater(disp, 1e-3)



if __name__ == "__main__":
    unittest.main()
