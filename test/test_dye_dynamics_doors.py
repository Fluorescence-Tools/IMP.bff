"""The dye roads by file path: `attach_dye_to_pdb` and `run_dye_langevin`.

Placing a dye on a residue and integrating its motion is built on IMP --
its hierarchies, its decorators, its integrators. These two functions are
that computation behind file paths and arrays, which is what lets them be
wrapped without IMP's own SWIG interfaces and reached from a build that
carries IMP as a private library and has no IMP in Python at all (PRD-139).

The tests therefore assert two different things. That the doors *work*:
a dye is placed, a trajectory comes back with the shapes it promises, and
the dye actually moves. And that they are reachable *the way they are meant
to be*: where `IMP.bff.get_build()` says `core+imp`, running one of them
must not import any IMP module beyond `IMP.bff` itself -- if it ever does,
the wheel has quietly grown a dependency on IMP's Python.

Skipped where the layer is not built: the IMP-free core has no dye
dynamics, and says so by not having the names.
"""

import os
import subprocess
import sys
import unittest

import numpy as np

import IMP.bff

_HAS_DOORS = hasattr(IMP.bff, "run_dye_langevin")


def _inputs():
    return (str(IMP.bff.get_structure_dir("1DG3.pdb")),
            str(IMP.bff.get_structure_dir("alexa488_r48.pdb")),
            str(IMP.bff.get_structure_dir("alexa488_r48.mol2")))


@unittest.skipUnless(_HAS_DOORS, "this build has no connection layer")
class TestAttach(unittest.TestCase):

    def test_a_dye_is_placed_and_the_side_chain_goes(self):
        protein, dye, _ = _inputs()
        stripped = IMP.bff.attach_dye_to_pdb(protein, dye, "A", 481)
        self.assertGreater(stripped, 0, "the site's side chain was left in place")

    def test_it_writes_the_labelled_structure(self):
        import tempfile
        protein, dye, _ = _inputs()
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "labelled.pdb")
            IMP.bff.attach_dye_to_pdb(protein, dye, "A", 481, out)
            self.assertTrue(os.path.isfile(out))
            with open(out) as fh:
                atoms = [l for l in fh if l.startswith("ATOM") or l.startswith("HETATM")]
            self.assertGreater(len(atoms), 100)

    def test_a_site_that_is_not_there_is_refused(self):
        protein, dye, _ = _inputs()
        with self.assertRaises(Exception):
            IMP.bff.attach_dye_to_pdb(protein, dye, "Z", 99999)


@unittest.skipUnless(_HAS_DOORS, "this build has no connection layer")
class TestLangevin(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        protein, dye, mol2 = _inputs()
        cls.traj = IMP.bff.run_dye_langevin(
            protein, dye, mol2, "A", 481, n_steps=500, write_every=100,
            minimize_steps=50, seed=1)

    def test_the_trajectory_has_the_shape_it_promises(self):
        t = self.traj
        self.assertGreater(t.n_frames, 0)
        self.assertGreater(t.n_atoms, 0)
        c = np.asarray(t.coordinates)
        self.assertEqual(c.shape, (t.n_frames, t.n_atoms, 3))
        self.assertTrue(np.all(np.isfinite(c)))
        for name in ("times_fs", "potential_energy", "kinetic_energy"):
            self.assertEqual(len(np.asarray(getattr(t, name))), t.n_frames, name)

    def test_the_dye_moves_and_stays_finite(self):
        c = np.asarray(self.traj.coordinates)
        spread = np.sqrt(((c - c.mean(axis=0)) ** 2).sum(-1).mean())
        self.assertGreater(spread, 0.0, "the dye did not move at all")
        self.assertLess(spread, 50.0, "the dye flew apart")

    def test_the_energies_are_finite(self):
        for name in ("potential_energy", "kinetic_energy"):
            v = np.asarray(getattr(self.traj, name))
            self.assertTrue(np.all(np.isfinite(v)), name)

    def test_an_unknown_integrator_is_refused(self):
        protein, dye, mol2 = _inputs()
        with self.assertRaises(Exception):
            IMP.bff.run_dye_langevin(protein, dye, mol2, "A", 481, n_steps=10,
                                     integrator="verlet-ish")


@unittest.skipUnless(_HAS_DOORS and IMP.bff.get_build() == "core+imp",
                     "only the IMP-linked wheel build makes this promise")
class TestNoImpInPython(unittest.TestCase):
    """The whole point: IMP is a linked library, not an import."""

    def test_running_a_dye_imports_no_imp_module(self):
        script = (
            "import sys, IMP.bff\n"
            "p = str(IMP.bff.get_structure_dir('1DG3.pdb'))\n"
            "d = str(IMP.bff.get_structure_dir('alexa488_r48.pdb'))\n"
            "m = str(IMP.bff.get_structure_dir('alexa488_r48.mol2'))\n"
            "IMP.bff.run_dye_langevin(p, d, m, 'A', 481, n_steps=20,"
            " write_every=10, minimize_steps=0, seed=1)\n"
            # marked, because the C++ side writes progress to stdout too
            "print('MODULES=' + repr(sorted(n for n in sys.modules"
            " if n.split('.')[0] == 'IMP')))\n")
        out = subprocess.run([sys.executable, "-c", script], capture_output=True,
                             text=True)
        self.assertEqual(out.returncode, 0, out.stderr[-2000:])
        marked = [l for l in out.stdout.splitlines() if l.startswith("MODULES=")]
        self.assertTrue(marked, out.stdout[-2000:])
        loaded = eval(marked[-1][len("MODULES="):])
        self.assertTrue(all(n.startswith("IMP.bff") or n == "IMP" for n in loaded),
                        "an IMP module was imported: %r" % (loaded,))


if __name__ == "__main__":
    unittest.main()
