"""A dye on a structure, as an ordinary simulation object.

Placing a dye on a residue and integrating its motion is built on IMP --
its hierarchies, its decorators, its integrators. `MolecularProbeSimulation` is that
computation behind arrays and `IMP::bff::ProbeSimulation`, which is what
lets it be wrapped without IMP's own SWIG interfaces and reached from a
build that carries IMP as a private library and has no IMP in Python at all
(PRD-139).

The tests assert three things. That it *works*: a dye is placed, it moves,
and a trajectory comes back with the shapes it promises. That it keeps the
*shared shape*: parameters are JSON, an unknown one is refused, and the
object is a `ProbeSimulation` like the grid walk. And that it is reachable
*the way it is meant to be*: where `IMP.bff.get_build()` says `core+imp`,
running it must not import any IMP module beyond `IMP.bff` itself -- if it
ever does, the package has quietly grown a dependency on IMP's Python.

Skipped where the layer is not built: the IMP-free core has no dye
dynamics, and says so by not having the name.
"""

import json
import subprocess
import sys
import unittest

import numpy as np

import IMP.bff

_HAS_DOORS = hasattr(IMP.bff, "MolecularProbeSimulation")


def _inputs():
    """The structure and the dye as arrays -- a frame, not a file path."""
    protein = IMP.bff.load_protein_frames(
        str(IMP.bff.get_structure_dir("1DG3.pdb")))[0]
    dye = IMP.bff.load_protein_frames(
        str(IMP.bff.get_structure_dir("alexa488_r48.pdb")))[0]
    return protein, dye, str(IMP.bff.get_structure_dir("alexa488_r48.mol2"))


def _simulation(**parameters):
    protein, dye, mol2 = _inputs()
    parameters.setdefault("seed", 1)
    return IMP.bff.MolecularProbeSimulation(protein, dye, mol2, "A", 481,
                                 json.dumps(parameters))


@unittest.skipUnless(_HAS_DOORS, "this build has no connection layer")
class TestSetup(unittest.TestCase):

    def test_the_dye_is_placed_and_the_side_chain_goes(self):
        sim = _simulation()
        self.assertGreater(sim.n_stripped, 0, "the side chain was left in place")
        self.assertGreater(sim.n_atoms, 10)
        self.assertEqual(len(sim.atom_names), sim.n_atoms)

    def test_the_labelled_structure_holds_the_dye(self):
        """The frame that comes back is the protein *with* the label on it."""
        sim = _simulation()
        protein, dye, _ = _inputs()
        labelled = sim.get_labelled_frame()
        # every atom still there, less the stripped side chain, plus the dye
        self.assertEqual(labelled.n_atoms,
                         protein.n_atoms - sim.n_stripped + dye.n_atoms)

    def test_a_site_that_is_not_there_is_refused(self):
        protein, dye, mol2 = _inputs()
        with self.assertRaises(Exception):
            IMP.bff.MolecularProbeSimulation(protein, dye, mol2, "Z", 99999)


@unittest.skipUnless(_HAS_DOORS, "this build has no connection layer")
class TestTheSharedShape(unittest.TestCase):
    """What every simulation in the module answers to."""

    def test_it_is_a_probe_simulation(self):
        sim = _simulation()
        self.assertIsInstance(sim, IMP.bff.ProbeSimulation)
        self.assertEqual(sim.get_type(), "dye-langevin")
        self.assertTrue(sim.has_energy())

    def test_parameters_are_json_and_report_what_is_in_force(self):
        sim = _simulation(temperature=310.0, integrator="bd")
        p = json.loads(sim.get_parameters())
        self.assertEqual(p["temperature"], 310.0)
        self.assertEqual(p["integrator"], "bd")

    def test_a_parameter_it_does_not_have_is_refused(self):
        protein, dye, mol2 = _inputs()
        with self.assertRaises(Exception):
            IMP.bff.MolecularProbeSimulation(protein, dye, mol2, "A", 481,
                                  json.dumps({"viscosity": 1.0}))

    def test_an_unknown_integrator_is_refused(self):
        protein, dye, mol2 = _inputs()
        with self.assertRaises(Exception):
            IMP.bff.MolecularProbeSimulation(protein, dye, mol2, "A", 481,
                                  json.dumps({"integrator": "verlet-ish"}))

    def test_parameters_are_fixed_once_it_is_built(self):
        sim = _simulation()
        with self.assertRaises(Exception):
            sim.set_parameters(json.dumps({"temperature": 400.0}))


@unittest.skipUnless(_HAS_DOORS, "this build has no connection layer")
class TestRunning(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.sim = _simulation()
        cls.sim.minimize(50)
        cls.traj = cls.sim.run(300, write_every=100)

    def test_minimize_lowers_the_energy(self):
        sim = _simulation()
        before = sim.get_potential_energy()
        after = sim.minimize(50)
        self.assertLessEqual(after, before)

    def test_step_moves_the_dye(self):
        sim = _simulation()
        before = np.asarray(sim.positions).copy()
        sim.step(100)
        after = np.asarray(sim.positions)
        self.assertEqual(before.shape, (sim.n_atoms, 3))
        self.assertGreater(float(np.abs(after - before).max()), 0.0)

    def test_positions_can_be_put_back(self):
        sim = _simulation()
        start = np.asarray(sim.positions).copy()
        sim.step(50)
        sim.set_positions(start.ravel().tolist())
        np.testing.assert_allclose(np.asarray(sim.positions), start, atol=1e-9)

    def test_a_wrong_number_of_coordinates_is_refused(self):
        sim = _simulation()
        with self.assertRaises(Exception):
            sim.set_positions([0.0, 0.0, 0.0])

    def test_the_trajectory_has_the_shape_it_promises(self):
        t = self.traj
        self.assertGreater(t.n_frames, 0)
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


@unittest.skipUnless(_HAS_DOORS and IMP.bff.get_build() == "core+imp",
                     "only the IMP-linked wheel build makes this promise")
class TestNoImpInPython(unittest.TestCase):
    """The whole point: IMP is a linked library, not an import."""

    def test_running_a_dye_imports_no_imp_module(self):
        script = (
            "import sys, IMP.bff\n"
            "p = IMP.bff.load_protein_frames("
            "str(IMP.bff.get_structure_dir('1DG3.pdb')))[0]\n"
            "d = IMP.bff.load_protein_frames("
            "str(IMP.bff.get_structure_dir('alexa488_r48.pdb')))[0]\n"
            "m = str(IMP.bff.get_structure_dir('alexa488_r48.mol2'))\n"
            "IMP.bff.MolecularProbeSimulation(p, d, m, 'A', 481).run(20, 10)\n"
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
