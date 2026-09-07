from __future__ import division
import unittest
import os

import tempfile

import math
import json
import numpy as np
import numpy.testing

import IMP
import IMP.core
import IMP.atom
import IMP.em
import IMP.bff


create_references = False

mdl = IMP.Model()
hier = IMP.atom.read_pdb(
    IMP.bff.get_example_path('structure/T4L/3GUN.pdb'), 
    mdl
)

av_parameter = {
    "linker_length": 20.0,
    "radii": (3.5, 0.0, 0.0),
    "linker_width": 0.5,
    "allowed_sphere_radius": 2.0,
    "contact_volume_thickness": 0.0,
    "contact_volume_trapped_fraction": -1,
    "simulation_grid_resolution": 0.5
}

def get_av(
    hier,
    residue_index: int = 132,
    atom_name: str = "CB",
    av_parameter: dict = av_parameter,
):
    av_p = IMP.Particle(mdl)
    sel = IMP.atom.Selection(hier)
    sel.set_atom_type(IMP.atom.AtomType(atom_name))
    sel.set_residue_index(residue_index)        
    source = sel.get_selected_particles()[0]
    IMP.bff.AV.do_setup_particle(mdl, av_p, source, **av_parameter)
    av = IMP.bff.AV(mdl, av_p)
    return av


class Tests(unittest.TestCase):
    """
    Tests for the bff.ProbeNetworkRestraint class.
    """

    def test_decorate_particle(self):
        # %%
        # The AV network restraint takes a fps.json file as an input.
        # fps.json files can contain multiple sets of distances for scoring
        # structures.
        fps_json_path = IMP.bff.get_example_path("structure/T4L/fret.fps.json")
        score_set_c1 = "chi2_C2_33p"
        fret_restraint = IMP.bff.ProbeNetworkRestraint(
            hier, str(fps_json_path), 
            score_set=score_set_c1,
            n_samples=500000
        )

        # Default mode (PRD-105): lattice-anchored grids, quadrature
        # distances -- deterministic, so the pin is tight. The legacy value
        # (11.918, MC) is pinned in test_av_lattice.py.
        v = fret_restraint.unprotected_evaluate(None)
        # Changed 2026-08-24 twice: the asymmetric chi2 became one function
        # (this restraint had the error bars swapped), and the restraint is
        # now worth 0.5 * chi2 -- what a Gaussian is worth, and what
        # IMP::core::Harmonic scores -- rather than the 0.25 a doubled
        # halving left it at.
        # this restraint had the error bars swapped (a model that was too
        # large was judged against error_neg), so an asymmetric pair now
        # divides by the other error. The AV distances are untouched -- the
        # model_ref assertions below are unchanged -- only the score they are
        # weighed with. See okf/validation/two_chi2_conventions.md.
        # 22.699... until 2026-08-27, when the shipped `strip_mask` of each
        # position began to be honoured (PRD-106): every T4L position asks for
        # its own side chain to be removed, and the volumes computed against
        # the unstripped structure were smaller. `test_strip_acceptance.py`
        # pins the mask behaviour itself; `test_av_lattice.py` keeps measuring
        # the lattice on a mask-free copy of this fixture, where the old number
        # still holds.
        # 22.5277 until 2026-09-01, when the accessible contact volume every
        # position of this file asks for stopped being ignored (PRD-121 G9).
        # The weighting pulls each cloud towards the surface, the 33 <R_DA>
        # shorten by 3.04 A and the chi-square halves.
        # Olga's radii were the default for part of 2026-09-01 and gave
        # 11.286617190950983 here (the 33 <R_DA> a further 0.40 A shorter,
        # 43.88 -> 43.48, and a better +0.44 A / 3.91 A against the file's own
        # experimental distances). IMP's own radii are the default again
        # (AV::set_radii_source) and this returns exactly to what it was
        # before that day -- the flip moved no other input.
        self.assertAlmostEqual(11.326777201821047, v, places=6)
        self.assertEqual(v, fret_restraint.unprotected_evaluate(None))

        # Re-measured 2026-09-01, when the accessible contact volume this file
        # asks for at every site stopped being ignored (PRD-121 G9). The
        # previous row was the same distances with the weighting inert; every
        # one of them shortened, and the agreement with `experiment_ref`
        # improved from a bias of +3.93 A (rmsd 5.78) to **+0.84 A** (rmsd
        # 3.93) -- which is the check that this is the model the file's fitted
        # trapped fractions were meant to produce.
        # Re-measured again 2026-09-01 with Olga's radii (see the score
        # above): every one of the 33 shortens once more, by 0.40 A on
        # average, and the bias against `experiment_ref` falls to +0.44 A.
        model_ref = np.array(
            [
                45.9, 36.7, 38.2, 48.7, 36.6, 38.7, 42.9, 50.5, 55.4, 48.3, 50.7,
                27.2, 47.5, 41.5, 42.2, 58.1, 53.3, 50.2, 33.5, 51.9, 48.8, 40.4,
                47.9, 40.6, 43.4, 44.7, 33.7, 37.2, 39.3, 45.5, 42.3, 30.5, 42.5
            ]
        )
        experiment_ref = np.array(
            [
                44.7, 39.7, 39.1, 47.2, 36.8, 37.5, 42.4, 50.6, 55.2, 48.1, 48.1,
                30.7, 45.8, 45.5, 42.2, 56.9, 46.1, 47.6, 37.1, 54.5, 49.3, 48.6,
                45.1, 40.3, 47.8, 42.2, 36.6, 34.2, 30.6, 39.8, 33.8, 37.8, 38.2
            ]
        )

        experimental_distances = fret_restraint.get_used_distances()
        model = list()
        experiment = list()
        for key in experimental_distances:
            e_dist = experimental_distances[key]
            d = json.loads(e_dist.get_json())
            model_distance = fret_restraint.get_model_distance(
                d["position1_name"],
                d["position2_name"],
                d["Forster_radius"],
                d["distance_type"]
            )
            e_dist.score_model(model_distance)
            model.append(model_distance)
            experiment.append(d['distance'])
        np.testing.assert_almost_equal(model_ref, model, decimal=0)
        np.testing.assert_almost_equal(experiment_ref, experiment_ref, decimal=0)


class R1PositionsTests(unittest.TestCase):
    """fps.json rotamer-ensemble positions (R1, PRD-108) and the C++ scorer."""

    def test_r1_position_warns_and_the_docking_filter_removes_it(self):
        import IMP.bff as fio
        payload = json.load(open(IMP.bff.get_example_path("structure/T4L/fret.fps.json")))
        positions, distances = payload["Positions"], payload["Distances"]
        first_pos = next(iter(positions))
        # an R1 position next to the AV ones, and a distance to it
        positions["r1_probe"] = {
            "chain_identifier": positions[first_pos]["chain_identifier"],
            "residue_seq_number": positions[first_pos]["residue_seq_number"],
            "atom_name": "CA", "simulation_type": "R1",
            "rotamer_library": "AlexaFluor 488 C1R cutoff30",
        }
        distances["r1_probe_dist"] = {
            "position1_name": "r1_probe", "position2_name": first_pos,
            "distance": 50.0, "error_neg": 5.0, "error_pos": 5.0,
            "Forster_radius": 52.0, "distance_type": "RDAMeanE",
        }
        score_sets = json.loads(json.dumps(payload.get("χ²", {})))
        # put the R1 distance into the score set so the restraint decorates it
        score_sets["chi2_C2_33p"]["distances"].append("r1_probe_dist")

        def _run(fps_path):
            log_path = fps_path + ".log"
            IMP.set_log_level(IMP.WARNING)
            target = IMP.SetLogTarget(IMP.TextOutput(log_path))
            r = IMP.bff.ProbeNetworkRestraint(hier, fps_path, score_set="chi2_C2_33p", n_samples=1000)
            r.unprotected_evaluate(None)
            del target
            return open(log_path).read() if os.path.exists(log_path) else ""

        with tempfile.TemporaryDirectory() as tmp:
            mixed = os.path.join(tmp, "mixed.fps.json")
            json.dump({"Positions": positions, "Distances": distances, "χ²": score_sets}, open(mixed, "w"))
            # the C++ side warns loudly when it meets a non-AV simulation_type
            self.assertIn("simulation_type 'R1'", _run(mixed))
            # the same filter a docking run must apply first
            kept_doc = fio.fps_positions_for_docking(
                json.dumps(positions), json.dumps(distances))
            kept_p = json.loads(kept_doc.positions)
            kept_d = json.loads(kept_doc.distances)
            self.assertNotIn("r1_probe", kept_p)
            self.assertNotIn("r1_probe_dist", kept_d)
            self.assertEqual(set(kept_p), set(positions) - {"r1_probe"})
            clean = os.path.join(tmp, "clean.fps.json")
            fio.write_fps_json(clean, json.dumps(kept_p), json.dumps(kept_d),
                               json.dumps(payload.get("χ²", {})))
            self.assertNotIn("simulation_type", _run(clean))




class TestProbeNetworkRestraintSet(unittest.TestCase):
    """`probe_network_restraint_set`: what `AVNetworkRestraintWrapper` did.

    The wrapper was a `%pythoncode` class subclassing
    `IMP.pmi.restraints.RestraintBase`, and nothing tested it -- it was named
    in `test_public_api_names.py` as a lazy name that resolves, and that was
    all. The two branches are tested here.
    """

    fps_json = str(IMP.bff.get_example_path("structure/T4L/fret.fps.json"))
    score_set = "chi2_C2_33p"

    def _hierarchy(self):
        m = IMP.Model()
        h = IMP.atom.read_pdb(
            IMP.bff.get_example_path("structure/T4L/3GUN.pdb"), m)
        return m, h

    def test_network_branch_scores_like_the_restraint(self):
        """Without `mean_position_restraint` the set holds the network itself,
        so it scores what `ProbeNetworkRestraint` scores."""
        m, h = self._hierarchy()
        rs = IMP.bff.probe_network_restraint_set(
            h, self.fps_json, score_set=self.score_set)
        direct = IMP.bff.ProbeNetworkRestraint(h, self.fps_json,
                                            score_set=self.score_set)
        self.assertAlmostEqual(rs.unprotected_evaluate(None),
                               direct.unprotected_evaluate(None), places=6)

    def test_mean_position_branch_is_one_restraint_per_distance(self):
        """With it, the set holds an `AVMeanDistanceRestraint` per distance --
        an approximation, so it is a different number, not the same one."""
        m, h = self._hierarchy()
        n_distances = len(
            IMP.bff.ProbeNetworkRestraint(
                h, self.fps_json,
                score_set=self.score_set).get_used_distances())
        rs = IMP.bff.probe_network_restraint_set(
            h, self.fps_json, score_set=self.score_set,
            mean_position_restraint=True)
        self.assertEqual(len(rs.get_restraints()), n_distances)
        self.assertGreater(rs.unprotected_evaluate(None), 0.0)

    def test_occupy_volume_gives_the_volumes_a_radius_and_a_mass(self):
        m, h = self._hierarchy()
        rs = IMP.bff.probe_network_restraint_set(
            h, self.fps_json, score_set=self.score_set, occupy_volume=True)
        avs = IMP.bff.ProbeNetworkRestraint(
            h, self.fps_json, score_set=self.score_set).get_used_avs()
        self.assertGreater(len(avs), 0)
        m2, h2 = self._hierarchy()
        net = IMP.bff.ProbeNetworkRestraint(h2, self.fps_json,
                                         score_set=self.score_set)
        used = net.get_used_avs()
        IMP.bff.set_av_xyzr_mass(used)
        for av in used:
            p = av.get_particle()
            self.assertTrue(IMP.core.XYZR.get_is_setup(p))
            self.assertTrue(IMP.atom.Mass.get_is_setup(p))
            self.assertGreater(IMP.core.XYZR(p).get_radius(), 0.0)

    def test_missing_fps_json_raises(self):
        m, h = self._hierarchy()
        with self.assertRaises(IOError):
            IMP.bff.probe_network_restraint_set(h, "/no/such/file.fps.json")
