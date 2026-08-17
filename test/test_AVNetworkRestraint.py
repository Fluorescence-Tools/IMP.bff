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
    Tests for the bff.AVNetworkRestraint class.
    """

    def test_decorate_particle(self):
        # %%
        # The AV network restraint takes a fps.json file as an input.
        # fps.json files can contain multiple sets of distances for scoring
        # structures.
        fps_json_path = IMP.bff.get_example_path("structure/T4L/fret.fps.json")
        score_set_c1 = "chi2_C2_33p"
        fret_restraint = IMP.bff.AVNetworkRestraint(
            hier, str(fps_json_path), 
            score_set=score_set_c1,
            n_samples=500000
        )

        # Default mode (PRD-105): lattice-anchored grids, quadrature
        # distances -- deterministic, so the pin is tight. The legacy value
        # (11.918, MC) is pinned in test_av_lattice.py.
        v = fret_restraint.unprotected_evaluate(None)
        self.assertAlmostEqual(13.505098587465483, v, places=6)
        self.assertEqual(v, fret_restraint.unprotected_evaluate(None))

        model_ref = np.array(
            [
                49.6, 39.3, 42.0, 53.0, 39.9, 41.1, 45.7, 53.3, 61.6, 52.3, 53.7,
                29.3, 50.7, 47.2, 46.6, 61.7, 57.2, 51.7, 35.8, 57.2, 54.4, 42.6,
                52.7, 42.8, 47.4, 47.2, 35.5, 37.5, 40.9, 51.0, 47.4, 33.8, 47.6
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
        import IMP.bff.fret.io as fio
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
            r = IMP.bff.AVNetworkRestraint(hier, fps_path, score_set="chi2_C2_33p", n_samples=1000)
            r.unprotected_evaluate(None)
            del target
            return open(log_path).read() if os.path.exists(log_path) else ""

        with tempfile.TemporaryDirectory() as tmp:
            mixed = os.path.join(tmp, "mixed.fps.json")
            json.dump({"Positions": positions, "Distances": distances, "χ²": score_sets}, open(mixed, "w"))
            # the C++ side warns loudly when it meets a non-AV simulation_type
            self.assertIn("simulation_type 'R1'", _run(mixed))
            # the Python filter is what a docking run must apply first
            kept_p, kept_d = fio.fps_positions_for_docking(positions, distances)
            self.assertNotIn("r1_probe", kept_p)
            self.assertNotIn("r1_probe_dist", kept_d)
            self.assertEqual(set(kept_p), set(positions) - {"r1_probe"})
            clean = os.path.join(tmp, "clean.fps.json")
            fio.write_fps_json(clean, kept_p, kept_d, payload.get("χ²", {}), validate=True)
            self.assertNotIn("simulation_type", _run(clean))


