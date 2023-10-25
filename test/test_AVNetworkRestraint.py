from __future__ import division
import unittest

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
        with open(fps_json_path) as fp:
            fps_json = json.load(fp)
        score_set_c1 = "chi2_C2_33p"
        fret_restraint = IMP.bff.AVNetworkRestraint(
            hier, str(fps_json_path), 
            score_set=score_set_c1,
            n_samples=500000
        )

        v = fret_restraint.unprotected_evaluate(None)
        self.assertAlmostEqual(11.917975852594935, v, places=0)

        model_ref = np.array(
            [
                49. , 38. , 41.7, 52.6, 38.8, 38.7, 44.3, 52. , 60.3, 50.3, 52.7,
                30. , 49.7, 46. , 46.7, 61.5, 56.2, 51.6, 35.4, 57.2, 53.5, 43.2,
                52.2, 43.4, 47.5, 47.4, 35.3, 37.9, 41.3, 50.8, 45.9, 33.1, 47.4
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
