from __future__ import division
import unittest

import tempfile

import math
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
    Tests for the bff.AV class.
    """

    def test_decorate_particle(self):
        atom_name = "CB"
        residue_index = 55
        mdl = IMP.Model()
        hier = IMP.atom.read_pdb(
            IMP.bff.get_example_path('structure/T4L/3GUN.pdb'), 
            mdl
        )        
        av_p = IMP.Particle(mdl)
        sel = IMP.atom.Selection(hier)
        sel.set_atom_type(IMP.atom.AtomType(atom_name))
        sel.set_residue_index(residue_index)        
        source = sel.get_selected_particles()[0]
        IMP.bff.AV.do_setup_particle(mdl, av_p, source, **av_parameter)
        av = IMP.bff.AV(mdl, av_p)
        np.testing.assert_almost_equal(av.get_source_coordinates(), (-11.589, 16.405, 17.556), decimal=3)
        np.testing.assert_almost_equal(av.get_mean_position(), (-15.9244, 19.2183, 20.1207), decimal=3)
        np.testing.assert_almost_equal(av.get_radii(), (3.5, 0, 0), decimal=3)
        self.assertEqual(av.get_parameters_are_optimized(), False)
        self.assertEqual(str(av.get_source()), '"Atom CB of residue 55"')

    def test_AccessibleVolumeDecorator(self):
        """
        Test the AccessibleVolumeDecorator class.
        """
        av1 = get_av(hier)
        av_mp = IMP.core.XYZ(av1)
        ref = (0., 0., 0.)
        np.testing.assert_allclose(av_mp.get_coordinates(), ref)
                
        av1.resample()  # Updates the AV
        ref = (0.31708, -25.513668, -1.132486)
        np.testing.assert_allclose(av_mp.get_coordinates(), ref, rtol=0.1)

    def test_access_av_feature(self):
        av1 = get_av(hier)
        av1_map = av1.get_map()
        bounds = 0.0, 20

        # PathMaps derive from IMP.em.DensityMap
        # write OBSTACLES to map
        with tempfile.NamedTemporaryFile(suffix=".mrc") as temp_file:
            IMP.em.write_map(av1_map, temp_file.name)

        pm_features = [
            IMP.bff.PM_TILE_PENALTY,             # Penality of visiting a tile
            IMP.bff.PM_TILE_COST,                # Cost of a path to the tile
            IMP.bff.PM_TILE_DENSITY,             # Density of tile
            IMP.bff.PM_TILE_COST_DENSITY,        # Cost * Density of tile
            IMP.bff.PM_TILE_PATH_LENGTH,         # Path length to tile (cost * grid spacing)
            IMP.bff.PM_TILE_PATH_LENGTH_DENSITY, # Path length to tile * density
            IMP.bff.PM_TILE_ACCESSIBLE_DENSITY,  # Density of tiles with path length in bounds
            # IMP.bff.PM_TILE_FEATURE,             # Additional feature of tile (accessed by name)
            # IMP.bff.PM_TILE_ACCESSIBLE_FEATURE   # Feature of tile with path length in bounds
        ]
        erw = IMP.em.MRCReaderWriter()
        
        for feature in pm_features:
            with tempfile.NamedTemporaryFile(suffix=".mrc") as temp_file:
                fn = temp_file.name
                fn_ref =  "./references/av_reference_%s.mrc" % feature
                if create_references:
                    fn = fn_ref
                IMP.bff.write_path_map(av1_map, fn, feature, bounds)
                
                em_map_ref = IMP.em.DensityMap()
                em_map_ref = IMP.em.read_map(fn_ref, erw)

                em_map = IMP.em.DensityMap()
                em_map = IMP.em.read_map(fn, erw)

    def test_av_random_points(self):
        n_samples = 10
        # create an AV in an inaccessible region
        av_parameter = {
            "linker_length": 20.0,
            "radii": (3.5, 0.0, 0.0),
            "linker_width": 0.5,
            "allowed_sphere_radius": 1.0,
            "contact_volume_thickness": 0.0,
            "contact_volume_trapped_fraction": -1,
            "simulation_grid_resolution": 0.5
        }
        av1 = get_av(hier, 99, av_parameter=av_parameter)
        
        # There should be no points in the map
        m1 = av1.get_map()
        k1 = m1.get_xyz_density()
        self.assertEqual(len(k1), 0)

        # No random samples on the map should be drawen
        p1 = IMP.bff.av_random_points(av1,  n_samples)
        self.assertEqual(len(p1), 0)

        # Increase the allowed sphere radius
        av_parameter['allowed_sphere_radius'] = 2.0
        av2 = get_av(hier, 99, av_parameter=av_parameter)
        m2 = IMP.bff.av_random_points(av2,  n_samples)
        self.assertEqual(len(m2), n_samples * 4)

    def test_av_random_distances(self):
        av1 = get_av(hier)
        av2 = get_av(hier, residue_index=55)
        distances = IMP.bff.av_random_distances(av1, av2, 500000)
        self.assertAlmostEqual(np.mean(distances), 54.38254912351925, places=1)

    def test_av_av_distance(self):
        av1 = get_av(hier)
        av2 = get_av(hier, residue_index=55)
        forster_radius = 52.0
        n_samples = 500000
        distance_types = [
            IMP.bff.DYE_PAIR_DISTANCE_E,             # Mean FRET averaged distance R_E
            IMP.bff.DYE_PAIR_DISTANCE_MEAN,          # Mean distance <R_DA>
            IMP.bff.DYE_PAIR_DISTANCE_MP,            # Distance between AV mean positions
            IMP.bff.DYE_PAIR_EFFICIENCY,             # Mean FRET efficiency
            # IMP.bff.DYE_PAIR_DISTANCE_DISTRIBUTION,  # (reserved for Distance distributions)
            # IMP.bff.DYE_PAIR_XYZ_DISTANCE            # Distance between XYZ of dye particles
        ]
        refs_distances = [53.614223677143364, 54.379459680068706, 52.11946290280115, 0.4563226915509543]
        for t, ref in zip(distance_types, refs_distances):
            v = IMP.bff.av_distance(
                av1, av2,
                forster_radius=forster_radius,
                distance_type=t,
                n_samples=n_samples
            )
            self.assertAlmostEqual(v, ref, places=1)
        
        # Test distance between AV and empty AV
        # create an AV in an inaccessible region
        av_parameter = {
            "linker_length": 20.0,
            "radii": (3.5, 0.0, 0.0),
            "linker_width": 0.5,
            "allowed_sphere_radius": 1.0,
            "contact_volume_thickness": 0.0,
            "contact_volume_trapped_fraction": -1,
            "simulation_grid_resolution": 0.5
        }
        av3 = get_av(hier, 99, av_parameter=av_parameter)
        v = IMP.bff.av_distance(
                av1, av3,
                forster_radius=forster_radius,
                distance_type=t,
                n_samples=n_samples
        )
        # If a distance cannot be computed returns a nan
        self.assertEqual(math.isnan(v), True)

    def test_distance_distributions(self):
        av1 = get_av(hier)
        av2 = get_av(hier, residue_index=55)
        rda_start, rda_stop, n_bins = 0, 100, 32
        n_samples = 10000
        rda = np.linspace(rda_start, rda_stop, n_bins)
        p_rda = IMP.bff.av_distance_distribution(av1, av2, rda, n_samples=n_samples)
        p_rda_ref = np.array(
            [   0.,    0.,    0.,    0.,    0.,    0.,    0.,    7.,   16.,
                62.,  144.,  235.,  419.,  646.,  919., 1131., 1348., 1425.,
                1305., 1085.,  716.,  392.,  126.,   24.,    0.,    0.,    0.,
                0.,    0.,    0.,    0.,    0.])
        ssdev = np.sum((p_rda_ref - p_rda)**2.)
        self.assertEqual(ssdev < 30000, True)

