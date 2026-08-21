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
        # Lattice-anchored (PRD-105) value; the legacy source-anchored value
        # (-15.9244, 19.2183, 20.1207) is pinned in test_av_lattice.py.
        # Values moved 2026-08-19 with the search stencil: the default is now 74
        # (the LabelLib reference metric) rather than 26. The 26 stencil's
        # isopath surface is cubic, so it under-reached along the diagonals and
        # returned ~15 % less volume -- 26146 voxels against an analytic 33510
        # for an obstacle-free 20 A linker, where 74 gives 30682. The AV is
        # therefore both larger and differently shaped, and every quantity
        # derived from it moves. See AV::get_search_stencil.
        np.testing.assert_almost_equal(av.get_mean_position(), (-15.0462, 19.6699, 20.228), decimal=3)
        np.testing.assert_almost_equal(av.get_radii(), (3.5, 0, 0), decimal=3)
        self.assertEqual(av.get_parameters_are_optimized(), False)
        self.assertEqual(str(av.get_source()), '"Atom CB of residue 55"')

    def test_get_radii_returns_all_three(self):
        """`get_radii` must return radius1, radius2, radius3 -- in that order.

        It returned ``{radius1, radius3, radius3}``: radius2 was never returned,
        contradicting the method's own docstring. The existing assertion above
        could not catch it, because it uses an AV1 dye where radius2 and radius3
        are both 0 and the duplication is invisible. Three *distinct* radii are
        what make the bug observable, so that is what this pins.

        The only consumer is ``restraints/network.py``'s ``max(av.get_radii())``,
        so the effect was latent -- wrong only when radius2 is the largest, which
        needs an AV3 dye. Nothing in the repo uses AV3 yet, which is exactly why
        it survived.
        """
        mdl = IMP.Model()
        hier = IMP.atom.read_pdb(
            IMP.bff.get_example_path('structure/T4L/3GUN.pdb'), mdl)
        av_p = IMP.Particle(mdl)
        sel = IMP.atom.Selection(hier)
        sel.set_atom_type(IMP.atom.AtomType("CB"))
        sel.set_residue_index(55)
        source = sel.get_selected_particles()[0]
        IMP.bff.AV.do_setup_particle(mdl, av_p, source, **av_parameter)
        av = IMP.bff.AV(mdl, av_p)

        av.set_radius1(5.0)
        av.set_radius2(4.5)
        av.set_radius3(1.5)
        np.testing.assert_almost_equal(av.get_radii(), (5.0, 4.5, 1.5), decimal=6)

        # radius2 largest: the case the old code got wrong and max() propagated
        av.set_radius1(1.0)
        av.set_radius2(9.0)
        av.set_radius3(2.0)
        np.testing.assert_almost_equal(av.get_radii(), (1.0, 9.0, 2.0), decimal=6)
        self.assertAlmostEqual(max(av.get_radii()), 9.0, places=6)

    def test_av3_matches_labellib_rule(self):
        """AV3 must be Olga's AV3, i.e. LabelLib's.

        LabelLib grades by **all** the radii it is given: `excludeConcentricSpheres`
        (`FlexLabel/src/FlexLabel.cxx:235`) sorts them, builds
        `rhos = LinSpaced(n + 1, 0, 1)`, and writes `min(ref, rhos[i])` over the
        shell between consecutive `atom_vdW + radius[i]`. The density is
        therefore **the fraction of probe radii that fit**: {1/3, 2/3, 1}.

        bff ignored `radius2` and `radius3` entirely -- its "AV3" was AV1 at
        `radius1`, and order-dependent with it. The four properties below are what
        the rule implies, and each one failed before:

        1. density is the mean of the per-radius indicators;
        2. the level set at k/n equals the AV1 volume of the k-th largest radius;
        3. `AV3(r, r, r) == AV1(r)`;
        4. the result does not depend on the order of the radii (LabelLib sorts).
        """
        import numpy as np
        from IMP.bff import compute_av

        rng = np.random.default_rng(0)
        xyz = rng.normal(scale=8.0, size=(150, 3))
        vdw = np.full(len(xyz), 1.7)
        source = np.zeros(3)
        kw = dict(linker_length=20.0, linker_width=1.5,
                  grid_resolution=1.0, allowed_sphere_radius=2.55)

        def dens(radii):
            av = compute_av(xyz, vdw, source, kw["linker_length"],
                            kw["linker_width"], radii, kw["grid_resolution"],
                            kw["allowed_sphere_radius"])
            ng = av.get_ng()
            return np.asarray(av.get_density(), dtype=float).reshape(ng, ng, ng)

        R = (3.5, 2.5, 1.0)
        a3 = dens(R)
        a1 = {r: dens((r, 0.0, 0.0)) for r in R}

        # 1. the mean-of-indicators rule
        pred = sum((a1[r] > 0).astype(float) for r in R) / 3.0
        np.testing.assert_allclose(a3, pred, atol=1e-6)

        # 2. level sets are the per-radius AV1 volumes
        for lvl, r in ((1.0, 3.5), (2.0 / 3.0, 2.5), (1.0 / 3.0, 1.0)):
            self.assertEqual(int((a3 >= lvl - 1e-6).sum()),
                             int((a1[r] > 0).sum()))

        # 3. equal radii collapse to AV1
        for r in (2.0, 3.0):
            np.testing.assert_allclose(dens((r, r, r)), dens((r, 0.0, 0.0)),
                                       atol=1e-6)

        # 4. order independence
        for perm in ((1.0, 2.5, 3.5), (2.5, 3.5, 1.0)):
            np.testing.assert_allclose(dens(perm), a3, atol=1e-6)

        # and the grading must survive the front door, not be binarised away
        self.assertGreater(len(np.unique(np.round(a3[a3 > 0], 6))), 1)

    def test_stencil_74_matches_the_reference_metric(self):
        """The 74 stencil reproduces LabelLib's path metric; 26 does not.

        With **no obstacles at all** the accessible volume must be a sphere of
        the linker length, so this needs no reference implementation to check
        against -- the analytic answer is 4/3 pi L^3.

        The 26 stencil ({1,2,3} squared offsets) has a cubic isopath surface and
        under-reaches along the diagonals, losing ~15 % of the volume in every
        shell. The 74 stencil is LabelLib's `essentialNeighbours()` shell set
        ({1,2,3,5,6}); bff reaches it with a sqrt(6) neighbour radius, which also
        admits the six squared-length-4 axis jumps -- redundant, since their cost
        equals two unit steps, so the path lengths are the same.

        Not a default: 74's longest jump is ~2.45 voxels and tunnels through
        thinner walls. See `AV::get_search_stencil`.
        """
        import numpy as np

        L, H = 20.0, 1.0
        mdl = IMP.Model()
        root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(mdl))
        # one atom far away: obstacle-free without an empty particle list
        far = IMP.Particle(mdl)
        fd = IMP.core.XYZR.setup_particle(far)
        fd.set_coordinates(IMP.algebra.Vector3D(1e6, 1e6, 1e6))
        fd.set_radius(0.01)
        root.add_child(IMP.atom.Hierarchy.setup_particle(far))
        src = IMP.Particle(mdl)
        sd = IMP.core.XYZR.setup_particle(src)
        sd.set_coordinates(IMP.algebra.Vector3D(0, 0, 0))
        sd.set_radius(0.0)
        root.add_child(IMP.atom.Hierarchy.setup_particle(src))

        def volume(stencil):
            ap = IMP.Particle(mdl)
            IMP.bff.AV.do_setup_particle(
                mdl, ap, src, linker_length=L, linker_width=1.5,
                radii=IMP.algebra.Vector3D(1.0, 0.0, 0.0),
                allowed_sphere_radius=0.0, simulation_grid_resolution=H)
            av = IMP.bff.AV(mdl, ap)
            av.set_search_stencil(stencil)
            av.resample()
            pts = np.asarray(av.get_map().get_xyz_density(), dtype=float)
            return int((pts[:, 3] > 0).sum()) if pts.size else 0

        analytic = 4.0 / 3.0 * np.pi * L ** 3          # 33510
        n26, n74 = volume(26), volume(74)

        # 74 is within 10 % of the analytic sphere; 26 is not
        self.assertGreater(n74 / analytic, 0.90)
        self.assertLess(n26 / analytic, 0.85)
        # and 74 recovers the volume 26 loses
        self.assertGreater(n74, n26 * 1.10)

    def test_stencil_compensation_recovers_the_reference_volume(self):
        """The speed stencil, compensated, reproduces the reference volume.

        A coarse stencil overestimates path length, so its accessible volume is
        that of a *shorter* linker. The bias is a property of the stencil, not
        of the structure, so scaling the linker length removes it: measured
        obstacle-free over L = 12-25 A and dye radii 1.0-3.5 A, the scale is
        1.0551 +/- 0.0021.

        That makes the speed option usable: stencil 26 is ~1.9x faster than 74
        but returns ~84 % of its volume; compensated it returns ~99.5 % and is
        still faster (measured 2.1 ms against 3.5 ms on T4L site 22).

        Compensation is **opt-in**, so asking for a stencil gives that stencil's
        own answer unless it is requested -- that is what keeps the historical
        stencil-30 pins valid.
        """
        import numpy as np

        mdl = IMP.Model()
        root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(mdl))
        far = IMP.Particle(mdl)
        fd = IMP.core.XYZR.setup_particle(far)
        fd.set_coordinates(IMP.algebra.Vector3D(1e6, 1e6, 1e6))
        fd.set_radius(0.01)
        root.add_child(IMP.atom.Hierarchy.setup_particle(far))
        src = IMP.Particle(mdl)
        sd = IMP.core.XYZR.setup_particle(src)
        sd.set_coordinates(IMP.algebra.Vector3D(0, 0, 0))
        sd.set_radius(0.0)
        root.add_child(IMP.atom.Hierarchy.setup_particle(src))

        def volume(stencil, compensate):
            ap = IMP.Particle(mdl)
            IMP.bff.AV.do_setup_particle(
                mdl, ap, src, linker_length=20.0, linker_width=1.5,
                radii=IMP.algebra.Vector3D(1.0, 0.0, 0.0),
                allowed_sphere_radius=0.0, simulation_grid_resolution=1.0)
            av = IMP.bff.AV(mdl, ap)
            av.set_search_stencil(stencil)
            av.set_compensate_stencil(compensate)
            av.resample()
            pts = np.asarray(av.get_map().get_xyz_density(), dtype=float)
            return int((pts[:, 3] > 0).sum()) if pts.size else 0

        ref = volume(74, False)
        raw = volume(26, False)
        comp = volume(26, True)

        self.assertLess(raw / ref, 0.90)          # the bias is real
        self.assertGreater(comp / ref, 0.97)      # and compensation removes it
        self.assertLess(comp / ref, 1.03)

        # opt-in: compensation must not move the reference stencil at all
        self.assertEqual(volume(74, True), ref)

    def test_AccessibleVolumeDecorator(self):
        """
        Test the AccessibleVolumeDecorator class.
        """
        av1 = get_av(hier)
        av_mp = IMP.core.XYZ(av1)
        ref = (0., 0., 0.)
        np.testing.assert_allclose(av_mp.get_coordinates(), ref)
                
        av1.resample()  # Updates the AV
        # lattice-anchored (PRD-105), 74-stencil default since 2026-08-19
        ref = (-0.323894, -25.204008, -3.823595)
        np.testing.assert_allclose(av_mp.get_coordinates(), ref, atol=1e-3)

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
        self.assertAlmostEqual(np.mean(distances), 55.8469, delta=0.15)

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
        # lattice-anchored (PRD-105); MC with 500k samples
        refs_distances = [54.8151, 55.8598, 52.9990, 0.4211]
        for t, ref in zip(distance_types, refs_distances):
            v = IMP.bff.av_distance(
                av1, av2,
                forster_radius=forster_radius,
                distance_type=t,
                n_samples=n_samples
            )
            self.assertAlmostEqual(v, ref, delta=0.12)   # MC with 500k samples: ~3.5 sd
        
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
        # Regenerated 2026-08-19 with the 74-stencil default (mean of 15 x 10k
        # runs, so the reference itself is not one noisy draw). The distribution
        # broadens and shifts out because the AV is larger under the reference
        # metric; see AV::get_search_stencil.
        p_rda_ref = np.array(
            [   0.,    0.,    0.,    0.,    0.,    0.,    1.,    6.,   27.,
               58.,  124.,  220.,  357.,  547.,  778.,  991., 1176., 1319.,
             1334., 1204.,  928.,  588.,  267.,   71.,    5.,    0.,    0.,
                0.,    0.,    0.,    0.,    0.])
        # two independent 10k-sample histograms differ by ~2N = 20000 in
        # summed squared deviation on average; 30000 failed about one run
        # in ten
        ssdev = np.sum((p_rda_ref - p_rda)**2.)
        self.assertEqual(ssdev < 60000, True)

