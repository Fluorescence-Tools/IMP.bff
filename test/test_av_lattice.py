"""PRD-105: AV evaluation on the semi space-fixed lattice.

The battery runs over the six valid flag modes of ``AVNetworkRestraint``:

    mode                space_fixed  shared_map  distance
    default             yes          yes         quad
    shared-mc           yes          yes         mc
    lattice-private     yes          no          quad
    lattice-private-mc  yes          no          mc
    legacy-quad         no           no          quad
    legacy              no           no          mc

Tier 1 (bit-exact): occupancy, maps and mean position of every space-fixed
mode equal a forced full recompute; shared == private on the lattice; the
legacy mode equals the pre-PRD-105 code (pins in
``references/prd105_legacy_pins.json``).

Tier 2 (deterministic-accurate): quadrature distances are the same number
every call and independent of evaluation order; their error against the
exact double-sum oracle is within the stated budget. MC-distance modes carry
the pre-existing RNG order-dependence defect and are expected-failure on the
order test.
"""
from __future__ import division
import hashlib
import json
import os
import unittest

import numpy as np

import IMP
import IMP.atom
import IMP.core
import IMP.bff

try:
    import RMF
    import IMP.rmf
    HAVE_RMF = True
except ImportError:
    HAVE_RMF = False


MODES = {
    "default":            (True,  True,  "quad"),
    "shared-mc":          (True,  True,  "mc"),
    "lattice-private":    (True,  False, "quad"),
    "lattice-private-mc": (True,  False, "mc"),
    "legacy-quad":        (False, False, "quad"),
    "legacy":             (False, False, "mc"),
}
SPACE_FIXED_MODES = [m for m, (sf, _, _) in MODES.items() if sf]
QUAD_MODES = [m for m, (_, _, d) in MODES.items() if d == "quad"]
MC_MODES = [m for m, (_, _, d) in MODES.items() if d == "mc"]

FPS_JSON = IMP.bff.get_example_path("structure/T4L/fret.fps.json")
PDB = IMP.bff.get_example_path("structure/T4L/3GUN.pdb")
RMF_FN = IMP.bff.get_example_path("structure/T4L/t4l_docking.rmf3")
HERE = os.path.dirname(os.path.abspath(__file__))
LEGACY_PINS = os.path.join(HERE, "references", "prd105_legacy_pins.json")

# Random samples for the mc modes: enough that a scoring comparison is not
# dominated by noise, few enough to keep the suite fast.
N_SAMPLES = 20000
QUAD_K = 100


def make_restraint(mode, hier, score_set="chi2_C1_33p", n_samples=N_SAMPLES,
                   quad_k=QUAD_K):
    space_fixed, shared_map, distance = MODES[mode]
    return IMP.bff.AVNetworkRestraint(
        hier, FPS_JSON, "AVNetworkRestraint_" + mode, score_set, n_samples,
        space_fixed=space_fixed, shared_map=shared_map, distance=distance,
        quad_k=quad_k)


def open_trajectory():
    m = IMP.Model()
    f = RMF.open_rmf_file_read_only(RMF_FN)
    hier = IMP.rmf.create_hierarchies(f, m)[0]
    IMP.rmf.load_frame(f, RMF.FrameID(0))
    return m, f, hier


def av_names(restraint):
    return sorted(json.loads(restraint.get_diagnostics_json())["avs"])


def densities(restraint):
    """{name: xyz-density array} of the restraint's own AV handles."""
    out = {}
    for name in av_names(restraint):
        av = restraint.get_used_av(name)
        out[name] = np.array(av.get_map().get_xyz_density(), dtype=np.float64)
    return out


def sha(arr):
    return hashlib.sha256(np.asarray(arr, dtype=np.float64).tobytes()).hexdigest()


def exact_mean_distance(p1, p2):
    """Exact weighted mean distance between two (x, y, z, w) clouds."""
    p1 = np.asarray(p1, dtype=np.float64)
    p2 = np.asarray(p2, dtype=np.float64)
    w1 = p1[:, 3] / p1[:, 3].sum()
    w2 = p2[:, 3] / p2[:, 3].sum()
    acc = 0.0
    step = 2000
    for i in range(0, len(p1), step):
        d = np.sqrt(((p1[i:i + step, None, :3] - p2[None, :, :3]) ** 2).sum(-1))
        acc += (w1[i:i + step, None] * w2[None, :] * d).sum()
    return acc


def exact_mean_efficiency(p1, p2, r0):
    p1 = np.asarray(p1, dtype=np.float64)
    p2 = np.asarray(p2, dtype=np.float64)
    w1 = p1[:, 3] / p1[:, 3].sum()
    w2 = p2[:, 3] / p2[:, 3].sum()
    acc = 0.0
    step = 2000
    for i in range(0, len(p1), step):
        d = np.sqrt(((p1[i:i + step, None, :3] - p2[None, :, :3]) ** 2).sum(-1))
        acc += (w1[i:i + step, None] * w2[None, :] / (1.0 + (d / r0) ** 6)).sum()
    return acc


AV_PARAMETER = {
    "linker_length": 20.0,
    "radii": (3.5, 0.0, 0.0),
    "linker_width": 0.5,
    "allowed_sphere_radius": 2.0,
    "contact_volume_thickness": 0.0,
    "contact_volume_trapped_fraction": -1,
    "simulation_grid_resolution": 0.5,
}


def make_av(mdl, hier, residue_index, av_parameter=AV_PARAMETER,
            space_fixed=True):
    av_p = IMP.Particle(mdl)
    sel = IMP.atom.Selection(hier)
    sel.set_atom_type(IMP.atom.AtomType("CB"))
    sel.set_residue_index(residue_index)
    source = sel.get_selected_particles()[0]
    IMP.bff.AV.do_setup_particle(mdl, av_p, source, **av_parameter)
    av = IMP.bff.AV(mdl, av_p)
    if not space_fixed:
        av.set_space_fixed(False)
    return av


class TestModeMatrix(unittest.TestCase):
    """Construction and one evaluation in every mode; invalid modes refuse."""

    @classmethod
    def setUpClass(cls):
        cls.mdl = IMP.Model()
        cls.hier = IMP.atom.read_pdb(PDB, cls.mdl)

    def test_every_mode_evaluates(self):
        for mode in MODES:
            with self.subTest(mode=mode):
                r = make_restraint(mode, self.hier, score_set="chi2_C2_33p")
                sf, sm, dist = MODES[mode]
                self.assertEqual(r.get_space_fixed(), sf)
                self.assertEqual(r.get_shared_map(), sm)
                self.assertEqual(r.get_distance_method(), dist)
                self.assertEqual(r.get_quad_k(), QUAD_K)
                v = r.unprotected_evaluate(None)
                self.assertTrue(np.isfinite(v))
                d = json.loads(r.get_diagnostics_json())
                self.assertEqual(d["evaluations"], 1)
                self.assertEqual(len(d["avs"]), 17)
                if sm:
                    # T4L: one spacing, four inflation radii -> four classes
                    self.assertEqual(d["shared_map_classes"], 4)
                else:
                    self.assertEqual(d["shared_map_classes"], 0)

    def test_defaults_select_the_default_mode(self):
        r = IMP.bff.AVNetworkRestraint(self.hier, FPS_JSON,
                                       score_set="chi2_C2_33p")
        self.assertTrue(r.get_space_fixed())
        self.assertTrue(r.get_shared_map())
        self.assertEqual(r.get_distance_method(), "quad")
        self.assertEqual(r.get_quad_k(), 100)

    def test_keyword_arguments(self):
        # positional and keyword spellings agree (the constructor is
        # overloaded for deserialization; kwargs come from a SWIG shadow)
        r1 = IMP.bff.AVNetworkRestraint(self.hier, FPS_JSON, "n", "chi2_C2_33p",
                                        123, False, False, "mc", 7)
        r2 = IMP.bff.AVNetworkRestraint(self.hier, FPS_JSON, name="n",
                                        score_set="chi2_C2_33p", n_samples=123,
                                        space_fixed=False, shared_map=False,
                                        distance="mc", quad_k=7)
        for r in (r1, r2):
            self.assertEqual(r.get_n_samples(), 123)
            self.assertFalse(r.get_space_fixed())
            self.assertEqual(r.get_distance_method(), "mc")
            self.assertEqual(r.get_quad_k(), 7)
        with self.assertRaises(TypeError):
            IMP.bff.AVNetworkRestraint(self.hier, FPS_JSON, bogus=1)

    def test_shared_map_requires_space_fixed(self):
        with self.assertRaises(Exception):
            IMP.bff.AVNetworkRestraint(self.hier, FPS_JSON,
                                       score_set="chi2_C2_33p",
                                       space_fixed=False, shared_map=True)

    def test_distance_method_is_checked(self):
        with self.assertRaises(Exception):
            IMP.bff.AVNetworkRestraint(self.hier, FPS_JSON,
                                       score_set="chi2_C2_33p",
                                       distance="fft")

    def test_av_registry_requires_space_fixed(self):
        av = make_av(self.mdl, self.hier, 55, space_fixed=False)
        reg = IMP.bff.AVOccupancyRegistry(IMP.atom.get_leaves(self.hier))
        with self.assertRaises(Exception):
            av.set_occupancy_registry(reg)


@unittest.skipUnless(HAVE_RMF, "needs RMF")
class TestTier1Exactness(unittest.TestCase):
    """Bit-exact contract of the space-fixed modes."""

    def test_fast_path_equals_full_recompute(self):
        for mode in SPACE_FIXED_MODES:
            with self.subTest(mode=mode):
                m, f, hier = open_trajectory()
                r = make_restraint(mode, hier)
                frames = list(f.get_root_frames())
                for fr in frames[:4]:
                    IMP.rmf.load_frame(f, fr)
                    r.unprotected_evaluate(None)
                    for name in av_names(r):
                        av = r.get_used_av(name)
                        d_fast = np.array(av.get_map().get_xyz_density())
                        mp_fast = np.array(av.get_mean_position())
                        av.resample(True, True)  # force_full
                        d_full = np.array(av.get_map().get_xyz_density())
                        mp_full = np.array(av.get_mean_position())
                        self.assertEqual(d_fast.shape, d_full.shape)
                        self.assertTrue(np.array_equal(d_fast, d_full),
                                        "%s/%s differs" % (mode, name))
                        self.assertTrue(np.array_equal(mp_fast, mp_full))
                # the forced recomputes are counted as full updates
                d = json.loads(r.get_diagnostics_json())
                self.assertGreater(d["av_totals"]["full"], 0)

    def test_shared_equals_private(self):
        m1, f1, h1 = open_trajectory()
        m2, f2, h2 = open_trajectory()
        r_shared = make_restraint("default", h1)
        r_private = make_restraint("lattice-private", h2)
        for fr in list(f1.get_root_frames())[:5]:
            IMP.rmf.load_frame(f1, fr)
            IMP.rmf.load_frame(f2, fr)
            v1 = r_shared.unprotected_evaluate(None)
            v2 = r_private.unprotected_evaluate(None)
            self.assertEqual(v1, v2)
            d1 = densities(r_shared)
            d2 = densities(r_private)
            for name in d1:
                self.assertTrue(np.array_equal(d1[name], d2[name]), name)
                self.assertEqual(r_shared.get_used_av(name).get_lattice_window(),
                                 r_private.get_used_av(name).get_lattice_window())
        d = json.loads(r_shared.get_diagnostics_json())
        self.assertEqual(d["shared_map_classes"], 4)

    def test_lattice_is_absolute(self):
        # The window of an AV is a function of the lattice-quantised source
        # only: voxel 0 sits at k0 * spacing, edge = 2*floor(ll/h + 1/2) + 1.
        m, f, hier = open_trajectory()
        r = make_restraint("default", hier)
        r.unprotected_evaluate(None)
        for name in av_names(r):
            av = r.get_used_av(name)
            h = av.get_simulation_grid_resolution()
            ll = av.get_linker_length()
            kx, ky, kz, n = av.get_lattice_window()
            half = int(np.floor(ll / h + 0.5))
            self.assertEqual(n, 2 * half + 1)
            src = np.array(av.get_source_coordinates())
            q = np.floor(src / h + 0.5).astype(int)
            self.assertEqual([kx, ky, kz], list(q - half))
            origin = np.array(av.get_map().get_origin())
            np.testing.assert_allclose(origin, np.array([kx, ky, kz]) * h,
                                       atol=1e-5)
            # the source sits in the central voxel
            self.assertEqual(av.get_map().get_header().get_nx(), n)


class TestLegacyByteIdentity(unittest.TestCase):
    """space_fixed=False reproduces the pre-PRD-105 code exactly."""

    @classmethod
    def setUpClass(cls):
        with open(LEGACY_PINS) as fh:
            cls.pins = json.load(fh)
        cls.mdl = IMP.Model()
        cls.hier = IMP.atom.read_pdb(PDB, cls.mdl)

    def test_maps_and_means(self):
        for grid in (0.5, 1.5):
            for asr in (1.0, 2.0):
                par = dict(AV_PARAMETER, allowed_sphere_radius=asr,
                           simulation_grid_resolution=grid)
                for ri in (132, 55, 99):
                    key = "%s_%s_%s" % (grid, asr, ri)
                    pin = self.pins[key]
                    av = make_av(self.mdl, self.hier, ri, par, space_fixed=False)
                    self.assertFalse(av.get_space_fixed())
                    d = np.array(av.get_map().get_xyz_density())
                    self.assertEqual(len(d), pin["n"], key)
                    self.assertEqual(sha(d), pin["sha256_d"], key)
                    av.resample()
                    d2 = np.array(av.get_map().get_xyz_density())
                    self.assertEqual(len(d2), pin["n2"], key)
                    self.assertEqual(sha(d2), pin["sha256_d2"], key)
                    self.assertEqual(list(av.get_mean_position()), pin["mean"], key)

    def test_pinned_reference_values(self):
        # the pre-PRD-105 pins of test_AccessibleVolume / test_AVNetworkRestraint
        av = make_av(self.mdl, self.hier, 55, space_fixed=False)
        np.testing.assert_almost_equal(av.get_mean_position(),
                                       (-15.9244, 19.2183, 20.1207), decimal=3)
        av1 = make_av(self.mdl, self.hier, 132, space_fixed=False)
        av1.resample()
        np.testing.assert_allclose(IMP.core.XYZ(av1).get_coordinates(),
                                   (0.31708, -25.513668, -1.132486), rtol=0.1)
        r = make_restraint("legacy", self.hier, score_set="chi2_C2_33p",
                           n_samples=500000)
        self.assertAlmostEqual(11.917975852594935, r.unprotected_evaluate(None),
                               places=0)

    def test_restraint_mp_distances(self):
        r = make_restraint("legacy", self.hier, score_set="chi2_C2_33p",
                           n_samples=1000)
        r.unprotected_evaluate(None)
        for k, dd in r.get_used_distances().items():
            v = r.get_model_distance(dd.position_1, dd.position_2, 52.0,
                                     IMP.bff.DYE_PAIR_DISTANCE_MP)
            self.assertEqual(v, self.pins["restraint_mp"][k], k)

    @unittest.skipUnless(HAVE_RMF, "needs RMF")
    def test_trajectory_mp_distances(self):
        m, f, hier = open_trajectory()
        r = make_restraint("legacy", hier, n_samples=1000)
        for i, fr in enumerate(list(f.get_root_frames())[:5]):
            IMP.rmf.load_frame(f, fr)
            r.unprotected_evaluate(None)
            for k, dd in r.get_used_distances().items():
                v = r.get_model_distance(dd.position_1, dd.position_2, 52.0,
                                         IMP.bff.DYE_PAIR_DISTANCE_MP)
                self.assertEqual(v, self.pins["traj_mp"][i][k], (i, k))


@unittest.skipUnless(HAVE_RMF, "needs RMF")
class TestDeterminism(unittest.TestCase):
    """Tier 2: quadrature scores are call-order independent; MC is not."""

    def _scores_in_order(self, mode, order):
        m, f, hier = open_trajectory()
        r = make_restraint(mode, hier)
        frames = list(f.get_root_frames())
        out = {}
        for i in order:
            IMP.rmf.load_frame(f, frames[i])
            out[i] = r.unprotected_evaluate(None)
        return out

    def _check_order_independent(self, mode):
        a = self._scores_in_order(mode, [0, 1, 2, 3])
        b = self._scores_in_order(mode, [3, 1, 0, 2])
        for i in a:
            self.assertEqual(a[i], b[i], (mode, i))

    def test_quad_modes_are_order_independent(self):
        for mode in QUAD_MODES:
            with self.subTest(mode=mode):
                self._check_order_independent(mode)

    def test_repeat_evaluation_is_identical_in_quad_modes(self):
        for mode in QUAD_MODES:
            with self.subTest(mode=mode):
                m, f, hier = open_trajectory()
                r = make_restraint(mode, hier)
                v = [r.unprotected_evaluate(None) for _ in range(3)]
                self.assertEqual(v[0], v[1])
                self.assertEqual(v[1], v[2])


def _make_mc_order_test(mode):
    @unittest.expectedFailure
    def test(self):
        # Documents the pre-existing defect: MC distances draw from a global,
        # unseeded generator, so scores depend on the evaluation order.
        # Replaced by "quad", not patched (PRD-105).
        TestDeterminism._check_order_independent(self, mode)
    test.__name__ = "test_mc_order_independence_" + mode.replace("-", "_")
    test.__doc__ = "xfail: %s carries the MC RNG order-dependence defect" % mode
    return test


for _mode in MC_MODES:
    _t = _make_mc_order_test(_mode)
    setattr(TestDeterminism, _t.__name__, _t)


class TestQuadratureAccuracy(unittest.TestCase):
    """Tier 2 error budget of the lattice quadrature against the exact oracle."""

    # 1.5 A grid: clouds of a few thousand points, so the exact O(N*M)
    # oracle stays cheap; the budget is stated for the 2.0 A screening grid.
    PAR = dict(AV_PARAMETER, simulation_grid_resolution=1.5)

    @classmethod
    def setUpClass(cls):
        cls.mdl = IMP.Model()
        cls.hier = IMP.atom.read_pdb(PDB, cls.mdl)
        cls.av1 = make_av(cls.mdl, cls.hier, 132, cls.PAR)
        cls.av2 = make_av(cls.mdl, cls.hier, 55, cls.PAR)
        cls.p1 = np.array(cls.av1.get_map().get_xyz_density())
        cls.p2 = np.array(cls.av2.get_map().get_xyz_density())
        cls.exact_mean = exact_mean_distance(cls.p1, cls.p2)
        cls.exact_eff = exact_mean_efficiency(cls.p1, cls.p2, 52.0)

    def test_quadrature_points(self):
        for k in (10, 100, 1000):
            pts = np.array(self.av1.get_quadrature_points(k)).reshape(-1, 4)
            self.assertLessEqual(len(pts), k)
            self.assertGreater(len(pts), 0)
            # block centroids preserve the weighted mean exactly (to fp)
            mean_cloud = (self.p1[:, :3] * self.p1[:, 3:]).sum(0) / self.p1[:, 3].sum()
            mean_quad = (pts[:, :3] * pts[:, 3:]).sum(0) / pts[:, 3].sum()
            np.testing.assert_allclose(mean_quad, mean_cloud, atol=1e-9)
            self.assertAlmostEqual(pts[:, 3].sum(), self.p1[:, 3].sum(), places=6)
        # a small cloud is returned as is
        n = len(self.p1)
        pts = np.array(self.av1.get_quadrature_points(n + 1)).reshape(-1, 4)
        self.assertEqual(len(pts), n)

    def test_error_budget_k100(self):
        # Budget (PRD-105): MC-50k class, i.e. <= 0.16 A on the mean
        # distance. With the second-order block correction the measured
        # error is ~1e-4 A here and <= 0.005 A on T4L @ 2.0 A; assert an
        # order of magnitude inside the budget.
        q = IMP.bff.av_distance_quadrature(self.av1, self.av2, 52.0,
                                           IMP.bff.DYE_PAIR_DISTANCE_MEAN, 100)
        err = abs(q - self.exact_mean)
        print("\nquad K=100 mean-distance error vs exact: %.5f A" % err)
        self.assertLess(err, 0.016)
        e = IMP.bff.av_distance_quadrature(self.av1, self.av2, 52.0,
                                           IMP.bff.DYE_PAIR_EFFICIENCY, 100)
        print("quad K=100 efficiency error vs exact: %.6f" % abs(e - self.exact_eff))
        self.assertLess(abs(e - self.exact_eff), 0.0005)
        # R_E is the FRET-averaged distance of the same efficiency
        re = IMP.bff.av_distance_quadrature(self.av1, self.av2, 52.0,
                                            IMP.bff.DYE_PAIR_DISTANCE_E, 100)
        self.assertAlmostEqual(re, 52.0 * (1.0 / e - 1.0) ** (1.0 / 6.0), places=9)

    def test_k_curve_converges(self):
        errs = []
        ks = (10, 30, 100, 300, 1000)
        for k in ks:
            q = IMP.bff.av_distance_quadrature(self.av1, self.av2, 52.0,
                                               IMP.bff.DYE_PAIR_DISTANCE_MEAN, k)
            errs.append(abs(q - self.exact_mean))
        print("\nquad K-curve (K: err/A):",
              ", ".join("%d: %.4f" % (k, e) for k, e in zip(ks, errs)))
        self.assertLess(errs[-1], 0.002)
        self.assertLess(errs[-1], errs[0])
        # exact structural quantities are exact
        mp = IMP.bff.av_distance_quadrature(self.av1, self.av2, 52.0,
                                            IMP.bff.DYE_PAIR_DISTANCE_MP, 100)
        self.assertEqual(mp, IMP.bff.av_distance(self.av1, self.av2, 52.0,
                                                 IMP.bff.DYE_PAIR_DISTANCE_MP, 1))

    def test_empty_av_gives_nan(self):
        par = dict(self.PAR, allowed_sphere_radius=1.0)
        av3 = make_av(self.mdl, self.hier, 99, par)
        self.assertEqual(len(av3.get_map().get_xyz_density()), 0)
        v = IMP.bff.av_distance_quadrature(self.av1, av3, 52.0,
                                           IMP.bff.DYE_PAIR_DISTANCE_MEAN, 100)
        self.assertTrue(np.isnan(v))

    def test_quad_error_estimate(self):
        r = make_restraint("default", self.hier, score_set="chi2_C2_33p")
        r.unprotected_evaluate(None)
        est = r.get_quad_error_estimate(1000)
        print("\nrestraint quad error estimate (K=100 vs 1000): %.5f A" % est)
        self.assertLess(est, 0.05)
        self.assertGreater(est, 0.0)


@unittest.skipUnless(HAVE_RMF, "needs RMF")
class TestDiagnostics(unittest.TestCase):
    """Skip / local / full, rolls and moved-bead counts are reported."""

    def test_repeat_evaluation_skips(self):
        for mode in SPACE_FIXED_MODES:
            with self.subTest(mode=mode):
                m, f, hier = open_trajectory()
                r = make_restraint(mode, hier)
                r.unprotected_evaluate(None)
                d0 = json.loads(r.get_diagnostics_json())
                self.assertEqual(d0["av_totals"]["skip"], 0)
                self.assertEqual(d0["av_totals"]["full"], 17)
                for _ in range(5):
                    r.unprotected_evaluate(None)
                d = json.loads(r.get_diagnostics_json())
                self.assertEqual(d["av_totals"]["skip"], 5 * 17)
                self.assertEqual(d["av_totals"]["full"], 17)
                self.assertEqual(d["av_totals"]["local"], 0)
                self.assertEqual(d["evaluations"], 6)
                for sm in d["shared_maps"]:
                    # one full raster, then every update was a skip
                    self.assertEqual(sm["full"], 1)
                    self.assertEqual(sm["local"], 0)
                    self.assertGreaterEqual(sm["skip"], 5)
                    self.assertEqual(sm["moved_last"], 159)

    def test_trajectory_rolls_and_moves(self):
        m, f, hier = open_trajectory()
        r = make_restraint("default", hier)
        for fr in list(f.get_root_frames())[:6]:
            IMP.rmf.load_frame(f, fr)
            r.unprotected_evaluate(None)
        d = json.loads(r.get_diagnostics_json())
        # every bead moves in every frame of the T4L docking trajectory
        self.assertEqual(d["shared_map_classes"], 4)
        for sm in d["shared_maps"]:
            self.assertEqual(sm["moved_last"], 159)
            self.assertGreater(sm["voxels"], 0)
        self.assertGreater(d["av_totals"]["rolls"], 0)
        self.assertGreater(d["av_totals"]["full"], 17)
        for name, a in d["avs"].items():
            self.assertEqual(len(a["window"]), 4)
            self.assertEqual(a["skip"] + a["local"] + a["full"], 6)


class TestOccupancyMap(unittest.TestCase):
    """The lattice raster: deltas, rolls and shared reads are exact."""

    def _system(self, n=40, seed=1):
        rng = np.random.RandomState(seed)
        m = IMP.Model()
        ps = []
        for i in range(n):
            p = IMP.Particle(m)
            IMP.core.XYZR.setup_particle(
                p, IMP.algebra.Sphere3D(IMP.algebra.Vector3D(*rng.uniform(-8, 8, 3)),
                                        rng.uniform(1.0, 2.5)))
            ps.append(p)
        return m, ps

    def test_shared_read_equals_private_window(self):
        m, ps = self._system()
        shared = IMP.bff.AVOccupancyMap(1.0, 1.25, ps)
        shared.request_window(-6, -5, -4, 11, 11, 11)
        shared.update()
        priv = IMP.bff.AVOccupancyMap(1.0, 1.25, ps)
        priv.set_window(-6, -5, -4, 11, 11, 11)
        priv.update()
        a = np.array(shared.get_window(-6, -5, -4, 11, 11, 11))
        b = np.array(priv.get_window(-6, -5, -4, 11, 11, 11))
        self.assertTrue(np.array_equal(a, b))
        self.assertGreater(a.max(), 0)
        # outside the extent reads as zero
        far = np.array(shared.get_window(500, 500, 500, 3, 3, 3))
        self.assertEqual(far.sum(), 0)
        # windows are windows of the same lattice: a shifted window equals
        # the shifted read
        c = np.array(shared.get_window(-5, -5, -4, 11, 11, 11)).reshape(11, 11, 11)
        self.assertTrue(np.array_equal(c[:, :, :10], a.reshape(11, 11, 11)[:, :, 1:]))

    def test_local_delta_equals_full_raster(self):
        m, ps = self._system()
        occ = IMP.bff.AVOccupancyMap(1.0, 1.75, ps)
        occ.set_window(-6, -6, -6, 13, 13, 13)
        occ.update()
        self.assertEqual(occ.get_number_of_full_updates(), 1)
        # nothing moved: skip
        self.assertFalse(occ.update())
        self.assertEqual(occ.get_number_of_skips(), 1)
        # move a few particles: local delta
        for p in ps[:3]:
            x = IMP.core.XYZ(p)
            x.set_coordinates(x.get_coordinates() + IMP.algebra.Vector3D(0.7, -0.4, 1.1))
        self.assertTrue(occ.update())
        self.assertEqual(occ.get_number_of_local_updates(), 1)
        self.assertEqual(occ.get_number_of_moved_last(), 3)
        a = np.array(occ.get_window(-6, -6, -6, 13, 13, 13))
        ref = IMP.bff.AVOccupancyMap(1.0, 1.75, ps)
        ref.set_window(-6, -6, -6, 13, 13, 13)
        ref.update()
        b = np.array(ref.get_window(-6, -6, -6, 13, 13, 13))
        self.assertTrue(np.array_equal(a, b))
        # many moved: falls back to a full raster, still exact
        for p in ps:
            x = IMP.core.XYZ(p)
            x.set_coordinates(x.get_coordinates() + IMP.algebra.Vector3D(-0.3, 0.2, 0.1))
        occ.update()
        self.assertEqual(occ.get_number_of_full_updates(), 2)
        ref.update(True)
        self.assertTrue(np.array_equal(
            np.array(occ.get_window(-6, -6, -6, 13, 13, 13)),
            np.array(ref.get_window(-6, -6, -6, 13, 13, 13))))

    def test_roll_slab_plus_delta_equals_full_raster(self):
        m, ps = self._system(60, 2)
        occ = IMP.bff.AVOccupancyMap(1.0, 1.25, ps)
        occ.set_window(-6, -6, -6, 13, 13, 13)
        occ.update()
        for shift in ((1, 0, 0), (0, -2, 1), (3, 3, -3), (-1, 1, 0)):
            # a few atoms move, and the window rolls by whole voxels
            for p in ps[:2]:
                x = IMP.core.XYZ(p)
                x.set_coordinates(x.get_coordinates() + IMP.algebra.Vector3D(0.5, 0.5, -0.5))
            k = np.array(occ.get_extent()[:3]) + np.array(shift)
            occ.set_window(int(k[0]), int(k[1]), int(k[2]), 13, 13, 13)
            occ.update()
            ref = IMP.bff.AVOccupancyMap(1.0, 1.25, ps)
            ref.set_window(int(k[0]), int(k[1]), int(k[2]), 13, 13, 13)
            ref.update()
            a = np.array(occ.get_window(int(k[0]), int(k[1]), int(k[2]), 13, 13, 13))
            b = np.array(ref.get_window(int(k[0]), int(k[1]), int(k[2]), 13, 13, 13))
            self.assertTrue(np.array_equal(a, b), shift)
        self.assertEqual(occ.get_number_of_rolls(), 4)
        self.assertEqual(occ.get_number_of_local_updates(), 4)
        self.assertEqual(occ.get_number_of_full_updates(), 1)

    def test_grow_on_demand(self):
        m, ps = self._system()
        occ = IMP.bff.AVOccupancyMap(1.0, 1.25, ps)
        occ.update()
        e0 = occ.get_extent()
        self.assertEqual(occ.get_number_of_grows(), 1)
        occ.request_window(20, 20, 20, 5, 5, 5)
        occ.update()
        e1 = occ.get_extent()
        self.assertEqual(occ.get_number_of_grows(), 2)
        self.assertGreaterEqual(e1[0] + e1[3], 25)
        self.assertLessEqual(e1[0], e0[0])
        # covered request: no regrow, and nothing moved: skip
        occ.request_window(0, 0, 0, 3, 3, 3)
        occ.update()
        self.assertEqual(occ.get_number_of_grows(), 2)
        self.assertEqual(occ.get_number_of_skips(), 1)


class TestAVHandle(unittest.TestCase):
    """The AV decorator's lattice surface."""

    @classmethod
    def setUpClass(cls):
        cls.mdl = IMP.Model()
        cls.hier = IMP.atom.read_pdb(PDB, cls.mdl)

    def test_space_fixed_flag_lives_on_the_particle(self):
        av = make_av(self.mdl, self.hier, 55)
        self.assertTrue(av.get_space_fixed())
        av.set_space_fixed(False)
        self.assertFalse(IMP.bff.AV(self.mdl, av.get_particle()).get_space_fixed())
        av.set_space_fixed(True)
        self.assertTrue(av.get_space_fixed())

    def test_private_lattice_av_skips_and_recomputes(self):
        av = make_av(self.mdl, self.hier, 55)
        av.resample()
        self.assertEqual(av.get_number_of_full_updates(), 1)
        g = av.get_result_generation()
        av.resample()
        self.assertEqual(av.get_number_of_skips(), 1)
        self.assertEqual(av.get_result_generation(), g)
        d0 = np.array(av.get_map().get_xyz_density())
        av.resample(True, True)
        self.assertEqual(av.get_result_generation(), g + 1)
        self.assertTrue(np.array_equal(d0, np.array(av.get_map().get_xyz_density())))
        # a fresh handle on the same particle rebuilds and agrees
        av2 = IMP.bff.AV(self.mdl, av.get_particle())
        self.assertTrue(np.array_equal(d0, np.array(av2.get_map().get_xyz_density())))

    def test_used_av_shares_the_restraint_map(self):
        r = make_restraint("default", self.hier, score_set="chi2_C2_33p")
        r.unprotected_evaluate(None)
        name = av_names(r)[0]
        av = r.get_used_av(name)
        self.assertEqual(av.get_number_of_full_updates(), 1)
        self.assertIsNotNone(av.get_occupancy_registry())
        self.assertEqual(len(r.get_occupancy_registry().get_maps()), 4)


if __name__ == '__main__':
    unittest.main()
