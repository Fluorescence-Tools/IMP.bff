"""Four things an accessible volume could not do, and one it claimed to.

* **`minimum_distance`** -- the closest approach of two volumes. Every other
  distance here is an average; this is the bound a crosslink-style upper limit
  is written against.
* **`rmp_from_model_distance`** -- the inverse of the model distance in the
  separation of the two volumes. An experiment reports a distance; a restraint
  between two points needs an `Rmp`, and the two differ by several angstrom in
  a way that depends on the shape of both clouds.
* **`cloud_overlap` / `av_overlap`** -- how much of a volume touches something.
* **`write_av`** -- the volume in a format a viewer opens.
* **`chain_weighting`** -- an fps.json field that this module declared, and
  documented, and did not read. It reads it now.
"""

import os

import numpy as np
import pytest

import IMP
import IMP.atom
import IMP.core
import IMP.bff

MEAN = IMP.bff.PROBE_PAIR_DISTANCE_MEAN
MEAN_E = IMP.bff.PROBE_PAIR_DISTANCE_E
EFF = IMP.bff.PROBE_PAIR_EFFICIENCY
MP = IMP.bff.PROBE_PAIR_DISTANCE_MP
MIN = IMP.bff.PROBE_PAIR_DISTANCE_MIN


@pytest.fixture(scope="module")
def clouds():
    """Two gaussian blobs, the second 40 A along x from the first."""
    rng = np.random.RandomState(0)
    a = np.column_stack([rng.normal(0, 5, (2000, 3)), np.ones(2000)])
    b = np.column_stack([rng.normal(0, 5, (2000, 3)) + [40, 0, 0], np.ones(2000)])
    return a, b


# ---------------------------------------------------------------------------
# the vocabulary
# ---------------------------------------------------------------------------

def test_rmin_joins_the_fps_json_distance_vocabulary():
    assert IMP.bff.probe_pair_distance_type_name(MIN) == "Rmin"
    assert IMP.bff.probe_pair_distance_type("Rmin") == MIN
    # and the names that were already there still answer
    for name in ("RDAMean", "RDAMeanE", "Rmp", "Efficiency", "pRDA"):
        assert IMP.bff.probe_pair_distance_type_name(
            IMP.bff.probe_pair_distance_type(name)) == name


# ---------------------------------------------------------------------------
# minimum distance
# ---------------------------------------------------------------------------

def test_the_minimum_is_a_bound_and_not_an_average(clouds):
    a, b = clouds
    A, B = a.ravel(), b.ravel()
    lo = IMP.bff.minimum_distance(A, B)
    assert lo < IMP.bff.cloud_model_distance(A, B, MEAN)
    assert lo < IMP.bff.cloud_model_distance(A, B, MP)
    # it is exactly the smallest pair distance, which a sample would over-state
    exact = np.sqrt(((a[:, None, :3] - b[None, :, :3]) ** 2).sum(-1)).min()
    assert lo == pytest.approx(exact, abs=1e-9)
    assert IMP.bff.cloud_model_distance(A, B, MIN) == pytest.approx(lo)


def test_a_zero_weight_point_is_not_in_the_volume():
    """Radius-zero masking leaves zero-weight points in the array; a minimum
    that counted them would report a contact that is not there."""
    near = np.array([[0., 0., 0., 1.], [10., 0., 0., 0.]])
    far = np.array([[20., 0., 0., 1.]])
    # 20, from the weighted point -- not 10, from the zero-weight one
    assert IMP.bff.minimum_distance(near.ravel(), far.ravel()) == \
        pytest.approx(20.0)


def test_empty_clouds_give_an_infinite_minimum():
    assert np.isinf(IMP.bff.minimum_distance([], [1., 2., 3., 1.]))


# ---------------------------------------------------------------------------
# the inversion
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("dtype,target", [(MEAN, 45.0), (MEAN_E, 45.0),
                                          (EFF, 0.35)])
def test_the_inversion_round_trips(clouds, dtype, target):
    a, b = clouds
    A, B = a.ravel(), b.ravel()
    r = IMP.bff.rmp_from_model_distance(A, B, target, dtype, 52.0)

    m1 = np.asarray(IMP.bff.points_weighted_mean(A))
    m2 = np.asarray(IMP.bff.points_weighted_mean(B))
    u = m2 - m1
    r0 = np.linalg.norm(u)
    shifted = b.copy()
    shifted[:, :3] += (r - r0) * (u / r0)

    assert IMP.bff.cloud_model_distance(A, shifted.ravel(), MP) == \
        pytest.approx(r, abs=1e-6), "the answer is a mean-position separation"
    assert IMP.bff.cloud_model_distance(A, shifted.ravel(), dtype, 52.0) == \
        pytest.approx(target, abs=0.01), "and it reproduces the target there"


def test_rmp_and_xyz_conventions_invert_to_themselves(clouds):
    a, b = clouds
    for dtype in (MP, IMP.bff.PROBE_PAIR_XYZ_DISTANCE):
        assert IMP.bff.rmp_from_model_distance(
            a.ravel(), b.ravel(), 37.0, dtype) == pytest.approx(37.0)


def test_a_distance_the_volumes_cannot_reach_is_refused(clouds):
    """Silence would be worse: a restraint built on an unreachable target is a
    restraint that pulls forever."""
    a, b = clouds
    with pytest.raises(Exception) as e:
        IMP.bff.rmp_from_model_distance(a.ravel(), b.ravel(), 2.0, MEAN, 52.0)
    assert "no separation reproduces" in str(e.value)


def test_the_mean_distance_exceeds_the_mean_position_distance(clouds):
    """Why the inversion is needed at all: the two are not the same number and
    the gap is not a constant."""
    a, b = clouds
    A, B = a.ravel(), b.ravel()
    assert IMP.bff.cloud_model_distance(A, B, MEAN) > \
        IMP.bff.cloud_model_distance(A, B, MP)


# ---------------------------------------------------------------------------
# overlap
# ---------------------------------------------------------------------------

def test_overlap_is_a_weight_fraction(clouds):
    a, b = clouds
    A = a.ravel()
    assert IMP.bff.cloud_overlap(A, a[:60, :3].ravel(), 6.0) > 0.5
    assert IMP.bff.cloud_overlap(A, b[:60, :3].ravel(), 6.0) == 0.0
    assert IMP.bff.cloud_overlap(A, a[:, :3].ravel(), 1000.0) == \
        pytest.approx(1.0)
    assert IMP.bff.cloud_overlap(A, [], 5.0) == 0.0


def test_overlap_rejects_a_reference_that_is_not_three_per_point(clouds):
    a, _ = clouds
    with pytest.raises(Exception):
        IMP.bff.cloud_overlap(a.ravel(), [1.0, 2.0], 5.0)


# ---------------------------------------------------------------------------
# a real volume: export, overlap, chain weighting
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def t4l_av():
    IMP.set_log_level(IMP.SILENT)
    pdb = IMP.bff.get_example_path("structure/T4L/3GUN.pdb")
    return IMP.bff.compute_av_from_structure(pdb, "A", 132, "CB", 20.5, 1.5,
                                             3.5, 0.0, 0.0, 1.5)


@pytest.mark.parametrize("ext", ["xyz", "pqr", "dx"])
def test_a_volume_writes_the_format_its_extension_names(t4l_av, tmp_path, ext):
    path = tmp_path / f"av.{ext}"
    IMP.bff.write_av(t4l_av, str(path))
    assert path.is_file() and path.stat().st_size > 0
    text = path.read_text()
    if ext == "xyz":
        assert int(text.splitlines()[0]) == t4l_av.get_n_points()
        assert len(text.splitlines()[2].split()) == 5   # element x y z weight
    elif ext == "pqr":
        assert text.startswith("ATOM  ") and text.rstrip().endswith("END")
        assert len(text.splitlines()) == t4l_av.get_n_points() + 1
    else:
        n = t4l_av.get_ng()
        assert f"counts {n} {n} {n}" in text
        assert "object 3 class array" in text
        assert "component \"data\" value 3" in text


def test_an_unknown_extension_is_refused(t4l_av, tmp_path):
    with pytest.raises(Exception) as e:
        IMP.bff.write_av(t4l_av, str(tmp_path / "av.bogus"))
    assert "extension" in str(e.value)


def test_a_volume_overlaps_its_own_site_and_not_the_far_side(tmp_path):
    IMP.set_log_level(IMP.SILENT)
    m = IMP.Model()
    h = IMP.atom.read_pdb(IMP.bff.get_example_path("structure/T4L/3GUN.pdb"), m,
                          IMP.atom.NonWaterNonHydrogenPDBSelector())
    src = IMP.atom.Selection(h, chain_id="A", residue_index=132,
                             atom_type=IMP.atom.AtomType("CB")
                             ).get_selected_particles()[0]
    av = IMP.bff.AV.setup_particle(m, IMP.Particle(m), src)
    av.set_av_parameter('{"linker_length": 20.5, "linker_width": 1.5,'
                        ' "radius1": 3.5, "simulation_grid_resolution": 1.5}')
    av.resample()

    # 8 A, not 5: the volume already excludes the protein's van der Waals
    # envelope inflated by the dye radius, so nothing is within ~5 A of an atom
    # centre and every overlap at that radius is zero
    near = IMP.bff.av_overlap(av, h, "resi 130-134", 8.0)
    far = IMP.bff.av_overlap(av, h, "resi 10", 8.0)
    whole = IMP.bff.av_overlap(av, h, "", 8.0)
    assert 0.0 < near <= 1.0
    assert far < near
    assert whole >= near
    # the selection language is the one this module already speaks, in either
    # dialect
    assert IMP.bff.av_overlap(av, h, "resid 130 to 134", 8.0) == \
        pytest.approx(near)
    # and the clearance itself: nothing touches at a radius the dye cannot
    # reach. The threshold is exactly the smallest heavy-atom radius of the
    # obstacle set plus the dye radius, because the carve inflates every
    # obstacle by the dye radius and the closest a voxel can sit to an atom
    # centre is that sum. Under the default radii -- IMP's, united-atom
    # (AV::set_radii_source) -- that is oxygen 1.70 + 3.5 = **5.20 A**:
    # measured 0.0 at 5.20 and 0.000366 at 5.21. Under Olga's, which were the
    # default for part of 2026-09-01, it is 1.49 + 3.5 = 4.99. Pinned from
    # both sides so a change of radii set shows up here as a moved edge
    # rather than as an assertion that still passes for the wrong reason.
    assert IMP.bff.av_overlap(av, h, "", 5.20) == 0.0
    assert IMP.bff.av_overlap(av, h, "", 5.21) > 0.0


# ---------------------------------------------------------------------------
# chain weighting
# ---------------------------------------------------------------------------

def test_the_uniform_weighting_is_the_unweighted_volume():
    w = IMP.bff.LinkerWeighting()
    assert w.get_is_uniform()
    for x in (0.0, 5.0, 500.0):
        assert w.get_weight(x) == 1.0
    assert w.get_supported_fraction(1.0) == 1.0


def test_a_table_interpolates_and_clamps():
    w = IMP.bff.LinkerWeighting(0.0, 10.0, [0.0, 1.0, 0.0])
    assert not w.get_is_uniform()
    assert w.get_weight(0.0) == pytest.approx(0.0)
    assert w.get_weight(5.0) == pytest.approx(1.0)
    assert w.get_weight(2.5) == pytest.approx(0.5)
    assert w.get_weight(7.5) == pytest.approx(0.5)
    # clamped, not NaN: one NaN weight makes every mean and distance NaN
    assert w.get_weight(-100.0) == pytest.approx(0.0)
    assert w.get_weight(1e6) == pytest.approx(0.0)
    assert np.isfinite(w.get_weight(1e6))


def test_a_table_needs_at_least_two_points():
    with pytest.raises(Exception):
        IMP.bff.LinkerWeighting(0.0, 10.0, [1.0])
    with pytest.raises(Exception):
        IMP.bff.LinkerWeighting(10.0, 0.0, [1.0, 1.0])


def test_the_shipped_table_loads():
    w = IMP.bff.linker_weighting(20.5)
    assert not w.get_is_uniform()
    assert w.get_x_min() > 0.0 and w.get_x_max() > w.get_x_min()
    assert len(w.get_weights()) > 100
    assert max(w.get_weights()) == pytest.approx(1.0)


def test_a_linker_longer_than_the_table_gets_no_weighting():
    assert IMP.bff.linker_weighting(1e6).get_is_uniform()


def test_chain_weighting_is_read_from_the_fps_json_field():
    """The field was declared in the schema, documented, and never read: a file
    that asked for it got the unweighted volume and no complaint."""
    IMP.set_log_level(IMP.SILENT)
    m = IMP.Model()
    h = IMP.atom.read_pdb(IMP.bff.get_example_path("structure/T4L/3GUN.pdb"), m,
                          IMP.atom.NonWaterNonHydrogenPDBSelector())
    src = IMP.atom.Selection(h, chain_id="A", residue_index=132,
                             atom_type=IMP.atom.AtomType("CB")
                             ).get_selected_particles()[0]

    def build(flag):
        av = IMP.bff.AV.setup_particle(m, IMP.Particle(m), src)
        av.set_av_parameter(
            '{"linker_length": 20.5, "linker_width": 1.5, "radius1": 3.5,'
            ' "simulation_grid_resolution": 1.5, "chain_weighting": %s}'
            % ("true" if flag else "false"))
        av.resample()
        return av

    off, on = build(False), build(True)
    assert off.get_chain_weighting() is False
    assert on.get_chain_weighting() is True
    d = IMP.algebra.get_distance(IMP.core.XYZ(off).get_coordinates(),
                                 IMP.core.XYZ(on).get_coordinates())
    assert d > 1.0, "the flag changed nothing, so it is still not being read"
    # and it can be turned back off
    on.set_chain_weighting(False)
    assert on.get_chain_weighting() is False


# ---------------------------------------------------------------------------
# the conventions that are not a single distance, and the ones that are exact
# ---------------------------------------------------------------------------

def test_the_distribution_is_available_as_a_distribution(clouds):
    """`pRDA` is refused as a scalar and answered as a histogram: `p(R_DA)`,
    which is what a fluorescence decay is a function of."""
    a, b = clouds
    axis = np.arange(0.0, 90.0, 2.5)
    hist = np.asarray(IMP.bff.cloud_distance_distribution(a.ravel(), b.ravel(),
                                                          axis))
    assert hist.shape == (len(axis) - 1,)
    assert hist.sum() == pytest.approx(1.0)
    assert (hist >= 0.0).all()

    # its mean is the mean distance the scalar reduction reports
    centres = 0.5 * (axis[1:] + axis[:-1])
    assert float((hist * centres).sum()) == pytest.approx(
        IMP.bff.cloud_model_distance(a.ravel(), b.ravel(), MEAN), abs=1.5)

    raw = np.asarray(IMP.bff.cloud_distance_distribution(
        a.ravel(), b.ravel(), axis, 50000, 0, False))
    assert raw.sum() > 1.0, "unnormalised should carry the weight"
    np.testing.assert_allclose(raw / raw.sum(), hist, atol=1e-12)


def test_the_distribution_rejects_an_axis_that_is_not_an_axis(clouds):
    a, b = clouds
    for bad in ([5.0], [10.0, 5.0, 20.0], []):
        with pytest.raises(Exception):
            IMP.bff.cloud_distance_distribution(a.ravel(), b.ravel(), bad)


def test_a_distribution_is_not_returned_as_a_distance(clouds):
    """`pRDA` reduced to a number would be R_E under another name -- a value
    that looks like an answer and is a different quantity."""
    a, b = clouds
    PR = IMP.bff.PROBE_PAIR_DISTANCE_DISTRIBUTION
    with pytest.raises(Exception) as e:
        IMP.bff.cloud_model_distance(a.ravel(), b.ravel(), PR, 52.0)
    assert "distribution" in str(e.value)
    with pytest.raises(Exception):
        IMP.bff.rmp_from_model_distance(a.ravel(), b.ravel(), 40.0, PR, 52.0)


def test_a_bound_is_not_inverted_from_a_sample(clouds):
    """The inversion samples pairs once and translates them. A *minimum* over
    that sample is biased high, and a restraint built on a biased bound is
    worse than none, so it refuses rather than answering."""
    a, b = clouds
    with pytest.raises(Exception) as e:
        IMP.bff.rmp_from_model_distance(a.ravel(), b.ravel(), 15.0, MIN, 52.0)
    assert "biased" in str(e.value)


@pytest.fixture(scope="module")
def two_avs():
    IMP.set_log_level(IMP.SILENT)
    m = IMP.Model()
    h = IMP.atom.read_pdb(IMP.bff.get_example_path("structure/T4L/3GUN.pdb"), m,
                          IMP.atom.NonWaterNonHydrogenPDBSelector())

    def build(residue, chain_weighting=False):
        src = IMP.atom.Selection(h, chain_id="A", residue_index=residue,
                                 atom_type=IMP.atom.AtomType("CB")
                                 ).get_selected_particles()[0]
        av = IMP.bff.AV.setup_particle(m, IMP.Particle(m), src)
        av.set_av_parameter(
            '{"linker_length": 20.5, "linker_width": 1.5, "radius1": 3.5,'
            ' "simulation_grid_resolution": 1.5, "chain_weighting": %s}'
            % ("true" if chain_weighting else "false"))
        av.resample()
        return av

    return m, h, build


def test_rmin_reaches_the_av_level_distances(two_avs):
    """The type was added to the fps.json vocabulary, so a file may declare it.
    Both AV-level distance functions have to know it, or such a file silently
    gets the mean distance instead."""
    m, h, build = two_avs
    a, b = build(132), build(44)
    exact = IMP.bff.av_distance(a, b, 52.0, MIN)
    assert 0.0 < exact < IMP.bff.av_distance(a, b, 52.0, MP)
    # both routes are exact here, so the quadrature must not approximate it
    assert IMP.bff.av_distance_quadrature(a, b, 52.0, MIN) == \
        pytest.approx(exact, abs=1e-9)
    # and it is not silently the mean distance
    assert exact != pytest.approx(IMP.bff.av_distance(a, b, 52.0, MEAN))


def test_the_decorator_writes_the_same_formats_as_the_value(two_avs, tmp_path):
    """`write_av` has two overloads -- the AV decorator and the value
    `compute_av` returns -- and the decorator's one also does the grid formats
    IMP.em writes."""
    m, h, build = two_avs
    av = build(132)
    for ext in ("xyz", "pqr", "dx", "mrc"):
        path = tmp_path / f"decorated.{ext}"
        IMP.bff.write_av(av, str(path))
        assert path.is_file() and path.stat().st_size > 0
    with pytest.raises(Exception):
        IMP.bff.write_av(av, str(tmp_path / "decorated.bogus"))


def test_chain_weighting_reaches_the_exported_grid(two_avs, tmp_path):
    """The weighting is folded into the map's density once per raster, so
    everything downstream sees it -- including a grid written for a viewer,
    which goes out through a different path than the point cloud."""
    m, h, build = two_avs
    plain = tmp_path / "plain.dx"
    weighted = tmp_path / "weighted.dx"
    IMP.bff.write_av(build(132), str(plain))
    IMP.bff.write_av(build(132, chain_weighting=True), str(weighted))

    def values(path):
        text = path.read_text()
        body = text.split("data follows\n", 1)[1].split("attribute")[0]
        return np.array([float(x) for x in body.split()])

    a, b = values(plain), values(weighted)
    assert a.shape == b.shape
    assert not np.allclose(a, b), "the exported grid ignored the weighting"
    assert b.max() < a.max(), "weighting should scale voxels down, not up"
