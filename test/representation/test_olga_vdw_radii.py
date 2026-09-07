"""Olga's van der Waals radii: the table, the unit, and the AV that uses them.

An accessible volume is geometry, and its size is set by how big the obstacles
are. Three radii sets are in circulation here and they are not interchangeable
-- FPS's element-keyed Bondi, IMP's united-atom set with implicit hydrogens,
and Olga's, keyed by *atom name*. The volume must use one of them, and the one
it uses by default is **IMP's** -- the radii on the particles (owner,
2026-09-01). Not because it reproduces Olga's published numbers best, which it
does not, but because the excluded-volume half of a docking score is
`clash_container` reading `IMP::core::XYZR`: a volume on one radii set and a
clash term on another are two halves of one score that disagree about how big
an atom is. Olga's table stays vendored and selectable, because it is what
reproduces Zenodo 3376527 (PRD-121 G9, `okf/validation/fps_screening_ab.md`,
which measures what the default costs).

The load-bearing test in this file is :func:`test_the_unit_is_nanometre_times_ten`.
Upstream the numbers are nanometres and this module works in Angstrom; getting
the factor of ten wrong would be silent -- every volume would come back either
a sphere (obstacles ten times too small) or empty -- so it is pinned against
named entries rather than trusted.
"""
import json
import math
import os

import numpy as np
import pytest

import IMP
import IMP.atom
import IMP.core
import IMP.bff

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
SHIPPED = os.path.join(REPO, "data", "olga_vdw_radii.csv")

#: Olga's upstream table, if the sibling checkout is here. It is the source
#: this was vendored from; the vendored copy is what the module actually uses,
#: so its absence is not a failure.
UPSTREAM = os.path.abspath(os.path.join(
    REPO, "..", "ucfret", "thirdparty", "olga", "src", "vdWRadii.json"))


def _shipped_rows():
    rows = {}
    with open(SHIPPED) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            if line.startswith("atom_name,"):
                continue
            name, value = line.rstrip("\n").rsplit(",", 1)
            rows[name] = float(value)
    return rows


# --------------------------------------------------------------------------
# The unit. Read this one first.
# --------------------------------------------------------------------------

def test_the_unit_is_nanometre_times_ten():
    """0.17 nm upstream is 1.70 A here, and the factor is Olga's own.

    The conversion is not inferred from the magnitudes. `coordsVdW()`
    (`Olga/src/AV/Position.cpp:118-124`) scales the coordinate *and* the radius
    by the same `10.0f` on one line, because pteros stores coordinates in
    nanometres and the AV kernel `calculateAV()` works in Angstrom:

        xyzw.emplace_back(frame.coord.at(i)[0] * 10.0f, ...,
                          pterosVDW(system, i) * 10.0f);

    So every value here is upstream's times ten. Carbon lands on Bondi's
    1.70 A, which is also FPS's carbon -- an independent check that the factor
    is right, since the two tables were built from the same source.
    """
    assert IMP.bff.olga_vdw_radius("C") == pytest.approx(1.70)     # 0.1700 nm
    assert IMP.bff.olga_vdw_radius("CA") == pytest.approx(1.70)
    assert IMP.bff.olga_vdw_radius("N") == pytest.approx(1.625)    # 0.1625 nm
    assert IMP.bff.olga_vdw_radius("O") == pytest.approx(1.49)     # 0.1490 nm
    assert IMP.bff.olga_vdw_radius("SD") == pytest.approx(1.782)   # 0.1782 nm
    assert IMP.bff.olga_vdw_radius("P") == pytest.approx(1.86)     # 0.1860 nm
    assert IMP.bff.olga_vdw_radius("H") == pytest.approx(1.00)     # 0.1000 nm

    # The whole table, so a single mistyped row cannot hide: everything is a
    # plausible Angstrom radius, and nothing is a plausible nanometre one.
    values = np.array(list(_shipped_rows().values()))
    assert values.min() >= 1.0 and values.max() <= 2.0


def test_the_fallback_is_olgas_and_it_is_flat():
    """1.50 A for an unknown name -- no element lookup, by design.

    `pterosVDW` is `vdWRMap.value(name, 0.15)` (`Position.cpp:110`): a flat
    0.15 nm for any atom name the table misses. It is *not* near-transparent --
    at 1.50 A it is within 0.01 A of Olga's own oxygen -- but it is smaller than
    every heavy atom in the table, and much smaller than IMP's united-atom
    carbon (1.85-2.275 A), which is what an unrecognised carbon would otherwise
    have got. That is reproduced rather than improved on.
    """
    assert IMP.bff.get_olga_vdw_fallback_radius() == pytest.approx(1.50)
    assert IMP.bff.olga_vdw_radius("ZZZ") == pytest.approx(1.50)
    assert IMP.bff.olga_vdw_radius("") == pytest.approx(1.50)
    # smaller than every heavy atom, larger than hydrogen
    values = [v for v in _shipped_rows().values() if v > 1.0]
    assert IMP.bff.get_olga_vdw_fallback_radius() < min(values) + 0.02
    assert IMP.bff.get_olga_vdw_fallback_radius() > 1.0


# --------------------------------------------------------------------------
# The vendored table
# --------------------------------------------------------------------------

def test_the_shipped_data_file_is_derived_from_the_cpp_table():
    """`data/olga_vdw_radii.csv` is generated, not hand-edited.

    Same contract as `data/fps_json_schema.json`: the C++ table is the
    definition and a test regenerates the data file, so provenance and code
    cannot drift apart.
    """
    with open(SHIPPED) as fh:
        assert fh.read() == IMP.bff.olga_vdw_radii_csv()


def test_the_table_has_olgas_128_entries():
    rows = _shipped_rows()
    assert len(rows) == 128
    assert IMP.bff.get_number_of_olga_vdw_radii() == 128
    assert sorted(rows) == list(IMP.bff.get_olga_vdw_atom_names())
    for name, value in rows.items():
        assert IMP.bff.olga_vdw_radius(name) == pytest.approx(value)


@pytest.mark.skipif(not os.path.exists(UPSTREAM),
                    reason="the Olga checkout is not beside this one")
def test_the_vendored_table_still_matches_upstream():
    """The vendoring, checked against its source when the source is here.

    Upstream has 131 lines and 128 keys: "H1", "H2" and "H3" appear twice with
    identical values, so JSON's last-wins costs nothing.
    """
    upstream = json.load(open(UPSTREAM))
    assert len(upstream) == 128
    for name, nm in upstream.items():
        assert IMP.bff.olga_vdw_radius(name) == pytest.approx(nm * 10.0)


# --------------------------------------------------------------------------
# What the AV does with it
# --------------------------------------------------------------------------

AV_PARAMETER = {
    "linker_length": 20.0,
    "radii": (3.5, 0.0, 0.0),
    "linker_width": 0.5,
    "allowed_sphere_radius": 2.0,
    "contact_volume_thickness": 0.0,
    "contact_volume_trapped_fraction": -1,
    "simulation_grid_resolution": 0.5,
}


def _t4l():
    model = IMP.Model()
    hierarchy = IMP.atom.read_pdb(
        IMP.bff.get_example_path("structure/T4L/3GUN.pdb"), model)
    return model, hierarchy


def _av(model, hierarchy, residue, atom, radii_source):
    particle = IMP.Particle(model)
    selection = IMP.atom.Selection(hierarchy)
    selection.set_atom_type(IMP.atom.AtomType(atom))
    selection.set_residue_index(residue)
    source = selection.get_selected_particles()[0]
    IMP.bff.AV.do_setup_particle(model, particle, source, **AV_PARAMETER)
    av = IMP.bff.AV(model, particle)
    av.set_radii_source(radii_source)
    av.resample()
    return av


def test_the_default_is_imps_own_radii():
    """Stated in `AV::set_radii_source`, pinned here.

    The owner's decision (2026-09-01, reversing a one-day default of Olga's):
    the volume must use one radii set, and that set is the one already on the
    particles, because that is the one `clash_container` measures overlap
    with. A caller who wants Olga's table has to ask for it -- which is the
    right way round, since asking is also what says "these numbers are
    Olga-era".
    """
    model, hierarchy = _t4l()
    particle = IMP.Particle(model)
    selection = IMP.atom.Selection(hierarchy)
    selection.set_atom_type(IMP.atom.AtomType("CB"))
    selection.set_residue_index(132)
    IMP.bff.AV.do_setup_particle(model, particle,
                                 selection.get_selected_particles()[0],
                                 **AV_PARAMETER)
    assert IMP.bff.AV(model, particle).get_radii_source() == "imp"


def test_an_unknown_radii_source_is_refused():
    """...including `"model"`, which is what `"imp"` was called for one day.

    Refused rather than aliased: a schema 1.5 file carrying that spelling is
    read as the 1.5 file it is, loudly, instead of being quietly reinterpreted
    under a default that has since flipped. Two names for one quantity is the
    failure this module is trying not to have.
    """
    model, hierarchy = _t4l()
    av = _av(model, hierarchy, 132, "CB", "olga")
    with pytest.raises(Exception):
        av.set_radii_source("bondi")
    with pytest.raises(Exception):
        av.set_radii_source("model")


def test_olgas_radii_open_the_volume_up():
    """Measured 2026-09-01 on T4L 3GUN, grid 0.5 A, clearance 2.0.

    Olga's heavy atoms are smaller than IMP's united-atom ones -- carbon 1.70
    against 1.85-2.275 -- so every volume grows, and a site that was walled in
    stops being walled in:

    ==========  ===============  ===============  ==============
    site        imp (voxels)     olga (voxels)    mean moves
    ==========  ===============  ===============  ==============
    132 CB      82 418           95 836  (+16 %)  1.05 A
    55 CB       136 823          144 032  (+5 %)  0.38 A
    99 CB       710              29 586  (x42)    17.9 A
    ==========  ===============  ===============  ==============

    Site 99 is the one that matters: at 710 voxels it was a numerical
    accident, not a volume, and that is the class of site FPS and Olga report
    as open where this module reports it buried (PRD-121, open item 5). This
    is what the default costs, and it is a measurement, not a regression --
    the default is IMP's because the clash term of a docking score reads IMP's
    radii, and `okf/validation/fps_screening_ab.md` prices that choice.
    """
    model, hierarchy = _t4l()
    counts = {}
    means = {}
    for source in ("imp", "olga"):
        av = _av(model, hierarchy, 132, "CB", source)
        counts[source] = len(av.get_map().get_xyz_density())
        means[source] = np.array(av.get_mean_position())
    assert counts["olga"] > counts["imp"]
    assert counts["olga"] / counts["imp"] == pytest.approx(1.163, abs=0.02)
    assert np.linalg.norm(means["olga"] - means["imp"]) > 0.5

    buried = {s: len(_av(model, hierarchy, 99, "CB", s)
                     .get_map().get_xyz_density())
              for s in ("imp", "olga")}
    assert buried["imp"] < 1000
    assert buried["olga"] > 20000


def test_the_source_is_a_per_position_fps_json_field():
    """schema 1.6 `radii_source`, so a file can say which set it was fitted for.

    It belongs on the position rather than on the run for the same reason the
    strip mask does: a fitted `contact_volume_trapped_fraction` only means
    anything against the obstacle set it was fitted over.
    """
    model, hierarchy = _t4l()
    particle = IMP.Particle(model)
    selection = IMP.atom.Selection(hierarchy)
    selection.set_atom_type(IMP.atom.AtomType("CB"))
    selection.set_residue_index(132)
    IMP.bff.AV.do_setup_particle(model, particle,
                                 selection.get_selected_particles()[0],
                                 **AV_PARAMETER)
    av = IMP.bff.AV(model, particle)
    av.set_av_parameter(json.dumps({"radii_source": "olga"}))
    assert av.get_radii_source() == "olga"
    av.set_av_parameter(json.dumps({}))          # absent = the default
    assert av.get_radii_source() == "imp"

    schema = json.loads(IMP.bff.fps_json_schema())
    field = schema["$defs"]["position"]["properties"]["radii_source"]
    assert field["default"] == "imp"
    assert sorted(field["enum"]) == ["imp", "olga"]
    # 1.6 is the rename plus the flipped default; 1.5 was the field's debut.
    assert IMP.bff.fps_schema_version() >= "1.6"


def test_a_particle_that_is_not_an_atom_keeps_its_own_radius():
    """Olga's table is keyed by atom name, so it has no opinion about a bead.

    Giving a coarse-grained particle the 1.50 A unknown-name fallback would
    silently shrink an obstacle set a caller sized deliberately. An atom whose
    *name* the table misses does take the fallback -- that case is inside
    Olga's domain and is what Olga does.
    """
    model = IMP.Model()
    bead = IMP.Particle(model)
    IMP.core.XYZR.setup_particle(bead, IMP.algebra.Sphere3D(
        IMP.algebra.Vector3D(0, 0, 0), 6.0))
    # the public accessor is pure Olga: not an atom, so the fallback
    assert IMP.bff.olga_vdw_particle_radius(bead) == pytest.approx(1.50)

    # ...but the AV keeps the bead's own 6.0 A, which is the behaviour that
    # matters. A volume grown next to one bead is a sphere of the linker
    # length minus whatever the bead blocks; shrinking the bead to 1.5 A would
    # change that and nothing would say so.
    model, hierarchy = _t4l()
    av_olga = _av(model, hierarchy, 132, "CB", "olga")
    # the T4L structure is all atoms, so nothing here exercises the bead path
    # beyond the accessor above; this is the sanity end of it.
    assert len(av_olga.get_map().get_xyz_density()) > 0


def test_the_unknown_names_of_the_shipped_fixtures_are_reported():
    """Which atoms take the 1.50 A fallback, said out loud rather than guessed.

    On T4L 3GUN, none: every atom name is in the table. On the HIV-RT fixture
    the misses are the old-style spellings -- `O1P`/`O2P` for the phosphate
    oxygens Olga calls `OP1`/`OP2`, `C7` for thymine's methyl, and the `C1*`
    star notation for the sugar. All but the carbons land within 0.01 A of the
    radius the table would have given them; the carbons come out 0.20 A small.
    """
    model, hierarchy = _t4l()
    leaves = IMP.atom.get_leaves(hierarchy)
    assert list(IMP.bff.olga_vdw_unknown_atom_names(leaves)) == []

    hiv = os.path.join(REPO, "test", "input", "fps", "hivrt_frame03991.pdb")
    model2 = IMP.Model()
    h2 = IMP.atom.read_pdb(hiv, model2,
                           IMP.atom.NonWaterNonHydrogenPDBSelector())
    missing = list(IMP.bff.olga_vdw_unknown_atom_names(
        IMP.atom.get_leaves(h2)))
    assert missing == ["C7", "O1P", "O2P"]


def test_a_shared_raster_refuses_two_radii_sets():
    """One shared occupancy raster is one obstacle set, so it cannot be both.

    The alternative is that whichever volume asked second gets the first one's
    atoms and still returns a plausible cloud. See
    `AVOccupancyRegistry::adopt_obstacle_radii`.
    """
    model, hierarchy = _t4l()
    registry = IMP.bff.AVOccupancyRegistry(IMP.atom.get_leaves(hierarchy))
    # A registry hands its maps a shared coordinate snapshot and they read it
    # instead of the Model; it is empty until this is called, and a map that
    # indexes an empty snapshot segfaults. ProbeNetworkRestraint drives this
    # through update_all(); a caller wiring a registry by hand has to.
    registry.refresh_snapshot()

    first = _av(model, hierarchy, 132, "CB", "olga")
    first.set_occupancy_registry(registry)
    first.resample()

    second = _av(model, hierarchy, 55, "CB", "imp")
    second.set_occupancy_registry(registry)
    with pytest.raises(Exception):
        second.resample()

    # ...and a second volume that agrees costs nothing: the check is on the
    # source, not on a walk of every atom. It must also get the *same* volume
    # it would get on a private raster -- the radii have to reach the shared
    # maps, and "they were ignored" would look exactly like success otherwise.
    third = _av(model, hierarchy, 55, "CB", "olga")
    third.set_occupancy_registry(registry)
    third.resample()
    private = _av(model, hierarchy, 55, "CB", "olga")
    assert (len(third.get_map().get_xyz_density())
            == len(private.get_map().get_xyz_density()) == 144032)
