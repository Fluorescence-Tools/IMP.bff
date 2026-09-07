"""The probe library as a `.mmfdb.pto`, and the units it hands out.

The bundled library is a CIF with two categories nobody could look up:
`_bff_probe` and `_bff_probe_spectrum`, invented because no dictionary in the stack
declared an item for a quantum yield, an extinction coefficient or a
fluorescence spectrum. Dictionary 1.8 declares them -- and declares them onto
the `probes`, `optical_properties` and `spectra` tables this stack's databases
already have, so a file and a database row now say the same thing in the same
words.

What is pinned here:

* the container round-trips **exactly** -- the same numbers come back, and a
  property the library does not carry stays absent rather than becoming zero;
* a Forster radius derived through the container equals the one derived from
  the CIF, so the conversion added no error;
* **what a probe is**, which survives the round trip. A spin label that came
  back as a dye would claim a Forster radius it cannot have;
* **the unit**. R0 is Angstrom everywhere in this package, and since
  2026-08-25 it is Angstrom at the source: the spectra stay nanometres and
  `forster_radius` converts once, so `probe_pto_forster_radius` no longer has a
  factor of ten to apply. The magnitude is still pinned, because getting the
  unit wrong scores every pair against an R0 ten times too small and that does
  not look like an error -- it looks like a protein where no pair is worth
  measuring.
"""

import json
import math
import os

import pytest

import IMP.bff as bff

from test_labelizer_pto import _objects  # the from-scratch EBML walker

HERE = os.path.dirname(os.path.abspath(__file__))

#: Where the dictionary that declares these items lives.
_DICT = os.path.join(HERE, "..", "..", "..", "mmfdb", "src", "mmfdb", "data",
                     "mmfdb_flr_ext.dic")


def _item_is_declared(item):
    """Does `mmfdb_flr_ext.dic` carry a save frame for this data item?

    Cheap on purpose: a save frame is `save_<item>` at the start of a line, so
    the whole check is a substring test. Skips rather than fails when ../mmfdb
    is not checked out -- the vocabulary cannot be verified from this repo
    alone, and a missing sibling is not a defect in this one.
    """
    if not os.path.exists(_DICT):
        pytest.skip("../mmfdb is not checked out; the vocabulary "
                    "cannot be verified")
    with open(_DICT) as handle:
        return ("\nsave_%s\n" % item) in handle.read()


@pytest.fixture(scope="module")
def container(tmp_path_factory):
    path = str(tmp_path_factory.mktemp("dyes") / "dyes.mmfdb.pto")
    n = bff.probe_library_to_pto(path)
    return path, n


def test_every_dye_in_the_library_is_written(container):
    path, n = container
    assert n == len(dict(bff.read_probe_library()))
    assert n > 30, "the bundled library should not have shrunk"


def test_the_scalars_and_spectra_round_trip_exactly(container):
    path, _ = container
    source = dict(bff.read_probe_library())
    back = dict(bff.probe_read_pto(path))
    assert set(back) == set(source)

    for name, dye in source.items():
        got = back[name]
        for field in ("quantum_yield", "extinction_coefficient", "lifetime",
                      "radius", "hydrodynamic_radius"):
            a, b = getattr(dye, field), getattr(got, field)
            # NaN is how the library spells "not carried", and it has to stay
            # NaN rather than arrive as zero.
            assert math.isnan(a) == math.isnan(b), "%s.%s" % (name, field)
            if not math.isnan(a):
                assert a == b, "%s.%s" % (name, field)
        assert dye.spectrum.size() == got.spectrum.size(), name


def test_a_property_the_library_lacks_is_absent_not_zero(container):
    """The distinction the whole profile exists to keep.

    A quantum yield of 0 is a real and unusual claim about a dye -- the
    bundled library carries several -- so an unmeasured property cannot be
    written as 0 and must simply not appear.
    """
    path, _ = container
    rows = json.loads(_objects(path)["probes.json"][2])["rows"]
    source = dict(bff.read_probe_library())

    saw_absent = saw_real_zero = False
    for row in rows:
        name = row["_flr_probe_list.chromophore_name"]
        present = {p["_mmfdb_optical_property.property_name"]
                   for p in row["optical_properties"]}
        dye = source[name]
        if math.isnan(dye.lifetime):
            assert "fluorescence_lifetime" not in present, name
            saw_absent = True
        if not math.isnan(dye.quantum_yield) and dye.quantum_yield == 0.0:
            assert "quantum_yield" in present, name
            saw_real_zero = True
    assert saw_absent, "no absent property in the library to check"
    assert saw_real_zero, "the library used to carry a genuine zero quantum yield"


def test_the_forster_radius_is_unchanged_by_the_round_trip(container):
    """The container must not cost accuracy: R0 goes as the sixth root of the
    overlap, so a truncated spectrum would show up here."""
    path, _ = container
    for donor, acceptor in [("AlexaFluor488", "AlexaFluor647"),
                            ("AlexaFluor488", "AlexaFluor594"),
                            ("LumiprobeCy3", "LumiprobeCy5")]:
        from_cif = bff.forster_radius(bff.get_probe(donor),
                                      bff.get_probe(acceptor))
        from_pto = bff.probe_pto_forster_radius(path, donor, acceptor)
        assert from_pto == pytest.approx(from_cif, abs=1e-9)


def test_the_container_hands_out_angstrom(container):
    """Both routes are Angstrom, and they agree.

    This used to pin a factor of ten: `forster_radius` returned nanometres and
    `probe_pto_forster_radius` was where the package converted. Since 2026-08-25
    R0 is Angstrom at the source -- the spectra stay nanometres, the
    conversion happens once inside `forster_radius` -- so what is left to
    check is that reading the dyes out of a container gives the same R0 as
    reading them out of the CIF, which is this function's actual claim.
    """
    path, _ = container
    direct = bff.forster_radius(bff.get_probe("AlexaFluor488"),
                                bff.get_probe("AlexaFluor647"))
    from_container = bff.probe_pto_forster_radius(path, "AlexaFluor488",
                                                "AlexaFluor647")
    assert from_container == pytest.approx(direct)
    # Angstrom, not nanometres: a dye-pair R0 is tens, not units. This is the
    # assertion that would catch a stray factor of ten reappearing at either
    # end -- the consumers (LlFretOptions, av_distance) all default to 52.
    assert 30.0 < from_container < 90.0


def test_a_missing_dye_names_the_container_it_is_missing_from(container):
    path, _ = container
    with pytest.raises(Exception) as excinfo:
        bff.probe_pto_forster_radius(path, "NotADye", "AlexaFluor647")
    assert "NotADye" in str(excinfo.value)


# ---------------------------------------------------------------------------
# Conformance
# ---------------------------------------------------------------------------

def test_the_framing_parses_without_any_of_our_code(container):
    path, _ = container
    objects = _objects(path)
    assert set(objects) == {"README", "probes.json", "spectra.json",
                            "probe_model.json"}
    kinds = {k for k, _, _ in objects.values()}
    assert kinds == {"readme", "mfdb.probes", "mfdb.spectra", "mfdb.model"}
    # `label.*` is the label-score namespace and `drot.*` the rotamer one;
    # generic probe data belongs to neither.
    assert not any(k.startswith(("label.", "drot.", "rot.")) for k in kinds)


def test_the_two_tables_declare_what_their_rows_are(container):
    """A probe row and a spectrum row count different things, and the profile
    forbids joining them by position."""
    path, _ = container
    model = json.loads(_objects(path)["probe_model.json"][2])
    artifacts = model["_artifacts"]
    assert artifacts["probes.json"]["_mmfdb_artifact.row_grain"] == "species"
    assert artifacts["spectra.json"]["_mmfdb_artifact.row_grain"] == "spectrum"

    edges = {(e["_mmfdb_edge.source_node_id"], e["_mmfdb_edge.target_node_id"]):
             e for e in model["_edges"]}
    edge = edges[("probes.json", "spectra.json")]
    assert edge["_mmfdb_edge.relationship_type"] == "maps_rows_of"
    # The join is named, not positional.
    assert edge["_mmfdb_edge.source_row_column"] == \
        "_flr_probe_list.chromophore_name"
    assert edge["_mmfdb_edge.target_row_column"] == "_mmfdb_spectrum.probe_id"


def test_every_column_name_is_a_dictionary_item(container):
    """The reason the container exists rather than another CIF of our own."""
    path, _ = container
    rows = json.loads(_objects(path)["probes.json"][2])["rows"]
    for row in rows[:5]:
        for key in row:
            if key == "optical_properties":
                continue
            # Looked up in the dictionary, not matched against a prefix list.
            # A prefix check passes for `_mmfdb_probe.dipole_atoms` whether or
            # not that item was ever declared, which is the whole failure this
            # container exists to avoid -- a name that *looks* like a
            # dictionary item and resolves to nothing.
            #
            # `_mmfdb_probe.*` is the local extension: what a probe is, who
            # sells it, its dipole and its formal charges have no flrCIF item,
            # and squatting on `_flr_probe_list.` for them would collide with a
            # future IHM-FLR one.
            assert _item_is_declared(key), (
                "%s is not declared in mmfdb_flr_ext.dic" % key)
        assert bff.mfdb_is_term("_mmfdb_probe.probe_type",
                                row["_mmfdb_probe.probe_type"])
        for prop in row["optical_properties"]:
            for key in prop:
                assert key.startswith("_mmfdb_optical_property."), key
            assert bff.mfdb_is_term("_mmfdb_optical_property.unit",
                                    prop.get("_mmfdb_optical_property.unit",
                                             "dimensionless"))

    spectra = json.loads(_objects(path)["spectra.json"][2])["rows"]
    for row in spectra[:5]:
        for key in row:
            assert key.startswith("_mmfdb_spectrum."), key
        assert bff.mfdb_is_term("_mmfdb_spectrum.spectrum_type",
                                row["_mmfdb_spectrum.spectrum_type"])


def test_a_spectrum_is_one_row_carrying_arrays(container):
    """`wavelengths` and `intensity_values` are declared as arrays, and they
    have to be arrays.

    They were first written a point per row, with a scalar sitting in a field
    named `wavelengths` -- wrong against the dictionary that declares them,
    wrong against the `spectra` table that stores them as blobs, and ten times
    the size: 9.4 MB for a library whose CIF source is 0.9 MB, almost all of it
    the same four JSON keys repeated forty-five thousand times.
    """
    path, _ = container
    rows = json.loads(_objects(path)["spectra.json"][2])["rows"]
    source = dict(bff.read_probe_library())

    # Two rows per dye that has a spectrum, not two per point.
    with_spectrum = [d for d in source.values() if d.spectrum.size() > 0]
    assert len(rows) == 2 * len(with_spectrum)

    for row in rows:
        wl = row["_mmfdb_spectrum.wavelengths"]
        iv = row["_mmfdb_spectrum.intensity_values"]
        assert isinstance(wl, list) and isinstance(iv, list)
        assert len(wl) == len(iv) > 100, "a spectrum is a curve, not a point"
        assert wl == sorted(wl), "the wavelength axis must be ordered"
        assert row["_mmfdb_spectrum.wavelength_unit"] == "nm"
        assert row["_mmfdb_spectrum.intensity_unit"] == "normalized"

    # And the file is the size that shape implies, not the other one.
    assert os.path.getsize(path) < 3 * 1024 * 1024, (
        "%.1f MB -- the point-per-row shape is back"
        % (os.path.getsize(path) / 1048576.0))


def test_the_container_says_what_it_conforms_to(container):
    path, _ = container
    model = json.loads(_objects(path)["probe_model.json"][2])
    assert model["_container"]["_mmfdb_container.profile"] == "PTO.MFDB"
    assert model["_container"]["_mmfdb_container.dictionary_version"] == \
        bff.MFDB_DICTIONARY_VERSION
    assert model["_operation"]["_mmfdb_operation.algorithm"] == "probe_library"
    # Attribution travels with the data.
    assert "_mmfdb_artifact.source" in model["_attribution"]


def test_the_readme_explains_the_normalisation(container):
    """"normalized" means the peak is 1, not that the area is -- and a spectral
    overlap integral that assumed unit area would be silently wrong."""
    path, _ = container
    text = _objects(path)["README"][2].decode("ascii")
    assert "RFC 8794" in text
    assert "peak is 1" in text and "unit area" in text


# ---------------------------------------------------------------------------
# One unit, everywhere
# ---------------------------------------------------------------------------

#: Every route in the package that hands back a Forster radius.
def _all_r0_routes(container, donor, acceptor):
    return {
        "forster_radius": bff.forster_radius(bff.get_probe(donor),
                                             bff.get_probe(acceptor)),
        "forster_radius_from_spectra": bff.forster_radius_from_spectra(
            donor, acceptor, 2.0 / 3.0),
        "probe_pto_forster_radius": bff.probe_pto_forster_radius(
            container, donor, acceptor),
    }


@pytest.mark.parametrize("donor,acceptor", [
    ("Alexa488", "Alexa647"),          # the spelling everyone types
    ("AlexaFluor488", "AlexaFluor647"),  # the library key
    ("AlexaFluor 488", "AlexaFluor 647"),  # the vendor's, with a space
    ("Atto532", "Atto643"),            # imported dyes
    ("Cy3", "Cy5"),                    # a vendor-prefixed key
])
def test_every_route_gives_the_same_angstrom(container, donor, acceptor):
    """**Every Forster radius in this package is Angstrom.**

    There are three ways to ask for one, and they must not disagree -- not by a
    factor of ten, and not by a rounding. They did once: `forster_radius`
    returned nanometres while everything consuming an R0 was Angstrom, and
    separately `forster_radius_from_spectra` matched names on whitespace alone
    so `Alexa488` worked through one door and not the other.

    The band check is the one that catches a returning factor of ten. A dye
    pair R0 is tens of Angstrom; a value near 5 is nanometres.
    """
    path, _ = container
    routes = _all_r0_routes(path, donor, acceptor)
    values = list(routes.values())
    for name, value in routes.items():
        assert value == pytest.approx(values[0], abs=1e-9), name
        assert 20.0 < value < 120.0, "%s returned %.3f -- nanometres?" % (
            name, value)


def test_the_consumers_default_to_the_same_unit():
    """A default is a claim about a unit too. If these drifted to nanometres
    the port would score every pair against an R0 ten times too small, which
    reads as a protein where no pair is worth measuring rather than as a bug."""
    assert 20.0 < bff.LlFretOptions().forster_radius < 120.0
    # And the efficiency curve crosses one half at R0, whatever R0 is.
    r0 = bff.LlFretOptions().forster_radius
    assert bff.fret_efficiency(r0, r0) == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# What a probe is, as opposed to how it is modelled
# ---------------------------------------------------------------------------

def test_the_bundled_library_says_every_row_is_a_dye():
    """`_bff_probe.probe_type` is read, not assumed.

    The column was added when `ProbeLibrary` became `ProbeLibrary`; before that
    the second column of the CIF was called `probe_type` and held the *vendor*,
    which is a different question and was never read at all.
    """
    library = dict(bff.read_probe_library())
    assert len(library) >= 40
    for name, probe in library.items():
        assert probe.probe_type == bff.PROBE_DYE, name
        assert bff.probe_is_fluorescent(probe.probe_type)
    # The vendor is kept now rather than dropped on the floor.
    assert library["ATTO532"].vendor == "ATTO"
    assert library["LumiprobeCy3"].vendor == "Lumiprobe"


def test_the_type_names_round_trip_and_an_unknown_one_is_unspecified():
    for kind in (bff.PROBE_DYE, bff.PROBE_FLUORESCENT_PROTEIN,
                 bff.PROBE_SPIN_LABEL, bff.PROBE_UNSPECIFIED):
        name = bff.probe_type_to_string(kind)
        assert bff.probe_type_from_string(name) == kind
        # Every spelling is a dictionary term, because the container writes it.
        assert bff.mfdb_is_term("_mmfdb_probe.probe_type", name)
    # A library that predates the column, or one with a typo, is not an error:
    # a probe whose kind was not stated is still a usable species.
    for junk in ("", "Dye", "nitroxide", "fluorophore"):
        assert bff.probe_type_from_string(junk) == bff.PROBE_UNSPECIFIED


def test_a_spin_label_is_refused_a_forster_radius_by_name():
    """A nitroxide absorbs and emits nothing, so R0 does not exist for it.

    The point of the message is that it names the *reason*. Before the probe
    type existed the only thing the code could say was "both dyes need a
    spectrum", which reads as a missing file rather than a category error.
    """
    donor = bff.get_probe("AlexaFluor488")
    label = bff.get_probe("AlexaFluor647")
    label.name = "MTSSL"
    label.probe_type = bff.PROBE_SPIN_LABEL

    with pytest.raises(Exception) as excinfo:
        bff.forster_radius(donor, label)
    message = str(excinfo.value)
    assert "MTSSL" in message and "spin label" in message

    # And it is refused on either side of the pair.
    with pytest.raises(Exception):
        bff.forster_radius(label, donor)
    # While the two real dyes still work.
    assert bff.forster_radius(donor, bff.get_probe("AlexaFluor647")) > 40.0


def test_the_probe_type_survives_a_container_round_trip(tmp_path):
    """Losing it would silently turn a spin label back into something that
    claims a Forster radius."""
    probes = {}
    probes["AlexaFluor488"] = bff.get_probe("AlexaFluor488")
    label = bff.get_probe("AlexaFluor647")
    label.name = "MTSSL"
    label.probe_type = bff.PROBE_SPIN_LABEL
    label.vendor = "Toronto Research Chemicals"
    probes["MTSSL"] = label

    out = str(tmp_path / "mixed.mmfdb.pto")
    bff.probe_write_pto(out, probes, "a dye and a spin label")
    back = dict(bff.probe_read_pto(out))

    assert back["AlexaFluor488"].probe_type == bff.PROBE_DYE
    assert back["MTSSL"].probe_type == bff.PROBE_SPIN_LABEL
    assert back["MTSSL"].vendor == "Toronto Research Chemicals"
    with pytest.raises(Exception):
        bff.forster_radius(back["AlexaFluor488"], back["MTSSL"])


def test_the_shipped_container_still_matches_the_library():
    """`data/dyes.mmfdb.pto` is derived from the CIF, so it can go stale.

    It did: it was written before `probe_type` and `vendor` existed and came
    back `unspecified` with no vendor for every one of its 40 probes, while the
    CIF beside it had both. Nothing caught it, because every other test in this
    file writes its own container into `tmp_path` and reads that back -- which
    proves the round trip and says nothing about the file that actually ships.
    """
    shipped = os.path.join(HERE, "..", "..", "data", "dyes.mmfdb.pto")
    if not os.path.exists(shipped):
        pytest.skip("the container is generated, not tracked")

    container = dict(bff.probe_read_pto(shipped))
    library = dict(bff.read_probe_library())
    assert set(container) == set(library)

    for name, probe in library.items():
        got = container[name]
        assert got.probe_type == probe.probe_type, name
        assert got.vendor == probe.vendor, name
        assert got.quantum_yield == pytest.approx(probe.quantum_yield), name
        assert got.extinction_coefficient == pytest.approx(
            probe.extinction_coefficient), name
        assert len(got.spectrum.wavelength) == len(probe.spectrum.wavelength), name

    # And the number it exists to hand out agrees with the CIF's, in Angstrom.
    direct = bff.forster_radius(library["AlexaFluor488"],
                                library["AlexaFluor647"])
    assert bff.probe_pto_forster_radius(shipped, "Alexa488", "Alexa647") == \
        pytest.approx(direct)
    assert 40.0 < direct < 70.0, "an R0 near 5 would be nanometres"


def test_no_probe_field_is_silently_dropped(tmp_path):
    """Every field of a `Probe` survives the container, or this fails.

    Written after three of them did not. `dipole_atoms`, `positive_atoms` and
    `negative_atoms` had no flrCIF item, so the writer simply left them out and
    the reader had nothing to look for -- and because a `Probe` with no dipole
    is *legal* (it means isotropic), the loss produced a valid-looking probe
    rather than an error.

    That is the damaging shape: an anisotropic probe comes back isotropic,
    kappa^2 falls back to 2/3, and the R0 that follows is wrong by a factor
    nothing in the output reveals.

    The test enumerates the fields rather than checking the three, so a field
    added later is covered without anyone remembering to come back here.
    """
    probe = bff.get_probe("AlexaFluor488")
    # Every field set to something distinguishable from its default, so a
    # dropped one cannot coincide with what a fresh Probe already holds.
    probe.vendor = "TestVendor"
    probe.probe_type = bff.PROBE_FLUORESCENT_PROTEIN
    probe.lifetime = 4.1
    probe.radius = 3.5
    probe.hydrodynamic_radius = 7.0
    probe.dipole_atoms = ["C1", "C2"]
    probe.chromophore_center_atom = "C7"
    probe.reactive_probe_name = "Alexa488-maleimide"
    probe.probe_origin = "extrinsic"
    probe.probe_link_type = "covalent"
    probe.positive_atoms = ["N1"]
    probe.negative_atoms = ["O1", "O2"]

    out = str(tmp_path / "every_field.mmfdb.pto")
    bff.probe_write_pto(out, {"AlexaFluor488": probe}, "every field set")
    back = dict(bff.probe_read_pto(out))["AlexaFluor488"]

    dropped = []
    for field in sorted(n for n in dir(probe)
                        if not n.startswith("_")
                        and not callable(getattr(probe, n))):
        if field in ("spectrum", "this", "thisown"):
            continue          # `spectrum` is checked elementwise by the
                              # round-trip test; `this`/`thisown` are SWIG's
        before, after = getattr(probe, field), getattr(back, field)
        if hasattr(before, "__len__") and not isinstance(before, str):
            before, after = list(before), list(after)
        if before != after:
            dropped.append("%s: %r -> %r" % (field, before, after))
    assert not dropped, "the container dropped:\n  " + "\n  ".join(dropped)
