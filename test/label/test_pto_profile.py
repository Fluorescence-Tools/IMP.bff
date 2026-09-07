"""The compiled-in vocabulary against the dictionary it was copied from.

`PtoProfile.h` carries the MMFDB controlled vocabularies as C++ constants
rather than parsing `mmfdb_flr_ext.dic` at run time, because that dictionary
lives in a sibling repository which is not a dependency of this one and writing
a file must not require a checkout to be present. That is the same choice
`ProbeLibrary.h` and `FPS.h` already made for flrCIF item names.

The choice comes with an obligation, and this file is it: a term is **checked,
never recalled**. Every value the C++ claims is re-read from the dictionary
here, so a term that is renamed or removed upstream fails the suite instead of
being written into containers for another year.

Skipped, not failed, when the sibling checkout is absent -- a contributor
without `../mmfdb` can still run the suite, they just cannot verify this.
"""

import os
import re

import pytest

import IMP.bff as bff

MMFDB_DATA = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "..", "mmfdb", "src", "mmfdb", "data",
)
#: The extension dictionaries this package's terms come from. `_mmfdb_column.*`
#: and `_mmfdb_container.*` are declared in the workflow extension and the rest
#: in the fluorescence one, so both have to be searched.
DICTS = [os.path.join(MMFDB_DATA, "mmfdb_flr_ext.dic"),
         os.path.join(MMFDB_DATA, "mmfdb_workflow_ext.dic")]
#: The one that carries `_dictionary.version` for the terms this package pins.
DICT = DICTS[0]

pytestmark = pytest.mark.skipif(
    not all(os.path.exists(d) for d in DICTS),
    reason="../mmfdb is not checked out; the vocabulary cannot be verified",
)


def _enumeration(item):
    """The `_item_enumeration.value` list an item declares, from the .dic files.

    The dictionaries are CIF save-frame files: the frame for an item runs from
    `save_<item>` to the next `save_`, and its enumeration is the block after
    `_item_enumeration.value` -- one value per line, optionally followed by an
    `_item_enumeration.detail` column, in which case the value is the line's
    first token.
    """
    for path in DICTS:
        text = open(path).read()
        marker = "save_%s\n" % item
        if marker not in text:
            continue
        start = text.index(marker)
        nxt = text.find("\nsave_", start + 1)
        frame = text[start: nxt if nxt != -1 else len(text)]
        if "_item_enumeration.value" not in frame:
            continue
        body = frame.split("_item_enumeration.value", 1)[1]
        values = []
        for raw in body.splitlines():
            line = raw.strip()
            if not line:
                if values:
                    break
                continue
            if line.startswith("_item_enumeration."):
                continue
            if line.startswith(("_", "save_", ";", "loop_")):
                break
            values.append(line.split()[0])
        if values:
            return values
    return None


#: The C++ accessor -> the dictionary item it was copied from.
VOCABULARIES = {
    "_mmfdb_artifact.artifact_kind": bff.mfdb_artifact_kinds,
    "_mmfdb_artifact.data_format": bff.mfdb_data_formats,
    "_mmfdb_artifact.row_grain": bff.mfdb_row_grains,
    "_mmfdb_edge.relationship_type": bff.mfdb_relationship_types,
    "_mmfdb_column.units": bff.mfdb_units,
    "_mmfdb_operation.operation_type": bff.mfdb_operation_types,
    "_mmfdb_label_score.score_type": bff.mfdb_label_score_types,
    "_mmfdb_label_score.status": bff.mfdb_label_score_statuses,
    "_mmfdb_label_score.definition": bff.mfdb_label_score_definitions,
    "_mmfdb_optical_property.unit": bff.mfdb_optical_property_units,
    "_mmfdb_probe.probe_type": bff.mfdb_probe_types,
    "_mmfdb_spectrum.spectrum_type": bff.mfdb_spectrum_types,
    "_mmfdb_spectrum.wavelength_unit": bff.mfdb_spectrum_wavelength_units,
    "_mmfdb_spectrum.intensity_unit": bff.mfdb_spectrum_intensity_units,
}


@pytest.mark.parametrize("item", sorted(VOCABULARIES))
def test_every_compiled_term_is_still_in_the_dictionary(item):
    """A subset is fine; an invented value is not.

    The C++ deliberately lists only the values this package can legitimately
    write -- `artifact_kind` has seventy terms upstream and nothing here
    produces a photon stream -- so the test is containment, not equality.
    """
    declared = _enumeration(item)
    assert declared, "no enumeration found for %s" % item
    ours = list(VOCABULARIES[item]())
    assert ours, "%s has no compiled values" % item
    unknown = [v for v in ours if v not in declared]
    assert not unknown, (
        "%s: %r are not in mmfdb_flr_ext.dic. A term that is genuinely "
        "missing belongs in the dictionary first." % (item, unknown)
    )


def test_the_recorded_dictionary_version_matches_the_file():
    """A file records the dictionary revision its terms came from, so a later
    reader can tell a renamed term from a typo. That only works if it is
    true."""
    text = open(DICT).read()
    match = re.search(r"_dictionary\.version\s+(\S+)", text)
    assert match, "the dictionary declares no version"
    assert bff.MFDB_DICTIONARY_VERSION == match.group(1), (
        "PtoProfile.h says dictionary %s, mmfdb_flr_ext.dic says %s -- "
        "re-check the vocabularies and bump the constant"
        % (bff.MFDB_DICTIONARY_VERSION, match.group(1))
    )


def test_the_label_site_grain_is_declared_upstream():
    """This package needed a grain the dictionary did not have and added it.

    If it disappears upstream, every container this package writes becomes
    non-conformant, and that should be loud.
    """
    assert "label_site" in _enumeration("_mmfdb_artifact.row_grain")


def test_the_labelizer_score_types_cover_every_parameter():
    """Each of the reference's seven terms maps onto a dictionary value."""
    declared = _enumeration("_mmfdb_label_score.score_type")
    for tag in ["cs", "se", "ss", "ce", "tp", "cr", "me"]:
        assert bff.ll_score_type(tag) in declared, tag
    assert bff.ll_score_type("combined") in declared


def test_units_are_a_positive_claim_or_absent():
    """`dimensionless` is a claim; an empty unit means unknown, and the two
    must not be conflated."""
    assert "dimensionless" in bff.mfdb_units()
    assert "angstroms" in bff.mfdb_units()
    columns = [bff.MfdbColumn("x", "", "", "no unit stated")]
    assert "units" not in bff.mfdb_columns_json(columns)
    columns = [bff.MfdbColumn("x", "dimensionless", "", "a ratio")]
    assert "dimensionless" in bff.mfdb_columns_json(columns)


# ---------------------------------------------------------------------------
# Attribution
# ---------------------------------------------------------------------------

#: The attribution items, which are plain values rather than enumerations and
#: so are checked for existence rather than for their contents.
ATTRIBUTION_ITEMS = [
    "_mmfdb_artifact.author",
    "_mmfdb_artifact.citation",
    "_mmfdb_artifact.license",
    "_mmfdb_artifact.terms",
    "_mmfdb_artifact.terms_url",
    "_mmfdb_artifact.source",
    "_mmfdb_artifact.redistributed_via",
]


@pytest.mark.parametrize("item", ATTRIBUTION_ITEMS)
def test_the_attribution_items_are_declared_upstream(item):
    """Added in dictionary 1.8 so a licence notice can travel with the data.

    A converted rotamer library is somebody else's measurement in this
    package's format, and its terms have to ride along with it; before this
    they lived in a writer's private header JSON.
    """
    text = "".join(open(d).read() for d in DICTS)
    assert "save_%s\n" % item in text, "%s is not in the dictionary" % item


def test_terms_may_be_words_when_there_is_no_spdx_identifier():
    """Not every licence has an SPDX id, and forcing one would state something
    the upstream did not. "free for academic use" is the common case."""
    a = bff.MfdbAttribution()
    a.author = "Shapovalov MV, Dunbrack RL"
    a.citation = "10.1016/j.str.2011.03.019"
    a.terms = "free for academic use"
    a.source = "dun2010bbdep.bin"
    a.redistributed_via = "FASPR (MIT), Bioinformatics 2020;36:3758-3765"
    tags = {t.item: t.value for t in bff.mfdb_attribution_tags(a)}

    assert tags["_mmfdb_artifact.terms"] == "free for academic use"
    assert "_mmfdb_artifact.license" not in tags, "an empty field is omitted"
    # The chain matters: the licences differ along it.
    assert "FASPR" in tags["_mmfdb_artifact.redistributed_via"]
    assert tags["_mmfdb_artifact.source"] == "dun2010bbdep.bin"


def test_an_spdx_identifier_is_equally_acceptable():
    a = bff.MfdbAttribution()
    a.license = "GPL-3.0-only"
    tags = {t.item: t.value for t in bff.mfdb_attribution_tags(a)}
    assert tags["_mmfdb_artifact.license"] == "GPL-3.0-only"
    assert "_mmfdb_artifact.terms" not in tags


def test_an_artifact_with_neither_a_license_nor_terms_is_refused():
    """Unknown terms and unrestricted terms are different facts, and the
    difference matters to whoever redistributes the file."""
    a = bff.MfdbAttribution()
    a.author = "somebody"
    with pytest.raises(Exception) as excinfo:
        bff.mfdb_attribution_tags(a)
    assert "license" in str(excinfo.value) and "terms" in str(excinfo.value)


# ---------------------------------------------------------------------------
# The fluorophore vocabulary
# ---------------------------------------------------------------------------

#: Items declared in dictionary 1.8 for dye and optical-component properties.
#: Nothing in the loaded dictionaries described a quantum yield, an extinction
#: coefficient or a fluorescence spectrum before -- the only near matches were
#: `_em_detector.detective_quantum_efficiency` and the crystallographic
#: `_refine.ls_extinction_*` and NMR `_pdbx_nmr_spectral_*` families, none of
#: which is about a dye. Both categories describe tables that already exist in
#: this stack's databases; the dictionary is catching up to them.
OPTICAL_ITEMS = [
    "_mmfdb_optical_property.id",
    "_mmfdb_optical_property.probe_id",
    "_mmfdb_optical_property.property_name",
    "_mmfdb_optical_property.property_value",
    "_mmfdb_optical_property.unit",
    "_mmfdb_spectrum.id",
    "_mmfdb_spectrum.probe_id",
    "_mmfdb_spectrum.spectrum_type",
    "_mmfdb_spectrum.wavelengths",
    "_mmfdb_spectrum.intensity_values",
    "_mmfdb_spectrum.wavelength_unit",
    "_mmfdb_spectrum.intensity_unit",
]


@pytest.mark.parametrize("item", OPTICAL_ITEMS)
def test_the_dye_property_items_are_declared(item):
    text = "".join(open(d).read() for d in DICTS)
    assert "save_%s\n" % item in text, "%s is not in the dictionary" % item


def test_the_declared_items_map_onto_the_tables_that_exist():
    """The point of declaring them: a `.cif` and a database row say the same
    thing in the same words. If the declared column names drifted from the
    schema they describe, they would say different things."""
    text = "".join(open(d).read() for d in DICTS)
    for item, table, column in [
        ("_mmfdb_optical_property.property_name", "optical_properties",
         "property_name"),
        ("_mmfdb_optical_property.property_value", "optical_properties",
         "property_value"),
        ("_mmfdb_spectrum.spectrum_type", "spectra", "spectrum_type"),
        ("_mmfdb_spectrum.wavelengths", "spectra", "wavelengths"),
        ("_mmfdb_spectrum.intensity_values", "spectra", "intensity_values"),
    ]:
        frame = text[text.index("save_%s\n" % item):]
        frame = frame[:frame.index("\nsave_")]
        assert "_mmfdb_schema.table_name  %s" % table in frame, item
        assert "_mmfdb_schema.column_name %s" % column in frame, item


def test_the_spectrum_types_cover_the_optical_components_too():
    """The spectra table is shared with filters and detectors -- a filter's
    transmission is the same shape as a dye's emission. Declaring only the
    three dye types would have made every stored filter row non-conformant."""
    declared = _enumeration("_mmfdb_spectrum.spectrum_type")
    for value in ("absorption", "excitation", "emission"):
        assert value in declared
    for value in ("transmission", "reflectance", "quantum_efficiency",
                  "responsivity"):
        assert value in declared, "%s is stored in this stack" % value


def test_the_wavelength_unit_is_the_one_the_data_uses():
    """`nm`, not `nanometres`. The stored spectra are labelled `nm`, and a
    dictionary that renamed them would invalidate every row for nothing."""
    declared = _enumeration("_mmfdb_spectrum.wavelength_unit")
    assert "nm" in declared
    assert "nanometres" not in declared


def test_the_optical_property_names_are_guidance_not_a_closed_set():
    """`property_name` deliberately carries no `_item_enumeration`.

    An enumeration in a DDL2 dictionary is a closed set, and
    `schema_from_dictionary` turns it into a SQL CHECK constraint. The primary
    dataset in this stack would violate one: 19882 optical properties are
    stored under 76 names written by importers that never agreed. Declaring
    the set closed would make a fresh database reject the data it exists to
    hold -- a worse failure than an uncontrolled column.

    So the canonical names live in the item's *description*, and this checks
    they are all still named there. `mfdb_check_term` still refuses anything
    else **on write**, which is the half that can be enforced without lying
    about what is stored.
    """
    assert _enumeration("_mmfdb_optical_property.property_name") is None, (
        "an enumeration here would generate a CHECK the stored data violates")

    text = "".join(open(d).read() for d in DICTS)
    frame = text[text.index("save__mmfdb_optical_property.property_name\n"):]
    frame = frame[:frame.index("\nsave_")]
    for name in bff.mfdb_optical_property_names():
        assert name in frame, "%s is not documented in the item" % name

    # And the writer side is still closed.
    assert bff.mfdb_is_term("_mmfdb_optical_property.property_name",
                            "quantum_yield")
    assert not bff.mfdb_is_term("_mmfdb_optical_property.property_name", "qy")


def test_the_spectrum_type_enumeration_is_closed_because_the_data_is():
    """The contrast that makes the case above a judgement rather than a dodge:
    every `spectrum_type` actually stored is one of the seven declared, so a
    closed set is a true statement and the CHECK it generates is wanted."""
    declared = _enumeration("_mmfdb_spectrum.spectrum_type")
    assert declared is not None
    assert set(bff.mfdb_spectrum_types()) == set(declared)
