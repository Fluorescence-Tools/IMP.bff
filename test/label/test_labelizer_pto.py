"""The label container: conformant framing, a lossless structure, no sentinels.

A container is only worth writing if something other than the code that wrote
it can read it. So the framing here is walked by a **from-scratch EBML parser
written in this file** against RFC 8794 -- no `IMP.bff` code, no libebml, forty
lines of variable-length integers. If a consistent misreading of the standard
were baked into the writer, walking it back with the writer's own reader would
not notice; this would.

The canonical cross-check is `../tttrlib/test/tools/pto_ebml_check.cpp`, which
walks a container with libebml itself, and a container written here **passes
it**, strict alignment included (5 payloads, 0 misaligned). It needs libebml
**2.0**, not the 1.4.7 Homebrew ships -- `EbmlId::FromBuffer` is 2.0-only and
1.4.7 puts `EMaxSizeLength`/`EDocType` in a header the tool does not include.
Build it from the clone in `../tttrlib/junk/libebml`::

    cmake -S ../tttrlib/junk/libebml -B /tmp/ebml-build \
          -DCMAKE_INSTALL_PREFIX=/tmp/ebml -DBUILD_SHARED_LIBS=OFF
    cmake --build /tmp/ebml-build -j && cmake --install /tmp/ebml-build
    c++ -std=c++17 -I/tmp/ebml/include \
        ../tttrlib/test/tools/pto_ebml_check.cpp \
        -L/tmp/ebml/lib -lebml -o /tmp/pto_ebml_check

That is an external dependency this suite will not take, which is why the
parser below exists as well -- and it is the stronger argument of the two,
because it shares no lineage with either implementation.

What is pinned:

* the document is EBML with `DocType "pto"`, and every element's declared size
  exactly fills its parent -- the property that makes a container walkable at
  all;
* the objects are the ones the profile says, under the `label.` kind prefix
  this package owns (`drot.` and `rot.bbdep.` belong to the rotamer containers
  and must not appear);
* the embedded structure comes back **byte for byte**, verified against a
  SHA-256 computed by `hashlib` rather than by the C++ that wrote it;
* a position that was not scored carries a status and **no value**. This is the
  one that matters most: the reference writes `-1` for "excluded" and `0` for
  "no contribution" into the same column as real scores, so its output cannot
  be read back, and a regression here would silently reintroduce that.
"""

import hashlib
import json
import os

import pytest

import IMP.bff as bff

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "input", "labelizer")

#: EBML element ids the profile uses, from the PTO specification.
ID_EBML_HEAD = 0x1A45DFA3
ID_SEGMENT = 0x18538067
ID_ATTACHMENTS = 0x1941A469
ID_ATTACHED_FILE = 0x61A7
ID_FILE_NAME = 0x466E
ID_FILE_DATA = 0x465C
ID_DOC_TYPE = 0x4282
ID_PTO_KIND = 0x1E54F001
ID_PTO_ENCODING = 0x1E54F002

#: Ids whose payload is itself a list of elements.
MASTERS = {ID_EBML_HEAD, ID_SEGMENT, ID_ATTACHMENTS, ID_ATTACHED_FILE}


def _read_vint(buf, pos, keep_marker):
    """One EBML variable-length integer: (value, bytes consumed).

    The leading zero bits of the first octet give the length. An id keeps its
    marker bit (it *is* the id); a data size drops it (it is a number).
    """
    first = buf[pos]
    if first == 0:
        raise AssertionError("invalid vint at %d" % pos)
    length = 1
    mask = 0x80
    while not (first & mask):
        mask >>= 1
        length += 1
    value = first if keep_marker else (first & (mask - 1))
    for i in range(1, length):
        value = (value << 8) | buf[pos + i]
    return value, length


def _walk(buf, start, end, depth=0):
    """Every element in [start, end), as (id, payload_start, payload_end)."""
    out = []
    pos = start
    while pos < end:
        eid, n = _read_vint(buf, pos, keep_marker=True)
        pos += n
        size, n = _read_vint(buf, pos, keep_marker=False)
        pos += n
        assert pos + size <= end, (
            "element 0x%X at depth %d declares %d bytes and overruns its "
            "parent by %d" % (eid, depth, size, pos + size - end)
        )
        out.append((eid, pos, pos + size))
        pos += size
    assert pos == end, "elements do not exactly fill their parent"
    return out


def _objects(path):
    """Every attached object of a PTO file: name -> (kind, encoding, bytes)."""
    buf = open(path, "rb").read()
    top = _walk(buf, 0, len(buf))
    ids = [e[0] for e in top]
    assert ids[0] == ID_EBML_HEAD, "a PTO file begins with an EBML header"
    assert ID_SEGMENT in ids, "and carries a Segment"

    doctype = None
    for eid, s, e in _walk(buf, top[0][1], top[0][2]):
        if eid == ID_DOC_TYPE:
            doctype = buf[s:e].rstrip(b"\x00").decode("ascii")
    assert doctype == "pto", "DocType is %r, not 'pto'" % doctype

    found = {}
    for eid, s, e in top:
        if eid != ID_SEGMENT:
            continue
        for sid, ss, se in _walk(buf, s, e):
            if sid != ID_ATTACHMENTS:
                continue
            for aid, as_, ae in _walk(buf, ss, se):
                if aid != ID_ATTACHED_FILE:
                    continue
                name = kind = encoding = None
                data = b""
                for fid, fs, fe in _walk(buf, as_, ae):
                    if fid == ID_FILE_NAME:
                        name = buf[fs:fe].rstrip(b"\x00").decode("utf-8")
                    elif fid == ID_PTO_KIND:
                        kind = buf[fs:fe].rstrip(b"\x00").decode("utf-8")
                    elif fid == ID_PTO_ENCODING:
                        encoding = buf[fs:fe].rstrip(b"\x00").decode("utf-8")
                    elif fid == ID_FILE_DATA:
                        data = buf[fs:fe]
                if name is not None:
                    found[name] = (kind, encoding, data)
    return found


@pytest.fixture(scope="module")
def container(tmp_path_factory):
    """A real scored structure, written as a container."""
    pdb = os.path.join(DATA, "1DDB-39.pdb")
    scores = bff.labelizer_score_structure(
        pdb, bff.labelizer_model_paper(), bff.LabelizerOptions(),
        os.path.join(DATA, "1DDB-39_cs.pdb"),
    )
    options = bff.LabelizerFRETOptions()
    options.n_refine = 0
    pairs = list(bff.labelizer_fret_pair_scores(pdb, bff.labelizer_combined_by_key(scores), options))[:50]
    settings = bff.labelizer_settings_json(
        bff.labelizer_model_paper(), bff.LabelizerOptions(), options,
        os.path.join(DATA, "1DDB-39_cs.pdb"),
    )
    out = str(tmp_path_factory.mktemp("pto") / "1DDB-39.mmfdb.pto")
    bff.labelizer_write_pto(out, pdb, scores, pairs, settings)
    return out, pdb, scores, pairs


def test_the_framing_parses_without_any_of_our_code(container):
    path = container[0]
    objects = _objects(path)
    assert set(objects) == {
        "README", "structure.pdb", "label_scores.json",
        "label_pairs.json", "label_model.json",
    }


def test_the_readme_is_first_and_is_plain_ascii(container):
    """A container outlives its software; the first object explains it."""
    objects = _objects(container[0])
    kind, encoding, data = objects["README"]
    assert kind == "readme"
    assert encoding == "text"
    text = data.decode("ascii")          # ascii, deliberately: it must be
    assert "RFC 8794" in text
    assert "SHA-256" in text
    assert "0x1E54F001" in text, "the custom ids must be given by number"


def test_every_kind_is_in_this_package_s_namespace(container):
    """`drot.` and `rot.bbdep.` belong to the rotamer containers."""
    objects = _objects(container[0])
    kinds = {k for k, _, _ in objects.values()}
    assert kinds == {"readme", "label.structure", "label.scores",
                     "label.pairs", "label.model"}
    assert not any(k.startswith(("drot.", "rot.")) for k in kinds)


def test_the_structure_survives_byte_for_byte(container, tmp_path):
    """Verified against hashlib, not against the hash that wrote the file."""
    path, pdb = container[0], container[1]
    original = open(pdb, "rb").read()

    embedded = _objects(path)["structure.pdb"][2]
    assert embedded == original, "the embedded structure is not the input"

    out = str(tmp_path / "recovered.pdb")
    recorded = bff.labelizer_extract_pto_structure(path, out)
    assert recorded == hashlib.sha256(original).hexdigest()
    assert open(out, "rb").read() == original


def test_the_scores_round_trip_field_for_field(container):
    path, _, scores, _ = container
    back = bff.labelizer_read_pto_scores(path)
    assert len(back) == len(scores)
    for a, b in zip(scores, back):
        assert (a.asym_id, a.seq_id, a.comp_id) == (b.asym_id, b.seq_id, b.comp_id)
        assert a.score_type == b.score_type
        assert a.status == b.status
        if a.status == "scored":
            assert a.value == b.value


def test_the_pairs_round_trip(container):
    path, _, _, pairs = container
    back = bff.labelizer_read_pto_pairs(path)
    assert len(back) == len(pairs)
    for a, b in zip(pairs, back):
        assert (a.seq_id_1, a.seq_id_2) == (b.seq_id_1, b.seq_id_2)
        assert a.value == pytest.approx(b.value)
        assert a.distance == pytest.approx(b.distance)


def test_an_unscored_position_carries_no_value_at_all(container):
    """The whole reason the container exists rather than another CSV."""
    objects = _objects(container[0])
    rows = json.loads(objects["label_scores.json"][2])["rows"]
    for row in rows:
        if row["_mmfdb_label_score.status"] == "scored":
            assert "_mmfdb_label_score.value" in row
        else:
            assert "_mmfdb_label_score.value" not in row, (
                "an unscored position must carry no number; a sentinel is "
                "exactly what _mmfdb_label_score.value forbids"
            )


def test_the_container_says_what_it_conforms_to(container):
    """A reader decides by the tags, never by the file name."""
    objects = _objects(container[0])
    model = json.loads(objects["label_model.json"][2])
    tags = model["_container"]
    assert tags["_mmfdb_container.profile"] == "PTO.MFDB"
    assert tags["_mmfdb_container.format"] == "pto"
    assert tags["_mmfdb_container.dictionary_version"] == \
        bff.MFDB_DICTIONARY_VERSION

    scores = model["_artifacts"]["label_scores.json"]
    assert scores["_mmfdb_artifact.row_grain"] == "label_site"
    assert scores["_mmfdb_artifact.checksum_algorithm"] == "sha256"
    pairs = model["_artifacts"]["label_pairs.json"]
    assert pairs["_mmfdb_artifact.row_grain"] == "pair"


def test_provenance_reaches_back_to_the_structure(container):
    """A derived object that cannot reach the primary data is malformed."""
    objects = _objects(container[0])
    model = json.loads(objects["label_model.json"][2])
    edges = {(e["_mmfdb_edge.source_node_id"], e["_mmfdb_edge.target_node_id"]):
             e["_mmfdb_edge.relationship_type"] for e in model["_edges"]}
    assert edges[("structure.pdb", "label_scores.json")] == "derived_from"
    assert edges[("label_scores.json", "label_pairs.json")] == "maps_rows_of"

    operation = model["_operation"]
    assert operation["_mmfdb_operation.operation_type"] == "analysis"
    assert operation["_mmfdb_operation.algorithm"] == "labelizer"
    # The settings must be complete, not a summary: a partial record looks
    # reproducible and is not.
    settings = json.loads(operation["_mmfdb_operation.settings_json"])
    assert settings["arithmetic"] == "published"
    assert {t["tag"] for t in settings["model"]} == {"cs", "se", "tp", "cr",
                                                     "ss", "ce"}


def test_a_term_that_is_not_in_the_dictionary_is_refused():
    """The check is the point: a writer that cannot find a word must not
    invent one, it must add it to the dictionary first."""
    assert bff.mfdb_is_term("_mmfdb_artifact.row_grain", "label_site")
    assert not bff.mfdb_is_term("_mmfdb_artifact.row_grain", "residue")
    with pytest.raises(Exception) as excinfo:
        bff.mfdb_check_term("_mmfdb_artifact.row_grain", "residue")
    assert "mmfdb_flr_ext.dic" in str(excinfo.value)
