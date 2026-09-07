"""Session persistence: chinet's JSONL format, in bff.

Ported from chinet's chinet/session.py (phase 2 of removing chinet from
chisurf). A Session saves the graph it holds -- nodes with their ports,
plus the free-floating ports chisurf's fit parameters are -- in chinet's
exact on-disk format, so existing .csp projects load with no importer;
loading resolves the cross-references (session.nodes, node.ports,
port.link) by _id and rebuilds live objects.

What is deliberately not here (see include/Session.h):

- chinet's DB registry: a node or port reaches a file only through
  add_node()/add_port(), never by construction alone.
- the MMFDB backend, schema.py's chinet.session.v1 conversion, and the
  write_to_db/read_from_db object plumbing.
- the legacy monolithic {"session": ..., "objects": [...]} format is read
  (old projects on disk), not written.

test/session/chinet_fixture.jsonl is a real chinet write and is the spec:
the fixture tests reconstruct it field for field, and a re-saved session
must reproduce the same document per _id.
"""
import json
import uuid
from pathlib import Path

import pytest

from IMP.bff import Node, Port, Session, get_session

FIXTURE = Path(__file__).parent / "chinet_fixture.jsonl"

#: The fixture's ids, so the assertions name what they check.
SESSION_UID = "c6940af7-c4a8-46c6-a5ab-ee47537c2a15"
FIT_UID = "c763e84c-5cea-4ec3-ba14-ecf074976120"
GLOBAL_UID = "b8ec51b1-e58f-4550-86e1-de50619c6e1b"
X_UID = "c8d51fb7-c164-4daa-84c5-a734876e7bbe"
RATE_UID = "61295c1b-0c25-4341-92df-2c8e3f630428"
OUT_UID = "46ed8ddd-29ec-4b07-80e5-da9812528eb3"
KAPPA2_UID = "a705e6a7-34b3-4052-ad11-25693e9afbc9"

PORT_KEYS = ["_id", "name", "type", "precursor", "death", "fixed",
             "is_output", "is_reactive", "is_bounded", "value", "bounds",
             "link", "value_type", "prior"]
NODE_KEYS = ["_id", "name", "type", "precursor", "death", "callback",
             "callback_type", "valid", "ports"]
SESSION_KEYS = ["_id", "name", "type", "precursor", "death", "nodes"]


def build_session():
    """A graph shaped like the fixture, plus what the fixture does not
    cover: vectors, int ports, a reactive input, a prior, a free port."""
    s = Session()

    fit = Node("fit_0")
    x = Port(0.7)
    x.name = "x"
    x.fixed = True
    rate = Port(2.5, False, False, False, True, 0.0, 5.0)
    rate.name = "rate"
    out = Port(0.4, False, True)
    out.name = "out_00"
    fit.add_input_port("x", x)
    fit.add_input_port("rate", rate)
    fit.add_output_port("out_00", out)
    fit.set_callback("pass_on", "C")

    glob = Node("global_0")
    kappa2 = Port(0.4)
    kappa2.name = "kappa2"
    kappa2.prior = {"kind": "gaussian", "mu": 0.4, "sigma": 0.1}
    glob.add_input_port("kappa2", kappa2)
    glob.set_valid(True)

    extra = Node("extra_0")
    series = Port()
    series.name = "series"
    series.value = [1.5, 2.5, 3.5]
    count = Port(3)
    count.name = "count"
    trigger = Port(1.0)
    trigger.name = "trigger"
    trigger.reactive = True
    extra.add_input_port("series", series)
    extra.add_input_port("count", count)
    extra.add_input_port("trigger", trigger)

    free = Port(9.9, True)
    free.name = "free_param"

    out.link = kappa2  # out_00 follows kappa2 (glob stays valid)

    s.add_node("fit_0", fit)
    s.add_node("global_0", glob)
    s.add_node("extra_0", extra)
    s.add_port(free)
    return s


def docs_by_id(path):
    """A JSONL file's documents keyed by _id (field order lost, as a
    reader does; field order is asserted separately)."""
    out = {}
    for line in Path(path).read_text().splitlines():
        if line.strip():
            d = json.loads(line)
            out[d["_id"]] = d
    return out


# ---------------------------------------------------------------------------
# a. Round-trip of a built graph: every field
# ---------------------------------------------------------------------------


def test_round_trip_restores_every_port_field(tmp_path):
    s = build_session()
    path = str(tmp_path / "graph.jsonl")
    s.save(path)
    loaded = Session.load(path)
    assert loaded is not None

    fit = loaded.get_node("fit_0")
    x, rate, out = fit.get_port("x"), fit.get_port("rate"), fit.get_port("out_00")

    assert x.value == pytest.approx(0.7)
    assert x.fixed is True
    assert x.is_output is False
    assert x.bounded is False
    assert x.get_value_type() == 1
    assert x.get_is_vector() is False
    assert x.prior is None

    assert rate.value == pytest.approx(2.5)
    assert rate.bounded is True
    assert (rate.get_lower_bound(), rate.get_upper_bound()) == (0.0, 5.0)

    assert out.is_output is True
    assert out.get_is_reactive() is False


def test_round_trip_resolves_links(tmp_path):
    s = build_session()
    path = str(tmp_path / "graph.jsonl")
    s.save(path)
    loaded = Session.load(path)

    out = loaded.get_node("fit_0").get_port("out_00")
    kappa2 = loaded.get_node("global_0").get_port("kappa2")
    assert out.is_linked()
    assert out.link.uid == kappa2.uid  # the loaded kappa2 itself
    assert out.value == pytest.approx(0.4)  # read through the link
    kappa2.value = 0.25
    assert out.value == pytest.approx(0.25)  # still live, not copied


def test_round_trip_restores_vector_int_and_free_ports(tmp_path):
    s = build_session()
    path = str(tmp_path / "graph.jsonl")
    s.save(path)
    loaded = Session.load(path)

    series = loaded.get_node("extra_0").get_port("series")
    assert series.get_is_vector() is True
    assert series.get_value_type() == 3
    assert list(series.value) == pytest.approx([1.5, 2.5, 3.5])

    count = loaded.get_node("extra_0").get_port("count")
    assert count.get_value_type() == 0  # int scalar stays int
    assert count.get_is_vector() is False
    assert count.value == pytest.approx(3)
    assert count.fixed is False

    trigger = loaded.get_node("extra_0").get_port("trigger")
    assert trigger.get_is_reactive() is True

    free = loaded.get_port("free_param")
    assert free is not None
    assert free.get_node() is None  # free-floating, not adopted by a node
    assert free.fixed is True
    assert free.value == pytest.approx(9.9)


def test_round_trip_restores_node_state_and_keys(tmp_path):
    s = build_session()
    path = str(tmp_path / "graph.jsonl")
    s.save(path)
    loaded = Session.load(path)

    assert set(loaded.get_nodes().keys()) == {"fit_0", "global_0", "extra_0"}
    assert loaded.get_number_of_nodes() == 3
    assert loaded.get_number_of_ports() == 8  # 7 node ports + 1 free

    fit = loaded.get_node("fit_0")
    assert fit.get_callback() == "pass_on"
    assert fit.get_callback_type_string() == "C"
    assert set(fit.get_ports().keys()) == {"x", "rate", "out_00"}
    assert fit.get_node_valid() is False

    glob = loaded.get_node("global_0")
    assert glob.get_callback() == ""
    assert glob.get_callback_type_string() == ""
    assert glob.get_node_valid() is True  # raw flag, not is_valid()


def test_round_trip_preserves_prior(tmp_path):
    s = build_session()
    path = str(tmp_path / "graph.jsonl")
    s.save(path)
    loaded = Session.load(path)
    kappa2 = loaded.get_node("global_0").get_port("kappa2")
    assert kappa2.prior == {"kind": "gaussian", "mu": 0.4, "sigma": 0.1}


def test_round_trip_preserves_ids_and_precursors(tmp_path):
    s = build_session()
    ids = {"session": (s.uid, s.precursor)}
    for key in s.get_nodes().keys():
        n = s.get_node(key)
        ids[key] = (n.uid, n.precursor)
        for pk in n.get_ports().keys():
            p = n.get_port(pk)
            ids[pk] = (p.uid, p.precursor)
    for p in s.get_ports():
        ids[p.name] = (p.uid, p.precursor)

    path = str(tmp_path / "graph.jsonl")
    s.save(path)
    loaded = Session.load(path)

    assert (loaded.uid, loaded.precursor) == ids["session"]
    assert loaded.death == 0
    for key in loaded.get_nodes().keys():
        n = loaded.get_node(key)
        assert (n.uid, n.precursor) == ids[key]
        for pk in n.get_ports().keys():
            p = n.get_port(pk)
            assert (p.uid, p.precursor) == ids[pk]
    for p in loaded.get_ports():
        assert (p.uid, p.precursor) == ids[p.name]


def test_save_load_save_is_identical_per_id(tmp_path):
    s = build_session()
    p1 = str(tmp_path / "one.jsonl")
    p2 = str(tmp_path / "two.jsonl")
    p3 = str(tmp_path / "three.jsonl")
    s.save(p1)
    Session.load(p1).save(p2)
    Session.load(p2).save(p3)
    assert docs_by_id(p2) == docs_by_id(p1)
    assert docs_by_id(p3) == docs_by_id(p1)


def test_saved_field_order_matches_chinet(tmp_path):
    s = build_session()
    path = tmp_path / "graph.jsonl"
    s.save(str(path))
    lines = [json.loads(line) for line in path.read_text().splitlines()]

    assert list(lines[0].keys()) == SESSION_KEYS
    assert lines[0]["type"] == "session"
    assert lines[0]["name"] == "session"
    assert lines[0]["nodes"] == {
        key: s.get_node(key).uid for key in ("fit_0", "global_0", "extra_0")}

    node_lines = [d for d in lines if d.get("type") == "node"]
    port_lines = [d for d in lines if d.get("type") == "port"]
    assert len(node_lines) == 3 and len(port_lines) == 8
    for d in node_lines:
        assert list(d.keys()) == NODE_KEYS
    for d in port_lines:
        assert list(d.keys()) == PORT_KEYS

    fit_doc = next(d for d in node_lines if d["name"] == "fit_0")
    # the node's ports in the order they were added, chinet's dict order
    assert list(fit_doc["ports"].keys()) == ["x", "rate", "out_00"]
    assert fit_doc["ports"]["out_00"] == s.get_node("fit_0").get_port("out_00").uid


# ---------------------------------------------------------------------------
# b. The committed chinet fixture
# ---------------------------------------------------------------------------


def test_load_chinet_fixture():
    s = Session.load(str(FIXTURE))
    assert s is not None
    assert s.get_number_of_nodes() == 2
    # 4 ports in the fixture: x, rate, out_00 and kappa2 (the phase brief
    # said 5; the committed chinet write -- the spec -- holds 4)
    assert s.get_number_of_ports() == 4
    assert s.uid == SESSION_UID
    assert s.precursor == SESSION_UID
    assert s.name == "session"

    assert set(s.get_nodes().keys()) == {"fit_0", "global_0"}
    fit, glob = s.get_node("fit_0"), s.get_node("global_0")
    assert fit.uid == FIT_UID and glob.uid == GLOBAL_UID

    assert fit.get_callback() == "pass_on"
    assert fit.get_callback_type_string() == "C"
    assert fit.get_node_valid() is False

    x, rate, out = fit.get_port("x"), fit.get_port("rate"), fit.get_port("out_00")
    kappa2 = glob.get_port("kappa2")

    assert x.uid == X_UID
    assert x.fixed is True
    assert x.value == pytest.approx(0.7)
    assert x.get_value_type() == 1
    assert x.bounded is False
    assert x.prior is None

    assert rate.uid == RATE_UID
    assert rate.bounded is True
    assert (rate.get_lower_bound(), rate.get_upper_bound()) == (0.0, 5.0)
    assert rate.value == pytest.approx(2.5)
    assert rate.fixed is False

    assert out.uid == OUT_UID
    assert out.is_output is True
    assert out.link.uid == KAPPA2_UID
    assert out.link.uid == kappa2.uid  # link resolved to the loaded target
    assert out.value == pytest.approx(0.4)  # 0.4 via the link, not a copy

    assert kappa2.uid == KAPPA2_UID
    assert kappa2.value == pytest.approx(0.4)
    assert kappa2.get_is_output() is False
    assert glob.get_callback() == "" and glob.get_callback_type_string() == ""
    assert glob.get_node_valid() is False


def test_fixture_resaves_identically(tmp_path):
    """Load the real chinet write, save it twice: every document per _id
    must come back exactly as chinet wrote it."""
    p1 = str(tmp_path / "resave1.jsonl")
    p2 = str(tmp_path / "resave2.jsonl")
    Session.load(str(FIXTURE)).save(p1)
    Session.load(p1).save(p2)
    original = docs_by_id(FIXTURE)
    assert docs_by_id(p1) == original
    assert docs_by_id(p2) == original


# ---------------------------------------------------------------------------
# c. The legacy monolithic format
# ---------------------------------------------------------------------------


def monolithic_text(indent=None):
    docs = [json.loads(line) for line in FIXTURE.read_text().splitlines()
            if line.strip()]
    return json.dumps({"session": docs[0], "objects": docs[1:]}, indent=indent)


def test_load_legacy_monolithic(tmp_path):
    path = tmp_path / "legacy.jsonl"
    path.write_text(monolithic_text())
    s = Session.load(str(path))
    assert s is not None
    assert s.get_number_of_nodes() == 2
    assert s.get_number_of_ports() == 4  # as the fixture: x, rate, out_00, kappa2
    assert s.uid == SESSION_UID
    fit = s.get_node("fit_0")
    assert fit.get_callback() == "pass_on"
    out = fit.get_port("out_00")
    kappa2 = s.get_node("global_0").get_port("kappa2")
    assert out.link.uid == kappa2.uid
    assert out.value == pytest.approx(0.4)


def test_load_legacy_monolithic_pretty_printed(tmp_path):
    """The legacy writer could pretty-print across many lines; the first
    line does not parse as an object, which is how it is recognised."""
    path = tmp_path / "legacy_pretty.json"
    path.write_text(monolithic_text(indent=2))
    s = Session.load(str(path))
    assert s is not None
    assert s.get_number_of_nodes() == 2
    assert s.get_node("fit_0").get_port("rate").bounded is True


def test_int_type_codes_round_trip(tmp_path):
    """Codes 0 and 2 (int scalar, int vector) cannot be built from
    Python's value writes (lists go float); documents carry them, and
    persistence must keep them -- including writing JSON integers."""
    int_vector = {"_id": "vec", "name": "v", "type": "port",
                  "precursor": "vec", "death": 0, "fixed": False,
                  "is_output": False, "is_reactive": False,
                  "is_bounded": False, "value": [1, 2, 3],
                  "bounds": [0.0, 0.0], "link": None, "value_type": 2,
                  "prior": None}
    int_scalar = {"_id": "cnt", "name": "c", "type": "port",
                  "precursor": "cnt", "death": 0, "fixed": True,
                  "is_output": False, "is_reactive": False,
                  "is_bounded": False, "value": 7,
                  "bounds": [0.0, 0.0], "link": None, "value_type": 0,
                  "prior": None}
    node = {"_id": "n", "name": "n", "type": "node", "precursor": "n",
            "death": 0, "callback": "", "callback_type": "",
            "valid": False, "ports": {"v": "vec", "c": "cnt"}}
    session = {"_id": "s", "name": "session", "type": "session",
               "precursor": "s", "death": 0, "nodes": {"n": "n"}}
    path = tmp_path / "legacy_ints.json"
    path.write_text(json.dumps({"session": session,
                                "objects": [node, int_vector, int_scalar]}))

    s = Session.load(str(path))
    n = s.get_node("n")
    vec, cnt = n.get_port("v"), n.get_port("c")
    assert vec.get_value_type() == 2 and vec.get_is_vector() is True
    assert list(vec.value) == pytest.approx([1.0, 2.0, 3.0])
    assert cnt.get_value_type() == 0 and cnt.fixed is True
    assert cnt.value == pytest.approx(7)

    out = str(tmp_path / "ints.jsonl")
    s.save(out)
    saved = docs_by_id(out)
    assert saved["vec"]["value"] == [1, 2, 3]  # integers, not 1.0, 2.0, 3.0
    assert saved["vec"]["value_type"] == 2
    assert saved["cnt"]["value"] == 7
    assert saved["cnt"]["value_type"] == 0
    assert Session.load(out).get_node("n").get_port("v").get_value_type() == 2


# ---------------------------------------------------------------------------
# d. Empty sessions, empty files, errors
# ---------------------------------------------------------------------------


def test_empty_session_round_trip(tmp_path):
    s = Session()
    assert s.get_number_of_nodes() == 0
    assert s.get_number_of_ports() == 0
    path = tmp_path / "empty.jsonl"
    s.save(str(path))

    # byte-level: the one line a chinet save of an empty session wrote
    assert path.read_text() == (
        '{"_id": "%s", "name": "session", "type": "session", '
        '"precursor": "%s", "death": 0, "nodes": {}}\n' % (s.uid, s.uid))

    loaded = Session.load(str(path))
    assert loaded is not None
    assert loaded.get_number_of_nodes() == 0
    assert loaded.get_number_of_ports() == 0
    assert loaded.uid == s.uid
    assert loaded.precursor == s.precursor
    assert loaded.name == "session"


def test_load_empty_file_is_none(tmp_path):
    path = tmp_path / "nothing.jsonl"
    path.write_text("")
    assert Session.load(str(path)) is None


def test_jsonl_without_session_line_is_none(tmp_path):
    """chinet returns None when no document is a session; so does this."""
    docs = [json.loads(line) for line in FIXTURE.read_text().splitlines()
            if line.strip()]
    path = tmp_path / "headless.jsonl"
    path.write_text("\n".join(json.dumps(d) for d in docs[1:]) + "\n")
    assert Session.load(str(path)) is None


def test_bad_files_raise(tmp_path):
    garbage = tmp_path / "garbage.jsonl"
    garbage.write_text("{not json\n")
    with pytest.raises(RuntimeError):
        Session.load(str(garbage))
    with pytest.raises(RuntimeError):
        Session.load(str(tmp_path / "missing.jsonl"))
    with pytest.raises(RuntimeError):
        Session().save(str(tmp_path / "no" / "such" / "dir" / "f.jsonl"))


def test_clear_lets_go_of_everything(tmp_path):
    s = build_session()
    uid = s.uid
    s.clear()
    assert s.get_number_of_nodes() == 0
    assert s.get_number_of_ports() == 0
    assert s.uid == uid  # the session's own identity is untouched
    s.save(str(tmp_path / "cleared.jsonl"))
    assert len(docs_by_id(tmp_path / "cleared.jsonl")) == 1


# ---------------------------------------------------------------------------
# The module-level default session
# ---------------------------------------------------------------------------


def test_get_session_is_a_lazy_singleton():
    assert get_session() is get_session()
    assert isinstance(get_session(), Session)
    assert get_session().name == "session"


def test_no_silent_registration():
    """Ports and nodes exist only where they are held: constructing them
    puts nothing in a session (chinet's DB registered everything)."""
    stray = Port(1.0)
    stray.name = "stray"
    Node("stray")
    assert Session().get_number_of_ports() == 0
    assert Session().get_number_of_nodes() == 0
    # the default session starts empty too, whenever it is first made
    assert isinstance(get_session(), Session)


def test_fresh_uids_are_uuid_shaped():
    for obj in (Session(), Port(1.0), Node("n")):
        parsed = uuid.UUID(obj.uid)
        assert str(parsed) == obj.uid  # canonical dashed hex, 8-4-4-4-12
    s = Session()
    assert s.precursor == s.uid  # chinet: precursor defaults to the oid


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v", "-p", "no:cacheprovider"]))
