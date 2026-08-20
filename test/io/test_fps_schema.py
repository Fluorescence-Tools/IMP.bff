"""fps.json schema — drift guardrail and validation behaviour (PRD-97)."""

import json
from pathlib import Path

import IMP.bff
import IMP.bff as fps_schema

REPO = Path(__file__).resolve().parents[2]


def test_shipped_json_schema_matches_authored_definition():
    """data/fps_json_schema.json is derived, not hand-edited (drift guard)."""
    shipped = json.loads((REPO / "data" / "fps_json_schema.json").read_text())
    assert shipped == fps_schema.to_json_schema()


def test_all_intree_network_fps_json_positions_and_distances_conform():
    """Positions/Distances of every shipped network-dialect file validate."""
    for name in ("T4L/fret.fps.json", "GBP/mGBP2_FP.fps.json",
                 "GBP/hGBP1.fps.json", "TG2/flex.fps.json"):
        payload = json.loads(
            (REPO / "examples" / "structure" / name).read_text())
        for pname, pos in payload.get("Positions", {}).items():
            errors, _ = fps_schema.validate_position(pos, pname)
            assert not errors, (name, errors)
        positions = set(payload.get("Positions", {}))
        for dname, dist in payload.get("Distances", {}).items():
            errors, _ = fps_schema.validate_distance(
                dist, dname, position_names=positions)
            assert not errors, (name, errors)


def test_flat_dialect_templates_conform():
    pos = json.loads(
        (REPO / "data" / "template_av_position.fps.json").read_text())
    errors, warnings = fps_schema.validate_position(pos, "template")
    assert not errors, errors
    assert not warnings, warnings
    dist = json.loads(
        (REPO / "data" / "template_av_pair.fps.json").read_text())
    errors, warnings = fps_schema.validate_distance(dist, "template")
    assert not errors, errors
    assert not warnings, warnings


def test_validate_rejects_wrong_types_and_unknown_enum():
    errors, _ = fps_schema.validate_position(
        {"linker_length": "long", "simulation_type": "AV9"}, "p")
    assert any("linker_length" in e for e in errors)
    assert any("simulation_type" in e for e in errors)


def test_validate_requires_distance_fields_and_position_refs():
    payload = {
        "Positions": {"a": {"residue_seq_number": 1}},
        "Distances": {"d": {"position1_name": "a", "position2_name": "ghost",
                            "distance": 40.0, "error_neg": 2.0,
                            "error_pos": 2.0}},
        "χ²": {"s": {"distances": ["d", "missing"]}},
    }
    errors, _ = fps_schema.fps_schema_validate(payload)
    assert any("ghost" in e for e in errors)
    assert any("missing" in e for e in errors)
    # a required distance field left out
    del payload["Distances"]["d"]["distance"]
    errors, _ = fps_schema.fps_schema_validate(payload)
    assert any("missing required field 'distance'" in e for e in errors)


def test_xyz_position_requires_coordinates():
    errors, _ = fps_schema.validate_position(
        {"simulation_type": "XYZ", "x": 1.0, "y": 2.0}, "p")
    assert any("'z'" in e for e in errors)


def test_unknown_fields_warn_but_do_not_fail():
    errors, warnings = fps_schema.validate_position(
        {"residue_seq_number": 5, "my_custom_key": 1}, "p")
    assert not errors
    assert any("my_custom_key" in w for w in warnings)
