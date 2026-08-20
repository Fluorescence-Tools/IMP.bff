/*
 * fps.json: the definition, the reader, the writer.
 *
 * The C++ owns the format -- the field tables, the JSON-Schema derivation, the
 * validator, both doors and both legacy readers. Everything below is one
 * `json.loads`/`json.dumps` per boundary crossing, because an fps.json object
 * carries an open set of keys whose types differ per key and whose unknown
 * members survive the round trip: that is a JSON object, and a dict is what
 * Python calls one.
 *
 * `read_evaluators_json`'s `factory` stays here for the same reason
 * `radial_diffusion_map` stayed in `IMP_bff.quenching.i`: it is a Python
 * callable, and this module defines the *format* -- what an evaluator *is*
 * belongs to the application that passes the factory.
 */

IMP_SWIG_VALUE(IMP::bff, FPSField, FPSFields);
IMP_SWIG_VALUE(IMP::bff, FPSValidation, FPSValidations);
IMP_SWIG_VALUE(IMP::bff, FPSDocument, FPSDocuments);

%rename(_fps_json_schema) IMP::bff::fps_json_schema;
%rename(_fps_distance_types) IMP::bff::fps_distance_types;
%rename(_validate_position) IMP::bff::validate_position;
%rename(_validate_distance) IMP::bff::validate_distance;
%rename(_fps_schema_validate) IMP::bff::fps_schema_validate;
%rename(_read_old_lps_txt) IMP::bff::read_old_lps_txt;
%rename(_read_old_distances_txt) IMP::bff::read_old_distances_txt;
%rename(_read_fps_json) IMP::bff::read_fps_json;
%rename(_write_fps_json) IMP::bff::write_fps_json;
%rename(_fps_positions_for_docking) IMP::bff::fps_positions_for_docking;
%rename(_read_evaluators_json) IMP::bff::read_evaluators_json;
%rename(_write_evaluators_json) IMP::bff::write_evaluators_json;

%include "IMP/bff/FPSSchema.h"
%include "IMP/bff/FPSIO.h"

%template(FPSFieldList) std::vector<IMP::bff::FPSField>;

%pythoncode %{
import json as _json

#: ``simulation_type`` values a file may carry. AV1: one dye radius; AV3:
#: three; XYZ: a fixed mean position (no volume simulation), carried by legacy
#: C# FPS conversions; R1: a rotamer ensemble -- a FRETpredict-style library
#: placed in the residue's backbone frame and Boltzmann-screened against the
#: structure (R2 = sampled library, R3 = mixture are reserved).
SIMULATION_TYPES = tuple(_IMP_bff.fps_simulation_types())

#: The subset ``IMP::bff::AVNetworkRestraint`` understands. R1 is missing, and
#: that is the point: the C++ never reads ``simulation_type`` and would score an
#: R1 position as an AV1 with its AV parameters. Filter with
#: :func:`fps_positions_for_docking` first.
AV_SIMULATION_TYPES = tuple(_IMP_bff.fps_av_simulation_types())

#: ``distance_type`` values and their flrCIF spellings. RDAMean = mean inter-dye
#: distance <R_DA>; RDAMeanE = FRET-averaged distance <R_DA>_E; Rmp = distance
#: between mean dye positions R_mp.
DISTANCE_TYPES = dict(_IMP_bff._fps_distance_types())

SCHEMA_VERSION = _IMP_bff.fps_schema_version()


def _fields_as_dict(fields):
    """A C++ field table in the shape the Python tables had.

    ``{name: {"type": <python type>, "default": ..., "flrcif": ..., ...}}``.
    The types are the Python ones a validator would have compared against, so a
    caller reading the table still reads ``spec["type"] is float``.
    """
    types = {_IMP_bff.FPS_STRING: str, _IMP_bff.FPS_NUMBER: float,
             _IMP_bff.FPS_INTEGER: int, _IMP_bff.FPS_BOOLEAN: bool,
             _IMP_bff.FPS_ARRAY: list}
    out = {}
    for f in fields:
        spec = {
            "type": types[f.type],
            "default": (None if not f.default_json
                        else _json.loads(f.default_json)),
            "flrcif": f.flrcif or None,
            "dialects": tuple(f.dialects),
        }
        if f.authored:
            spec["authored"] = f.authored
        if f.enum_values:
            spec["enum"] = tuple(f.enum_values)
        if f.required:
            spec["required"] = True
        out[f.name] = spec
    return out


POSITION_FIELDS = _fields_as_dict(_IMP_bff.fps_position_fields())
DISTANCE_FIELDS = _fields_as_dict(_IMP_bff.fps_distance_fields())
SCORE_SET_FIELDS = _fields_as_dict(_IMP_bff.fps_score_set_fields())


def to_json_schema():
    """The JSON-Schema document for the network-dialect fps.json.

    The shipped ``data/fps_json_schema.json`` is this function's output; a test
    regenerates it and fails on drift, so the C++ field tables stay the single
    authored definition.
    """
    return _json.loads(_IMP_bff._fps_json_schema())


def validate_position(position, name="position"):
    """Validate one position object. Returns ``(errors, warnings)``."""
    report = _IMP_bff._validate_position(_json.dumps(position), str(name))
    return list(report.errors), list(report.warnings)


def validate_distance(distance, name="distance", position_names=None):
    """Validate one distance object. Returns ``(errors, warnings)``.

    *position_names* cross-checks the two ends against the defined positions;
    ``None`` skips the check, which is what a single-object check wants.
    """
    report = _IMP_bff._validate_distance(
        _json.dumps(distance), str(name),
        [] if position_names is None else sorted(str(n) for n in position_names))
    return list(report.errors), list(report.warnings)


def fps_schema_validate(payload):
    """Validate a parsed network-dialect fps.json payload.

    ``errors`` violate the schema; ``warnings`` flag unknown fields and
    suspicious-but-legal values. An empty ``errors`` list means it conforms.
    """
    report = _IMP_bff._fps_schema_validate(_json.dumps(payload))
    return list(report.errors), list(report.warnings)


#: The name the schema section carried before it was one module.
validate = fps_schema_validate


def read_old_lps_txt(path, pdb_paths=None):
    """Read a C# FPS labeling-positions ``.txt`` file.

    :returns: ``(positions, molecules)`` -- the positions dict, and the unique
        molecule names in order of appearance.
    """
    doc = _IMP_bff._read_old_lps_txt(
        str(path), [] if pdb_paths is None else [str(p) for p in pdb_paths])
    return _json.loads(doc.positions), list(doc.molecules)


def read_old_distances_txt(path):
    """Read a C# FPS experimental-distances ``.txt`` file."""
    return _json.loads(_IMP_bff._read_old_distances_txt(str(path)))


def read_fps_json(path, pdb_paths=None, validate=False):
    """Load an fps.json labeling file, or old C# ``.txt`` format files.

    A path that does not end in ``.json`` is read as a legacy positions file,
    and a ``Distances.txt`` beside it is read for the distances.

    :returns: ``(positions, distances, score_sets, extra)``.
    """
    doc = _IMP_bff._read_fps_json(
        str(path), [] if pdb_paths is None else [str(p) for p in pdb_paths],
        bool(validate))
    return (_json.loads(doc.positions), _json.loads(doc.distances),
            _json.loads(doc.score_sets), _json.loads(doc.extra))


def write_fps_json(path, positions, distances, score_sets=None, extra=None,
                   validate=False, **kwargs):
    """Write an fps.json file.

    :param score_sets: the ``"χ²"`` section; falsy omits the key.
    :param validate: check the assembled payload and raise rather than write a
        non-conforming file.
    """
    _IMP_bff._write_fps_json(
        str(path), _json.dumps(positions), _json.dumps(distances),
        _json.dumps(score_sets or {}), _json.dumps(extra or {}), bool(validate))


def fps_positions_for_docking(positions, distances=None):
    """Keep only the positions the C++ AV scorer understands.

    Rotamer-ensemble positions (``simulation_type == "R1"``) are Python-only;
    ``IMP::bff::AVNetworkRestraint`` would score them as AV1 with their AV
    parameters. Returns ``(positions, distances)`` restricted to AV1/AV3/XYZ
    positions and to distances whose two ends survive.
    """
    doc = _IMP_bff._fps_positions_for_docking(
        _json.dumps(positions), _json.dumps({} if distances is None else distances))
    kept = _json.loads(doc.positions)
    return (kept, {}) if distances is None else (kept, _json.loads(doc.distances))


def write_evaluators_json(path, evaluators):
    """Replace the ``Evaluators`` key of an fps.json file, keeping the rest."""
    _IMP_bff._write_evaluators_json(str(path), _json.dumps(
        [ev.to_dict() if hasattr(ev, "to_dict") else dict(ev)
         for ev in evaluators]))


def read_evaluators_json(path, factory=None):
    """Read the ``Evaluators`` list from an fps.json file.

    :param factory: called on each evaluator dict to instantiate an application
        object (ChiSurf passes its ``evaluators.from_dict``). Entries the
        factory raises on are skipped. Without a factory the raw dicts come
        back -- this module defines the *format*; what an evaluator *is*
        belongs to the application.
    """
    raw = _json.loads(_IMP_bff._read_evaluators_json(str(path)))
    if factory is None:
        return raw
    out = []
    for d in raw:
        try:
            out.append(factory(d))
        except Exception:
            pass
    return out
%}
