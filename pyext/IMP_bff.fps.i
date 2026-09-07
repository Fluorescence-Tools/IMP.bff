/*
 * fps.json: the definition, the reader, the writer.
 *
 * The C++ owns the format -- the field tables, the JSON-Schema derivation, the
 * validator, both doors and both legacy readers. It also owns the dict bridges
 * now: every function takes or returns JSON **text**, because an fps.json
 * object carries an open set of keys whose types differ per key and whose
 * unknown members survive the round trip -- that is a JSON object, and JSON
 * text is the one representation SWIG can marshal. A caller passes a dict as
 * `json.dumps(...)` and parses the text back with the same call; the round trip
 * is no longer hidden in an `.i` wrapper.
 *
 * The module constants are calls now: `fps_simulation_types()`,
 * `fps_av_simulation_types()`, `fps_distance_types()` (a dict proxy),
 * `fps_schema_version()` and `fps_position_fields()` / `fps_distance_fields()` /
 * `fps_score_set_fields()`.
 */

IMP_SWIG_VALUE(IMP::bff, FPSField, FPSFields);
IMP_SWIG_VALUE(IMP::bff, FPSValidation, FPSValidations);
IMP_SWIG_VALUE(IMP::bff, FPSDocument, FPSDocuments);
// The distance-measurement record the network restraint consumes; it is
// serialisable, so its value semantics carry the cereal hook too.
IMP_SWIG_VALUE(IMP::bff, AVPairDistanceMeasurement, AVPairDistanceMeasurements)
IMP_SWIG_VALUE_SERIALIZE_IMPL(IMP::bff, AVPairDistanceMeasurement)

%attribute(IMP::bff::FPSValidation, bool, is_valid, get_is_valid);

%include "IMP/bff/FPS.h"

%template(FPSFieldList) std::vector<IMP::bff::FPSField>;
%template(MapStringAVPairDistanceMeasurement) std::map<std::string, IMP::bff::AVPairDistanceMeasurement>;
