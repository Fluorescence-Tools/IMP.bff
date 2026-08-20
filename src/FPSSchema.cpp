/**
 * \file FPSSchema.cpp
 * \brief What an fps.json file may say, and whether one says it.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/FPSSchema.h>

#include <IMP/bff/internal/json.h>

#include <IMP/exception.h>

#include <algorithm>
#include <set>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

std::string fps_schema_version() { return "1.0"; }

// Named, not anonymous: IMP compiles this module as one translation unit.
namespace fps_schema {

//! One row of a field table, spelled the way the Python `dict(...)` calls were.
FPSField field(const std::string& name, FPSFieldType type,
               const std::string& default_json, const std::string& flrcif,
               const std::string& authored,
               const std::vector<std::string>& dialects,
               const std::vector<std::string>& enum_values = {},
               bool required = false) {
    FPSField f;
    f.name = name;
    f.type = type;
    f.default_json = default_json;
    f.flrcif = flrcif;
    f.authored = authored;
    f.dialects = dialects;
    f.enum_values = enum_values;
    f.required = required;
    return f;
}

const std::vector<std::string> NETWORK_FLAT = {"network", "flat"};
const std::vector<std::string> FLAT = {"flat"};
const std::vector<std::string> CSFPS = {"csfps"};

std::string json_type_name(FPSFieldType t) {
    switch (t) {
        case FPS_STRING: return "string";
        case FPS_NUMBER: return "number";
        case FPS_INTEGER: return "integer";
        case FPS_BOOLEAN: return "boolean";
        case FPS_ARRAY: return "array";
    }
    return "string";
}

//! The Python type name the old validator put in its message.
std::string python_type_name(FPSFieldType t) {
    switch (t) {
        case FPS_STRING: return "str";
        case FPS_NUMBER: return "float";
        case FPS_INTEGER: return "int";
        case FPS_BOOLEAN: return "bool";
        case FPS_ARRAY: return "list";
    }
    return "str";
}

//! What `type(value).__name__` would say, so the messages read as they did.
std::string value_type_name(const nlohmann::json& v) {
    if (v.is_string()) return "str";
    if (v.is_boolean()) return "bool";
    if (v.is_number_integer()) return "int";
    if (v.is_number_float()) return "float";
    if (v.is_array()) return "list";
    if (v.is_object()) return "dict";
    if (v.is_null()) return "NoneType";
    return "object";
}

//! `repr()` of a JSON scalar, for the message. Strings get single quotes.
std::string repr(const nlohmann::json& v) {
    if (v.is_string()) return "'" + v.get<std::string>() + "'";
    if (v.is_boolean()) return v.get<bool>() ? "True" : "False";
    if (v.is_null()) return "None";
    return v.dump();
}

bool type_matches(const nlohmann::json& v, FPSFieldType t) {
    switch (t) {
        case FPS_STRING: return v.is_string();
        // bool is an int subclass in Python; a bool where a number is expected
        // is a bug, and the old validator said so explicitly.
        case FPS_BOOLEAN: return v.is_boolean();
        case FPS_INTEGER: return v.is_number_integer() && !v.is_boolean();
        // An int is acceptable where a float is expected.
        case FPS_NUMBER: return v.is_number() && !v.is_boolean();
        case FPS_ARRAY: return v.is_array();
    }
    return false;
}

void check_fields(const nlohmann::json& obj, const std::vector<FPSField>& fields,
                  const std::string& where, FPSValidation& out) {
    if (!obj.is_object()) {
        out.errors.push_back(where + ": expected an object, got " +
                             value_type_name(obj));
        return;
    }
    for (std::size_t i = 0; i < fields.size(); ++i) {
        const FPSField& spec = fields[i];
        if (!obj.contains(spec.name)) {
            if (spec.required) {
                out.errors.push_back(where + ": missing required field '" +
                                     spec.name + "'");
            }
            continue;
        }
        const nlohmann::json& value = obj[spec.name];
        if (!type_matches(value, spec.type)) {
            out.errors.push_back(where + "." + spec.name + ": expected " +
                                 python_type_name(spec.type) + ", got " +
                                 value_type_name(value) + " (" + repr(value) +
                                 ")");
            continue;
        }
        if (!spec.enum_values.empty()) {
            const std::string as_text =
                    value.is_string() ? value.get<std::string>() : value.dump();
            if (std::find(spec.enum_values.begin(), spec.enum_values.end(),
                          as_text) == spec.enum_values.end()) {
                std::vector<std::string> sorted = spec.enum_values;
                std::sort(sorted.begin(), sorted.end());
                std::ostringstream listed;
                listed << "[";
                for (std::size_t k = 0; k < sorted.size(); ++k) {
                    if (k) listed << ", ";
                    listed << "'" << sorted[k] << "'";
                }
                listed << "]";
                out.errors.push_back(where + "." + spec.name + ": " +
                                     repr(value) + " is not one of " +
                                     listed.str());
            }
        }
    }
    for (nlohmann::json::const_iterator it = obj.begin(); it != obj.end(); ++it) {
        bool known = false;
        for (std::size_t i = 0; i < fields.size() && !known; ++i) {
            known = fields[i].name == it.key();
        }
        if (!known) {
            out.warnings.push_back(where + ": unknown field '" + it.key() + "'");
        }
    }
}

nlohmann::json parse(const std::string& text, const std::string& what) {
    try {
        return nlohmann::json::parse(text);
    } catch (const std::exception& e) {
        IMP_THROW("Cannot parse " << what << " as JSON: " << e.what(),
                  ValueException);
    }
}

nlohmann::json field_property(const FPSField& spec) {
    nlohmann::json prop;
    prop["type"] = json_type_name(spec.type);
    if (!spec.enum_values.empty()) prop["enum"] = spec.enum_values;
    if (!spec.default_json.empty()) {
        prop["default"] = nlohmann::json::parse(spec.default_json);
    }
    if (!spec.flrcif.empty()) prop["x-flrcif-item"] = spec.flrcif;
    if (!spec.authored.empty()) prop["description"] = spec.authored;
    prop["x-dialects"] = spec.dialects;
    return prop;
}

nlohmann::json object_schema(const std::vector<FPSField>& fields) {
    nlohmann::json properties = nlohmann::json::object();
    std::vector<std::string> required;
    for (std::size_t i = 0; i < fields.size(); ++i) {
        properties[fields[i].name] = field_property(fields[i]);
        if (fields[i].required) required.push_back(fields[i].name);
    }
    std::sort(required.begin(), required.end());
    nlohmann::json schema;
    schema["type"] = "object";
    schema["properties"] = properties;
    // Unknown keys are preserved by readers and writers but are not silently
    // part of the format: this stays true so old files validate, and the
    // validator reports them as warnings.
    schema["additionalProperties"] = true;
    if (!required.empty()) schema["required"] = required;
    return schema;
}

}  // namespace fps_schema

std::vector<std::string> fps_simulation_types() {
    return {"AV1", "AV3", "XYZ", "R1"};
}

std::vector<std::string> fps_av_simulation_types() {
    return {"AV1", "AV3", "XYZ"};
}

std::map<std::string, std::string> fps_distance_types() {
    std::map<std::string, std::string> out;
    out["RDAMean"] = "<R_DA>";
    out["RDAMeanE"] = "<R_DA>_E";
    out["Rmp"] = "R_mp";
    return out;
}

std::vector<FPSField> fps_position_fields() {
    using namespace fps_schema;
    std::vector<FPSField> f;
    // --- labelling site ---------------------------------------------------
    f.push_back(field("chain_identifier", FPS_STRING, "\"\"", "",
                      "Chain (asym) id of the attachment atom; empty matches any.",
                      NETWORK_FLAT));
    f.push_back(field("residue_seq_number", FPS_INTEGER, "0", "",
                      "Author residue sequence number of the attachment atom.",
                      NETWORK_FLAT));
    f.push_back(field("residue_name", FPS_STRING, "\"\"", "",
                      "Residue name of the attachment atom (informational).",
                      FLAT));
    f.push_back(field("atom_name", FPS_STRING, "\"CA\"", "",
                      "Atom name the dye linker attaches to (e.g. CA, CB).",
                      NETWORK_FLAT));
    // --- AV geometry (flrCIF FPS AV parameters) ---------------------------
    f.push_back(field("linker_length", FPS_NUMBER, "20.0",
                      "_flr_FPS_AV_parameter.linker_length", "", NETWORK_FLAT));
    f.push_back(field("linker_width", FPS_NUMBER, "0.5",
                      "_flr_FPS_AV_parameter.linker_width", "", NETWORK_FLAT));
    f.push_back(field("radius1", FPS_NUMBER, "3.5",
                      "_flr_FPS_AV_parameter.probe_radius_1", "", NETWORK_FLAT));
    f.push_back(field("radius2", FPS_NUMBER, "0.0",
                      "_flr_FPS_AV_parameter.probe_radius_2", "", NETWORK_FLAT));
    f.push_back(field("radius3", FPS_NUMBER, "0.0",
                      "_flr_FPS_AV_parameter.probe_radius_3", "", NETWORK_FLAT));
    f.push_back(field("simulation_type", FPS_STRING, "\"AV1\"", "",
                      "Label model: AV1 (one radius), AV3 (three radii), XYZ "
                      "(fixed mean position, no simulation), or R1 (rotamer "
                      "ensemble: screened rotamer library, Python-only).",
                      {"network", "flat", "csfps"}, fps_simulation_types()));
    // --- rotamer-ensemble positions (R1, imp.bff dialect, PRD-108) --------
    f.push_back(field("rotamer_library", FPS_STRING, "\"\"", "",
                      "Rotamer library name (IMP.bff registry / FRETpredict "
                      "spelling, e.g. 'AlexaFluor 488 C1R cutoff30'). Required "
                      "for simulation_type R1.",
                      FLAT));
    f.push_back(field("dye_name", FPS_STRING, "\"\"", "",
                      "Dye name for R0 from spectra (e.g. 'AlexaFluor 488').",
                      FLAT));
    f.push_back(field("temperature", FPS_NUMBER, "298.15", "",
                      "Screening temperature (K) of an R1 ensemble.", FLAT));
    f.push_back(field("electrostatic", FPS_BOOLEAN, "false", "",
                      "Add Debye-Hueckel electrostatics to the R1 screening.",
                      FLAT));
    f.push_back(field("potential", FPS_STRING, "\"lj\"", "",
                      "Screening potential of an R1 ensemble: 'lj' or 'gauss'.",
                      FLAT, {"lj", "gauss"}));
    f.push_back(field("simulation_grid_resolution", FPS_NUMBER, "1.5", "",
                      "AV grid spacing in Angstrom. Related to (but not the same "
                      "item as) _flr_FPS_global_parameter.AV_min_grid_A, which is "
                      "the *minimum* grid constant of FPS's relative-grid scheme.",
                      NETWORK_FLAT));
    // --- accessible-contact-volume extensions (imp.bff dialect) -----------
    f.push_back(field("allowed_sphere_radius", FPS_NUMBER, "1.5",
                      "_flr_FPS_global_parameter.AV_allowed_sphere", "", FLAT));
    f.push_back(field("contact_volume_thickness", FPS_NUMBER, "0.0", "",
                      "Thickness of the contact layer above the molecular surface "
                      "for accessible-contact-volume (ACV) weighting; 0 disables.",
                      FLAT));
    f.push_back(field("contact_volume_trapped_fraction", FPS_NUMBER, "-1.0", "",
                      "Fraction of dye density trapped in the contact volume; "
                      "negative disables ACV re-weighting.",
                      FLAT));
    f.push_back(field("min_sphere_volume_fraction", FPS_NUMBER, "0.0", "",
                      "Minimum accessible fraction of the free-dye sphere volume "
                      "below which a position is flagged as buried.",
                      FLAT));
    f.push_back(field("anchor_atoms", FPS_STRING, "\"\"", "",
                      "Explicit anchor atom names overriding the default linker "
                      "anchor search (comma separated).",
                      FLAT));
    f.push_back(field("strip_mask", FPS_STRING, "\"\"", "",
                      "Selection mask of atoms removed as obstacles before the AV "
                      "simulation, in the PyMOL dialect 'chain <id> and resid <n> "
                      "and [not] name A+B+...' ('+'-separated lists, as in PyMOL "
                      "itself). Empty means the default strip: the attachment "
                      "residue's side chain minus the attachment atom.",
                      FLAT));
    f.push_back(field("chain_weighting", FPS_BOOLEAN, "false", "",
                      "Weight AV grid points by linker-chain statistics instead "
                      "of uniformly.",
                      FLAT));
    // --- rigid-body / multi-body assignment (ChiSurf dialect) -------------
    f.push_back(field("body_id", FPS_INTEGER, "0", "",
                      "Index of the rigid body (input PDB) this position labels; "
                      "docking moves bodies relative to each other.",
                      {"network", "csfps"}));
    // --- fixed positions (XYZ simulation_type, legacy C# conversions) -----
    f.push_back(field("x", FPS_NUMBER, "",
                      "_flr_FPS_mean_probe_position.mpp_xcoord", "", CSFPS));
    f.push_back(field("y", FPS_NUMBER, "",
                      "_flr_FPS_mean_probe_position.mpp_ycoord", "", CSFPS));
    f.push_back(field("z", FPS_NUMBER, "",
                      "_flr_FPS_mean_probe_position.mpp_zcoord", "", CSFPS));
    return f;
}

std::vector<FPSField> fps_distance_fields() {
    using namespace fps_schema;
    std::vector<std::string> types;
    const std::map<std::string, std::string> dt = fps_distance_types();
    // Insertion order in the Python was RDAMean, RDAMeanE, Rmp; a std::map
    // sorts, and those two orders happen to agree.
    for (std::map<std::string, std::string>::const_iterator it = dt.begin();
         it != dt.end(); ++it) {
        types.push_back(it->first);
    }
    std::vector<FPSField> f;
    f.push_back(field("position1_name", FPS_STRING, "",
                      "_flr_fret_distance_restraint.sample_probe_id_1", "",
                      NETWORK_FLAT, {}, true));
    f.push_back(field("position2_name", FPS_STRING, "",
                      "_flr_fret_distance_restraint.sample_probe_id_2", "",
                      NETWORK_FLAT, {}, true));
    f.push_back(field("distance", FPS_NUMBER, "",
                      "_flr_fret_distance_restraint.distance", "", NETWORK_FLAT,
                      {}, true));
    f.push_back(field("error_neg", FPS_NUMBER, "",
                      "_flr_fret_distance_restraint.distance_error_minus", "",
                      NETWORK_FLAT, {}, true));
    f.push_back(field("error_pos", FPS_NUMBER, "",
                      "_flr_fret_distance_restraint.distance_error_plus", "",
                      NETWORK_FLAT, {}, true));
    f.push_back(field("distance_type", FPS_STRING, "\"RDAMean\"",
                      "_flr_fret_distance_restraint.distance_type", "",
                      NETWORK_FLAT, types));
    f.push_back(field("Forster_radius", FPS_NUMBER, "52.0",
                      "_flr_fret_forster_radius.forster_radius", "",
                      NETWORK_FLAT));
    return f;
}

std::vector<FPSField> fps_score_set_fields() {
    using namespace fps_schema;
    std::vector<FPSField> f;
    // A score set (the "chi2" top-level section) selects a named subset of the
    // distances to score together. flrCIF groups restraints via
    // `_flr_fret_distance_restraint.group_id`; the named-set spelling here is
    // the authored fps.json form of the same concept.
    f.push_back(field("distances", FPS_ARRAY, "",
                      "_flr_fret_distance_restraint.group_id",
                      "Names of the Distances entries belonging to this set.",
                      {}, {}, true));
    f.push_back(field("maximum_NaNs_allowed", FPS_INTEGER, "", "",
                      "Olga screening: maximum number of NaN model distances "
                      "tolerated before a structure's score set is rejected.",
                      {"network"}));
    f.push_back(field("penalty_NaN", FPS_NUMBER, "", "",
                      "Olga screening: chi-square penalty added per NaN model "
                      "distance.",
                      {"network"}));
    return f;
}

std::string fps_json_schema() {
    const nlohmann::json position =
            fps_schema::object_schema(fps_position_fields());
    const nlohmann::json distance =
            fps_schema::object_schema(fps_distance_fields());
    const nlohmann::json score_set =
            fps_schema::object_schema(fps_score_set_fields());

    nlohmann::json out;
    out["$schema"] = "https://json-schema.org/draft/2020-12/schema";
    out["$id"] = "https://integrativemodeling.org/schemas/bff/fps.json";
    out["title"] = "fps.json — FRET labelling positions and distances";
    out["x-schema-version"] = fps_schema_version();
    out["type"] = "object";
    nlohmann::json properties;
    properties["Positions"] = {{"type", "object"},
                               {"additionalProperties", position}};
    properties["Distances"] = {{"type", "object"},
                               {"additionalProperties", distance}};
    properties["χ²"] = {{"type", "object"},
                        {"additionalProperties", score_set}};
    properties["Evaluators"] = {{"type", "array"},
                                {"items", {{"type", "object"}}}};
    out["properties"] = properties;
    out["additionalProperties"] = true;
    nlohmann::json defs;
    defs["position"] = position;
    defs["distance"] = distance;
    defs["score_set"] = score_set;
    out["$defs"] = defs;
    return out.dump(2);
}

namespace fps_schema {

void validate_position_json(const nlohmann::json& position,
                            const std::string& name, FPSValidation& out) {
    check_fields(position, fps_position_fields(), name, out);
    if (!position.is_object()) return;

    std::string stype = "AV1";
    if (position.contains("simulation_type") &&
        position["simulation_type"].is_string()) {
        stype = position["simulation_type"].get<std::string>();
    }
    if (stype == "R1") {
        std::string lib;
        if (position.contains("rotamer_library") &&
            position["rotamer_library"].is_string()) {
            lib = position["rotamer_library"].get<std::string>();
        }
        const std::size_t a = lib.find_first_not_of(" \t\r\n");
        if (a == std::string::npos) {
            out.errors.push_back(
                    name + ": simulation_type R1 requires 'rotamer_library'");
        }
        const char* av_keys[] = {"linker_length", "linker_width", "radius1"};
        for (int i = 0; i < 3; ++i) {
            if (position.contains(av_keys[i])) {
                out.warnings.push_back(
                        name + ": '" + av_keys[i] +
                        "' is an AV parameter; an R1 position ignores it "
                        "(Python side) or is scored as AV1 by "
                        "AVNetworkRestraint -- filter with "
                        "fps_positions_for_docking()");
            }
        }
    }
    if (stype == "XYZ") {
        const char* keys[] = {"x", "y", "z"};
        for (int i = 0; i < 3; ++i) {
            if (!position.contains(keys[i])) {
                out.errors.push_back(name + ": simulation_type XYZ requires '" +
                                     keys[i] + "'");
            }
        }
    } else if (stype == "AV3") {
        const char* keys[] = {"radius2", "radius3"};
        for (int i = 0; i < 2; ++i) {
            double value = 0.0;
            if (position.contains(keys[i]) && position[keys[i]].is_number()) {
                value = position[keys[i]].get<double>();
            }
            if (!(value > 0.0)) {
                out.warnings.push_back(
                        name + ": simulation_type AV3 with non-positive '" +
                        keys[i] + "' behaves like AV1");
            }
        }
    }
}

void validate_distance_json(const nlohmann::json& distance,
                            const std::string& name,
                            const std::set<std::string>* position_names,
                            FPSValidation& out) {
    check_fields(distance, fps_distance_fields(), name, out);
    if (position_names == nullptr || !distance.is_object()) return;
    const char* keys[] = {"position1_name", "position2_name"};
    for (int i = 0; i < 2; ++i) {
        if (!distance.contains(keys[i]) || !distance[keys[i]].is_string()) {
            continue;
        }
        const std::string pname = distance[keys[i]].get<std::string>();
        if (!pname.empty() && position_names->count(pname) == 0) {
            out.errors.push_back(name + "." + keys[i] + ": '" + pname +
                                 "' is not a defined position");
        }
    }
}

}  // namespace fps_schema

FPSValidation validate_position(const std::string& position_json,
                                const std::string& name) {
    FPSValidation out;
    fps_schema::validate_position_json(
            fps_schema::parse(position_json, "a position"), name, out);
    return out;
}

FPSValidation validate_distance(const std::string& distance_json,
                                const std::string& name,
                                const std::vector<std::string>& position_names) {
    FPSValidation out;
    std::set<std::string> names(position_names.begin(), position_names.end());
    fps_schema::validate_distance_json(
            fps_schema::parse(distance_json, "a distance"), name,
            position_names.empty() ? nullptr : &names, out);
    return out;
}

FPSValidation fps_schema_validate(const std::string& payload_json) {
    FPSValidation out;
    const nlohmann::json payload = fps_schema::parse(payload_json, "the payload");
    if (!payload.is_object()) {
        out.errors.push_back("payload: expected an object, got " +
                             fps_schema::value_type_name(payload));
        return out;
    }

    nlohmann::json positions = nlohmann::json::object();
    if (payload.contains("Positions")) {
        if (payload["Positions"].is_object()) {
            positions = payload["Positions"];
        } else {
            out.errors.push_back("Positions: expected an object");
        }
    }
    for (nlohmann::json::const_iterator it = positions.begin();
         it != positions.end(); ++it) {
        fps_schema::validate_position_json(
                it.value(), "Positions['" + it.key() + "']", out);
    }

    std::set<std::string> position_names;
    for (nlohmann::json::const_iterator it = positions.begin();
         it != positions.end(); ++it) {
        position_names.insert(it.key());
    }

    nlohmann::json distances = nlohmann::json::object();
    if (payload.contains("Distances")) {
        if (payload["Distances"].is_object()) {
            distances = payload["Distances"];
        } else {
            out.errors.push_back("Distances: expected an object");
        }
    }
    for (nlohmann::json::const_iterator it = distances.begin();
         it != distances.end(); ++it) {
        fps_schema::validate_distance_json(
                it.value(), "Distances['" + it.key() + "']", &position_names,
                out);
    }

    nlohmann::json score_sets = nlohmann::json::object();
    if (payload.contains("χ²")) {
        if (payload["χ²"].is_object()) {
            score_sets = payload["χ²"];
        } else {
            out.errors.push_back("χ²: expected an object");
        }
    }
    for (nlohmann::json::const_iterator it = score_sets.begin();
         it != score_sets.end(); ++it) {
        const std::string where = "χ²['" + it.key() + "']";
        fps_schema::check_fields(it.value(), fps_score_set_fields(), where, out);
        if (!it.value().is_object()) continue;
        if (!it.value().contains("distances") ||
            !it.value()["distances"].is_array()) {
            continue;
        }
        const nlohmann::json& refs = it.value()["distances"];
        for (std::size_t k = 0; k < refs.size(); ++k) {
            if (!refs[k].is_string()) continue;
            const std::string ref = refs[k].get<std::string>();
            if (!distances.contains(ref)) {
                out.errors.push_back(where + ": '" + ref +
                                     "' is not a defined distance");
            }
        }
    }

    if (payload.contains("Evaluators") && !payload["Evaluators"].is_null() &&
        !payload["Evaluators"].is_array()) {
        out.errors.push_back("Evaluators: expected an array");
    }
    return out;
}

IMPBFF_END_NAMESPACE
