/**
 * \file FPSSchema.cpp
 * \brief What an fps.json file may say, and whether one says it.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/FPSSchema.h>

#include <IMP/bff/internal/Text.h>

#include <IMP/bff/internal/json.h>

#include <IMP/bff/Base.h>

#include <algorithm>
#include <set>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

//! 1.7 (2026-09-01) added `coarse_clash`, `clash_radii_source` and
//! `clash_radii_scale` to the `Project` section: which van der Waals radii the
//! **excluded-volume** term measures overlap by, and by what factor. 1.6 gave
//! the accessible volume a radii source and left the other half of the docking
//! score reading `IMP::core::XYZR`; these three are that half. They are stored
//! together because `coarse_clash` decides whether the other two can be
//! honoured at all -- a coarse bead is one sphere per residue, not an atom --
//! so a project carrying one without the other resumes into an exception.
//! 1.6 (2026-09-01) renamed the `radii_source` value `"model"` to `"imp"` and
//! made it the **default**, reversing 1.5's one-day default of `"olga"`. The
//! excluded-volume half of a docking score is `clash_container`, which reads
//! `IMP::core::XYZR` -- the particles' radii -- so a volume built on Olga's
//! table would leave the two halves of one score disagreeing about how big an
//! atom is. `"olga"` is unchanged and still selectable, and is what an
//! Olga-era reproduction asks for; `okf/validation/fps_screening_ab.md`
//! measures what the default costs against the published numbers.
//! 1.5 (2026-09-01) added `radii_source` on a position: which van der Waals
//! set the accessible volume inflates obstacles by. It is on the position
//! because the rest of a position's parameters -- above all a fitted
//! `contact_volume_trapped_fraction` -- only mean anything against the
//! obstacle set they were fitted over.
//! 1.4 (2026-09-01) added the top-level `Project` section: the structures in
//! body order, the labelling source, the selected distances, the five per-mode
//! parameter blocks, the AV and conversion globals and the current poses
//! (PRD-121 G1). It is what makes a run repeatable, and it is additive like
//! every step before it -- a file without one is still valid. 1.3 (2026-08-31)
//! added the `ATOM` simulation type: a labelling position
//! that is a plain atom of the structure with no accessible volume. It is what
//! FPS's LP files spell `ATOM` (`LabelingPositions.cs:203`), and a distance
//! whose two ends are both of them is FPS's **bond** (PRD-121 G4). 1.2
//! (2026-08-31) added `reference_atoms` on a position: the frame a fixed
//! (`XYZ`) coordinate was measured in, which screening fits onto each library
//! structure (PRD-121 G3). 1.1 (2026-08-27) added the document-level
//! `strip_mask` and `schema_version` fields. Every step is additive: every
//! 1.0 file is a valid 1.7 file -- except for a 1.5 file that spells this
//! field `"model"`, which is refused rather than reinterpreted (that spelling
//! existed for one day and no file outside this repository carries it).
std::string fps_schema_version() { return "1.7"; }

// Named, not anonymous: IMP compiles this module as one translation unit.
namespace fps_schema {

//! One row of a field table.
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
// A fourth dialect, and the only one that is not historical: the `Project`
// section is this module's own extension (PRD-121 G1, schema 1.4). It is named
// rather than folded into `flat` so a reader can tell "imp.bff wrote this
// template" from "imp.bff wrote this run".
const std::vector<std::string> PROJECT = {"project"};

std::string json_type_name(FPSFieldType t) {
    switch (t) {
        case FPS_STRING: return "string";
        case FPS_NUMBER: return "number";
        case FPS_INTEGER: return "integer";
        case FPS_BOOLEAN: return "boolean";
        case FPS_ARRAY: return "array";
        case FPS_OBJECT: return "object";
    }
    return "string";
}

//! The type name a validation message names the expected type by.
std::string python_type_name(FPSFieldType t) {
    switch (t) {
        case FPS_STRING: return "str";
        case FPS_NUMBER: return "float";
        case FPS_INTEGER: return "int";
        case FPS_BOOLEAN: return "bool";
        case FPS_ARRAY: return "list";
        case FPS_OBJECT: return "dict";
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
        // A bool where a number is expected is a bug, not a zero or a one,
        // and is reported as one.
        case FPS_BOOLEAN: return v.is_boolean();
        case FPS_INTEGER: return v.is_number_integer() && !v.is_boolean();
        // An int is acceptable where a float is expected.
        case FPS_NUMBER: return v.is_number() && !v.is_boolean();
        case FPS_ARRAY: return v.is_array();
        case FPS_OBJECT: return v.is_object();
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
    return {"AV1", "AV3", "XYZ", "R1", "ATOM"};
}

std::vector<std::string> fps_av_simulation_types() {
    // What ProbeNetworkRestraint can build a network out of. `XYZ` and `ATOM`
    // are points rather than volumes -- they carry no cloud -- but they are
    // scorable positions, which is what this list answers.
    return {"AV1", "AV3", "XYZ", "ATOM"};
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
                      "(fixed mean position, no simulation), ATOM (a plain "
                      "atom of the structure, no volume -- a crosslink or "
                      "anchor endpoint; a distance between two of them is a "
                      "bond), or R1 (rotamer ensemble: screened rotamer "
                      "library, Python-only). A distance with an XYZ or ATOM "
                      "end is scored as R_mp whatever its distance_type says, "
                      "because that end is a point and has no cloud to "
                      "average over (FPS: FilterEngine.cs:305-307).",
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
    // Negative means *derive it* from linker_width and the grid (see
    // IMP::bff::AV::get_effective_allowed_sphere_radius), and that is what an
    // absent key has always done. The "1.5" declared here before was never the
    // value an omitted key took, so the schema described a behaviour nothing
    // had. It cannot be spelled as "no default" -- that makes the key
    // *required* in this table -- so the sentinel is the default.
    f.push_back(field("allowed_sphere_radius", FPS_NUMBER, "-1.0",
                      "_flr_FPS_global_parameter.AV_allowed_sphere", "", FLAT));
    f.push_back(field("contact_volume_thickness", FPS_NUMBER, "0.0", "",
                      "Depth of the contact layer above the dye-excluded surface, "
                      "A: a cloud voxel is in contact when excluded volume lies "
                      "within it. Quantised to whole grid steps, so a thickness "
                      "below simulation_grid_resolution is no layer at all; "
                      "0 disables.",
                      FLAT));
    f.push_back(field("contact_volume_trapped_fraction", FPS_NUMBER, "-1.0", "",
                      "Share of the cloud's total weight carried by the contact "
                      "layer; must be below 1, and negative disables the "
                      "re-weighting. Fitted per site -- it does not come from "
                      "the structure.",
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
                      "Atoms removed as obstacles before this position's volume "
                      "is computed, as a selection expression: "
                      "'chain A and resid 132 and not name CA+CB+C+N+O', or the "
                      "same in the other spelling, 'chain A and resid 132 and "
                      "not name CA CB C N O'. The document's own strip_mask, if "
                      "it has one, is removed as well -- the two are a union, "
                      "not an override. Empty means the default strip: the "
                      "attachment residue's side chain minus the attachment "
                      "atom.",
                      FLAT));
    f.push_back(field("chain_weighting", FPS_BOOLEAN, "false", "",
                      "Weight AV grid points by linker-chain statistics instead "
                      "of uniformly.",
                      FLAT));
    // 1.5: the radii set the position's other parameters were calibrated
    // against. Named on the position, not on the run, for the same reason the
    // strip mask is -- a fitted contact_volume_trapped_fraction only means
    // anything against the obstacle set it was fitted over.
    // 1.6: "model" respelled "imp" and made the default, because the clash
    // term of a docking score reads the particles' radii and the volume must
    // agree with it. See AV::set_radii_source.
    f.push_back(field("radii_source", FPS_STRING, "\"imp\"", "",
                      "Van der Waals radii the accessible volume inflates "
                      "obstacles by. 'imp' (the default) is the radius each "
                      "particle carries, which after IMP's read_pdb is the "
                      "CHARMM-derived united-atom set with implicit hydrogens "
                      "(carbon 1.85-2.275 A) -- the same radii the clash term "
                      "of a docking score reads. 'olga' is Olga's name-keyed "
                      "table -- C 1.70, N 1.625, O 1.49, S 1.782, P 1.86, "
                      "H 1.00 A, and 1.50 A for an atom name it does not carry "
                      "-- and is what reproduces Olga-era numbers.",
                      FLAT, {"imp", "olga"}));
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
    // The frame `x`/`y`/`z` were measured in. Without it a fixed position is
    // a coordinate in someone else's box; with it, screening superimposes
    // these atoms on each library structure and carries the coordinate
    // along. The default is a real value ("[]") and not the empty string on
    // purpose: an empty default in this table means *required*.
    f.push_back(field("reference_atoms", FPS_ARRAY, "[]", "",
                      "Atoms defining the frame the fixed (XYZ) coordinate "
                      "was measured in, as objects: {chain_identifier, "
                      "residue_seq_number, atom_name, x, y, z} and optionally "
                      "atom_serial, where x/y/z are the atom's coordinates in "
                      "that frame. Screening Kabsch-fits them onto each "
                      "structure, transports x/y/z of the position through "
                      "the fit and reports the fit RMSD (FPS: RefRMSD). "
                      "Three atoms are the minimum that fixes a rotation.",
                      {"flat", "csfps"}));
    return f;
}

std::vector<FPSField> fps_distance_fields() {
    using namespace fps_schema;
    std::vector<std::string> types;
    const std::map<std::string, std::string> dt = fps_distance_types();
    // Sorted, which for these three is also the order they are documented in:
    // RDAMean, RDAMeanE, Rmp.
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

std::string combined_strip_mask(const std::string& document_mask,
                                const std::string& position_mask) {
    const std::string a = internal::trimmed(document_mask);
    const std::string b = internal::trimmed(position_mask);
    if (a.empty()) return b;
    if (b.empty()) return a;
    // Parenthesised, because either half may be an `or`-chain of its own and
    // `x or y and z` does not mean `(x or y) and z`.
    return "(" + a + ") or (" + b + ")";
}

std::vector<FPSField> fps_document_fields() {
    using namespace fps_schema;
    std::vector<FPSField> f;
    f.push_back(field("strip_mask", FPS_STRING, "\"\"", "",
                      "Atoms removed as obstacles for **every** position in the "
                      "file, as a selection expression -- what the structure "
                      "carries that no volume should be blocked by: "
                      "'resname HOH SOL WAT NA CL', 'not polymer', a "
                      "crystallisation additive. A position's own strip_mask "
                      "adds to this rather than replacing it, because the two "
                      "answer different questions: this one is about the "
                      "structure, that one about the labelling site.",
                      FLAT));
    f.push_back(field("schema_version", FPS_STRING, "\"\"", "",
                      "The fps.json schema this file was written against. "
                      "Absent means a file written before the schema was "
                      "versioned.",
                      FLAT));
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

std::vector<FPSField> fps_mode_parameter_fields() {
    using namespace fps_schema;
    std::vector<FPSField> f;
    // FPS's `FPSParameters` (`SimulationBase.cs:11`), defaults from
    // `ProjectData()` Dock block (`ProjectData.cs:67-81`). Every field carries
    // a real default: an empty `default_json` is how a field is declared to
    // have none, and a project is meant to be writable one key at a time.
    f.push_back(field("viscosity_factor", FPS_NUMBER, "1.0", "",
                      "FPS integrator damping. Stored, not obeyed: this port "
                      "uses IMP's optimisers, not FPS's damped Verlet "
                      "integrator (PRD-121 D1).",
                      PROJECT));
    f.push_back(field("time_step_factor", FPS_NUMBER, "1.0", "",
                      "FPS integrator timestep scale. Stored, not obeyed.",
                      PROJECT));
    f.push_back(field("max_iterations", FPS_INTEGER, "200000", "",
                      "FPS's integrator step budget. NOT this port's budget -- "
                      "200000 damped Verlet steps and 500 conjugate-gradient "
                      "iterations are different quantities; see 'iterations'.",
                      PROJECT));
    f.push_back(field("max_force", FPS_NUMBER, "400.0", "",
                      "Past MaxForce*err^2/2 a restraint goes linear instead "
                      "of parabolic. Honoured. FPS ships 400 for "
                      "Dock/Sample/Screening and 10000 for "
                      "Refine/Error estimation.",
                      PROJECT));
    f.push_back(field("clash_tolerance", FPS_NUMBER, "1.0", "",
                      "The overlap in angstrom that costs one chi2 unit; "
                      "k_clash = 2/tolerance^2. Honoured.",
                      PROJECT));
    f.push_back(field("rkt", FPS_NUMBER, "10.0", "",
                      "FPS's rkT, i.e. 1/kT. Sample mode is the only reader.",
                      PROJECT));
    f.push_back(field("e_tolerance", FPS_NUMBER, "100.0", "",
                      "FPS's ETolerance, which is dead in all five of its own "
                      "modes. Stored for provenance.",
                      PROJECT));
    f.push_back(field("k_tolerance", FPS_NUMBER, "0.001", "",
                      "FPS convergence tolerance. Stored, not obeyed.",
                      PROJECT));
    f.push_back(field("f_tolerance", FPS_NUMBER, "0.001", "",
                      "FPS force tolerance. Stored, not obeyed.", PROJECT));
    f.push_back(field("t_tolerance", FPS_NUMBER, "0.02", "",
                      "FPS torque tolerance. Stored, not obeyed.", PROJECT));
    f.push_back(field("optimize_selected", FPS_STRING, "\"Selected\"", "",
                      "Which DISTANCES contribute: the score set, every "
                      "distance in the file, or the score set then everything. "
                      "Clashes are never gated.",
                      PROJECT, {"Selected", "All", "SelectedThenAll"}));
    f.push_back(field("iterations", FPS_INTEGER, "500", "",
                      "This port's own iteration budget for the mode.",
                      PROJECT));
    return f;
}

std::vector<FPSField> fps_conversion_fields() {
    using namespace fps_schema;
    std::vector<FPSField> f;
    f.push_back(field("forster_radius", FPS_NUMBER, "52.0", "",
                      "R0 the distance-conversion polynomial is fitted at, A "
                      "(FPS ConversionParameters.R0).",
                      PROJECT));
    f.push_back(field("polynomial_order", FPS_INTEGER, "3", "",
                      "Order of the conversion polynomial "
                      "(FPS ConversionParameters.PolynomOrder).",
                      PROJECT));
    return f;
}

std::vector<FPSField> fps_av_global_fields() {
    using namespace fps_schema;
    std::vector<FPSField> f;
    f.push_back(field("grid_size", FPS_NUMBER, "0.2", "",
                      "Relative grid spacing: dg = max(min(0.2L, 0.2W, "
                      "0.4R1, 0.4R2, 0.4R3), min_grid_size).",
                      PROJECT));
    f.push_back(field("min_grid_size", FPS_NUMBER, "0.4", "",
                      "Floor of the derived grid spacing, A.", PROJECT));
    f.push_back(field("linker_initial_sphere", FPS_NUMBER, "0.5", "",
                      "FPS's LinkerInitialSphere. NOT this module's derived "
                      "clearance -- the two are different quantities, which is "
                      "the trap that made a 4.5 A linker width return a "
                      "silently empty volume (PRD-121 G8). Stored for "
                      "provenance.",
                      PROJECT));
    f.push_back(field("link_search_nodes", FPS_INTEGER, "3", "",
                      "How many voxels FPS's link search hops at a time.",
                      PROJECT));
    f.push_back(field("e_samples", FPS_INTEGER, "200000", "",
                      "Samples FPS draws for <R_DA>_E.", PROJECT));
    return f;
}

std::vector<FPSField> fps_project_fields() {
    using namespace fps_schema;
    std::vector<FPSField> f;
    f.push_back(field("schema_version", FPS_STRING, "\"\"", "",
                      "The schema this project was written against.",
                      PROJECT));
    f.push_back(field("mode", FPS_STRING, "\"None\"", "",
                      "FPS's ProjectFPSMode: everything that moves bodies is "
                      "Dock, screening is Filter.",
                      PROJECT, {"None", "Dock", "Filter"}));
    f.push_back(field("source", FPS_STRING, "\"\"", "",
                      "Provenance: the banner line of the legacy text export "
                      "this project was converted from, when it was.",
                      PROJECT));
    f.push_back(field("structures", FPS_ARRAY, "[]", "",
                      "One structure path per rigid body, IN BODY ORDER. A "
                      "relative path is resolved against the project file.",
                      PROJECT));
    f.push_back(field("labelling_json", FPS_STRING, "\"\"", "",
                      "The labelling source as an fps.json path. Empty, with "
                      "positions_path also empty, means this file carries its "
                      "own Positions and Distances.",
                      PROJECT));
    f.push_back(field("positions_path", FPS_STRING, "\"\"", "",
                      "The labelling source as a legacy C# FPS positions .txt "
                      "(FPS's LabelingPositionsPath).",
                      PROJECT));
    f.push_back(field("distances_path", FPS_STRING, "\"\"", "",
                      "The legacy distances .txt beside it (FPS's "
                      "DistancesPath).",
                      PROJECT));
    f.push_back(field("score_set", FPS_STRING, "\"\"", "",
                      "The named chi2 score set this run scores; empty is "
                      "every distance.",
                      PROJECT));
    f.push_back(field("selected_distances", FPS_ARRAY, "[]", "",
                      "Selected distances BY NAME. Empty defers to score_set.",
                      PROJECT));
    f.push_back(field("selected_flags", FPS_ARRAY, "[]", "",
                      "FPS's SelectedDistances: a Boolean[] parallel to the "
                      "distances file, so positional rather than named. Only "
                      "the legacy reader fills it; kept so a converted project "
                      "is not silently re-interpreted.",
                      PROJECT));
    f.push_back(field("fixed_body", FPS_INTEGER, "0", "",
                      "The body_id held still while the others move.",
                      PROJECT));
    f.push_back(field("ev_weight", FPS_NUMBER, "1.0", "",
                      "Weight of the excluded-volume (clash) term.", PROJECT));
    f.push_back(field("sigma_da", FPS_NUMBER, "6.0", "",
                      "Width of the mean-position transfer function, A.",
                      PROJECT));
    f.push_back(field("shuffle", FPS_NUMBER, "10.0", "",
                      "Initial random displacement of the mobile bodies, A; "
                      "0 starts from the input pose.",
                      PROJECT));
    f.push_back(field("refine_av_cycles", FPS_INTEGER, "0", "",
                      "FPS-style refinement cycles after a dock: re-sample the "
                      "volumes in the docked context and minimise again.",
                      PROJECT));
    f.push_back(field("coarse_clash", FPS_BOOLEAN, "true", "",
                      "Clash detection on one coarse bead per residue rather "
                      "than all atoms. Minimisation only, about three times "
                      "cheaper per step. Must be false for clash_radii_source "
                      "or clash_radii_scale to be honoured.",
                      PROJECT));
    f.push_back(field("clash_radii_source", FPS_STRING, "\"imp\"", "",
                      "Which van der Waals radii the EXCLUDED-VOLUME term "
                      "measures overlap by: 'imp' the particles' own "
                      "united-atom radii, 'olga' the Bondi-scale table FPS's "
                      "ClashTolerance was calibrated against. The position-"
                      "level radii_source is the volume's answer to the same "
                      "question; this is the clash term's.",
                      PROJECT, {"imp", "olga"}));
    f.push_back(field("clash_radii_scale", FPS_NUMBER, "1.0", "",
                      "Multiplies whatever clash_radii_source gives. 1.0 "
                      "changes nothing. It cannot stand in for the source: on "
                      "HIV-RT the scale reproducing Olga's pair count is 0.82, "
                      "its total overlap 0.83 and its energy 0.85.",
                      PROJECT));
    f.push_back(field("poses", FPS_ARRAY, "[]", "",
                      "The current rigid-body poses, as capture_poses writes "
                      "them: one {body_id, t, q} per body. [] means the input "
                      "structures as they stand; anything else is a state a "
                      "run reached, so storing the project stores the run.",
                      PROJECT));
    f.push_back(field("parameters", FPS_OBJECT, "{}", "",
                      "The five per-mode parameter blocks, keyed by FPS's own "
                      "mode names: Dock, Refine, Error estimation, Sample, "
                      "Screening.",
                      PROJECT));
    f.push_back(field("conversion", FPS_OBJECT, "{}", "",
                      "FPS's ConversionParameters.", PROJECT));
    f.push_back(field("av", FPS_OBJECT, "{}", "",
                      "FPS's AVGlobalParameters.", PROJECT));
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
    // The Project section. Its three object-valued fields carry only
    // `"type": "object"` from the field table, so the real sub-schemas are
    // spliced in here -- one place, so `parameters` in `$defs` and
    // `parameters` under `properties` cannot describe different things.
    nlohmann::json project = fps_schema::object_schema(fps_project_fields());
    const nlohmann::json mode_parameters =
            fps_schema::object_schema(fps_mode_parameter_fields());
    project["properties"]["parameters"] = {
            {"type", "object"}, {"additionalProperties", mode_parameters}};
    project["properties"]["conversion"] =
            fps_schema::object_schema(fps_conversion_fields());
    project["properties"]["av"] =
            fps_schema::object_schema(fps_av_global_fields());
    properties["Project"] = project;
    // Document-level fields sit beside the three sections rather than inside
    // one of them: they are statements about the structure, not about a
    // position or a distance.
    const std::vector<FPSField> document = fps_document_fields();
    for (std::size_t i = 0; i < document.size(); ++i) {
        properties[document[i].name] = fps_schema::field_property(document[i]);
    }
    out["properties"] = properties;
    out["additionalProperties"] = true;
    nlohmann::json defs;
    defs["position"] = position;
    defs["distance"] = distance;
    defs["score_set"] = score_set;
    defs["project"] = project;
    defs["mode_parameters"] = mode_parameters;
    out["$defs"] = defs;
    return out.dump(2);
}

namespace fps_schema {

void validate_position_json(const nlohmann::json& position,
                            const std::string& name, FPSValidation& out) {
    check_fields(position, fps_position_fields(), name, out);
    if (!position.is_object()) return;

    // A reference atom is a coordinate plus a way of finding the same atom in
    // another structure; either half alone is useless, so both are required
    // of each element rather than of the array.
    if (position.contains("reference_atoms") &&
        position["reference_atoms"].is_array()) {
        const nlohmann::json& atoms = position["reference_atoms"];
        for (std::size_t i = 0; i < atoms.size(); ++i) {
            const std::string where =
                    name + ".reference_atoms[" + std::to_string(i) + "]";
            if (!atoms[i].is_object()) {
                out.errors.push_back(where + ": expected an object, got " +
                                     value_type_name(atoms[i]));
                continue;
            }
            const char* coords[] = {"x", "y", "z"};
            for (int k = 0; k < 3; ++k) {
                if (!atoms[i].contains(coords[k]) ||
                    !type_matches(atoms[i][coords[k]], FPS_NUMBER)) {
                    out.errors.push_back(where + ": missing or non-numeric '" +
                                         std::string(coords[k]) + "'");
                }
            }
            const bool by_site =
                    atoms[i].contains("residue_seq_number") &&
                    atoms[i].contains("atom_name") &&
                    atoms[i]["atom_name"].is_string() &&
                    !atoms[i]["atom_name"].get<std::string>().empty();
            const bool by_serial =
                    atoms[i].contains("atom_serial") &&
                    type_matches(atoms[i]["atom_serial"], FPS_INTEGER) &&
                    atoms[i]["atom_serial"].get<int>() > 0;
            if (!by_site && !by_serial) {
                out.errors.push_back(
                        where + ": needs 'residue_seq_number' with "
                                "'atom_name', or a positive 'atom_serial', to "
                                "be found in another structure");
            }
        }
    }

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
                        "ProbeNetworkRestraint -- filter with "
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
        // A frame of one or two atoms does not fix a rotation, so the point it
        // transports is arbitrary rather than merely uncertain. Warn rather
        // than error: the file is well-formed, the frame is not usable.
        if (position.contains("reference_atoms") &&
            position["reference_atoms"].is_array() &&
            position["reference_atoms"].size() > 0 &&
            position["reference_atoms"].size() < 3) {
            out.warnings.push_back(
                    name + ": 'reference_atoms' has " +
                    std::to_string(position["reference_atoms"].size()) +
                    " atoms; three are needed to fix a rotation, so this "
                    "frame will not be fitted");
        }
    } else if (stype == "ATOM") {
        // An ATOM position *is* an atom of the structure, so it needs a way of
        // naming one; unlike XYZ it carries no coordinate of its own.
        const bool named = position.contains("residue_seq_number") &&
                           type_matches(position["residue_seq_number"],
                                        FPS_INTEGER) &&
                           position["residue_seq_number"].get<int>() > 0;
        if (!named) {
            out.errors.push_back(
                    name + ": simulation_type ATOM requires a positive "
                           "'residue_seq_number' -- it names an atom of the "
                           "structure rather than carrying a coordinate");
        }
        const char* av_keys[] = {"linker_length", "linker_width", "radius1",
                                 "radius2", "radius3"};
        for (int i = 0; i < 5; ++i) {
            if (position.contains(av_keys[i])) {
                out.warnings.push_back(
                        name + ": '" + av_keys[i] +
                        "' is an AV parameter; an ATOM position has no volume "
                        "and ignores it");
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

//! One `Project` section. \p distance_names cross-checks selected_distances.
void validate_project_json(const nlohmann::json& project,
                           const std::string& name,
                           const std::set<std::string>* distance_names,
                           FPSValidation& out) {
    check_fields(project, fps_project_fields(), name, out);
    if (!project.is_object()) return;

    // The three nested blocks. `check_fields` above has already established
    // that each is an object where present, so a non-object was reported once
    // and is not reported again here.
    if (project.contains("parameters") && project["parameters"].is_object()) {
        const nlohmann::json& parameters = project["parameters"];
        const std::vector<FPSField> block = fps_mode_parameter_fields();
        for (nlohmann::json::const_iterator it = parameters.begin();
             it != parameters.end(); ++it) {
            check_fields(it.value(), block,
                         name + ".parameters['" + it.key() + "']", out);
        }
    }
    if (project.contains("conversion") && project["conversion"].is_object()) {
        check_fields(project["conversion"], fps_conversion_fields(),
                     name + ".conversion", out);
    }
    if (project.contains("av") && project["av"].is_object()) {
        check_fields(project["av"], fps_av_global_fields(), name + ".av", out);
    }

    if (project.contains("structures") && project["structures"].is_array()) {
        const nlohmann::json& structures = project["structures"];
        for (std::size_t i = 0; i < structures.size(); ++i) {
            if (!structures[i].is_string()) {
                out.errors.push_back(name + ".structures[" +
                                     std::to_string(i) + "]: expected a path, "
                                     "got " + value_type_name(structures[i]));
            }
        }
    }

    // A selection naming a distance the file does not define is the failure
    // mode that costs a whole run: the score set silently shrinks and the
    // number is over different data than the one it is compared against.
    if (distance_names != nullptr && !distance_names->empty() &&
        project.contains("selected_distances") &&
        project["selected_distances"].is_array()) {
        const nlohmann::json& selected = project["selected_distances"];
        for (std::size_t i = 0; i < selected.size(); ++i) {
            if (!selected[i].is_string()) continue;
            const std::string ref = selected[i].get<std::string>();
            if (distance_names->find(ref) == distance_names->end()) {
                out.errors.push_back(name + ".selected_distances: '" + ref +
                                     "' is not a defined distance");
            }
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

FPSValidation validate_project(const std::string& project_json,
                               const std::string& name) {
    FPSValidation out;
    fps_schema::validate_project_json(
            fps_schema::parse(project_json, "a project"), name, nullptr, out);
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

    if (payload.contains("Project")) {
        std::set<std::string> distance_names;
        for (nlohmann::json::const_iterator it = distances.begin();
             it != distances.end(); ++it) {
            distance_names.insert(it.key());
        }
        // A project that points at an external labelling file carries no
        // Distances of its own, and cross-checking a selection against an
        // empty set would reject every name in it.
        fps_schema::validate_project_json(
                payload["Project"], "Project",
                distance_names.empty() ? nullptr : &distance_names, out);
    }

    if (payload.contains("Evaluators") && !payload["Evaluators"].is_null() &&
        !payload["Evaluators"].is_array()) {
        out.errors.push_back("Evaluators: expected an array");
    }
    return out;
}

IMPBFF_END_NAMESPACE
