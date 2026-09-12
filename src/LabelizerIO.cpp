/**
 * \file LabelizerIO.cpp
 * \brief A scored structure as one `.mmfdb.pto` container.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/LabelizerIO.h>
#include <fstream>
#include <sstream>
#include <limits>
#include <IMP/bff/internal/json.h>
#include <IMP/bff/Pto.h>

IMPBFF_BEGIN_NAMESPACE

const char* const LABELIZER_PTO_README = "README";
const char* const LABELIZER_PTO_STRUCTURE = "structure.pdb";
const char* const LABELIZER_PTO_SCORES = "label_scores.json";
const char* const LABELIZER_PTO_PAIRS = "label_pairs.json";
const char* const LABELIZER_PTO_MODEL = "label_model.json";

namespace {

//! The README that goes in first, so the file explains itself.
/*! Not documentation *about* the format: the format telling a reader what it
    is, in the file, in language that survives this library not building. */
const char* const LABELIZER_README_TEXT =
    "PTO.MFDB label-site container\n"
    "=============================\n"
    "\n"
    "This is an EBML document (RFC 8794) with DocType \"pto\". Every payload is\n"
    "an AttachedFile carrying a FileName, a PtoKind saying what it is, a\n"
    "PtoEncoding saying how the bytes are coded, and the bytes. Nothing in it\n"
    "is compressed, encrypted, or stored outside the file.\n"
    "\n"
    "To walk it by hand: an EBML element is an id, then a Data Size, then the\n"
    "payload; both the id and the size are variable-length integers whose\n"
    "leading zero bits give the byte count. Attachments (0x1941A469) holds\n"
    "AttachedFile (0x61A7) elements; inside one, FileName is 0x466E, FileData\n"
    "is 0x465C, and the two custom ids PtoKind (0x1E54F001) and PtoEncoding\n"
    "(0x1E54F002) say what the payload means.\n"
    "\n"
    "Objects in this container:\n"
    "  README             this text\n"
    "  structure.pdb      the structure that was scored, byte for byte as it\n"
    "                     was read. Write its FileData to a file and check it\n"
    "                     against the SHA-256 in its checksum tag.\n"
    "  label_scores.json  one row per (position, score_type). A row with no\n"
    "                     \"value\" was NOT computed -- its \"status\" says why.\n"
    "                     Absence is information; there are no sentinels.\n"
    "  label_pairs.json   one row per pair of positions, when pairs were\n"
    "                     scored.\n"
    "  label_model.json   the complete settings the run used.\n"
    "\n"
    "Positions are named the way mmCIF names a position: asym_id (the chain,\n"
    "as the source file spells it -- an author chain id, not an assembly one)\n"
    "and seq_id. Every controlled value is a term from the MMFDB dictionary\n"
    "mmfdb_flr_ext.dic; the container tags record which version.\n"
    "\n"
    "Scores are likelihood ratios, not probabilities: they run from 0 to about\n"
    "4 and the combined score is unbounded above.\n";

std::string labelizer_read_file(const std::string& path) {
    std::ifstream in(path.c_str(), std::ios::binary);
    if (!in) {
        IMP_THROW("labelizer_write_pto: cannot read " << path, IOException);
    }
    std::ostringstream os;
    os << in.rdbuf();
    return os.str();
}

void labelizer_add(PtoWriter& w, const std::string& name, const std::string& kind,
            const std::string& encoding, const std::string& bytes) {
    w.add(name, kind, encoding, bytes.data(), bytes.size());
}

//! The tags of one object, flattened into the JSON its payload carries.
/*! PtoWriter frames bytes and carries a kind and an encoding; it has no tag
    element. Rather than add one -- which would change the container format,
    which this profile does not do -- the profile tags of an object ride in a
    `_tags` member of the object's own JSON. A payload that is not JSON keeps
    its tags in the model object instead. */
nlohmann::json labelizer_tags_json(const std::vector<MfdbTag>& tags) {
    nlohmann::json j = nlohmann::json::object();
    for (std::size_t i = 0; i < tags.size(); ++i) j[tags[i].item] = tags[i].value;
    return j;
}

nlohmann::json labelizer_score_row(const LabelizerScore& s) {
    nlohmann::json r;
    r["_mmfdb_label_score.asym_id"] = s.asym_id;
    r["_mmfdb_label_score.seq_id"] = s.seq_id;
    r["_mmfdb_label_score.comp_id"] = s.comp_id;
    r["_mmfdb_label_score.score_type"] = s.score_type;
    r["_mmfdb_label_score.definition"] = "labelizer";
    r["_mmfdb_label_score.status"] = s.status;
    // Absent means not computed. A sentinel here is exactly what the
    // dictionary forbids, and what makes the reference's CSVs unreadable.
    if (s.status == "scored") {
        r["_mmfdb_label_score.value"] = s.value;
        r["_mmfdb_label_score.units"] = "dimensionless";
    }
    return r;
}

std::vector<MfdbColumn> labelizer_score_columns() {
    std::vector<MfdbColumn> c;
    c.push_back(MfdbColumn("_mmfdb_label_score.asym_id", "",
                           "_mmfdb_label_score.asym_id",
                           "Chain, as the source file spells it (author id)."));
    c.push_back(MfdbColumn("_mmfdb_label_score.seq_id", "",
                           "_mmfdb_label_score.seq_id", "Residue number."));
    c.push_back(MfdbColumn("_mmfdb_label_score.comp_id", "",
                           "_mmfdb_label_score.comp_id", "Residue name."));
    c.push_back(MfdbColumn("_mmfdb_label_score.score_type", "",
                           "_mmfdb_label_score.score_type",
                           "What is scored; a dictionary term."));
    c.push_back(MfdbColumn("_mmfdb_label_score.definition", "",
                           "_mmfdb_label_score.definition",
                           "Whose definition of the score_type this follows."));
    c.push_back(MfdbColumn("_mmfdb_label_score.value", "dimensionless",
                           "_mmfdb_label_score.value",
                           "A likelihood ratio, unbounded above. Absent when "
                           "the score was not computed."));
    c.push_back(MfdbColumn("_mmfdb_label_score.status", "",
                           "_mmfdb_label_score.status",
                           "Why a position has no value, when it has none."));
    return c;
}

std::vector<MfdbColumn> labelizer_pair_columns() {
    std::vector<MfdbColumn> c;
    c.push_back(MfdbColumn("asym_id_1", "", "_flr_poly_probe_position.asym_id",
                           "Chain of the first position."));
    c.push_back(MfdbColumn("seq_id_1", "", "_flr_poly_probe_position.seq_id",
                           "Residue number of the first position."));
    c.push_back(MfdbColumn("asym_id_2", "", "_flr_poly_probe_position.asym_id",
                           "Chain of the second position."));
    c.push_back(MfdbColumn("seq_id_2", "", "_flr_poly_probe_position.seq_id",
                           "Residue number of the second position."));
    c.push_back(MfdbColumn("value", "dimensionless", "",
                           "The FRET pair score."));
    c.push_back(MfdbColumn("distance", "angstroms",
                           "_flr_fret_model_distance.distance",
                           "Probe-dye distance; the first conformation."));
    c.push_back(MfdbColumn("distance_2", "angstroms",
                           "_flr_fret_model_distance.distance",
                           "The second conformation's distance, when there "
                           "is one."));
    c.push_back(MfdbColumn("joined_label_score", "dimensionless", "",
                           "The combined label score of the two positions."));
    c.push_back(MfdbColumn("probe_model", "", "",
                           "How the dye position was obtained: cbeta, "
                           "alpha_cone or accessible_volume."));
    return c;
}

const char* labelizer_dye_model_name(ProbeModel m) {
    switch (m) {
        case PROBE_MODEL_CBETA: return "cbeta";
        case PROBE_MODEL_ALPHA_CONE: return "alpha_cone";
        case PROBE_MODEL_ACCESSIBLE_VOLUME: return "accessible_volume";
    }
    return "unknown";
}

std::string labelizer_object_text(const std::string& path, const std::string& name) {
    PtoReader r(path);
    const int i = r.find(name);
    if (i < 0) return "";
    const std::vector<unsigned char> d = r.data(r.objects()[i]);
    return std::string(d.begin(), d.end());
}

}  // namespace

// ---------------------------------------------------------------------------
// Writing
// ---------------------------------------------------------------------------

void labelizer_write_pto(const std::string& path, const std::string& pdb_path,
                  const std::vector<LabelizerScore>& scores,
                  const std::vector<LabelizerFRETPairScore>& pairs,
                  const std::string& settings_json) {
    const std::string structure = labelizer_read_file(pdb_path);
    const std::string structure_sum = mfdb_checksum(structure);

    // The scores, tidy: one row per (position, score_type).
    nlohmann::json score_rows = nlohmann::json::array();
    for (std::size_t i = 0; i < scores.size(); ++i) {
        score_rows.push_back(labelizer_score_row(scores[i]));
    }
    nlohmann::json score_doc;
    score_doc["rows"] = score_rows;

    nlohmann::json pair_rows = nlohmann::json::array();
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        const LabelizerFRETPairScore& p = pairs[i];
        nlohmann::json r;
        r["asym_id_1"] = p.asym_id_1;
        r["seq_id_1"] = p.seq_id_1;
        r["asym_id_2"] = p.asym_id_2;
        r["seq_id_2"] = p.seq_id_2;
        r["value"] = p.value;
        r["distance"] = p.distance;
        // NaN is not JSON; a single-conformation pair simply has no second
        // distance, and omitting it says that.
        if (p.distance_2 == p.distance_2) r["distance_2"] = p.distance_2;
        r["joined_label_score"] = p.joined_label_score;
        r["probe_model"] = labelizer_dye_model_name(p.probe_model);
        pair_rows.push_back(r);
    }
    nlohmann::json pair_doc;
    pair_doc["rows"] = pair_rows;

    // The model object carries the settings, the container-level conformance
    // tags, the operation, and the provenance edges -- everything that is
    // about the file rather than about a row.
    nlohmann::json model_doc;
    model_doc["_container"] = labelizer_tags_json(mfdb_container_tags());
    model_doc["_operation"] = labelizer_tags_json(mfdb_operation_tags(
            "analysis", "labelizer", settings_json, "IMP.bff",
            get_module_version()));
    model_doc["settings"] = nlohmann::json::parse(settings_json);

    nlohmann::json edges = nlohmann::json::array();
    edges.push_back(labelizer_tags_json(mfdb_edge_tags(
            LABELIZER_PTO_STRUCTURE, LABELIZER_PTO_SCORES, "derived_from")));
    if (!pairs.empty()) {
        // The pair table is coarser than the score table -- one row per two
        // positions -- so the edge names the columns that join them rather
        // than leaving it to position.
        edges.push_back(labelizer_tags_json(mfdb_edge_tags(
                LABELIZER_PTO_SCORES, LABELIZER_PTO_PAIRS, "maps_rows_of",
                "_mmfdb_label_score.seq_id", "seq_id_1")));
    }
    model_doc["_edges"] = edges;

    const std::string score_bytes = score_doc.dump(1);
    const std::string pair_bytes = pair_doc.dump(1);

    model_doc["_artifacts"] = nlohmann::json::object();
    model_doc["_artifacts"][LABELIZER_PTO_STRUCTURE] = labelizer_tags_json(mfdb_artifact_tags(
            LABELIZER_PTO_STRUCTURE, "processed_data", "text", "", -1, structure_sum));
    model_doc["_artifacts"][LABELIZER_PTO_SCORES] = labelizer_tags_json(mfdb_artifact_tags(
            LABELIZER_PTO_SCORES, "parameter_table", "json", "label_site",
            static_cast<long>(scores.size()), mfdb_checksum(score_bytes),
            labelizer_score_columns()));
    if (!pairs.empty()) {
        model_doc["_artifacts"][LABELIZER_PTO_PAIRS] = labelizer_tags_json(mfdb_artifact_tags(
                LABELIZER_PTO_PAIRS, "analysis_result", "json", "pair",
                static_cast<long>(pairs.size()), mfdb_checksum(pair_bytes),
                labelizer_pair_columns()));
    }
    const std::string model_bytes = model_doc.dump(1);

    PtoWriter w(path);
    labelizer_add(w, LABELIZER_PTO_README, "readme", "text", LABELIZER_README_TEXT);
    labelizer_add(w, LABELIZER_PTO_STRUCTURE, "label.structure", "text", structure);
    labelizer_add(w, LABELIZER_PTO_SCORES, "label.scores", "json", score_bytes);
    if (!pairs.empty()) {
        labelizer_add(w, LABELIZER_PTO_PAIRS, "label.pairs", "json", pair_bytes);
    }
    labelizer_add(w, LABELIZER_PTO_MODEL, "label.model", "json", model_bytes);
    w.close();
}

// ---------------------------------------------------------------------------
// Reading
// ---------------------------------------------------------------------------

std::vector<LabelizerScore> labelizer_read_pto_scores(const std::string& path) {
    const std::string text = labelizer_object_text(path, LABELIZER_PTO_SCORES);
    if (text.empty()) {
        IMP_THROW("labelizer_read_pto_scores: " << path << " carries no "
                  << LABELIZER_PTO_SCORES, IOException);
    }
    const nlohmann::json j = nlohmann::json::parse(text);
    std::vector<LabelizerScore> out;
    const nlohmann::json& rows = j.at("rows");
    for (nlohmann::json::const_iterator it = rows.begin(); it != rows.end();
         ++it) {
        LabelizerScore s;
        s.asym_id = it->at("_mmfdb_label_score.asym_id").get<std::string>();
        s.seq_id = it->at("_mmfdb_label_score.seq_id").get<int>();
        s.comp_id = it->at("_mmfdb_label_score.comp_id").get<std::string>();
        s.score_type = it->at("_mmfdb_label_score.score_type").get<std::string>();
        s.status = it->at("_mmfdb_label_score.status").get<std::string>();
        if (it->count("_mmfdb_label_score.value")) {
            s.value = it->at("_mmfdb_label_score.value").get<double>();
        }
        out.push_back(s);
    }
    return out;
}

std::vector<LabelizerFRETPairScore> labelizer_read_pto_pairs(const std::string& path) {
    const std::string text = labelizer_object_text(path, LABELIZER_PTO_PAIRS);
    std::vector<LabelizerFRETPairScore> out;
    if (text.empty()) return out;
    const nlohmann::json j = nlohmann::json::parse(text);
    const nlohmann::json& rows = j.at("rows");
    for (nlohmann::json::const_iterator it = rows.begin(); it != rows.end();
         ++it) {
        LabelizerFRETPairScore p;
        p.asym_id_1 = it->at("asym_id_1").get<std::string>();
        p.seq_id_1 = it->at("seq_id_1").get<int>();
        p.asym_id_2 = it->at("asym_id_2").get<std::string>();
        p.seq_id_2 = it->at("seq_id_2").get<int>();
        p.value = it->at("value").get<double>();
        p.distance = it->at("distance").get<double>();
        // The writer omits `distance_2` when there is no second conformation
        // (a NaN), so its absence is meaningful and must come back as a NaN.
        // Reading it at all is new: the field was written and never read, so
        // every two-state pair round-tripped through a container came back
        // with `distance_2 == 0`, which reads as a real distance of zero --
        // efficiency 1, and a plausible-looking score built on it.
        p.distance_2 = it->count("distance_2")
                               ? it->at("distance_2").get<double>()
                               : std::numeric_limits<double>::quiet_NaN();
        p.joined_label_score = it->at("joined_label_score").get<double>();
        const std::string m = it->at("probe_model").get<std::string>();
        p.probe_model = m == "accessible_volume" ? PROBE_MODEL_ACCESSIBLE_VOLUME
                    : m == "alpha_cone"        ? PROBE_MODEL_ALPHA_CONE
                                               : PROBE_MODEL_CBETA;
        out.push_back(p);
    }
    return out;
}

std::string labelizer_extract_pto_structure(const std::string& path,
                                     const std::string& out_pdb_path) {
    PtoReader r(path);
    const int i = r.find(LABELIZER_PTO_STRUCTURE);
    if (i < 0) {
        IMP_THROW("labelizer_extract_pto_structure: " << path
                  << " carries no structure", IOException);
    }
    const std::vector<unsigned char> d = r.data(r.objects()[i]);
    const std::string bytes(d.begin(), d.end());
    const std::string got = mfdb_checksum(bytes);

    // What the container says it should be. Verification is explicit: opening
    // a file never hashes a payload, extracting one always does.
    const std::string model = labelizer_object_text(path, LABELIZER_PTO_MODEL);
    if (!model.empty()) {
        const nlohmann::json j = nlohmann::json::parse(model);
        if (j.count("_artifacts") &&
            j["_artifacts"].count(LABELIZER_PTO_STRUCTURE)) {
            const nlohmann::json& a = j["_artifacts"][LABELIZER_PTO_STRUCTURE];
            if (a.count("_mmfdb_artifact.checksum")) {
                const std::string want =
                        a["_mmfdb_artifact.checksum"].get<std::string>();
                if (want != got) {
                    IMP_THROW("labelizer_extract_pto_structure: " << path
                              << " does not hold the bytes it says it does: "
                              << "recorded " << want << ", found " << got,
                              IOException);
                }
            }
        }
    }
    std::ofstream out(out_pdb_path.c_str(), std::ios::binary);
    if (!out) {
        IMP_THROW("labelizer_extract_pto_structure: cannot write " << out_pdb_path,
                  IOException);
    }
    out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    return got;
}

std::string labelizer_read_pto_settings(const std::string& path) {
    const std::string model = labelizer_object_text(path, LABELIZER_PTO_MODEL);
    if (model.empty()) return "";
    const nlohmann::json j = nlohmann::json::parse(model);
    return j.count("settings") ? j["settings"].dump() : std::string();
}

std::string labelizer_settings_json(const std::vector<LabelizerParameter>& model,
                             const LabelizerOptions& options,
                             const LabelizerFRETOptions& fret_options,
                             const std::string& conservation_path) {
    nlohmann::json j;
    nlohmann::json terms = nlohmann::json::array();
    for (std::size_t i = 0; i < model.size(); ++i) {
        nlohmann::json t;
        t["tag"] = model[i].tag;
        t["score_type"] = labelizer_score_type(model[i].tag);
        t["table"] = model[i].table;
        t["weight"] = model[i].weight;
        terms.push_back(t);
    }
    j["model"] = terms;
    j["arithmetic"] =
            options.model == LABELIZER_MODEL_PUBLISHED ? "published" : "corrected";
    j["probe_radius"] = options.probe_radius;
    j["n_sphere_points"] = options.n_sphere_points;
    j["exclusion_distance"] = options.exclusion_distance;
    j["exclusion_exposure"] = options.exclusion_exposure;
    j["exclusion_residue"] = options.exclusion_residue;
    j["hse_radius"] = options.hse_radius;
    j["conservation_path"] = conservation_path;

    nlohmann::json f;
    f["probe_model"] = labelizer_dye_model_name(fret_options.probe_model);
    f["refine_probe_model"] = labelizer_dye_model_name(fret_options.refine_probe_model);
    f["n_refine"] = fret_options.n_refine;
    f["forster_radius"] = fret_options.forster_radius;
    f["label_score_threshold"] = fret_options.label_score_threshold;
    f["linker_length"] = fret_options.linker_length;
    f["linker_width"] = fret_options.linker_width;
    f["r1"] = fret_options.r1;
    f["r2"] = fret_options.r2;
    f["r3"] = fret_options.r3;
    f["grid_resolution"] = fret_options.grid_resolution;
    f["alpha_cone_radius"] = fret_options.alpha_cone_radius;
    f["alpha_cone_offset"] = fret_options.alpha_cone_offset;
    j["fret"] = f;
    return j.dump();
}

IMPBFF_END_NAMESPACE
