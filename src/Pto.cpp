/**
 * \file Pto.cpp
 * \brief ptolib's implementation, compiled once for IMP.bff, and the thin
 *        PtoWriter/PtoReader face this package's writers and readers use.
 *
 * The container lives in `internal/ptolib.h` (a vendored copy of
 * https://github.com/tpeulen/ptolib). This is the one translation unit in the
 * module that defines PTOLIB_IMPLEMENTATION; the JSON dependency it needs for
 * `.dstore` column descriptions is the module's own vendored nlohmann copy.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#define PTOLIB_JSON_INCLUDE <IMP/bff/internal/json.h>
#define PTOLIB_IMPLEMENTATION
// ptolib directly, and *before* Pto.h, rather than letting Pto.h pull it in.
//
// ptolib is an stb-style single header: its declarations sit inside a
// `PTOLIB_H` guard and its implementation *outside* one, compiled only where
// `PTOLIB_IMPLEMENTATION` is defined. This file is that place. But `Pto.h`
// has its own include guard, and in IMP's module build every `.cpp` is
// concatenated into one translation unit (`bff_all.cpp`) in alphabetical
// order -- so `Faspr.cpp`, `Labelizer.cpp` and two others include `Pto.h`
// first, and by the time this line is reached that guard is already closed.
// The include would then do nothing at all, ptolib's implementation would
// never be compiled, and the module failed to link on every `pto::File`
// symbol. Including ptolib itself is immune: its declaration half is skipped
// by its own guard and its implementation half is not guarded against this.
//
// The standalone build never saw it, because there each `.cpp` is its own
// translation unit and this file gets a clean slate.
#include <IMP/bff/internal/ptolib.h>
#include <IMP/bff/Pto.h>

#include <IMP/bff/Base.h>

#include <cstdio>
#include <cstring>


IMPBFF_BEGIN_NAMESPACE

namespace {

PtoObject to_object(const pto::PtoObject& o) {
    PtoObject out;
    out.uid = o.uid;
    out.name = o.name;
    out.kind = o.kind;
    out.encoding = o.encoding;
    out.offset = static_cast<std::size_t>(o.offset);
    out.size = static_cast<std::size_t>(o.size);
    return out;
}

const std::string kBanner =
        std::string(pto::kDefaultBanner) +
        "\n"
        "This container was written by IMP.bff.\n"
        "https://github.com/tpeulen/IMP.bff\n";

}  // namespace

// --------------------------------------------------------------------------
// the writer
// --------------------------------------------------------------------------

PtoWriter::PtoWriter(const std::string& path, int /*doctype_version*/)
        : closed_(false) {
    if (!file_.create(path, "", kBanner)) {
        IMP_THROW("PTO: cannot open " << path << " for writing: " << file_.error(),
                  IOException);
    }
    file_.set_writing_app("IMP.bff");
}

PtoWriter::~PtoWriter() {
    try {
        close();
    } catch (...) {
        // a destructor is the wrong place to raise; close() reports properly
    }
}

PtoObject PtoWriter::add(const std::string& name, const std::string& kind,
                         const std::string& encoding, const void* data,
                         std::size_t size) {
    if (closed_) IMP_THROW("PTO: the writer is closed", IOException);
    const unsigned long long uid =
            file_.add(kind, encoding, name, static_cast<const unsigned char*>(data), size);
    if (uid == 0) {
        IMP_THROW("PTO: writing " << file_.filename() << " failed: " << file_.error(),
                  IOException);
    }
    PtoObject object = to_object(file_.object(uid));
    objects_.push_back(object);
    return object;
}

void PtoWriter::close() {
    if (closed_) return;
    closed_ = true;
    if (!file_.is_open()) return;
    const bool ok = file_.commit();
    const std::string why = file_.error();
    const std::string path = file_.filename();
    file_.close();
    if (!ok) IMP_THROW("PTO: writing " << path << " failed: " << why, IOException);
}

// --------------------------------------------------------------------------
// the reader
// --------------------------------------------------------------------------

bool PtoReader::looks_like_pto(const std::string& path) {
    return pto::is_pto_file(path);
}

PtoReader::PtoReader(const std::string& path) : doctype_version_(0) {
    if (!file_.open(path, false)) {
        // The message keeps "EBML" and "container" in it: bin/ programs match
        // on them to say something useful about a file that is not a .pto.
        IMP_THROW("PTO: cannot open container " << path << ": "
                  << (file_.error().empty() ? "not an EBML document with DocType pto"
                                            : file_.error()),
                  IOException);
    }
    doctype_version_ = static_cast<int>(file_.doctype_version());
    const std::vector<pto::PtoObject> objs = file_.objects();
    objects_.reserve(objs.size());
    for (std::size_t i = 0; i < objs.size(); ++i) objects_.push_back(to_object(objs[i]));
}

PtoReader::~PtoReader() {}

int PtoReader::find(const std::string& name) const {
    for (std::size_t i = 0; i < objects_.size(); ++i) {
        if (objects_[i].name == name) return static_cast<int>(i);
    }
    return -1;
}

std::vector<unsigned char> PtoReader::data(const PtoObject& object) const {
    return file_.read(object.uid);
}

IMPBFF_END_NAMESPACE

// -------- from PtoProfile.cpp --------
/**
 * (formerly PtoProfile.cpp, now a section of this file)
 * \brief PTO.MFDB — the vocabulary layer over the PTO container.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/bff_config.h>
#include <IMP/bff/internal/Sha256.h>
#include <IMP/bff/internal/json.h>


#include <algorithm>
#include <map>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

const char* const MFDB_PROFILE = "PTO.MFDB";
const char* const MFDB_PROFILE_VERSION = "1.1";
const char* const MFDB_PROFILE_READ_VERSION = "1";
const char* const MFDB_DICTIONARY_VERSION = "1.8";

namespace {

std::vector<std::string> pp_list(const char* const* v, std::size_t n) {
    return std::vector<std::string>(v, v + n);
}

}  // namespace

// ---------------------------------------------------------------------------
// The controlled vocabularies
// ---------------------------------------------------------------------------
//
// Transcribed from `../mmfdb/src/mmfdb/data/mmfdb_flr_ext.dic`, version 1.8.
// Only the values this package can legitimately write are listed: the
// dictionary's `artifact_kind` carries seventy terms for photon streams,
// imaging and project archives that nothing here produces, and listing them
// would invite a caller to pick one.

const std::vector<std::string>& mfdb_artifact_kinds() {
    static const char* const v[] = {
        "readme", "sample_metadata", "parameter_table", "analysis_result",
        "processed_data", "derived_product", "row_mapping",
        "external_reference", "json_summary", "table_export"};
    static const std::vector<std::string> l = pp_list(v, 10);
    return l;
}

const std::vector<std::string>& mfdb_data_formats() {
    static const char* const v[] = {"json", "csv",  "tsv",  "cif",
                                    "text", "bin",  "npy",  "dstore",
                                    "pto",  "unknown"};
    static const std::vector<std::string> l = pp_list(v, 10);
    return l;
}

const std::vector<std::string>& mfdb_row_grains() {
    static const char* const v[] = {
        "photon", "burst", "dwell",  "segment", "pixel", "voxel",
        "frame",  "line",  "molecule", "track", "spot",  "region",
        "curve_point", "state", "species", "channel", "pair", "file",
        "histogram_bin", "spectrum",
        // Added in dictionary 1.8 for this package: a position on a polymer
        // that a label could be attached to, the thing
        // `_flr_poly_probe_position` addresses. Nothing existing fitted --
        // `species` and `state` count what a molecule is, not a place on one.
        "label_site"};
    static const std::vector<std::string> l = pp_list(v, 21);
    return l;
}

const std::vector<std::string>& mfdb_relationship_types() {
    static const char* const v[] = {
        "included_in", "contains", "derived_from", "supersedes",
        "uses_external_reference", "parameter_depends_on", "parameter_of",
        "linked_to", "project_contains", "grouped_in", "measured_sample",
        "calibrated_by", "maps_rows_of"};
    static const std::vector<std::string> l = pp_list(v, 13);
    return l;
}

const std::vector<std::string>& mfdb_units() {
    static const char* const v[] = {
        "hours", "minutes", "seconds", "milliseconds", "microseconds",
        "nanoseconds", "picoseconds", "femtoseconds", "hertz", "kilohertz",
        "megahertz", "counts", "counts_per_second", "photons", "nanometres",
        "micrometres", "angstroms", "pixels", "degrees", "radians", "celsius",
        "kelvins", "molar", "millimolar", "micromolar", "nanomolar",
        "picomolar", "dimensionless"};
    static const std::vector<std::string> l = pp_list(v, 28);
    return l;
}

const std::vector<std::string>& mfdb_operation_types() {
    // A subset of the dictionary's enumeration: what a writer in *this*
    // package produces. A reader must expect the rest -- most of that
    // enumeration is photon-level and belongs to tttrlib and chisurf.
    //
    // `clustering` is not `ndxplorer_clustering`, which is one tool's
    // burst-selection step. This one covers grouping conformations: a rotamer
    // library built from a trajectory, MSM microstate construction,
    // `cluster_frames_leader`. Added for the .drot libraries, which have to be
    // able to say *which* protocol built them -- dihedral k-means and leader
    // RMSD do not produce comparable libraries, and recording `analysis` for
    // both would be a lie that nothing ever checks.
    static const char* const v[] = {"analysis", "analysis_run", "validation",
                                    "model_fitting", "external_tool",
                                    "import", "calibration", "clustering"};
    static const std::vector<std::string> l = pp_list(v, 8);
    return l;
}

const std::vector<std::string>& mfdb_optical_property_names() {
    static const char* const v[] = {
        "quantum_yield", "extinction_coefficient", "excitation_maximum",
        "emission_maximum", "fluorescence_lifetime", "anisotropy_fundamental",
        "steric_radius", "hydrodynamic_radius"};
    static const std::vector<std::string> l = pp_list(v, 8);
    return l;
}

const std::vector<std::string>& mfdb_optical_property_units() {
    static const char* const v[] = {"dimensionless", "nanometres", "angstroms",
                                    "nanoseconds", "per_molar_per_centimetre"};
    static const std::vector<std::string> l = pp_list(v, 5);
    return l;
}

const std::vector<std::string>& mfdb_spectrum_types() {
    static const char* const v[] = {"absorption", "excitation", "emission",
                                    "transmission", "reflectance",
                                    "quantum_efficiency", "responsivity"};
    static const std::vector<std::string> l = pp_list(v, 7);
    return l;
}

const std::vector<std::string>& mfdb_spectrum_wavelength_units() {
    // `nm`, deliberately: every stored spectrum in this stack is labelled that
    // way and the dictionary follows the data rather than renaming 2896 rows.
    static const char* const v[] = {"nm", "um"};
    static const std::vector<std::string> l = pp_list(v, 2);
    return l;
}

const std::vector<std::string>& mfdb_spectrum_intensity_units() {
    static const char* const v[] = {"normalized", "counts",
                                    "per_molar_per_centimetre",
                                    "dimensionless"};
    static const std::vector<std::string> l = pp_list(v, 4);
    return l;
}

const std::vector<std::string>& mfdb_label_score_types() {
    static const char* const v[] = {
        "conservation", "solvent_exposure", "secondary_structure",
        "charge_environment", "tryptophan_proximity", "cysteine_resemblance",
        "methionine_exclusion", "fret_sensitivity", "measurement", "combined"};
    static const std::vector<std::string> l = pp_list(v, 10);
    return l;
}

const std::vector<std::string>& mfdb_label_score_statuses() {
    static const char* const v[] = {"scored", "excluded", "unresolved",
                                    "unavailable"};
    static const std::vector<std::string> l = pp_list(v, 4);
    return l;
}

const std::vector<std::string>& mfdb_probe_types() {
    static const char* const v[] = {"dye", "fluorescent_protein", "spin_label",
                                    "unspecified"};
    static const std::vector<std::string> l = pp_list(v, 4);
    return l;
}

const std::vector<std::string>& mfdb_label_score_definitions() {
    static const char* const v[] = {"labelizer", "consurf", "dssp", "msms"};
    static const std::vector<std::string> l = pp_list(v, 4);
    return l;
}

namespace {

const std::vector<std::string>* pp_enumeration(const std::string& item) {
    if (item == "_mmfdb_artifact.artifact_kind") return &mfdb_artifact_kinds();
    if (item == "_mmfdb_artifact.data_format") return &mfdb_data_formats();
    if (item == "_mmfdb_artifact.row_grain") return &mfdb_row_grains();
    if (item == "_mmfdb_edge.relationship_type") return &mfdb_relationship_types();
    if (item == "_mmfdb_column.units") return &mfdb_units();
    if (item == "_mmfdb_operation.operation_type") return &mfdb_operation_types();
    if (item == "_mmfdb_label_score.score_type") return &mfdb_label_score_types();
    if (item == "_mmfdb_label_score.status") return &mfdb_label_score_statuses();
    if (item == "_mmfdb_label_score.definition")
        return &mfdb_label_score_definitions();
    if (item == "_mmfdb_optical_property.property_name")
        return &mfdb_optical_property_names();
    if (item == "_mmfdb_optical_property.unit")
        return &mfdb_optical_property_units();
    if (item == "_mmfdb_probe.probe_type") return &mfdb_probe_types();
    if (item == "_mmfdb_spectrum.spectrum_type") return &mfdb_spectrum_types();
    if (item == "_mmfdb_spectrum.wavelength_unit")
        return &mfdb_spectrum_wavelength_units();
    if (item == "_mmfdb_spectrum.intensity_unit")
        return &mfdb_spectrum_intensity_units();
    return 0;
}

}  // namespace

bool mfdb_is_term(const std::string& item, const std::string& value) {
    const std::vector<std::string>* e = pp_enumeration(item);
    if (!e) return false;
    return std::find(e->begin(), e->end(), value) != e->end();
}

void mfdb_check_term(const std::string& item, const std::string& value) {
    const std::vector<std::string>* e = pp_enumeration(item);
    if (!e) {
        IMP_THROW("mfdb_check_term: no enumeration is declared here for "
                  << item << "; add it from mmfdb_flr_ext.dic rather than "
                  << "writing an unchecked value", ValueException);
    }
    if (std::find(e->begin(), e->end(), value) != e->end()) return;
    std::ostringstream allowed;
    for (std::size_t i = 0; i < e->size(); ++i) {
        allowed << (i ? ", " : "") << (*e)[i];
    }
    IMP_THROW("mfdb_check_term: '" << value << "' is not a value of " << item
              << ". Allowed: " << allowed.str()
              << ". A term that is genuinely missing belongs in "
              << "mmfdb_flr_ext.dic first, not here.", ValueException);
}

// ---------------------------------------------------------------------------
// Columns
// ---------------------------------------------------------------------------

std::string mfdb_columns_json(const std::vector<MfdbColumn>& columns) {
    nlohmann::json j = nlohmann::json::array();
    for (std::size_t i = 0; i < columns.size(); ++i) {
        const MfdbColumn& c = columns[i];
        nlohmann::json e;
        e["name"] = c.name;
        // No unit means the unit is unknown, which is a different claim from
        // `dimensionless`; an empty string is therefore omitted, not written.
        if (!c.units.empty()) {
            mfdb_check_term("_mmfdb_column.units", c.units);
            e["units"] = c.units;
        }
        if (!c.item.empty()) e["item"] = c.item;
        if (!c.description.empty()) e["description"] = c.description;
        j.push_back(e);
    }
    return j.dump();
}

// ---------------------------------------------------------------------------
// The tag sets
// ---------------------------------------------------------------------------

std::vector<MfdbTag> mfdb_container_tags() {
    std::vector<MfdbTag> t;
    t.push_back(MfdbTag("_mmfdb_container.profile", MFDB_PROFILE));
    t.push_back(MfdbTag("_mmfdb_container.profile_version", MFDB_PROFILE_VERSION));
    t.push_back(MfdbTag("_mmfdb_container.profile_read_version",
                        MFDB_PROFILE_READ_VERSION));
    t.push_back(MfdbTag("_mmfdb_container.format", "pto"));
    t.push_back(MfdbTag("_mmfdb_container.dictionary_version",
                        MFDB_DICTIONARY_VERSION));
    return t;
}

std::vector<MfdbTag> mfdb_artifact_tags(const std::string& artifact_id,
                                        const std::string& artifact_kind,
                                        const std::string& data_format,
                                        const std::string& row_grain,
                                        long row_count,
                                        const std::string& checksum,
                                        const std::vector<MfdbColumn>& columns) {
    mfdb_check_term("_mmfdb_artifact.artifact_kind", artifact_kind);
    mfdb_check_term("_mmfdb_artifact.data_format", data_format);
    if (!row_grain.empty()) {
        mfdb_check_term("_mmfdb_artifact.row_grain", row_grain);
    }
    std::vector<MfdbTag> t;
    t.push_back(MfdbTag("_mmfdb_artifact.artifact_id", artifact_id));
    t.push_back(MfdbTag("_mmfdb_artifact.artifact_kind", artifact_kind));
    t.push_back(MfdbTag("_mmfdb_artifact.data_format", data_format));
    if (!row_grain.empty()) {
        t.push_back(MfdbTag("_mmfdb_artifact.row_grain", row_grain));
    }
    if (row_count >= 0) {
        std::ostringstream os;
        os << row_count;
        t.push_back(MfdbTag("_mmfdb_artifact.row_count", os.str()));
    }
    if (!checksum.empty()) {
        t.push_back(MfdbTag("_mmfdb_artifact.checksum", checksum));
        t.push_back(MfdbTag("_mmfdb_artifact.checksum_algorithm", "sha256"));
    }
    if (!columns.empty()) {
        t.push_back(MfdbTag("_mmfdb_column", mfdb_columns_json(columns)));
    }
    return t;
}

std::vector<MfdbTag> mfdb_operation_tags(const std::string& operation_type,
                                         const std::string& algorithm,
                                         const std::string& settings_json,
                                         const std::string& software_package,
                                         const std::string& software_version) {
    mfdb_check_term("_mmfdb_operation.operation_type", operation_type);
    std::vector<MfdbTag> t;
    t.push_back(MfdbTag("_mmfdb_operation.operation_type", operation_type));
    if (!algorithm.empty()) {
        t.push_back(MfdbTag("_mmfdb_operation.algorithm", algorithm));
    }
    t.push_back(MfdbTag("_mmfdb_operation.settings_json", settings_json));
    // The identity of a run: same settings, same artifact.
    t.push_back(MfdbTag("_mmfdb_operation.settings_hash",
                        mfdb_checksum(settings_json)));
    t.push_back(MfdbTag("_mmfdb_operation.software_package", software_package));
    t.push_back(MfdbTag("_mmfdb_operation.software_version", software_version));
    t.push_back(MfdbTag("_mmfdb_operation.dictionary_version",
                        MFDB_DICTIONARY_VERSION));
    return t;
}

std::vector<MfdbTag> mfdb_edge_tags(const std::string& source_node_id,
                                    const std::string& target_node_id,
                                    const std::string& relationship_type,
                                    const std::string& source_row_column,
                                    const std::string& target_row_column) {
    mfdb_check_term("_mmfdb_edge.relationship_type", relationship_type);
    std::vector<MfdbTag> t;
    t.push_back(MfdbTag("_mmfdb_edge.source_node_id", source_node_id));
    t.push_back(MfdbTag("_mmfdb_edge.target_node_id", target_node_id));
    t.push_back(MfdbTag("_mmfdb_edge.relationship_type", relationship_type));
    if (!source_row_column.empty()) {
        t.push_back(MfdbTag("_mmfdb_edge.source_row_column", source_row_column));
    }
    if (!target_row_column.empty()) {
        t.push_back(MfdbTag("_mmfdb_edge.target_row_column", target_row_column));
    }
    return t;
}

std::vector<MfdbTag> mfdb_attribution_tags(const MfdbAttribution& a) {
    if (a.license.empty() && a.terms.empty()) {
        IMP_THROW("mfdb_attribution_tags: an artifact needs either a license "
                  "(an SPDX identifier) or terms (the conditions in words). "
                  "Data whose terms are unknown is not the same as data that "
                  "is unrestricted, and the difference matters to whoever "
                  "redistributes it.", ValueException);
    }
    std::vector<MfdbTag> t;
    if (!a.author.empty())
        t.push_back(MfdbTag("_mmfdb_artifact.author", a.author));
    if (!a.citation.empty())
        t.push_back(MfdbTag("_mmfdb_artifact.citation", a.citation));
    if (!a.license.empty())
        t.push_back(MfdbTag("_mmfdb_artifact.license", a.license));
    if (!a.terms.empty())
        t.push_back(MfdbTag("_mmfdb_artifact.terms", a.terms));
    if (!a.terms_url.empty())
        t.push_back(MfdbTag("_mmfdb_artifact.terms_url", a.terms_url));
    if (!a.source.empty())
        t.push_back(MfdbTag("_mmfdb_artifact.source", a.source));
    if (!a.redistributed_via.empty())
        t.push_back(MfdbTag("_mmfdb_artifact.redistributed_via",
                            a.redistributed_via));
    return t;
}

std::string mfdb_checksum(const std::string& bytes) {
    return internal::sha256_hex(bytes);
}

IMPBFF_END_NAMESPACE
