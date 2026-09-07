/**
 * \file FPSProject.cpp
 * \brief The document that binds a FRET docking or screening run together.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/FPSProject.h>

#include <IMP/bff/FPSIO.h>
#include <IMP/bff/FPSSchema.h>
#include <IMP/bff/internal/Text.h>
#include <IMP/bff/internal/json.h>

#include <IMP/bff/Base.h>

#include <cstdlib>
#include <fstream>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

// Named, not anonymous: IMP compiles this module as one translation unit, so
// an anonymous namespace here would still collide with another .cpp's helpers.
namespace fps_project {
using IMP::bff::internal::ends_with;
using IMP::bff::internal::file_exists;
using IMP::bff::internal::trimmed;

std::string directory_of(const std::string& path) {
    const std::size_t at = path.find_last_of('/');
    return at == std::string::npos ? std::string() : path.substr(0, at);
}

//! Windows drive letters count: a legacy project's paths are `C:\...`.
/*! They will not exist on this machine, but treating one as *relative* and
    gluing the project's directory in front of it produces a path that is wrong
    in a way nobody can read. Left alone, it is obviously a Windows path. */
bool is_absolute(const std::string& path) {
    if (path.empty()) return false;
    if (path[0] == '/' || path[0] == '\\') return true;
    return path.size() > 1 && path[1] == ':';
}

//! A path as written, resolved against the project's own directory.
std::string resolved(const std::string& path, const std::string& project_path) {
    if (path.empty() || is_absolute(path)) return path;
    const std::string dir = directory_of(project_path);
    return dir.empty() ? path : dir + "/" + path;
}

//! `Double.ToString()` under an unknown culture.
/*! `strtod` reads `0.0005`; a German-locale FPS wrote `0,0005`, which `strtod`
    reads as **0** and stops at the comma. Retrying with the comma turned into
    a point is the difference between a tolerance and a silent zero. Only tried
    when the first parse consumed less than the whole token, so a value that
    parses cleanly is never rewritten. */
double to_double(const std::string& text, double fallback) {
    const std::string t = trimmed(text);
    if (t.empty()) return fallback;
    char* end = nullptr;
    const double value = std::strtod(t.c_str(), &end);
    if (end != nullptr && *end == '\0') return value;
    std::string swapped = t;
    bool changed = false;
    for (std::size_t i = 0; i < swapped.size(); ++i) {
        if (swapped[i] == ',') {
            swapped[i] = '.';
            changed = true;
        }
    }
    if (!changed) return fallback;
    end = nullptr;
    const double retried = std::strtod(swapped.c_str(), &end);
    return (end != nullptr && *end == '\0') ? retried : fallback;
}

int to_int(const std::string& text, int fallback) {
    const std::string t = trimmed(text);
    if (t.empty()) return fallback;
    char* end = nullptr;
    const long value = std::strtol(t.c_str(), &end, 10);
    if (end != nullptr && *end == '\0') return static_cast<int>(value);
    // A .NET `Int32` never prints a decimal point, but an option that was a
    // double in one build and an int in another does; take the integer part
    // rather than dropping the value.
    const double as_double = to_double(t, static_cast<double>(fallback));
    return static_cast<int>(as_double);
}

//! `Boolean.ToString()` is `True`/`False`; anything else is read leniently.
bool to_bool(const std::string& text) {
    std::string t = trimmed(text);
    for (std::size_t i = 0; i < t.size(); ++i) {
        t[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(t[i])));
    }
    return t == "true" || t == "1" || t == "yes";
}

std::vector<std::string> split_ws(const std::string& line) {
    std::vector<std::string> parts;
    std::istringstream in(line);
    std::string token;
    while (in >> token) parts.push_back(token);
    return parts;
}

double json_double(const nlohmann::json& obj, const char* key, double fallback) {
    if (!obj.is_object() || !obj.contains(key)) return fallback;
    const nlohmann::json& v = obj[key];
    return v.is_number() ? v.get<double>() : fallback;
}

int json_int(const nlohmann::json& obj, const char* key, int fallback) {
    if (!obj.is_object() || !obj.contains(key)) return fallback;
    const nlohmann::json& v = obj[key];
    if (v.is_number_integer()) return v.get<int>();
    if (v.is_number()) return static_cast<int>(v.get<double>());
    return fallback;
}

bool json_bool(const nlohmann::json& obj, const char* key, bool fallback) {
    if (!obj.is_object() || !obj.contains(key)) return fallback;
    const nlohmann::json& v = obj[key];
    if (v.is_boolean()) return v.get<bool>();
    // A 0/1 is accepted because a hand-written project file and a converted
    // one do not always agree about how a flag is spelled; anything else is
    // the default rather than a guess.
    if (v.is_number()) return v.get<double>() != 0.0;
    return fallback;
}

std::string json_string(const nlohmann::json& obj, const char* key,
                        const std::string& fallback) {
    if (!obj.is_object() || !obj.contains(key)) return fallback;
    const nlohmann::json& v = obj[key];
    return v.is_string() ? v.get<std::string>() : fallback;
}

std::vector<std::string> json_strings(const nlohmann::json& obj,
                                      const char* key) {
    std::vector<std::string> out;
    if (!obj.is_object() || !obj.contains(key) || !obj[key].is_array()) {
        return out;
    }
    const nlohmann::json& a = obj[key];
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].is_string()) out.push_back(a[i].get<std::string>());
    }
    return out;
}

nlohmann::json mode_parameters_json(const FPSModeParameters& p) {
    nlohmann::json out;
    out["viscosity_factor"] = p.viscosity_factor;
    out["time_step_factor"] = p.time_step_factor;
    out["max_iterations"] = p.max_iterations;
    out["max_force"] = p.max_force;
    out["clash_tolerance"] = p.clash_tolerance;
    out["rkt"] = p.rkt;
    out["e_tolerance"] = p.e_tolerance;
    out["k_tolerance"] = p.k_tolerance;
    out["f_tolerance"] = p.f_tolerance;
    out["t_tolerance"] = p.t_tolerance;
    out["optimize_selected"] = p.optimize_selected;
    out["iterations"] = p.iterations;
    return out;
}

FPSModeParameters mode_parameters_from_json(const nlohmann::json& obj,
                                            const FPSModeParameters& base) {
    FPSModeParameters p = base;
    p.viscosity_factor = json_double(obj, "viscosity_factor", p.viscosity_factor);
    p.time_step_factor = json_double(obj, "time_step_factor", p.time_step_factor);
    p.max_iterations = json_int(obj, "max_iterations", p.max_iterations);
    p.max_force = json_double(obj, "max_force", p.max_force);
    p.clash_tolerance = json_double(obj, "clash_tolerance", p.clash_tolerance);
    p.rkt = json_double(obj, "rkt", p.rkt);
    p.e_tolerance = json_double(obj, "e_tolerance", p.e_tolerance);
    p.k_tolerance = json_double(obj, "k_tolerance", p.k_tolerance);
    p.f_tolerance = json_double(obj, "f_tolerance", p.f_tolerance);
    p.t_tolerance = json_double(obj, "t_tolerance", p.t_tolerance);
    p.optimize_selected =
            json_string(obj, "optimize_selected", p.optimize_selected);
    p.iterations = json_int(obj, "iterations", p.iterations);
    return p;
}

//! One ` FieldName = value` line of an exported struct, applied to a block.
void apply_struct_field(FPSModeParameters& p, const std::string& field,
                        const std::string& value) {
    if (field == "ViscosityFactor") {
        p.viscosity_factor = to_double(value, p.viscosity_factor);
    } else if (field == "TimeStepFactor") {
        p.time_step_factor = to_double(value, p.time_step_factor);
    } else if (field == "MaxIterations") {
        p.max_iterations = to_int(value, p.max_iterations);
    } else if (field == "MaxForce") {
        p.max_force = to_double(value, p.max_force);
    } else if (field == "ClashTolerance") {
        p.clash_tolerance = to_double(value, p.clash_tolerance);
    } else if (field == "rkT") {
        p.rkt = to_double(value, p.rkt);
    } else if (field == "ETolerance") {
        p.e_tolerance = to_double(value, p.e_tolerance);
    } else if (field == "KTolerance") {
        p.k_tolerance = to_double(value, p.k_tolerance);
    } else if (field == "FTolerance") {
        p.f_tolerance = to_double(value, p.f_tolerance);
    } else if (field == "TTolerance") {
        p.t_tolerance = to_double(value, p.t_tolerance);
    } else if (field == "OptimizeSelected") {
        p.optimize_selected = trimmed(value);
    }
}

//! One option of the export: its key, its raw value, its parts.
struct Option {
    std::string key, value;
    //! Elements of an array value, split on whitespace.
    std::vector<std::string> elements;
    //! ` Field = value` lines of a struct value, in file order.
    std::vector<std::pair<std::string, std::string> > fields;
    bool is_array;

    Option() : is_array(false) {}
};

//! Split a `Key = value` line. False when the line is not one.
bool split_assignment(const std::string& line, std::string& key,
                      std::string& value) {
    const std::size_t at = line.find('=');
    if (at == std::string::npos) return false;
    key = trimmed(line.substr(0, at));
    value = trimmed(line.substr(at + 1));
    return !key.empty();
}

}  // namespace fps_project

std::vector<std::string> fps_mode_names() {
    return {"Dock", "Refine", "Error estimation", "Sample", "Screening"};
}

FPSModeParameters fps_mode_parameters(const std::string& mode) {
    // ProjectData.cs:61-149, and okf/references/fps-sampling-and-protocol.md
    // section 6. The four fields that differ between modes are
    // ViscosityFactor / TimeStepFactor / MaxIterations / MaxForce /
    // ClashTolerance / the three tolerances / OptimizeSelected; everything
    // else is the same in all five.
    FPSModeParameters p;  // the Dock block is the constructor's default
    if (mode == "Refine") {
        p.viscosity_factor = 0.7;
        p.time_step_factor = 0.5;
        p.max_iterations = 500000;
        p.max_force = 10000.0;
        p.clash_tolerance = 0.5;
        p.k_tolerance = 0.0005;
        p.f_tolerance = 0.0005;
        p.t_tolerance = 0.01;
        p.optimize_selected = "All";
    } else if (mode == "Error estimation") {
        p.viscosity_factor = 0.7;
        p.max_iterations = 100000;
        p.max_force = 10000.0;
        p.clash_tolerance = 0.5;
        p.optimize_selected = "All";
    } else if (mode == "Sample") {
        // MaxIterations counts *accepted* Metropolis moves per snapshot here,
        // not integrator steps, so 8000 is not 25x fewer of the same thing.
        p.max_iterations = 8000;
    } else if (mode == "Screening") {
        // Screening runs FilterEngine, which reads only OptimizeSelected and
        // the conversion R0; the rest of the block is stored and unused.
        p.iterations = 0;
    }
    return p;
}

std::string FPSModeParameters::get_json() const {
    return fps_project::mode_parameters_json(*this).dump(2);
}

FPSProject::FPSProject()
    : mode("None"), poses("[]"), fixed_body(0), ev_weight(1.0), sigma_da(6.0),
      shuffle(10.0), refine_av_cycles(0), coarse_clash(true),
      clash_radii_source("imp"), clash_radii_scale(1.0),
      dock(fps_mode_parameters("Dock")),
      refine(fps_mode_parameters("Refine")),
      error_estimation(fps_mode_parameters("Error estimation")),
      sample(fps_mode_parameters("Sample")),
      screening(fps_mode_parameters("Screening")) {}

FPSModeParameters FPSProject::get_parameters(const std::string& mode) const {
    if (mode == "Refine") return refine;
    if (mode == "Error estimation") return error_estimation;
    if (mode == "Sample") return sample;
    if (mode == "Screening") return screening;
    return dock;
}

void FPSProject::set_parameters(const std::string& mode,
                                const FPSModeParameters& parameters) {
    if (mode == "Dock") {
        dock = parameters;
    } else if (mode == "Refine") {
        refine = parameters;
    } else if (mode == "Error estimation") {
        error_estimation = parameters;
    } else if (mode == "Sample") {
        sample = parameters;
    } else if (mode == "Screening") {
        screening = parameters;
    }
}

std::string FPSProject::get_labelling_path() const {
    // Resolved here rather than at read time so the *stored* strings stay as
    // they were written: a project that says `hiv_rt.fps.json` keeps saying so
    // after a read and a write, and moving the directory keeps working. A read
    // that resolved eagerly would write absolute paths back out and quietly
    // pin the file to one machine.
    if (!labelling_json.empty()) {
        return fps_project::resolved(labelling_json, path);
    }
    if (!positions_path.empty()) {
        return fps_project::resolved(positions_path, path);
    }
    return path;
}

std::string FPSProject::get_distances_path() const {
    return fps_project::resolved(distances_path, path);
}

std::vector<std::string> FPSProject::get_structure_paths() const {
    std::vector<std::string> out;
    out.reserve(structures.size());
    for (std::size_t i = 0; i < structures.size(); ++i) {
        out.push_back(fps_project::resolved(structures[i], path));
    }
    return out;
}

DockingParameters FPSProject::get_docking_parameters(
        const std::string& mode) const {
    const FPSModeParameters p = get_parameters(mode);
    DockingParameters out;
    out.n_frames = p.iterations;
    out.max_force = p.max_force;
    out.clash_tolerance = p.clash_tolerance;
    out.optimize_selected = p.optimize_selected;
    out.score_set = score_set;
    out.fixed_body = fixed_body;
    out.ev_weight = ev_weight;
    out.sigma_da = sigma_da;
    out.shuffle_max_translation = shuffle;
    out.refine_av_cycles = refine_av_cycles;
    out.coarse_clash = coarse_clash;
    out.clash_radii_source = clash_radii_source;
    out.clash_radii_scale = clash_radii_scale;
    return out;
}

std::vector<std::string> FPSProject::resolve_selected_distances(
        const std::vector<std::string>& distance_names) const {
    std::vector<std::string> out;
    const std::size_t n =
            std::min(selected_flags.size(), distance_names.size());
    for (std::size_t i = 0; i < n; ++i) {
        if (selected_flags[i]) out.push_back(distance_names[i]);
    }
    return out;
}

std::vector<std::string> FPSProject::get_problems() const {
    std::vector<std::string> out;
    if (structures.empty()) {
        out.push_back("no structures: a project needs one per rigid body, in "
                      "body order");
    }
    const std::vector<std::string> paths = get_structure_paths();
    for (std::size_t i = 0; i < paths.size(); ++i) {
        if (!fps_project::file_exists(paths[i])) {
            out.push_back("structure " + std::to_string(i) + " is missing: " +
                          paths[i]);
        }
    }
    if (!labelling_json.empty() && !positions_path.empty()) {
        // Not fatal -- `get_labelling_path` has a rule and follows it -- but
        // a project that names two labelling sources is a project someone
        // half-converted, and the one that loses is silent otherwise.
        out.push_back("two labelling sources: labelling_json (" +
                      labelling_json + ") is used and positions_path (" +
                      positions_path + ") is ignored");
    }
    const std::string labelling = get_labelling_path();
    if (labelling.empty()) {
        out.push_back("no labelling source: set labelling_json, or "
                      "positions_path for the legacy pair");
    } else if (!fps_project::file_exists(labelling)) {
        out.push_back("the labelling file is missing: " + labelling);
    }
    if (!distances_path.empty() &&
        !fps_project::file_exists(get_distances_path())) {
        out.push_back("the distances file is missing: " + get_distances_path());
    }
    if (fixed_body < 0 ||
        (!structures.empty() &&
         fixed_body >= static_cast<int>(structures.size()))) {
        out.push_back("fixed_body " + std::to_string(fixed_body) +
                      " is not one of the " + std::to_string(structures.size()) +
                      " bodies");
    }
    if (ev_weight < 0) {
        out.push_back("ev_weight is negative (" + std::to_string(ev_weight) +
                      "): a clash would lower the score");
    }
    if (sigma_da <= 0) {
        out.push_back("sigma_da must be positive, not " +
                      std::to_string(sigma_da));
    }
    if (clash_radii_scale <= 0) {
        out.push_back("clash_radii_scale must be positive, not " +
                      std::to_string(clash_radii_scale));
    }
    if (clash_radii_source != "imp" && clash_radii_source != "olga") {
        out.push_back("clash_radii_source is '" + clash_radii_source +
                      "'; the radii sets are 'imp' and 'olga'");
    }
    // A coarse bead is one sphere per residue, so a per-atom radii table has
    // nothing to say about it; `dock_minimize` refuses the combination and
    // the project should say so before a run finds out.
    if (coarse_clash &&
        (clash_radii_source != "imp" || clash_radii_scale != 1.0)) {
        out.push_back("coarse_clash is on, so clash_radii_source / "
                      "clash_radii_scale cannot be honoured: a coarse bead is "
                      "one sphere per residue, not an atom");
    }
    return out;
}

std::string FPSProject::get_json() const {
    nlohmann::json out;
    out["schema_version"] = fps_schema_version();
    out["mode"] = mode;
    if (!source.empty()) out["source"] = source;
    out["structures"] = structures;
    out["labelling_json"] = labelling_json;
    out["positions_path"] = positions_path;
    out["distances_path"] = distances_path;
    out["score_set"] = score_set;
    out["selected_distances"] = selected_distances;
    out["selected_flags"] = selected_flags;
    out["fixed_body"] = fixed_body;
    out["ev_weight"] = ev_weight;
    out["sigma_da"] = sigma_da;
    out["shuffle"] = shuffle;
    out["refine_av_cycles"] = refine_av_cycles;
    out["coarse_clash"] = coarse_clash;
    out["clash_radii_source"] = clash_radii_source;
    out["clash_radii_scale"] = clash_radii_scale;
    // The poses are stored as the JSON array `capture_poses` produced, not as
    // its text: a run that is stored and continued must not depend on a string
    // surviving two escapings.
    nlohmann::json poses_value =
            nlohmann::json::parse(poses.empty() ? "[]" : poses, nullptr, false);
    if (poses_value.is_discarded() || !poses_value.is_array()) {
        poses_value = nlohmann::json::array();
    }
    out["poses"] = poses_value;

    nlohmann::json parameters = nlohmann::json::object();
    const std::vector<std::string> modes = fps_mode_names();
    for (std::size_t i = 0; i < modes.size(); ++i) {
        parameters[modes[i]] =
                fps_project::mode_parameters_json(get_parameters(modes[i]));
    }
    out["parameters"] = parameters;

    nlohmann::json conv;
    conv["forster_radius"] = conversion.forster_radius;
    conv["polynomial_order"] = conversion.polynomial_order;
    out["conversion"] = conv;

    nlohmann::json avg;
    avg["grid_size"] = av.grid_size;
    avg["min_grid_size"] = av.min_grid_size;
    avg["linker_initial_sphere"] = av.linker_initial_sphere;
    avg["link_search_nodes"] = av.link_search_nodes;
    avg["e_samples"] = av.e_samples;
    out["av"] = avg;

    return out.dump(2);
}

FPSProject fps_project_from_json(const std::string& project_json,
                                 const std::string& path) {
    using namespace fps_project;
    nlohmann::json obj;
    try {
        obj = nlohmann::json::parse(project_json);
    } catch (const std::exception& e) {
        IMP_THROW("Cannot parse the Project section as JSON: " << e.what(),
                  ValueException);
    }
    FPSProject p;
    p.path = path;
    if (!obj.is_object()) return p;

    p.mode = json_string(obj, "mode", p.mode);
    p.source = json_string(obj, "source", p.source);
    const std::vector<std::string> structures = json_strings(obj, "structures");
    p.structures = structures;
    p.labelling_json = json_string(obj, "labelling_json", "");
    p.positions_path = json_string(obj, "positions_path", "");
    p.distances_path = json_string(obj, "distances_path", "");
    p.score_set = json_string(obj, "score_set", "");
    p.selected_distances = json_strings(obj, "selected_distances");
    if (obj.contains("selected_flags") && obj["selected_flags"].is_array()) {
        const nlohmann::json& flags = obj["selected_flags"];
        for (std::size_t i = 0; i < flags.size(); ++i) {
            if (flags[i].is_boolean()) {
                p.selected_flags.push_back(flags[i].get<bool>() ? 1 : 0);
            } else if (flags[i].is_number()) {
                p.selected_flags.push_back(flags[i].get<double>() != 0 ? 1 : 0);
            }
        }
    }
    p.fixed_body = json_int(obj, "fixed_body", p.fixed_body);
    p.ev_weight = json_double(obj, "ev_weight", p.ev_weight);
    p.sigma_da = json_double(obj, "sigma_da", p.sigma_da);
    p.shuffle = json_double(obj, "shuffle", p.shuffle);
    p.refine_av_cycles = json_int(obj, "refine_av_cycles", p.refine_av_cycles);
    p.coarse_clash = json_bool(obj, "coarse_clash", p.coarse_clash);
    p.clash_radii_source =
            json_string(obj, "clash_radii_source", p.clash_radii_source);
    p.clash_radii_scale =
            json_double(obj, "clash_radii_scale", p.clash_radii_scale);
    if (obj.contains("poses") && obj["poses"].is_array()) {
        p.poses = obj["poses"].dump();
    }

    if (obj.contains("parameters") && obj["parameters"].is_object()) {
        const nlohmann::json& parameters = obj["parameters"];
        const std::vector<std::string> modes = fps_mode_names();
        for (std::size_t i = 0; i < modes.size(); ++i) {
            if (!parameters.contains(modes[i])) continue;
            p.set_parameters(modes[i],
                             mode_parameters_from_json(parameters[modes[i]],
                                                       p.get_parameters(modes[i])));
        }
    }
    if (obj.contains("conversion")) {
        const nlohmann::json& c = obj["conversion"];
        p.conversion.forster_radius =
                json_double(c, "forster_radius", p.conversion.forster_radius);
        p.conversion.polynomial_order =
                json_int(c, "polynomial_order", p.conversion.polynomial_order);
    }
    if (obj.contains("av")) {
        const nlohmann::json& a = obj["av"];
        p.av.grid_size = json_double(a, "grid_size", p.av.grid_size);
        p.av.min_grid_size = json_double(a, "min_grid_size", p.av.min_grid_size);
        p.av.linker_initial_sphere =
                json_double(a, "linker_initial_sphere",
                            p.av.linker_initial_sphere);
        p.av.link_search_nodes =
                json_int(a, "link_search_nodes", p.av.link_search_nodes);
        p.av.e_samples = json_int(a, "e_samples", p.av.e_samples);
    }
    return p;
}

FPSProject read_fps_project_txt(const std::string& path) {
    using namespace fps_project;
    std::ifstream in(path.c_str());
    if (!in) IMP_THROW("Cannot read " << path, IOException);

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[line.size() - 1] == '\r') {
            line.erase(line.size() - 1);
        }
        lines.push_back(line);
    }

    std::string banner;
    std::vector<Option> options;
    for (std::size_t i = 0; i < lines.size(); ++i) {
        const std::string& raw = lines[i];
        if (raw.empty()) continue;
        if (raw[0] == '#') {
            // Category header, or -- for the first two -- the application
            // banner and the save date. The banner is kept as provenance; the
            // categories carry no information a key does not, because every
            // key of the export is unique across categories.
            if (banner.empty()) banner = trimmed(raw.substr(1));
            continue;
        }
        if (std::isspace(static_cast<unsigned char>(raw[0]))) {
            // An indented line belongs to the option above it. Which kind it
            // is was decided by that option's own value, never by this line:
            // ` Field = value` is a legal array element.
            if (options.empty()) continue;
            Option& owner = options.back();
            if (owner.is_array) continue;  // consumed when the key was read
            std::string field, value;
            if (split_assignment(raw, field, value)) {
                owner.fields.push_back(std::make_pair(field, value));
            }
            continue;
        }
        Option option;
        if (!split_assignment(raw, option.key, option.value)) continue;
        // `System.String[]`, `System.Boolean[]`: an array prints its *type* on
        // the key line and its elements on the next one, which is written
        // even when the array is empty. Consume that line here so a blank one
        // is read as "no elements" rather than as a category separator.
        option.is_array = ends_with(option.value, "[]");
        if (option.is_array && i + 1 < lines.size()) {
            option.elements = split_ws(lines[i + 1]);
            ++i;
        }
        options.push_back(option);
    }

    FPSProject project;
    project.path = path;
    project.source = banner;
    for (std::size_t i = 0; i < options.size(); ++i) {
        const Option& o = options[i];
        if (o.key == "ProjectFPSMode") {
            project.mode = o.value;
        } else if (o.key == "MoleculesPaths") {
            project.structures = o.elements;
        } else if (o.key == "LabelingPositionsPath") {
            project.positions_path = o.value;
        } else if (o.key == "DistancesPath") {
            project.distances_path = o.value;
        } else if (o.key == "SelectedDistances") {
            for (std::size_t k = 0; k < o.elements.size(); ++k) {
                project.selected_flags.push_back(to_bool(o.elements[k]) ? 1 : 0);
            }
        } else if (o.key == "ConversionParameters") {
            for (std::size_t k = 0; k < o.fields.size(); ++k) {
                const std::string& f = o.fields[k].first;
                const std::string& v = o.fields[k].second;
                if (f == "R0") {
                    project.conversion.forster_radius =
                            to_double(v, project.conversion.forster_radius);
                } else if (f == "PolynomOrder") {
                    project.conversion.polynomial_order =
                            to_int(v, project.conversion.polynomial_order);
                }
            }
        } else if (o.key == "AVGlobalParameters") {
            for (std::size_t k = 0; k < o.fields.size(); ++k) {
                const std::string& f = o.fields[k].first;
                const std::string& v = o.fields[k].second;
                if (f == "GridSize") {
                    project.av.grid_size = to_double(v, project.av.grid_size);
                } else if (f == "MinGridSize") {
                    project.av.min_grid_size =
                            to_double(v, project.av.min_grid_size);
                } else if (f == "LinkerInitialSphere") {
                    project.av.linker_initial_sphere =
                            to_double(v, project.av.linker_initial_sphere);
                } else if (f == "LinkSearchNodes") {
                    project.av.link_search_nodes =
                            to_int(v, project.av.link_search_nodes);
                } else if (f == "ESamples") {
                    project.av.e_samples = to_int(v, project.av.e_samples);
                }
            }
        } else {
            // The five parameter blocks. FPS's property names, not its
            // category names: `DockParameters` sits in a category called
            // "Search parameters", and keying on the property is what makes
            // the mapping legible.
            std::string mode;
            if (o.key == "DockParameters") {
                mode = "Dock";
            } else if (o.key == "RefineParameters") {
                mode = "Refine";
            } else if (o.key == "ErrorEstimationParameters") {
                mode = "Error estimation";
            } else if (o.key == "SampleParameters") {
                mode = "Sample";
            } else if (o.key == "ScreeningParameters") {
                mode = "Screening";
            }
            if (mode.empty()) continue;  // a key from another build; skipped
            FPSModeParameters p = project.get_parameters(mode);
            for (std::size_t k = 0; k < o.fields.size(); ++k) {
                apply_struct_field(p, o.fields[k].first, o.fields[k].second);
            }
            project.set_parameters(mode, p);
        }
    }
    return project;
}

FPSProject read_fps_project(const std::string& path) {
    if (!fps_project::ends_with(path, ".json")) {
        return read_fps_project_txt(path);
    }
    std::ifstream in(path.c_str());
    if (!in) IMP_THROW("Cannot read " << path, IOException);
    std::ostringstream buffer;
    buffer << in.rdbuf();
    nlohmann::json payload;
    try {
        payload = nlohmann::json::parse(buffer.str());
    } catch (const std::exception& e) {
        IMP_THROW("Cannot parse " << path << " as JSON: " << e.what(),
                  ValueException);
    }
    if (payload.is_object() && payload.contains("Project")) {
        return fps_project_from_json(payload["Project"].dump(), path);
    }
    // A plain fps.json with no Project section: the defaults, pointed at it.
    // Not an error -- "run this fps.json" is a thing to want, and refusing it
    // would make the project a barrier rather than a convenience.
    FPSProject project;
    project.path = path;
    return project;
}

void write_fps_project(const std::string& path, const FPSProject& project,
                       const std::string& positions_json,
                       const std::string& distances_json,
                       const std::string& score_sets_json,
                       const std::string& extra_json, bool validate) {
    nlohmann::json extra;
    try {
        extra = nlohmann::json::parse(extra_json);
    } catch (const std::exception& e) {
        IMP_THROW("Cannot parse the extra keys as JSON: " << e.what(),
                  ValueException);
    }
    if (!extra.is_object()) extra = nlohmann::json::object();
    extra["Project"] = nlohmann::json::parse(project.get_json());
    // `write_fps_json` owns the byte layout of an fps.json, the validation and
    // the refusal to write a non-conforming file. A project file is an
    // fps.json with one more section, so it goes through that writer rather
    // than around it.
    write_fps_json(path, positions_json, distances_json, score_sets_json,
                   extra.dump(), validate);
}

IMPBFF_END_NAMESPACE
