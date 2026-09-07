/**
 * \file FPSIO.cpp
 * \brief fps.json, and the legacy C# FPS `.txt` files it replaced.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/internal/Text.h>
#include <IMP/bff/FPSIO.h>

#include <IMP/bff/FPSSchema.h>
#include <IMP/bff/internal/json.h>

#include <IMP/bff/Base.h>

#include <algorithm>
#include <dirent.h>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

IMPBFF_BEGIN_NAMESPACE

// Named, not anonymous: IMP compiles this module as one translation unit.
namespace fps_io {
using IMP::bff::internal::ends_with;
using IMP::bff::internal::file_exists;
using IMP::bff::internal::trimmed;



std::string directory_of(const std::string& path) {
    const std::size_t at = path.find_last_of('/');
    return at == std::string::npos ? std::string() : path.substr(0, at);
}

std::string basename_of(const std::string& path) {
    const std::size_t at = path.find_last_of('/');
    return at == std::string::npos ? path : path.substr(at + 1);
}

std::string stem_lowered(const std::string& path) {
    std::string base = basename_of(path);
    const std::size_t dot = base.find_last_of('.');
    if (dot != std::string::npos) base = base.substr(0, dot);
    for (std::size_t i = 0; i < base.size(); ++i) {
        base[i] = static_cast<char>(std::tolower(base[i]));
    }
    return base;
}


std::vector<std::string> pdbs_in(const std::string& dir) {
    std::vector<std::string> out;
    DIR* handle = opendir(dir.empty() ? "." : dir.c_str());
    if (handle == nullptr) return out;
    while (struct dirent* entry = readdir(handle)) {
        const std::string name(entry->d_name);
        if (!ends_with(name, ".pdb")) continue;
        out.push_back(dir.empty() ? name : dir + "/" + name);
    }
    closedir(handle);
    // `readdir` order is the filesystem's; sorted so the same directory gives
    // the same reading twice.
    std::sort(out.begin(), out.end());
    return out;
}

std::vector<std::string> split_ws(const std::string& line) {
    std::vector<std::string> parts;
    std::istringstream in(line);
    std::string token;
    while (in >> token) parts.push_back(token);
    return parts;
}


//! Chain, residue number and atom name of one atom serial in a PDB.
bool atom_by_serial(const std::string& pdb_path, int serial, std::string& chain,
                    int& res_seq, std::string& atom_name) {
    std::ifstream in(pdb_path.c_str());
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.compare(0, 6, "ATOM  ") != 0 &&
            line.compare(0, 6, "HETATM") != 0) {
            continue;
        }
        if (line.size() < 26) continue;
        const std::string serial_text = trimmed(line.substr(6, 5));
        char* end = nullptr;
        const long current = std::strtol(serial_text.c_str(), &end, 10);
        if (end == serial_text.c_str() || *end != '\0') continue;
        if (current != serial) continue;
        chain = trimmed(line.substr(21, 1));
        const std::string res_text = trimmed(line.substr(22, 4));
        res_seq = res_text.empty() ? 1 : std::atoi(res_text.c_str());
        atom_name = trimmed(line.substr(12, 4));
        return true;
    }
    return false;
}

nlohmann::json parse_or_throw(const std::string& text, const std::string& what) {
    try {
        return nlohmann::json::parse(text);
    } catch (const std::exception& e) {
        IMP_THROW("Cannot parse " << what << " as JSON: " << e.what(),
                  ValueException);
    }
}

void throw_if_invalid(const nlohmann::json& payload, const std::string& lead) {
    const FPSValidation report = fps_schema_validate(payload.dump());
    if (report.errors.empty()) return;
    std::ostringstream message;
    message << lead;
    for (std::size_t i = 0; i < report.errors.size(); ++i) {
        message << "\n  " << report.errors[i];
    }
    IMP_THROW(message.str(), ValueException);
}

}  // namespace fps_io

FPSDocument read_old_lps_txt(const std::string& path,
                             const std::vector<std::string>& pdb_paths) {
    FPSDocument out;
    std::vector<std::string> pdbs = pdb_paths;
    if (pdbs.empty()) pdbs = fps_io::pdbs_in(fps_io::directory_of(path));

    std::map<std::string, std::string> pdb_by_molecule;
    for (std::size_t i = 0; i < pdbs.size(); ++i) {
        pdb_by_molecule[fps_io::stem_lowered(pdbs[i])] = pdbs[i];
    }

    std::ifstream in(path.c_str());
    if (!in) IMP_THROW("Cannot read " << path, IOException);

    nlohmann::json positions = nlohmann::json::object();
    std::string line;
    while (std::getline(in, line)) {
        line = fps_io::trimmed(line);
        if (line.empty() || line[0] == '#') continue;
        const std::vector<std::string> parts = fps_io::split_ws(line);
        if (parts.size() < 5) continue;

        const std::string name = parts[0];
        const std::string molecule = parts[1];
        const std::string av_type = parts[3];

        if (std::find(out.molecules.begin(), out.molecules.end(), molecule) ==
            out.molecules.end()) {
            out.molecules.push_back(molecule);
        }
        const int body_id = static_cast<int>(
                std::find(out.molecules.begin(), out.molecules.end(), molecule) -
                out.molecules.begin());

        double linker_length = 0.0, linker_width = 0.0;
        double r1 = 0.0, r2 = 0.0, r3 = 0.0;
        int atom_serial = -1;
        bool av_is_atom = false;
        if (av_type == "AV1" && parts.size() >= 8) {
            linker_length = std::atof(parts[4].c_str());
            linker_width = std::atof(parts[5].c_str());
            r1 = std::atof(parts[6].c_str());
            atom_serial = std::atoi(parts[7].c_str());
        } else if (av_type == "AV3" && parts.size() >= 10) {
            linker_length = std::atof(parts[4].c_str());
            linker_width = std::atof(parts[5].c_str());
            r1 = std::atof(parts[6].c_str());
            r2 = std::atof(parts[7].c_str());
            r3 = std::atof(parts[8].c_str());
            atom_serial = std::atoi(parts[9].c_str());
        } else if (av_type == "ATOM" && parts.size() >= 5) {
            // FPS's fifth position type (`LabelingPositions.cs:203`): a plain
            // atom, `AVType = None` with `AtomID > 0`. No volume, no
            // coordinate of its own -- and a distance between two of them is
            // what FPS calls a *bond*. It was falling into the "line from a
            // dialect this does not know" branch below and being dropped, so a
            // legacy file's crosslinks silently disappeared.
            atom_serial = std::atoi(parts[4].c_str());
            av_is_atom = true;
        } else if (av_type == "XYZ" && parts.size() >= 7) {
            nlohmann::json fixed;
            fixed["simulation_type"] = "XYZ";
            fixed["x"] = std::atof(parts[4].c_str());
            fixed["y"] = std::atof(parts[5].c_str());
            fixed["z"] = std::atof(parts[6].c_str());
            fixed["body_id"] = body_id;
            positions[name] = fixed;
            continue;
        } else {
            // A line that does not parse is not an error: the format has no
            // version and no schema, so it is a line from a dialect this does
            // not know.
            continue;
        }

        std::string chain;
        int res_seq = 1;
        std::string atom_name = "CA";
        std::map<std::string, std::string>::const_iterator found =
                pdb_by_molecule.find(fps_io::stem_lowered(molecule));
        bool resolved = false;
        if (atom_serial >= 0 && found != pdb_by_molecule.end() &&
            fps_io::file_exists(found->second)) {
            resolved = fps_io::atom_by_serial(found->second, atom_serial, chain,
                                              res_seq, atom_name);
            if (!resolved) { chain.clear(); res_seq = 1; atom_name = "CA"; }
        }
        if (!resolved && atom_serial >= 0 && found == pdb_by_molecule.end()) {
            // No PDB for this molecule: carry the serial as a proxy, which is
            // what the caller gets and can still key on.
            res_seq = atom_serial;
        }

        nlohmann::json position;
        position["chain_identifier"] = chain;
        position["residue_seq_number"] = res_seq;
        position["atom_name"] = atom_name;
        position["simulation_type"] = av_type;
        position["body_id"] = body_id;
        if (av_is_atom) {
            // An ATOM position has no volume, so it carries no AV parameters:
            // writing zeroed ones would validate as an AV1 of radius 0.
            positions[name] = position;
            continue;
        }
        position["linker_length"] = linker_length;
        position["linker_width"] = linker_width;
        position["radius1"] = r1;
        position["radius2"] = r2;
        position["radius3"] = r3;
        position["simulation_grid_resolution"] = 1.5;
        positions[name] = position;
    }
    out.positions = positions.dump();
    return out;
}

std::string read_old_distances_txt(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) IMP_THROW("Cannot read " << path, IOException);

    nlohmann::json distances = nlohmann::json::object();
    std::string distance_type = "RDAMean";

    std::string first;
    std::getline(in, first);
    first = fps_io::trimmed(first);
    const std::vector<std::string> head = fps_io::split_ws(first);
    if (head.size() == 1) {
        distance_type = head[0];
    } else if (!head.empty()) {
        // Not a lone type token: the first line is already a distance.
        in.clear();
        in.seekg(0);
    }

    std::string line;
    while (std::getline(in, line)) {
        line = fps_io::trimmed(line);
        if (line.empty() || line[0] == '#') continue;
        const std::vector<std::string> parts = fps_io::split_ws(line);
        if (parts.size() < 5) continue;

        nlohmann::json d;
        d["position1_name"] = parts[0];
        d["position2_name"] = parts[1];
        d["distance"] = std::atof(parts[2].c_str());
        d["error_neg"] = std::atof(parts[3].c_str());
        d["error_pos"] = std::atof(parts[4].c_str());
        d["distance_type"] = distance_type;
        d["Forster_radius"] =
                parts.size() >= 6 ? std::atof(parts[5].c_str()) : 52.0;
        distances[parts[0] + "_" + parts[1]] = d;
    }
    return distances.dump();
}

std::vector<std::string> read_old_distances_order(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) IMP_THROW("Cannot read " << path, IOException);

    std::vector<std::string> out;
    std::string first;
    std::getline(in, first);
    first = fps_io::trimmed(first);
    const std::vector<std::string> head = fps_io::split_ws(first);
    if (head.size() != 1 && !head.empty()) {
        // Not a lone type token: the first line is already a distance.
        in.clear();
        in.seekg(0);
    }

    std::string line;
    while (std::getline(in, line)) {
        line = fps_io::trimmed(line);
        if (line.empty() || line[0] == '#') continue;
        const std::vector<std::string> parts = fps_io::split_ws(line);
        if (parts.size() < 5) continue;
        out.push_back(parts[0] + "_" + parts[1]);
    }
    return out;
}

FPSDocument read_fps_json(const std::string& path,
                          const std::vector<std::string>& pdb_paths,
                          bool validate) {
    if (!fps_io::ends_with(path, ".json")) {
        FPSDocument out = read_old_lps_txt(path, pdb_paths);
        const std::string dir = fps_io::directory_of(path);
        const std::string beside =
                dir.empty() ? std::string("Distances.txt")
                            : dir + "/Distances.txt";
        if (fps_io::file_exists(beside)) {
            out.distances = read_old_distances_txt(beside);
        }
        return out;
    }

    std::ifstream in(path.c_str());
    if (!in) IMP_THROW("Cannot read " << path, IOException);
    std::ostringstream buffer;
    buffer << in.rdbuf();
    nlohmann::json payload = fps_io::parse_or_throw(buffer.str(), path);

    if (validate) {
        fps_io::throw_if_invalid(
                payload, path + " does not conform to the fps.json schema:");
    }

    FPSDocument out;
    if (payload.contains("Positions")) {
        // A document-level `strip_mask` is folded into every position, so a
        // position taken from this document stands on its own. The document's
        // own key stays in `extra`, unchanged, for a writer to round-trip.
        const std::string document_mask =
                payload.contains("strip_mask") && payload["strip_mask"].is_string()
                        ? payload["strip_mask"].get<std::string>()
                        : std::string();
        nlohmann::json positions = payload["Positions"];
        if (!document_mask.empty()) {
            for (auto it = positions.begin(); it != positions.end(); ++it) {
                const std::string own =
                        it.value().contains("strip_mask") &&
                                        it.value()["strip_mask"].is_string()
                                ? it.value()["strip_mask"].get<std::string>()
                                : std::string();
                it.value()["strip_mask"] =
                        combined_strip_mask(document_mask, own);
            }
        }
        out.positions = positions.dump();
        payload.erase("Positions");
    }
    if (payload.contains("Distances")) {
        out.distances = payload["Distances"].dump();
        payload.erase("Distances");
    }
    if (payload.contains("χ²")) {
        out.score_sets = payload["χ²"].dump();
        payload.erase("χ²");
    }
    out.extra = payload.dump();
    return out;
}

void write_fps_json(const std::string& path, const std::string& positions_json,
                    const std::string& distances_json,
                    const std::string& score_sets_json,
                    const std::string& extra_json, bool validate) {
    nlohmann::json payload = fps_io::parse_or_throw(extra_json, "the extra keys");
    if (!payload.is_object()) payload = nlohmann::json::object();
    payload["Positions"] = fps_io::parse_or_throw(positions_json, "the positions");
    payload["Distances"] = fps_io::parse_or_throw(distances_json, "the distances");
    const nlohmann::json sets =
            fps_io::parse_or_throw(score_sets_json, "the score sets");
    if (!sets.is_null() && !(sets.is_object() && sets.empty())) {
        payload["χ²"] = sets;
    }
    if (validate) {
        fps_io::throw_if_invalid(payload,
                                 "refusing to write a non-conforming fps.json:");
    }
    std::ofstream out(path.c_str());
    if (!out) IMP_THROW("Cannot write " << path, IOException);
    out << payload.dump(2) << "\n";
}

FPSDocument fps_positions_for_docking(const std::string& positions_json,
                                      const std::string& distances_json) {
    const nlohmann::json positions =
            fps_io::parse_or_throw(positions_json, "the positions");
    const std::vector<std::string> allowed = fps_av_simulation_types();

    nlohmann::json kept = nlohmann::json::object();
    for (nlohmann::json::const_iterator it = positions.begin();
         it != positions.end(); ++it) {
        std::string stype = "AV1";
        if (it.value().is_object() && it.value().contains("simulation_type")) {
            const nlohmann::json& v = it.value()["simulation_type"];
            stype = v.is_string() ? v.get<std::string>() : v.dump();
        }
        if (std::find(allowed.begin(), allowed.end(), stype) != allowed.end()) {
            kept[it.key()] = it.value();
        }
    }

    FPSDocument out;
    out.positions = kept.dump();
    const nlohmann::json distances =
            fps_io::parse_or_throw(distances_json, "the distances");
    nlohmann::json kept_d = nlohmann::json::object();
    for (nlohmann::json::const_iterator it = distances.begin();
         it != distances.end(); ++it) {
        if (!it.value().is_object()) continue;
        const char* keys[] = {"position1_name", "position2_name"};
        bool both = true;
        for (int i = 0; i < 2 && both; ++i) {
            both = it.value().contains(keys[i]) &&
                   it.value()[keys[i]].is_string() &&
                   kept.contains(it.value()[keys[i]].get<std::string>());
        }
        if (both) kept_d[it.key()] = it.value();
    }
    out.distances = kept_d.dump();
    return out;
}

std::string read_evaluators_json(const std::string& path) {
    if (!fps_io::file_exists(path)) return "[]";
    std::ifstream in(path.c_str());
    if (!in) return "[]";
    std::ostringstream buffer;
    buffer << in.rdbuf();
    nlohmann::json payload;
    try {
        payload = nlohmann::json::parse(buffer.str());
    } catch (const std::exception&) {
        // An evaluator list is optional, and a caller asking for one is asking
        // whether there is one.
        return "[]";
    }
    if (!payload.is_object() || !payload.contains("Evaluators") ||
        !payload["Evaluators"].is_array()) {
        return "[]";
    }
    nlohmann::json out = nlohmann::json::array();
    const nlohmann::json& list = payload["Evaluators"];
    for (std::size_t i = 0; i < list.size(); ++i) {
        if (list[i].is_object()) out.push_back(list[i]);
    }
    return out.dump();
}

void write_evaluators_json(const std::string& path,
                           const std::string& evaluators_json) {
    nlohmann::json payload = nlohmann::json::object();
    if (fps_io::file_exists(path)) {
        std::ifstream in(path.c_str());
        std::ostringstream buffer;
        buffer << in.rdbuf();
        try {
            payload = nlohmann::json::parse(buffer.str());
        } catch (const std::exception&) {
            payload = nlohmann::json::object();
        }
        if (!payload.is_object()) payload = nlohmann::json::object();
    }
    payload["Evaluators"] =
            fps_io::parse_or_throw(evaluators_json, "the evaluators");
    std::ofstream out(path.c_str());
    if (!out) IMP_THROW("Cannot write " << path, IOException);
    out << payload.dump(2) << "\n";
}

IMPBFF_END_NAMESPACE
