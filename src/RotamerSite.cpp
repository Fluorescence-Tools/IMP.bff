/**
 * \file RotamerSite.cpp
 * \brief Placing a rotamer library on a residue: the backbone frame, the
 *        atom selectors, and the library registry.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/RotamerSite.h>

#include <IMP/bff/Scoring.h>
#include <IMP/bff/bff_config.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/json.h>

#include <IMP/exception.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

namespace {

std::string rs_upper(std::string s) {
    for (auto& c : s) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return s;
}

std::string rs_strip(const std::string& s) {
    const std::size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return std::string();
    const std::size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

//! `A48_C1R_cutoff30` -> `A48_C1R`: the stem without the cutoff suffix.
std::string strip_cutoff(const std::string& stem) {
    const std::size_t at = stem.find("_cutoff");
    return at == std::string::npos ? stem : stem.substr(0, at);
}

bool rs_exists(const std::string& path) {
    std::ifstream in(path.c_str());
    return in.good();
}

//! The registry, read once. Throws through the caller when unreadable.
const nlohmann::json& registry() {
    static const nlohmann::json cached = [] {
        const std::string path = get_data_path("rotamer_library") +
                                 "/libraries.json";
        std::ifstream in(path);
        if (!in) {
            IMP_THROW("cannot read the rotamer library registry " << path,
                      IOException);
        }
        return nlohmann::json::parse(in);
    }();
    return cached;
}

//! The registry's entry for a key, with the resolved spelling added.
nlohmann::json entry_for(const std::string& key, const std::string& asked) {
    const nlohmann::json& reg = registry();
    if (!reg.contains(key)) {
        IMP_THROW("Unknown rotamer library '" << asked << "'", ValueException);
    }
    nlohmann::json out = reg[key];
    out["name"] = key;
    out["library_name"] = asked;
    const int cutoff = library_name_cutoff(asked);
    out["cutoff"] = cutoff < 0 ? nlohmann::json() : nlohmann::json(cutoff);
    return out;
}

}  // namespace

// --------------------------------------------------------------------------
// The backbone frame
// --------------------------------------------------------------------------

void resolve_backbone_site(const std::vector<double>& coords,
                           const std::vector<std::string>& atom_names,
                           const std::vector<std::string>& chain_ids,
                           const std::vector<int>& residue_indices,
                           const std::string& chain, int residue,
                           double** out_view, int* n_out_view) {
    const std::string chain_id = rs_upper(rs_strip(chain));
    double found[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    bool have[3] = {false, false, false};
    const std::size_t n = atom_names.size();
    for (std::size_t i = 0; i < n && i * 3 + 2 < coords.size(); i++) {
        const std::string name = rs_upper(rs_strip(atom_names[i]));
        int which = -1;
        if (name == "CA") which = 0;
        else if (name == "N") which = 1;
        else if (name == "C") which = 2;
        else continue;
        if (have[which]) continue;
        if (!chain_id.empty() && i < chain_ids.size()) {
            const std::string frame_chain =
                    rs_upper(rs_strip(chain_ids[i]));
            if (!frame_chain.empty() && frame_chain != chain_id) continue;
        }
        if (i < residue_indices.size() && residue_indices[i] != -1 &&
            residue_indices[i] != residue) {
            continue;
        }
        found[which * 3 + 0] = coords[i * 3];
        found[which * 3 + 1] = coords[i * 3 + 1];
        found[which * 3 + 2] = coords[i * 3 + 2];
        have[which] = true;
    }
    const char* names[3] = {"CA", "N", "C"};
    for (int k = 0; k < 3; k++) {
        if (!have[k]) {
            IMP_THROW("Missing backbone atom " << names[k] << " for chain "
                      << (chain.empty() ? std::string("A") : chain)
                      << " residue " << residue, ValueException);
        }
    }
    double* out = internal::new_double_view(9, out_view, n_out_view);
    if (out == NULL) return;
    for (int k = 0; k < 9; k++) out[k] = found[k];
}

void backbone_rotation(const std::vector<double>& ca,
                       const std::vector<double>& n,
                       const std::vector<double>& c,
                       double** out_view, int* n_out_view) {
    double x[3] = {n[0] - ca[0], n[1] - ca[1], n[2] - ca[2]};
    double nx = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    for (int d = 0; d < 3; d++) x[d] /= nx;
    double yt[3] = {c[0] - ca[0], c[1] - ca[1], c[2] - ca[2]};
    double nt = std::sqrt(yt[0] * yt[0] + yt[1] * yt[1] + yt[2] * yt[2]);
    for (int d = 0; d < 3; d++) yt[d] /= nt;
    double z[3] = {x[1] * yt[2] - x[2] * yt[1],
                   x[2] * yt[0] - x[0] * yt[2],
                   x[0] * yt[1] - x[1] * yt[0]};
    double nz = std::sqrt(z[0] * z[0] + z[1] * z[1] + z[2] * z[2]);
    for (int d = 0; d < 3; d++) z[d] /= nz;
    double y[3] = {z[1] * x[2] - z[2] * x[1],
                   z[2] * x[0] - z[0] * x[2],
                   z[0] * x[1] - z[1] * x[0]};
    double* out = internal::new_double_view(9, out_view, n_out_view);
    if (out == NULL) return;
    for (int d = 0; d < 3; d++) out[d] = x[d];
    for (int d = 0; d < 3; d++) out[3 + d] = y[d];
    for (int d = 0; d < 3; d++) out[6 + d] = z[d];
}

void transform_library_to_site(const std::vector<double>& coords,
                               const std::vector<double>& ca,
                               const std::vector<double>& n,
                               const std::vector<double>& c,
                               double** out_view, int* n_out_view) {
    double* rot = NULL;
    int n_rot = 0;
    backbone_rotation(ca, n, c, &rot, &n_rot);
    double* out = internal::new_double_view(coords.size(), out_view,
                                            n_out_view);
    if (out == NULL || rot == NULL) {
        std::free(rot);
        return;
    }
    const std::size_t n3 = coords.size();
    for (std::size_t i = 0; i < n3; i += 3) {
        for (int d = 0; d < 3; d++) {
            out[i + d] = coords[i] * rot[d] + coords[i + 1] * rot[3 + d] +
                         coords[i + 2] * rot[6 + d] + ca[d];
        }
    }
    std::free(rot);
}

// --------------------------------------------------------------------------
// Atom selectors
// --------------------------------------------------------------------------

std::vector<int> selector_atom_indices(
        const std::vector<std::string>& atom_names,
        const std::vector<std::string>& selectors,
        const std::vector<std::string>& resnames) {
    std::vector<std::string> upper;
    upper.reserve(atom_names.size());
    for (const auto& name : atom_names) upper.push_back(rs_upper(rs_strip(name)));
    const bool check_res = !resnames.empty();
    std::vector<std::string> upper_res;
    if (check_res) {
        upper_res.reserve(resnames.size());
        for (const auto& r : resnames) upper_res.push_back(rs_upper(rs_strip(r)));
    }

    std::vector<int> out;
    out.reserve(selectors.size());
    for (const auto& raw : selectors) {
        // "NAME" or "NAME and resname RES"
        const std::string text = rs_strip(raw);
        std::string wanted = rs_upper(text);
        std::string want_res;
        const std::size_t at = text.find(" and ");
        if (at != std::string::npos) {
            wanted = rs_upper(rs_strip(text.substr(0, at)));
            const std::string tail = text.substr(at + 5);
            // the first `resname TOKEN` in the tail, case-insensitive
            std::istringstream tokens(tail);
            std::string head, value;
            while (tokens >> head >> value) {
                if (rs_upper(head) == "RESNAME") {
                    want_res = rs_upper(value);
                    break;
                }
            }
        }
        int found = -1;
        for (std::size_t i = 0; i < upper.size(); i++) {
            if (upper[i] != wanted) continue;
            if (!want_res.empty() && check_res) {
                if (i >= upper_res.size() || upper_res[i] != want_res) continue;
            }
            found = static_cast<int>(i);
            break;
        }
        if (found < 0) {
            std::string why;
            if (!want_res.empty() && check_res) {
                why = " (residue " + want_res +
                      " not found with that atom name)";
            }
            IMP_THROW("Atom selector '" << text
                      << "' matched no atom in the rotamer library" << why,
                      ValueException);
        }
        out.push_back(found);
    }
    return out;
}

// --------------------------------------------------------------------------
// The bundled library registry
// --------------------------------------------------------------------------

std::string normalize_library_name(const std::string& name) {
    std::string s = rs_strip(name);
    const std::string suffix_mark = " cutoff";
    // strip a trailing `cutoffNN`
    const std::size_t at = s.rfind(suffix_mark);
    if (at != std::string::npos) {
        const std::string tail = s.substr(at + suffix_mark.size());
        bool digits = !tail.empty();
        for (char ch : tail) {
            digits = digits && std::isdigit(static_cast<unsigned char>(ch));
        }
        if (digits) s = rs_strip(s.substr(0, at));
    }
    return s;
}

int library_name_cutoff(const std::string& name) {
    const std::size_t at = name.rfind("cutoff");
    if (at == std::string::npos) return -1;
    const std::string tail = name.substr(at + 6);
    if (tail.empty()) return -1;
    for (char ch : tail) {
        if (!std::isdigit(static_cast<unsigned char>(ch))) return -1;
    }
    return std::atoi(tail.c_str());
}

std::string rotamer_library_metadata(const std::string& name) {
    return entry_for(normalize_library_name(name), name).dump();
}

std::string library_filename(const std::string& metadata_json, int cutoff) {
    const nlohmann::json meta = nlohmann::json::parse(metadata_json);
    const std::string base = meta.value("filename", std::string());
    if (cutoff < 0) return base;
    return strip_cutoff(base) + "_cutoff" + std::to_string(cutoff);
}

std::string resolve_rotamer_library_path(const std::string& name,
                                         const std::string& lib_dir) {
    if (rs_exists(name)) return name;

    const nlohmann::json meta =
            nlohmann::json::parse(rotamer_library_metadata(name));
    const int cutoff = meta.contains("cutoff") && meta["cutoff"].is_number()
            ? meta["cutoff"].get<int>() : -1;
    const std::string filename = library_filename(meta.dump(), cutoff);
    const std::string stem = strip_cutoff(filename);

    // The canonical FRETpredict set first: the requested cutoff is the one
    // loaded here, not whatever the RMF templates happen to hold.
    const std::string data_dir = get_data_path("rotamer_library");
    const std::string canonical = data_dir + "/" + filename + ".bcif";
    if (rs_exists(canonical) &&
        rs_exists(data_dir + "/" + stem + ".pdb")) {
        return canonical;
    }

    const std::string template_dir =
            lib_dir.empty() ? get_template_dir("rotamer") : lib_dir;
    // the Python order: the cutoff-specific file first, then the RMF
    // templates, then the plain PDBs
    const std::string candidates[] = {
            template_dir + "/" + filename + ".bcif",
            template_dir + "/" + filename + ".rmf3",
            template_dir + "/" + stem + ".rmf3",
            template_dir + "/" + stem + ".pdb",
            template_dir + "/" + filename + ".pdb",
    };
    for (const std::string& path : candidates) {
        if (!rs_exists(path)) continue;
        const bool is_stem_rmf = path.size() > 5 &&
                path.substr(path.size() - 5) == ".rmf3" &&
                path.substr(template_dir.size() + 1,
                            path.size() - template_dir.size() - 6) == stem;
        if (is_stem_rmf && cutoff != -1 && cutoff != 30) {
            IMP_THROW("'" << name << "': only the cutoff-30 RMF template "
                      << path << " is available; the cutoff-" << cutoff
                      << " library needs " << filename
                      << ".bcif next to " << stem << ".pdb",
                      IOException);
        }
        return path;
    }
    IMP_THROW("No rotamer library found for '" << name << "'", IOException);
}

IMPBFF_END_NAMESPACE
