/**
 * \file Rotamer.cpp
 * \brief The rotamer layer: energy, site, ensemble, fps payload, FRET driver.
 *
 * Sections in the order of IMP/bff/Rotamer.h; each is marked with the file
 * it came from.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

// -------- from RotamerEnergy.cpp --------
/**
 * (formerly RotamerEnergy.cpp, now a section of this file)
 * \brief The interaction energy of each rotamer with its protein.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/Rotamer.h>
#include <IMP/bff/internal/OutputView.h>

#include <algorithm>
#include <cmath>

IMPBFF_BEGIN_NAMESPACE

std::vector<double> rotamer_interaction_energies(
        double* rotamer_coords, int n_rotamer_coords,
        double* protein_coords, int n_protein_coords,
        double* rmin_ij, int n_rmin_ij,
        double* eps_ij, int n_eps_ij,
        const std::vector<double>& q_probe,
        const std::vector<double>& q_protein,
        int n_rotamers, int n_probe_atoms, int n_protein_atoms,
        int potential,
        double steric_cutoff, double coulomb_cutoff,
        double debye_length, double coulomb_prefactor) {
    std::vector<double> out(static_cast<std::size_t>(std::max(0, n_rotamers)) * 2, 0.0);
    if (n_rotamers <= 0 || n_probe_atoms <= 0 || n_protein_atoms <= 0) return out;
    (void)n_rotamer_coords; (void)n_protein_coords;
    (void)n_rmin_ij; (void)n_eps_ij;

    const bool electrostatic =
            static_cast<int>(q_probe.size()) == n_probe_atoms &&
            static_cast<int>(q_protein.size()) == n_protein_atoms;
    const double steric_cut2 = steric_cutoff * steric_cutoff;
    const double coulomb_cut2 = coulomb_cutoff * coulomb_cutoff;

#pragma omp parallel for schedule(static)
    for (int r = 0; r < n_rotamers; ++r) {
        double steric = 0.0, coulomb = 0.0;
        const double* conf =
                rotamer_coords + static_cast<std::size_t>(r) * n_probe_atoms * 3;
        for (int a = 0; a < n_probe_atoms; ++a) {
            const double ax = conf[3 * a + 0];
            const double ay = conf[3 * a + 1];
            const double az = conf[3 * a + 2];
            const double qa = electrostatic ? q_probe[a] : 0.0;
            const double* rmin_row =
                    rmin_ij + static_cast<std::size_t>(a) * n_protein_atoms;
            const double* eps_row =
                    eps_ij + static_cast<std::size_t>(a) * n_protein_atoms;
            for (int b = 0; b < n_protein_atoms; ++b) {
                const double dx = ax - protein_coords[3 * b + 0];
                const double dy = ay - protein_coords[3 * b + 1];
                const double dz = az - protein_coords[3 * b + 2];
                const double d2 = dx * dx + dy * dy + dz * dz;

                // Squared distances until a cutoff passes: the square root is
                // the expensive part and most pairs never need it.
                if (d2 < steric_cut2) {
                    const double d = std::sqrt(d2);
                    const double ratio = rmin_row[b] / d;
                    if (potential == ROTAMER_POTENTIAL_GAUSS) {
                        steric += eps_row[b] * std::exp(-0.5 / (ratio * ratio));
                    } else {
                        const double r6 = ratio * ratio * ratio;
                        const double ratio6 = r6 * r6;
                        steric += eps_row[b] * (ratio6 * ratio6 - 2.0 * ratio6);
                    }
                }
                if (electrostatic && qa != 0.0 && q_protein[b] != 0.0 &&
                    d2 < coulomb_cut2) {
                    const double d = std::sqrt(d2);
                    coulomb += qa * q_protein[b] * coulomb_prefactor / d *
                               std::exp(-d / debye_length);
                }
            }
        }
        out[2 * r + 0] = steric;
        out[2 * r + 1] = coulomb;
    }
    return out;
}

namespace {
//! Padded axis-aligned bounding box of one conformer: xmin,ymin,zmin,xmax,...
inline void conformer_box(const double* xyz, int n_atoms, double pad, double* box) {
    for (int k = 0; k < 3; ++k) { box[k] = xyz[k]; box[3 + k] = xyz[k]; }
    for (int a = 1; a < n_atoms; ++a) {
        for (int k = 0; k < 3; ++k) {
            const double v = xyz[3 * a + k];
            if (v < box[k]) box[k] = v;
            if (v > box[3 + k]) box[3 + k] = v;
        }
    }
    for (int k = 0; k < 3; ++k) { box[k] -= pad; box[3 + k] += pad; }
}

inline bool boxes_overlap(const double* a, const double* b) {
    for (int k = 0; k < 3; ++k) {
        if (a[3 + k] < b[k] || b[3 + k] < a[k]) return false;
    }
    return true;
}
}  // namespace

namespace {
std::vector<double> rotamer_pair_energy_matrix_impl(
        const std::vector<double>& coords_a, const std::vector<double>& coords_b,
        const std::vector<double>& rmin, const std::vector<double>& eps,
        int n_a_conf, int n_a_atoms, int n_b_conf, int n_b_atoms,
        double r_cutoff, double aabb_pad, double r_floor) {
    std::vector<double> out(
            static_cast<std::size_t>(std::max(0, n_a_conf)) * std::max(0, n_b_conf), 0.0);
    if (n_a_conf <= 0 || n_b_conf <= 0 || n_a_atoms <= 0 || n_b_atoms <= 0) return out;
    if (rmin.size() != static_cast<std::size_t>(n_a_atoms) * n_b_atoms ||
        eps.size() != rmin.size()) {
        return out;   // no parameters means no interaction
    }

    std::vector<double> box_a(static_cast<std::size_t>(n_a_conf) * 6);
    std::vector<double> box_b(static_cast<std::size_t>(n_b_conf) * 6);
    for (int i = 0; i < n_a_conf; ++i) {
        conformer_box(&coords_a[static_cast<std::size_t>(i) * n_a_atoms * 3],
                      n_a_atoms, aabb_pad, &box_a[static_cast<std::size_t>(i) * 6]);
    }
    for (int j = 0; j < n_b_conf; ++j) {
        conformer_box(&coords_b[static_cast<std::size_t>(j) * n_b_atoms * 3],
                      n_b_atoms, aabb_pad, &box_b[static_cast<std::size_t>(j) * 6]);
    }

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n_a_conf; ++i) {
        const double* ca = &coords_a[static_cast<std::size_t>(i) * n_a_atoms * 3];
        for (int j = 0; j < n_b_conf; ++j) {
            if (!boxes_overlap(&box_a[static_cast<std::size_t>(i) * 6],
                               &box_b[static_cast<std::size_t>(j) * 6])) continue;
            const double* cb = &coords_b[static_cast<std::size_t>(j) * n_b_atoms * 3];
            double e = 0.0;
            for (int a = 0; a < n_a_atoms; ++a) {
                const double ax = ca[3 * a + 0], ay = ca[3 * a + 1], az = ca[3 * a + 2];
                const std::size_t row = static_cast<std::size_t>(a) * n_b_atoms;
                for (int b = 0; b < n_b_atoms; ++b) {
                    const double dx = ax - cb[3 * b + 0];
                    const double dy = ay - cb[3 * b + 1];
                    const double dz = az - cb[3 * b + 2];
                    const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
                    const double rm = rmin[row + b];
                    // Repulsive-only, and inside the cutoff. Both tests are on
                    // the unclamped distance.
                    if (!(d < rm) || !(d < r_cutoff)) continue;
                    const double safe = d > r_floor ? d : r_floor;
                    const double ratio = rm / safe;
                    const double r6 = ratio * ratio * ratio;
                    const double ratio6 = r6 * r6;
                    e += eps[row + b] * (ratio6 * ratio6 - 2.0 * ratio6);
                }
            }
            out[static_cast<std::size_t>(i) * n_b_conf + j] = e;
        }
    }
    return out;
}
}  // namespace

void pair_energy_matrix_kernel(const std::vector<double>& coords_a, const std::vector<double>& coords_b,
        const std::vector<double>& rmin, const std::vector<double>& eps,
        int n_a_conf, int n_a_atoms, int n_b_conf, int n_b_atoms,
        double** out_view, int* n_out_view,
        double r_cutoff, double aabb_pad, double r_floor) {
    internal::copy_to_view(rotamer_pair_energy_matrix_impl(coords_a, coords_b, rmin, eps, n_a_conf, n_a_atoms, n_b_conf, n_b_atoms, r_cutoff, aabb_pad, r_floor),
                           out_view, n_out_view);
}

std::vector<double> lj_pair_energies(
        double* coords, int n_coords,
        const std::vector<int>& index_a, const std::vector<int>& index_b,
        const std::vector<double>& rmin, const std::vector<double>& eps,
        int n_frames, int n_atoms, int n_pairs,
        bool repulsive_only, double r_floor) {
    std::vector<double> out(static_cast<std::size_t>(std::max(0, n_frames)), 0.0);
    if (n_frames <= 0 || n_pairs <= 0 || n_atoms <= 0) return out;
    (void)n_coords;

#pragma omp parallel for schedule(static)
    for (int f = 0; f < n_frames; ++f) {
        const double* xyz = coords + static_cast<std::size_t>(f) * n_atoms * 3;
        double e = 0.0;
        for (int p = 0; p < n_pairs; ++p) {
            const int ia = index_a[p], ib = index_b[p];
            const double dx = xyz[3 * ia + 0] - xyz[3 * ib + 0];
            const double dy = xyz[3 * ia + 1] - xyz[3 * ib + 1];
            const double dz = xyz[3 * ia + 2] - xyz[3 * ib + 2];
            const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
            // Both tests on the unclamped distance; only the ratio uses
            // the floor.
            if (repulsive_only && !(d < rmin[p])) continue;
            const double safe = d > r_floor ? d : r_floor;
            const double ratio = rmin[p] / safe;
            const double r6 = ratio * ratio * ratio;
            const double ratio6 = r6 * r6;
            e += eps[p] * (ratio6 * ratio6 - 2.0 * ratio6);
        }
        out[f] = e;
    }
    return out;
}

IMPBFF_END_NAMESPACE

// -------- from RotamerSite.cpp --------
/**
 * (formerly RotamerSite.cpp, now a section of this file)
 * \brief Placing a rotamer library on a residue: the backbone frame, the
 *        atom selectors, and the library registry.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/RotamerLibrary.h>
#include <IMP/bff/ProbeSampling.h>
#include <IMP/bff/internal/Text.h>

#include <IMP/bff/Scoring.h>
#include <IMP/bff/bff_config.h>
#include <IMP/bff/internal/json.h>

#include <IMP/bff/Base.h>

#include <cctype>
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

//! The family containers a rotamer library may be filed in, in the order
//! they are consulted. Named rather than globbed: a directory listing is the
//! one thing this lookup cannot do portably, and three names are the whole
//! shipped set. A container that is not there is simply skipped.
const char* const kFamilies[] = {"dyes", "spinlabels", "sidechains"};

//! `dyes.drot.pto::A48_C1R_cutoff10` when a family container holds
//! `library`, otherwise empty. Only each catalog is read.
std::string library_in_family(const std::string& data_dir,
                              const std::string& library) {
    for (std::size_t i = 0; i < sizeof(kFamilies) / sizeof(kFamilies[0]); ++i) {
        const std::string path =
                data_dir + "/" + kFamilies[i] + ".drot.pto";
        if (!rs_exists(path)) continue;
        const std::vector<std::string> listed = drot_catalog(path);
        for (std::size_t j = 0; j < listed.size(); ++j) {
            if (listed[j] == library) return path + "::" + library;
        }
    }
    return std::string();
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

}  // namespace

std::string rotamer_library_registry() { return registry().dump(); }

namespace {

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

std::string rotamer_library_metadata_for_path(const std::string& path) {
    // basename, extension off
    std::string stem = path;
    const std::size_t slash = stem.find_last_of("/\\");
    if (slash != std::string::npos) stem = stem.substr(slash + 1);
    const std::size_t dot = stem.find_last_of('.');
    if (dot != std::string::npos) stem = stem.substr(0, dot);

    const int cutoff = library_name_cutoff(stem);
    const std::string base = strip_cutoff(stem);
    const nlohmann::json& reg = registry();
    for (nlohmann::json::const_iterator it = reg.begin(); it != reg.end();
         ++it) {
        const std::string filename = it.value().value("filename",
                                                      std::string());
        if (strip_cutoff(filename) != base) continue;
        nlohmann::json out = it.value();
        out["name"] = it.key();
        out["library_name"] = it.key();
        out["cutoff"] = cutoff < 0 ? nlohmann::json() : nlohmann::json(cutoff);
        return out.dump();
    }
    return nlohmann::json::object().dump();
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
    // loaded here, not whatever the RMF templates happen to hold. Within it
    // `.drot` wins over `.bcif` -- it is the shipped form (PRD-118) and it
    // carries its own template, so unlike the frame store it needs no `.pdb`
    // beside it.
    //
    // The libraries ship one container per family, so the answer for a dye is
    // usually a locator into `dyes.drot.pto` rather than a file of its own.
    // A standalone container still wins when one is there: that is how a user
    // drops their own library in beside the shipped set.
    const std::string data_dir = get_data_path("rotamer_library");
    const std::string pto = data_dir + "/" + filename + ".drot.pto";
    if (rs_exists(pto)) return pto;
    const std::string drot = data_dir + "/" + filename + ".drot";
    if (rs_exists(drot)) return drot;
    const std::string in_family = library_in_family(data_dir, filename);
    if (!in_family.empty()) return in_family;
    const std::string canonical = data_dir + "/" + filename + ".bcif";
    if (rs_exists(canonical) &&
        rs_exists(data_dir + "/" + stem + ".pdb")) {
        return canonical;
    }

    const std::string template_dir =
            lib_dir.empty() ? get_template_dir("rotamer") : lib_dir;
    // most specific first: the cutoff-specific file, then the RMF
    // templates, then the plain PDBs
    const std::string candidates[] = {
            template_dir + "/" + filename + ".drot.pto",
            template_dir + "/" + filename + ".drot",
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

// --------------------------------------------------------------------------
// Reading a library
// --------------------------------------------------------------------------

namespace {

//! The directory part of a path, or "." when it has none.
std::string rs_dirname(const std::string& path) {
    const std::size_t slash = path.find_last_of('/');
    if (slash == std::string::npos) return std::string(".");
    return path.substr(0, slash);
}

//! The file part of a path.
std::string rs_basename(const std::string& path) {
    const std::size_t slash = path.find_last_of('/');
    if (slash == std::string::npos) return path;
    return path.substr(slash + 1);
}

//! A library file's stem: its name without a cutoff suffix or extensions.
/*! `A48_C1R_cutoff10.drot.pto` -> `A48_C1R`, which is what the shared `.pdb`
    template beside it is named. */
std::string rs_library_stem(const std::string& file_name) {
    std::string stem = file_name;
    const std::size_t dot = stem.find('.');
    if (dot != std::string::npos) stem = stem.substr(0, dot);
    const std::size_t cut = stem.find("_cutoff");
    if (cut != std::string::npos) stem = stem.substr(0, cut);
    return stem;
}

//! The residue names of a PDB's atom records, in file order.
std::vector<std::string> rs_resnames_from_pdb(const std::string& path) {
    std::vector<std::string> out;
    std::ifstream in(path);
    if (!in) return out;
    std::string line;
    while (std::getline(in, line)) {
        if (line.compare(0, 6, "ATOM  ") != 0 &&
            line.compare(0, 6, "HETATM") != 0) {
            continue;
        }
        if (line.size() < 20) continue;
        out.push_back(internal::trimmed(line.substr(17, 3)));
        if (out.size() >= 10000) break;
    }
    return out;
}

//! The bundled `<stem>.pdb` a registry entry's `filename` names, or "".
std::string rs_pdb_for_metadata(const nlohmann::json& meta) {
    if (!meta.is_object() || !meta.contains("filename")) return std::string();
    const std::string filename = meta["filename"].is_string()
                                         ? meta["filename"].get<std::string>()
                                         : std::string();
    const std::string stem = rs_library_stem(filename);
    if (stem.empty()) return std::string();
    const std::string dir = get_data_path("rotamer_library");
    const std::string candidates[2] = {dir + "/" + stem + ".pdb",
                                       dir + "/" + filename + ".pdb"};
    for (int i = 0; i < 2; ++i) {
        if (rs_exists(candidates[i])) return candidates[i];
    }
    return std::string();
}

}  // namespace

std::vector<std::string> infer_rotamer_resnames(
        const std::vector<std::string>& atom_names,
        const std::string& metadata_json) {
    std::vector<std::string> out;
    nlohmann::json meta = nlohmann::json::parse(metadata_json, NULL, false);
    if (meta.is_discarded() || !meta.is_object()) return out;

    // The dye's residue name is whatever its own selectors name.
    std::vector<std::string> selectors;
    const char* keys[4] = {"mu", "r", "positive", "negative"};
    for (int k = 0; k < 4; ++k) {
        if (!meta.contains(keys[k]) || meta[keys[k]].is_null()) continue;
        const nlohmann::json& v = meta[keys[k]];
        if (v.is_string()) {
            selectors.push_back(v.get<std::string>());
        } else if (v.is_array()) {
            for (nlohmann::json::const_iterator it = v.begin(); it != v.end();
                 ++it) {
                if (it->is_string()) selectors.push_back(it->get<std::string>());
            }
        }
    }
    const std::vector<std::string> resnames = selector_resnames(selectors);
    if (resnames.empty()) return out;
    const std::string dye_resname = resnames[0];

    // The linker's is the registry name's suffix: `A48_C1R` -> `C1R`.
    std::string linker_resname;
    const std::string name = meta.contains("name") && meta["name"].is_string()
                                     ? meta["name"].get<std::string>()
                                     : std::string();
    const std::size_t underscore = name.find_last_of('_');
    if (underscore != std::string::npos && underscore + 1 < name.size()) {
        const std::string tail = name.substr(underscore + 1);
        const bool shaped =
                (tail.size() == 3 || tail.size() == 4) &&
                std::isupper(static_cast<unsigned char>(tail[0])) &&
                tail[tail.size() - 1] == 'R' &&
                std::isupper(static_cast<unsigned char>(tail[tail.size() - 2]));
        if (shaped) linker_resname = tail;
    }
    if (linker_resname.empty()) {
        out.assign(atom_names.size(), dye_resname);
        return out;
    }

    // The linker's atoms, which the FRETpredict libraries name identically in
    // every dye. An atom outside the set is the dye's.
    static const char* kLinkerAtoms[] = {
            "CA", "HA", "C",  "O",  "C6",  "H10", "H11", "S1",  "C7",
            "C8", "H12", "C9", "O3", "N3",  "C10", "O4",  "C11", "H13",
            "H14", "C12", "H15", "H16", "C13", "H17", "H18", "C14", "H19",
            "H20", "C15", "H21", "H22", "N99", "H23", "N",   "H",   "HX2",
            "HX3"};
    const std::size_t n_linker = sizeof(kLinkerAtoms) / sizeof(kLinkerAtoms[0]);
    out.reserve(atom_names.size());
    for (std::size_t i = 0; i < atom_names.size(); ++i) {
        bool is_linker = false;
        for (std::size_t j = 0; j < n_linker; ++j) {
            if (atom_names[i] == kLinkerAtoms[j]) {
                is_linker = true;
                break;
            }
        }
        out.push_back(is_linker ? linker_resname : dye_resname);
    }
    return out;
}

RotamerLibrary load_rotamer_library(const std::string& name,
                                    const std::string& lib_dir) {
    const std::string locator = resolve_rotamer_library_path(name, lib_dir);
    // A locator addresses one library inside a family container
    // (`dyes.drot.pto::A48_C1R_cutoff10`); a plain path is a container of its
    // own. `read_drot` takes either, but everything here that reasons about
    // the file -- its suffix, its sidecars -- wants the container.
    std::string container = locator, inside;
    const std::size_t sep = locator.find("::");
    if (sep != std::string::npos) {
        container = locator.substr(0, sep);
        inside = locator.substr(sep + 2);
    }

    // A path names its dye only through the registry: the coordinates are in
    // the file, the transition-dipole and attachment selectors are not.
    const std::string metadata = rs_exists(name)
                                         ? rotamer_library_metadata_for_path(
                                                   container)
                                         : rotamer_library_metadata(name);
    nlohmann::json meta = nlohmann::json::parse(metadata, NULL, false);
    if (meta.is_discarded()) meta = nlohmann::json::object();

    const std::string dir = rs_dirname(container);
    const std::string file_name = rs_basename(container);
    RotamerLibrary lib;
    std::string sidecar_pdb;

    if (internal::ends_with(file_name, ".drot") ||
        internal::ends_with(file_name, ".drot.pto")) {
        // The shipped form (PRD-118): the template, the Z-matrix and the
        // per-conformer internal coordinates in one container, so nothing
        // beside it is read -- names, residues and elements are all in it.
        lib = read_drot(locator);
        // `.drot` keeps the weights the encoder found (cluster populations);
        // this loader's contract is the normalised form.
        normalize_weights_in_place(lib.weights);
        sidecar_pdb = dir + "/" +
                      rs_library_stem(inside.empty() ? file_name : inside) +
                      ".pdb";
    } else if (internal::ends_with(file_name, ".bcif") ||
               internal::ends_with(file_name, ".dcd")) {
        // FRETpredict's library set: `<stem>.pdb` for the names and residues,
        // frames beside it, and per-rotamer weights in a text file.
        sidecar_pdb = dir + "/" + rs_library_stem(file_name) + ".pdb";
        std::string stem_with_cutoff = file_name;
        const std::size_t dot = stem_with_cutoff.find_last_of('.');
        if (dot != std::string::npos) {
            stem_with_cutoff = stem_with_cutoff.substr(0, dot);
        }
        const std::string weights = dir + "/" + stem_with_cutoff +
                                    "_weights.txt";
        lib = load_rotamer_library_trajectory(
                sidecar_pdb, container, rs_exists(weights) ? weights : "");
    } else {
        IMP_THROW("Unsupported rotamer library file: "
                          << container
                          << " (an .rmf3 template is read through "
                             "read_rotamer_library_rmf, which needs IMP.rmf)",
                  ValueException);
    }

    lib.path = locator;
    lib.metadata = meta.dump();
    if (lib.resnames.empty() && rs_exists(sidecar_pdb)) {
        lib.resnames = rs_resnames_from_pdb(sidecar_pdb);
    }
    if (lib.resnames.empty()) {
        const std::string bundled = rs_pdb_for_metadata(meta);
        if (!bundled.empty()) lib.resnames = rs_resnames_from_pdb(bundled);
    }
    if (lib.resnames.empty()) {
        lib.resnames = infer_rotamer_resnames(lib.atom_names, lib.metadata);
    }
    return lib;
}

IMPBFF_END_NAMESPACE

// -------- from RotamerEnsemble.cpp --------
/**
 * (formerly RotamerEnsemble.cpp, now a section of this file)
 * \brief A rotamer library placed and screened at one labelling site.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/ProbeLibrary.h>
#include <IMP/bff/FPS.h>


#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

const char* const SIMULATION_TYPE_R1 = "R1";

namespace {

//! A metadata selector entry, which the registry writes as a string or a list.
std::vector<std::string> selector_list(const nlohmann::json& meta,
                                       const char* key) {
    std::vector<std::string> out;
    if (!meta.is_object() || !meta.contains(key) || meta[key].is_null()) {
        return out;
    }
    const nlohmann::json& v = meta[key];
    if (v.is_string()) {
        out.push_back(v.get<std::string>());
    } else if (v.is_array()) {
        for (nlohmann::json::const_iterator it = v.begin(); it != v.end();
             ++it) {
            if (it->is_string()) out.push_back(it->get<std::string>());
        }
    }
    return out;
}

//! `str()` of a number, spelled as Python spells it -- `298.15`, not
//! `298.150000`. The params map is provenance a reader reads back.
std::string number_text(double v) {
    std::ostringstream s;
    s << v;
    return s.str();
}

std::string bool_text(bool v) { return v ? "True" : "False"; }

//! The `(n, 4)` points and the `(n, 3)` dipoles of any states.
/*! The dipoles come back empty when the states carry none *or* carry a number
    that does not match their points -- an AV cloud has none, and a mismatched
    orientation array cannot be paired with a point by index, which is the only
    way \f$\kappa^2\f$ can use it. */
void states_arrays(const States& s, std::vector<double>& points,
                   std::vector<double>& mu) {
    internal::OwnedView p, o;
    s.get_points(&p.data, &p.size);
    s.get_orientations(&o.data, &o.size);
    points = p.vector();
    mu.clear();
    const std::size_t n = points.size() / 4;
    if (o.size > 0 && static_cast<std::size_t>(o.size) == n * 3) {
        mu = o.vector();
    }
}

}  // namespace

// --------------------------------------------------------------------------
// RotamerEnsemble
// --------------------------------------------------------------------------

RotamerEnsemble::RotamerEnsemble(
        const std::vector<double>& points,
        const std::vector<double>& attachment_point,
        const std::vector<double>& orientations,
        const std::string& position_name,
        const std::map<std::string, std::string>& params,
        const std::vector<double>& atoms,
        const std::vector<std::string>& atom_names,
        const std::vector<std::string>& resnames,
        const std::vector<double>& energies, double partition,
        const std::string& library, const std::string& chain, int residue)
    : States(points, attachment_point, orientations, position_name, params),
      atoms_(atoms), atom_names_(atom_names), resnames_(resnames),
      energies_(energies), partition_(partition), library_(library),
      chain_(chain), residue_(residue) {}

void RotamerEnsemble::get_atoms(double** out_view, int* n_out_view) const {
    internal::copy_to_view(atoms_, out_view, n_out_view);
}

void RotamerEnsemble::get_energies(double** out_view, int* n_out_view) const {
    internal::copy_to_view(energies_, out_view, n_out_view);
}

void RotamerEnsemble::get_centres(double** out_view, int* n_out_view) const {
    const int n = get_n_points();
    double* buffer = internal::new_double_view(static_cast<std::size_t>(n) * 3,
                                               out_view, n_out_view);
    if (buffer == NULL) return;
    for (int i = 0; i < n; ++i) {
        buffer[i * 3] = points_[i * 4];
        buffer[i * 3 + 1] = points_[i * 4 + 1];
        buffer[i * 3 + 2] = points_[i * 4 + 2];
    }
}

void RotamerEnsemble::get_weights(double** out_view, int* n_out_view) const {
    const int n = get_n_points();
    double* buffer = internal::new_double_view(static_cast<std::size_t>(n),
                                               out_view, n_out_view);
    if (buffer == NULL) return;
    for (int i = 0; i < n; ++i) buffer[i] = points_[i * 4 + 3];
}

int RotamerEnsemble::get_n_atoms() const {
    const int n = get_n_points();
    if (n == 0) return 0;
    return static_cast<int>(atoms_.size() / (static_cast<std::size_t>(n) * 3));
}

double RotamerEnsemble::get_effective_sample_size() const {
    double sum = 0.0, sum_sq = 0.0;
    for (int i = 0; i < get_n_points(); ++i) {
        const double w = points_[i * 4 + 3];
        if (!(w > 0.0) || !std::isfinite(w)) continue;
        sum += w;
        sum_sq += w * w;
    }
    if (sum_sq <= 0.0) return 0.0;
    return sum * sum / sum_sq;
}

FRETPairGeometry RotamerEnsemble::pair_geometry(const States& other,
                                                bool use_dipoles) const {
    std::vector<double> other_points, other_mu;
    states_arrays(other, other_points, other_mu);
    if (!use_dipoles) other_mu.clear();

    const int n1 = get_n_points();
    const std::size_t n2 = other_points.size() / 4;
    std::vector<double> c1(static_cast<std::size_t>(n1) * 3), w1(n1);
    for (int i = 0; i < n1; ++i) {
        c1[i * 3] = points_[i * 4];
        c1[i * 3 + 1] = points_[i * 4 + 1];
        c1[i * 3 + 2] = points_[i * 4 + 2];
        w1[i] = points_[i * 4 + 3];
    }
    std::vector<double> c2(n2 * 3), w2(n2);
    for (std::size_t i = 0; i < n2; ++i) {
        c2[i * 3] = other_points[i * 4];
        c2[i * 3 + 1] = other_points[i * 4 + 1];
        c2[i * 3 + 2] = other_points[i * 4 + 2];
        w2[i] = other_points[i * 4 + 3];
    }
    return fret_pair_geometry(
            c1, w1, c2, w2,
            use_dipoles ? orientations_ : std::vector<double>(), other_mu);
}

FRETPairEfficiencies RotamerEnsemble::pair_distribution(
        const States& other, double forster_radius, double tau0) const {
    return fret_pair_efficiencies(pair_geometry(other), forster_radius, tau0);
}

FRETPairEfficiencies RotamerEnsemble::pair_distribution_from_probes(
        const States& other, const std::string& donor,
        const std::string& acceptor, double tau0) const {
    const FRETPairGeometry geometry = pair_geometry(other);
    // Angstrom, like every other length here; the spectra are nanometres and
    // `forster_radius_from_spectra` converts once, inside.
    const double r0 = forster_radius_from_spectra(donor, acceptor,
                                                  geometry.kappa2_avg);
    return fret_pair_efficiencies(geometry, r0, tau0);
}

// --------------------------------------------------------------------------
// Placing a library at a site
// --------------------------------------------------------------------------

RotamerEnsemble RotamerEnsemble::from_frame(
        const ProteinFrame& frame, const std::string& chain, int residue,
        const RotamerLibrary& library, const RotamerSiteOptions& options,
        const std::string& position_name) {
    internal::OwnedView backbone;
    resolve_backbone_site(frame.coords, frame.atom_names, frame.chain_ids,
                          frame.residue_indices, chain, residue,
                          &backbone.data, &backbone.size);
    const std::vector<double> ca(backbone.data, backbone.data + 3);
    const std::vector<double> n_atom(backbone.data + 3, backbone.data + 6);
    const std::vector<double> c_atom(backbone.data + 6, backbone.data + 9);

    internal::OwnedView placed;
    transform_library_to_site(library.coords, ca, n_atom, c_atom, &placed.data,
                              &placed.size);
    const std::vector<double> rotamers = placed.vector();

    nlohmann::json meta = nlohmann::json::object();
    if (!library.metadata.empty()) {
        meta = nlohmann::json::parse(library.metadata, NULL, false);
        if (meta.is_discarded()) meta = nlohmann::json::object();
    }

    const RotamerScoreResult score = compute_rotamer_score(
            rotamers, frame.coords, frame.atom_names, frame.resnames,
            library.atom_names, selector_list(meta, "positive"),
            selector_list(meta, "negative"), library.resnames,
            frame.residue_indices, frame.chain_ids, residue, chain,
            library.weights, options.temperature, options.ignore_h,
            options.electrostatic, options.potential, options.sigma_scaling,
            options.epsilon_scaling);

    const int n_atoms = static_cast<int>(library.atom_names.size());
    const int n_rotamers =
            n_atoms > 0 ? static_cast<int>(rotamers.size() / (n_atoms * 3)) : 0;

    // Which atom is the chromophore centre, and which two span the transition
    // dipole, are the registry's (`r` and `mu`). A library outside the
    // registry has neither: its first atom stands for the centre and the
    // vector from its first atom to its second for the dipole -- a direction,
    // not the dye's, and an off-registry library should carry selectors.
    const std::vector<int> centre_idx = selector_atom_indices(
            library.atom_names, selector_list(meta, "r"), library.resnames);
    const int centre = centre_idx.empty() ? 0 : centre_idx[0];
    const std::vector<int> mu_idx = selector_atom_indices(
            library.atom_names, selector_list(meta, "mu"), library.resnames);
    const int mu_from = mu_idx.size() >= 2 ? mu_idx[0] : 0;
    const int mu_to = mu_idx.size() >= 2 ? mu_idx[1] : 1;

    std::vector<double> points(static_cast<std::size_t>(n_rotamers) * 4, 0.0);
    std::vector<double> mu(static_cast<std::size_t>(n_rotamers) * 3, 0.0);
    for (int k = 0; k < n_rotamers; ++k) {
        const std::size_t base = static_cast<std::size_t>(k) * n_atoms * 3;
        for (int d = 0; d < 3; ++d) {
            points[k * 4 + d] = rotamers[base + centre * 3 + d];
        }
        points[k * 4 + 3] =
                k < static_cast<int>(score.weights.size()) ? score.weights[k]
                                                           : 0.0;
        if (n_atoms > mu_to && n_atoms > mu_from) {
            double norm = 0.0;
            for (int d = 0; d < 3; ++d) {
                const double v = rotamers[base + mu_to * 3 + d] -
                                 rotamers[base + mu_from * 3 + d];
                mu[k * 3 + d] = v;
                norm += v * v;
            }
            norm = std::sqrt(norm);
            // Two atoms at one position have no direction between them. The
            // Python divided anyway and put three NaNs in the array, which
            // reach kappa2 and poison every pair the rotamer takes part in.
            if (norm > 0.0) {
                for (int d = 0; d < 3; ++d) mu[k * 3 + d] /= norm;
            }
        }
    }

    std::string library_name;
    if (meta.is_object() && meta.contains("library_name") &&
        meta["library_name"].is_string()) {
        library_name = meta["library_name"].get<std::string>();
    } else if (meta.is_object() && meta.contains("name") &&
               meta["name"].is_string()) {
        library_name = meta["name"].get<std::string>();
    } else {
        library_name = library.path;
    }

    std::map<std::string, std::string> params;
    params["simulation_type"] = SIMULATION_TYPE_R1;
    params["library"] = library_name;
    params["chain"] = chain;
    params["residue"] = number_text(residue);
    params["temperature"] = number_text(options.temperature);
    params["electrostatic"] = bool_text(options.electrostatic);
    params["potential"] = options.potential;
    params["ignore_h"] = bool_text(options.ignore_h);
    params["sigma_scaling"] = number_text(options.sigma_scaling);
    params["epsilon_scaling"] = number_text(options.epsilon_scaling);
    params["partition"] = number_text(score.partition);

    std::string name = position_name;
    if (name.empty()) name = chain + number_text(residue);

    return RotamerEnsemble(points, ca, mu, name, params, rotamers,
                           library.atom_names, library.resnames,
                           score.energies, score.partition, library_name,
                           chain, residue);
}

RotamerEnsemble RotamerEnsemble::from_site(const std::string& structure,
                                           const std::string& chain,
                                           int residue,
                                           const std::string& library,
                                           const RotamerSiteOptions& options,
                                           const std::string& position_name,
                                           int frame_index) {
    const std::vector<ProteinFrame> frames =
            load_protein_frames(structure, frame_index + 1);
    if (frame_index >= static_cast<int>(frames.size())) {
        IMP_THROW(structure << " has " << frames.size() << " frame(s), frame "
                            << frame_index << " requested",
                  IMP::ValueException);
    }
    return from_frame(frames[frame_index], chain, residue,
                      load_rotamer_library(library), options, position_name);
}

// --------------------------------------------------------------------------
// Weighted averaging over an ensemble
// --------------------------------------------------------------------------

void frame_weights_from_partitions(double* z_values, int n_frames, int n_pair,
                           double** out_view, int* n_out_view) {
    if (n_pair != 2) {
        throw std::invalid_argument("Z must have shape (n_frames, 2)");
    }
    double* out = internal::new_double_view(n_frames < 0 ? 0 : n_frames,
                                           out_view, n_out_view);
    if (out == nullptr) return;

    double total = 0.0;
    for (int i = 0; i < n_frames; ++i) {
        out[i] = z_values[2 * i] * z_values[2 * i + 1];
        total += out[i];
    }
    if (total == 0.0) {
        // uniform, not undefined: a frame where neither dye has an accessible
        // conformer says nothing about the others
        const double u = n_frames > 0 ? 1.0 / (double) n_frames : 0.0;
        for (int i = 0; i < n_frames; ++i) out[i] = u;
        return;
    }
    for (int i = 0; i < n_frames; ++i) out[i] /= total;
}

std::vector<double> weighted_average_sd_se(double* values, int n_values,
                                           double* weights, int n_weights) {
    if (n_values != n_weights) {
        throw std::invalid_argument("values and weights must be the same length");
    }
    std::vector<double> v, w;
    double total = 0.0;
    for (int i = 0; i < n_values; ++i) {
        // a frame where the dye could not be placed contributes nothing rather
        // than poisoning the mean
        if (!std::isfinite(values[i])) continue;
        v.push_back(values[i]);
        w.push_back(weights[i]);
        total += weights[i];
    }
    std::vector<double> out(3, std::nan(""));
    if (v.empty()) return out;

    double mean = 0.0;
    for (size_t i = 0; i < v.size(); ++i) mean += v[i] * (w[i] / total);
    double variance = 0.0;
    for (size_t i = 0; i < v.size(); ++i) {
        const double d = v[i] - mean;
        variance += d * d * (w[i] / total);
    }
    out[0] = mean;
    out[1] = std::sqrt(variance);
    // by the surviving frame count, not the effective count
    out[2] = std::sqrt(variance / (double) v.size());
    return out;
}

double effective_frame_fraction(double* weights, int n_weights) {
    std::vector<double> w;
    for (int i = 0; i < n_weights; ++i)
        if (weights[i] != 0.0) w.push_back(weights[i]);
    if (w.empty()) return 0.0;

    const double uniform = 1.0 / (double) w.size();
    double entropy = 0.0;
    for (size_t i = 0; i < w.size(); ++i)
        entropy -= w[i] * std::log(w[i] / uniform);
    return std::exp(entropy);
}

IMPBFF_END_NAMESPACE

// -------- from RotamerFps.cpp --------
/**
 * (formerly RotamerFps.cpp, now a section of this file)
 * \brief fps.json and rotamer ensembles.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/FRETPair.h>



IMPBFF_BEGIN_NAMESPACE

namespace {

//! The first of \p keys the object has as a non-empty value, or "".
/*! The fps.json aliases: a caller asks for every spelling at once and takes
    whichever the writer used. A key present but empty counts as absent, which
    is what an fps file written with blank fields means. */
std::string first_text(const nlohmann::json& payload,
                       const std::vector<std::string>& keys) {
    if (!payload.is_object()) return std::string();
    for (std::size_t i = 0; i < keys.size(); ++i) {
        if (!payload.contains(keys[i]) || payload[keys[i]].is_null()) continue;
        const nlohmann::json& v = payload[keys[i]];
        if (v.is_string()) {
            const std::string s = v.get<std::string>();
            if (!s.empty()) return s;
        } else if (v.is_number()) {
            std::ostringstream out;
            out << v.get<double>();
            return out.str();
        } else if (v.is_boolean()) {
            return v.get<bool>() ? "True" : "False";
        }
    }
    return std::string();
}

std::vector<std::string> keys_of(const char* a, const char* b = NULL,
                                 const char* c = NULL, const char* d = NULL) {
    std::vector<std::string> out;
    out.push_back(a);
    if (b) out.push_back(b);
    if (c) out.push_back(c);
    if (d) out.push_back(d);
    return out;
}

nlohmann::json parse_object(const std::string& text) {
    nlohmann::json out = nlohmann::json::parse(text, NULL, false);
    if (out.is_discarded() || !out.is_object()) return nlohmann::json::object();
    return out;
}

//! One position's entry, or an empty object when the document has no such key.
nlohmann::json position_entry(const nlohmann::json& positions,
                              const std::string& name) {
    if (name.empty() || !positions.contains(name)) {
        return nlohmann::json::object();
    }
    return positions[name];
}

//! A `params` value as a number, or NaN when it is missing or not one.
double param_number(const std::map<std::string, std::string>& params,
                    const std::string& key) {
    const std::map<std::string, std::string>::const_iterator it =
            params.find(key);
    if (it == params.end() || it->second.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    char* end = NULL;
    const double v = std::strtod(it->second.c_str(), &end);
    if (end == it->second.c_str() || *end != '\0') {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return v;
}

std::string param_text(const std::map<std::string, std::string>& params,
                       const std::string& key) {
    const std::map<std::string, std::string>::const_iterator it =
            params.find(key);
    return it == params.end() ? std::string() : it->second;
}

//! A `params` flag: 1 true, 0 false, -1 when it is not recorded.
int param_flag(const std::map<std::string, std::string>& params,
               const std::string& key) {
    const std::string text = param_text(params, key);
    if (text.empty()) return -1;
    return (text == "True" || text == "true" || text == "1") ? 1 : 0;
}

}  // namespace

RotamerPosition rotamer_position_from_payload(const std::string& name,
                                              const std::string& payload_json) {
    const nlohmann::json payload = parse_object(payload_json);
    RotamerPosition out;
    out.name = name;
    out.chain = first_text(payload, keys_of("chain_identifier", "chain",
                                            "segid"));
    const std::string residue = first_text(
            payload, keys_of("residue_seq_number", "residue", "resid"));
    out.residue = std::atoi(residue.c_str());
    const std::string atom = first_text(payload, keys_of("atom_name"));
    if (!atom.empty()) out.atom_name = atom;
    out.dye = first_text(payload, keys_of("dye_name", "dye", "chromophore"));
    out.library = first_text(payload, keys_of("rotamer_library",
                                              "library_name", "libname",
                                              "library"));
    out.role = first_text(payload, keys_of("role", "label_type"));
    return out;
}

RotamerDistance rotamer_distance_from_payload(
        const std::string& name, const std::string& payload_json,
        const std::string& positions_json) {
    const nlohmann::json payload = parse_object(payload_json);
    const nlohmann::json positions = parse_object(positions_json);
    RotamerDistance out;
    out.name = name;
    out.donor_position = first_text(
            payload, keys_of("position1_name", "donor_position",
                             "donor_position_name", "dye1_position"));
    out.acceptor_position = first_text(
            payload, keys_of("position2_name", "acceptor_position",
                             "acceptor_position_name", "dye2_position"));
    const nlohmann::json donor_entry =
            position_entry(positions, out.donor_position);
    const nlohmann::json acceptor_entry =
            position_entry(positions, out.acceptor_position);

    out.donor = first_text(payload, keys_of("donor", "donor_dye", "dye1",
                                            "dye_1"));
    if (out.donor.empty()) out.donor = first_text(donor_entry, keys_of("dye"));
    out.acceptor = first_text(payload, keys_of("acceptor", "acceptor_dye",
                                               "dye2", "dye_2"));
    if (out.acceptor.empty()) {
        out.acceptor = first_text(acceptor_entry, keys_of("dye"));
    }

    out.libname_1 = first_text(payload, keys_of("libname_1", "donor_library",
                                                "donor_rotamer_library",
                                                "library_1"));
    if (out.libname_1.empty()) {
        out.libname_1 = first_text(donor_entry,
                                   keys_of("rotamer_library", "library"));
    }
    out.libname_2 = first_text(payload, keys_of("libname_2",
                                                "acceptor_library",
                                                "acceptor_rotamer_library",
                                                "library_2"));
    if (out.libname_2.empty()) {
        out.libname_2 = first_text(acceptor_entry,
                                   keys_of("rotamer_library", "library"));
    }
    return out;
}

RotamerFpsSelection read_rotamer_fps(const std::string& path,
                                     const std::string& distance_name) {
    const FPSDocument document = read_fps_json(path);
    const nlohmann::json positions = parse_object(document.positions);
    const nlohmann::json distances = parse_object(document.distances);
    if (distances.empty()) {
        IMP_THROW("No distances found in " << path, ValueException);
    }
    std::string name = distance_name;
    if (name.empty()) name = distances.begin().key();
    if (!distances.contains(name)) {
        IMP_THROW("Distance '" << name << "' not found in " << path,
                  ValueException);
    }

    RotamerFpsSelection out;
    out.distance = rotamer_distance_from_payload(name, distances[name].dump(),
                                                 document.positions);
    if (!positions.contains(out.distance.donor_position)) {
        IMP_THROW("Donor position '" << out.distance.donor_position
                                     << "' not found in " << path,
                  ValueException);
    }
    if (!positions.contains(out.distance.acceptor_position)) {
        IMP_THROW("Acceptor position '" << out.distance.acceptor_position
                                        << "' not found in " << path,
                  ValueException);
    }
    out.donor = rotamer_position_from_payload(
            out.distance.donor_position,
            positions[out.distance.donor_position].dump());
    out.acceptor = rotamer_position_from_payload(
            out.distance.acceptor_position,
            positions[out.distance.acceptor_position].dump());
    out.positions = document.positions;
    out.distances = document.distances;
    // The two free-form sections, merged: what a caller wants of them is
    // "whatever else the file said", and which section a key sat in is the
    // writer's business.
    nlohmann::json extra = parse_object(document.score_sets);
    const nlohmann::json rest = parse_object(document.extra);
    for (nlohmann::json::const_iterator it = rest.begin(); it != rest.end();
         ++it) {
        extra[it.key()] = it.value();
    }
    out.extra = extra.dump();
    return out;
}

std::string rotamer_position_payload(const std::string& chain, int residue,
                                     const std::string& library,
                                     const std::string& atom_name,
                                     const std::string& dye,
                                     double temperature, int electrostatic,
                                     const std::string& potential) {
    nlohmann::json payload;
    payload["chain_identifier"] = chain;
    payload["residue_seq_number"] = residue;
    payload["atom_name"] = atom_name;
    payload["simulation_type"] = SIMULATION_TYPE_R1;
    payload["rotamer_library"] = library;
    if (!dye.empty()) payload["dye_name"] = dye;
    if (!std::isnan(temperature)) payload["temperature"] = temperature;
    if (electrostatic >= 0) payload["electrostatic"] = electrostatic != 0;
    if (!potential.empty()) payload["potential"] = potential;
    return payload.dump();
}

std::string rotamer_positions_payload(
        const std::map<std::string, RotamerEnsemble>& ensembles,
        const std::string& atom_name) {
    nlohmann::json out = nlohmann::json::object();
    for (std::map<std::string, RotamerEnsemble>::const_iterator it =
                 ensembles.begin();
         it != ensembles.end(); ++it) {
        const RotamerEnsemble& e = it->second;
        const std::map<std::string, std::string> params = e.get_params();
        out[it->first] = nlohmann::json::parse(rotamer_position_payload(
                e.get_chain(), e.get_residue(), e.get_library(), atom_name,
                param_text(params, "dye"), param_number(params, "temperature"),
                param_flag(params, "electrostatic"),
                param_text(params, "potential")));
    }
    return out.dump();
}

std::string distances_from_ensembles(
        const std::map<std::string, RotamerEnsemble>& ensembles,
        const std::vector<std::pair<std::string, std::string> >& pairs,
        double forster_radius, const std::string& distance_type, double error,
        double error_fraction, const std::string& kappa2) {
    if (kappa2 != "isotropic" && kappa2 != "dipoles") {
        IMP_THROW("kappa2 must be 'isotropic' or 'dipoles', not '" << kappa2
                                                                   << "'",
                  ValueException);
    }
    if (distance_type != "RDAMean" && distance_type != "RDAMeanE" &&
        distance_type != "Rmp") {
        IMP_THROW("unknown distance_type '" << distance_type << "'",
                  ValueException);
    }

    nlohmann::json out = nlohmann::json::object();
    for (std::size_t p = 0; p < pairs.size(); ++p) {
        const std::string& name1 = pairs[p].first;
        const std::string& name2 = pairs[p].second;
        std::map<std::string, RotamerEnsemble>::const_iterator i1 =
                ensembles.find(name1);
        std::map<std::string, RotamerEnsemble>::const_iterator i2 =
                ensembles.find(name2);
        if (i1 == ensembles.end() || i2 == ensembles.end()) {
            IMP_THROW("no ensemble for the pair (" << name1 << ", " << name2
                                                   << ")",
                      ValueException);
        }
        const RotamerEnsemble& e1 = i1->second;
        const RotamerEnsemble& e2 = i2->second;

        // The isotropic geometry drops the dipoles, which is what puts
        // kappa2 = 2/3 on every pair.
        const FRETPairGeometry geometry =
                e1.pair_geometry(e2, kappa2 == "dipoles");
        double value = 0.0;
        if (distance_type == "Rmp") {
            value = e1.dRmp(e2);
        } else if (distance_type == "RDAMean") {
            for (std::size_t k = 0; k < geometry.R.size(); ++k) {
                value += geometry.R[k] * geometry.weight[k];
            }
        } else {
            const FRETPairEfficiencies eff =
                    fret_pair_efficiencies(geometry, forster_radius);
            const double mean_e = eff.static_efficiency;
            if (mean_e <= 0.0) {
                for (std::size_t k = 0; k < geometry.R.size(); ++k) {
                    value += geometry.R[k] * geometry.weight[k];
                }
            } else if (mean_e >= 1.0) {
                value = 0.0;
            } else {
                value = forster_radius *
                        std::pow(1.0 / mean_e - 1.0, 1.0 / 6.0);
            }
        }
        const double err = error >= 0.0 ? error : error_fraction * value;
        nlohmann::json entry;
        entry["position1_name"] = name1;
        entry["position2_name"] = name2;
        entry["distance_type"] = distance_type;
        entry["distance"] = value;
        entry["error_neg"] = err;
        entry["error_pos"] = err;
        entry["Forster_radius"] = forster_radius;
        out[name1 + "_" + name2] = entry;
    }
    return out.dump();
}

std::map<std::string, RotamerEnsemble> rotamer_ensembles_from_fps(
        const std::string& fps_json, const std::string& structure,
        const std::map<std::string, std::string>& library_map,
        const RotamerSiteOptions& options, int frame_index) {
    const FPSDocument document = read_fps_json(fps_json);
    const nlohmann::json positions = parse_object(document.positions);
    std::map<std::string, RotamerEnsemble> out;
    if (positions.empty()) return out;

    const std::vector<ProteinFrame> frames =
            load_protein_frames(structure, frame_index + 1);
    if (frame_index >= static_cast<int>(frames.size())) {
        IMP_THROW(structure << " has " << frames.size() << " frame(s), frame "
                            << frame_index << " requested",
                  ValueException);
    }
    const ProteinFrame& frame = frames[frame_index];

    // One read per library and not one per position: a donor and an acceptor
    // labelled with the same dye are the common case.
    std::map<std::string, RotamerLibrary> libraries;
    for (nlohmann::json::const_iterator it = positions.begin();
         it != positions.end(); ++it) {
        const RotamerPosition position =
                rotamer_position_from_payload(it.key(), it.value().dump());
        std::string library_name = position.library;
        const std::map<std::string, std::string>::const_iterator mapped =
                library_map.find(it.key());
        if (mapped != library_map.end()) library_name = mapped->second;
        // A position that names no library is skipped and not guessed at.
        if (library_name.empty()) continue;
        if (libraries.find(library_name) == libraries.end()) {
            libraries[library_name] = load_rotamer_library(library_name);
        }
        out.insert(std::make_pair(
                it.key(),
                RotamerEnsemble::from_frame(frame, position.chain,
                                            position.residue,
                                            libraries[library_name], options,
                                            it.key())));
    }
    return out;
}

void write_rotamer_fps(const std::string& path,
                       const std::string& positions_json,
                       const std::string& distances_json,
                       const std::string& merge_into, bool validate) {
    nlohmann::json positions = nlohmann::json::object();
    nlohmann::json distances = nlohmann::json::object();
    nlohmann::json score_sets = nlohmann::json::object();
    nlohmann::json extra = nlohmann::json::object();
    if (!merge_into.empty()) {
        const FPSDocument document = read_fps_json(merge_into);
        positions = parse_object(document.positions);
        distances = parse_object(document.distances);
        score_sets = parse_object(document.score_sets);
        extra = parse_object(document.extra);
    }
    const nlohmann::json new_positions = parse_object(positions_json);
    for (nlohmann::json::const_iterator it = new_positions.begin();
         it != new_positions.end(); ++it) {
        positions[it.key()] = it.value();
    }
    const nlohmann::json new_distances = parse_object(distances_json);
    for (nlohmann::json::const_iterator it = new_distances.begin();
         it != new_distances.end(); ++it) {
        distances[it.key()] = it.value();
    }
    write_fps_json(path, positions.dump(), distances.dump(), score_sets.dump(),
                   extra.dump(), validate);
}

IMPBFF_END_NAMESPACE

// -------- from RotamerFret.cpp --------
/**
 * (formerly RotamerFret.cpp, now a section of this file)
 * \brief FRET over a trajectory from two screened rotamer libraries.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */



#include <iomanip>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! `np.savetxt`'s default: one row per line, `%.18e` per column.
void save_txt(const std::string& path, const std::vector<double>& values,
              int columns, const std::string& header = "") {
    std::ofstream out(path.c_str());
    if (!out) IMP_THROW("cannot write " << path, IOException);
    if (!header.empty()) {
        std::istringstream lines(header);
        std::string line;
        while (std::getline(lines, line)) out << "# " << line << "\n";
    }
    out << std::scientific << std::setprecision(18);
    for (std::size_t i = 0; i < values.size(); ++i) {
        out << values[i];
        const bool row_end = (columns <= 1) ||
                             ((i + 1) % static_cast<std::size_t>(columns) == 0);
        out << (row_end ? "\n" : " ");
    }
}

//! Every number of a whitespace-separated text file, in file order.
std::vector<double> load_txt(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) IMP_THROW("cannot read " << path, IOException);
    std::vector<double> out;
    std::string line;
    while (std::getline(in, line)) {
        const std::string trimmed = internal::trimmed(line);
        if (trimmed.empty() || trimmed[0] == '#') continue;
        std::istringstream fields(trimmed);
        double v;
        while (fields >> v) out.push_back(v);
    }
    return out;
}

//! The finite entries of \p values and their weights, both renormalised.
void finite_pairs(const std::vector<double>& values,
                  const std::vector<double>& weights,
                  std::vector<double>& out_values,
                  std::vector<double>& out_weights) {
    out_values.clear();
    out_weights.clear();
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (!std::isfinite(values[i])) continue;
        out_values.push_back(values[i]);
        out_weights.push_back(i < weights.size() ? weights[i] : 0.0);
    }
}

}  // namespace

void RotamerFRET::load_libraries() {
    if (residues_.size() != 2) {
        IMP_THROW("the residue list must contain exactly 2 residue numbers, "
                  "not " << residues_.size(),
                  ValueException);
    }
    while (chains_.size() < 2) chains_.push_back(std::string());
    lib_1_ = load_rotamer_library(libname_1_);
    lib_2_ = load_rotamer_library(libname_2_);
}

RotamerFRET::RotamerFRET(
        const std::string& protein, const std::vector<int>& residues,
        const std::vector<std::string>& chains, const std::string& donor,
        const std::string& acceptor, const std::string& libname_1,
        const std::string& libname_2, double temperature, bool electrostatic,
        const std::string& potential, bool ign_H, double sigma_scaling,
        double epsilon_scaling, bool fixed_R0, double r0,
        const std::string& r0lib, double z_cutoff,
        const std::string& output_prefix, int max_frames,
        const std::vector<double>& user_weights)
    : frames_(load_protein_frames(protein, max_frames)), residues_(residues),
      chains_(chains), donor_(donor), acceptor_(acceptor),
      libname_1_(libname_1), libname_2_(libname_2), r0lib_(r0lib),
      output_prefix_(output_prefix),
      site_(temperature, electrostatic, potential, ign_H, sigma_scaling,
            epsilon_scaling),
      fixed_r0_(fixed_R0), r0_(r0), z_cutoff_(z_cutoff),
      user_weights_(user_weights) {
    load_libraries();
}

RotamerFRET::RotamerFRET(const std::vector<ProteinFrame>& frames,
                         const std::vector<int>& residues,
                         const std::vector<std::string>& chains,
                         const std::string& donor,
                         const std::string& acceptor,
                         const std::string& libname_1,
                         const std::string& libname_2,
                         const RotamerSiteOptions& site, bool fixed_R0,
                         double r0, const std::string& r0lib, double z_cutoff,
                         const std::string& output_prefix,
                         const std::vector<double>& user_weights)
    : frames_(frames), residues_(residues), chains_(chains), donor_(donor),
      acceptor_(acceptor), libname_1_(libname_1), libname_2_(libname_2),
      r0lib_(r0lib), output_prefix_(output_prefix), site_(site),
      fixed_r0_(fixed_R0), r0_(r0), z_cutoff_(z_cutoff),
      user_weights_(user_weights) {
    load_libraries();
}

RotamerFRET RotamerFRET::from_frames(
        const std::vector<ProteinFrame>& frames,
        const std::vector<int>& residues,
        const std::vector<std::string>& chains, const std::string& donor,
        const std::string& acceptor, const std::string& libname_1,
        const std::string& libname_2, double temperature, bool electrostatic,
        const std::string& potential, bool ign_H, double sigma_scaling,
        double epsilon_scaling, bool fixed_R0, double r0,
        const std::string& r0lib, double z_cutoff,
        const std::string& output_prefix,
        const std::vector<double>& user_weights) {
    return RotamerFRET(frames, residues, chains, donor, acceptor, libname_1,
                       libname_2,
                       RotamerSiteOptions(temperature, electrostatic,
                                          potential, ign_H, sigma_scaling,
                                          epsilon_scaling),
                       fixed_R0, r0, r0lib, z_cutoff, output_prefix,
                       user_weights);
}

FRETFrameResult RotamerFRET::frame_fret(const ProteinFrame& frame) {
    const RotamerEnsemble donor = RotamerEnsemble::from_frame(
            frame, chains_[0], residues_[0], lib_1_, site_);
    const RotamerEnsemble acceptor = RotamerEnsemble::from_frame(
            frame, chains_[1], residues_[1], lib_2_, site_);

    FRETFrameResult out;
    out.z_donor = donor.get_partition();
    out.z_acceptor = acceptor.get_partition();

    const FRETPairGeometry geometry = donor.pair_geometry(acceptor);
    out.kappa2_avg = geometry.kappa2_avg;
    if (!fixed_r0_) {
        // FRETpredict's convention: R0 is recomputed per frame at that
        // frame's <kappa2> and then applied with the isotropic formula.
        r0_ = forster_radius_from_spectra(donor_, acceptor_,
                                          geometry.kappa2_avg, r0lib_);
        if (r0_ == 0.0) {
            const double nan = std::numeric_limits<double>::quiet_NaN();
            out.kappa2_avg = nan;
            out.static_efficiency = nan;
            out.dynamic1 = nan;
            out.dynamic2 = nan;
            return out;
        }
    }
    const FRETPairEfficiencies eff = fret_pair_efficiencies(geometry, r0_);
    out.static_efficiency = eff.static_efficiency;
    out.dynamic1 = eff.dynamic1;
    out.dynamic2 = eff.dynamic2;
    return out;
}

void RotamerFRET::trajectory_analysis() {
    const std::size_t n = frames_.size();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    z_values_.assign(n * 2, 0.0);
    k2_values_.assign(n, nan);
    estatic_.assign(n, nan);
    edynamic1_.assign(n, nan);
    edynamic2_.assign(n, nan);

    for (std::size_t i = 0; i < n; ++i) {
        const FRETFrameResult result = frame_fret(frames_[i]);
        z_values_[i * 2] = result.z_donor;
        z_values_[i * 2 + 1] = result.z_acceptor;
        if (result.z_donor <= z_cutoff_ || result.z_acceptor <= z_cutoff_) {
            continue;
        }
        k2_values_[i] = result.kappa2_avg;
        estatic_[i] = result.static_efficiency;
        edynamic1_[i] = result.dynamic1;
        edynamic2_[i] = result.dynamic2;
    }
}

void RotamerFRET::write_summary(const std::string& prefix,
                                const std::vector<double>& k2,
                                const std::vector<double>& es,
                                const std::vector<double>& ed1,
                                const std::vector<double>& ed2,
                                const std::vector<double>& weights) const {
    std::vector<double> rows;
    const std::vector<double>* quantities[4] = {&k2, &es, &ed1, &ed2};
    if (k2.size() == 1) {
        // One frame has no spread: the average is the value and the SD and SE
        // are undefined rather than zero.
        const double nan = std::numeric_limits<double>::quiet_NaN();
        for (int q = 0; q < 4; ++q) {
            rows.push_back((*quantities[q])[0]);
            rows.push_back(nan);
            rows.push_back(nan);
        }
    } else {
        for (int q = 0; q < 4; ++q) {
            std::vector<double> values, w;
            finite_pairs(*quantities[q], weights, values, w);
            const std::vector<double> stats = weighted_average_sd_se(
                    values.empty() ? NULL : &values[0],
                    static_cast<int>(values.size()),
                    w.empty() ? NULL : &w[0], static_cast<int>(w.size()));
            rows.insert(rows.end(), stats.begin(), stats.end());
        }
    }
    std::ostringstream path;
    path << prefix << "-data-" << residues_[0] << "-" << residues_[1] << ".dat";
    save_txt(path.str(), rows, 3,
             "quantity Average SD SE\nk2 Estatic Edynamic1 Edynamic2");
}

void RotamerFRET::save(const std::string& output_prefix) {
    if (z_values_.empty()) trajectory_analysis();
    const std::string prefix =
            output_prefix.empty() ? output_prefix_ : output_prefix;
    std::ostringstream tail;
    tail << "-" << residues_[0] << "-" << residues_[1] << ".dat";

    internal::OwnedView ws;
    frame_weights_from_partitions(z_values_.empty() ? NULL : &z_values_[0],
                                  static_cast<int>(z_values_.size() / 2), 2,
                                  &ws.data, &ws.size);

    save_txt(prefix + "-Z" + tail.str(), z_values_, 2);
    save_txt(prefix + "-w_s" + tail.str(), ws.vector(), 1);
    save_txt(prefix + "-k2" + tail.str(), k2_values_, 1);
    save_txt(prefix + "-Es" + tail.str(), estatic_, 1);
    save_txt(prefix + "-Ed1" + tail.str(), edynamic1_, 1);
    save_txt(prefix + "-Ed2" + tail.str(), edynamic2_, 1);

    std::vector<double> weights(k2_values_.size(), 1.0);
    if (!user_weights_.empty()) {
        if (user_weights_.size() != k2_values_.size()) {
            IMP_THROW("Weights array has size "
                              << user_weights_.size()
                              << " whereas the number of frames is "
                              << k2_values_.size(),
                      ValueException);
        }
        weights = user_weights_;
    }
    double total = 0.0;
    for (std::size_t i = 0; i < weights.size(); ++i) total += weights[i];
    if (total > 0.0) {
        for (std::size_t i = 0; i < weights.size(); ++i) weights[i] /= total;
    }
    write_summary(prefix, k2_values_, estatic_, edynamic1_, edynamic2_,
                  weights);
}

void RotamerFRET::reweight(bool boltzmann_weights,
                           const std::vector<double>& user_weights,
                           const std::string& output_prefix) {
    const std::string prefix =
            output_prefix.empty() ? output_prefix_ : output_prefix;
    std::ostringstream tail;
    tail << "-" << residues_[0] << "-" << residues_[1] << ".dat";

    if (boltzmann_weights) {
        std::vector<double> z = load_txt(output_prefix_ + "-Z" + tail.str());
        internal::OwnedView ws;
        frame_weights_from_partitions(z.empty() ? NULL : &z[0],
                                      static_cast<int>(z.size() / 2), 2,
                                      &ws.data, &ws.size);
        weights_ = ws.vector();
    }
    if (!user_weights.empty()) user_weights_ = user_weights;

    const std::vector<double> k2 = load_txt(output_prefix_ + "-k2" +
                                            tail.str());
    const std::vector<double> es = load_txt(output_prefix_ + "-Es" +
                                            tail.str());
    const std::vector<double> ed1 = load_txt(output_prefix_ + "-Ed1" +
                                             tail.str());
    const std::vector<double> ed2 = load_txt(output_prefix_ + "-Ed2" +
                                             tail.str());

    std::vector<double> weights(k2.size(), 1.0);
    if (weights_.size() == k2.size()) weights = weights_;
    if (!user_weights_.empty()) {
        if (user_weights_.size() != k2.size()) {
            IMP_THROW("Weights array has size "
                              << user_weights_.size()
                              << " whereas the number of frames is "
                              << k2.size(),
                      ValueException);
        }
        for (std::size_t i = 0; i < weights.size(); ++i) {
            weights[i] *= user_weights_[i];
        }
    }
    double total = 0.0;
    for (std::size_t i = 0; i < weights.size(); ++i) total += weights[i];
    if (total > 0.0) {
        for (std::size_t i = 0; i < weights.size(); ++i) weights[i] /= total;
    }
    write_summary(prefix, k2, es, ed1, ed2, weights);
}

void RotamerFRET::run() {
    trajectory_analysis();
    save();
}

void RotamerFRET::get_z_values(double** out_view, int* n_out_view) const {
    internal::copy_to_view(z_values_, out_view, n_out_view);
}
void RotamerFRET::get_k2_values(double** out_view, int* n_out_view) const {
    internal::copy_to_view(k2_values_, out_view, n_out_view);
}
void RotamerFRET::get_estatic_values(double** out_view, int* n_out_view) const {
    internal::copy_to_view(estatic_, out_view, n_out_view);
}
void RotamerFRET::get_edynamic1_values(double** out_view,
                                       int* n_out_view) const {
    internal::copy_to_view(edynamic1_, out_view, n_out_view);
}
void RotamerFRET::get_edynamic2_values(double** out_view,
                                       int* n_out_view) const {
    internal::copy_to_view(edynamic2_, out_view, n_out_view);
}

RotamerFRET rotamer_fret_from_fps(const std::string& fps_path,
                                  const std::string& protein,
                                  const std::string& distance_name,
                                  const std::string& output_prefix,
                                  bool fixed_R0, double r0, double temperature,
                                  bool electrostatic) {
    const RotamerFpsSelection selection = read_rotamer_fps(fps_path,
                                                           distance_name);
    std::vector<int> residues;
    residues.push_back(selection.donor.residue);
    residues.push_back(selection.acceptor.residue);
    std::vector<std::string> chains;
    chains.push_back(selection.donor.chain);
    chains.push_back(selection.acceptor.chain);

    const std::string donor = selection.distance.donor.empty()
                                      ? selection.donor.dye
                                      : selection.distance.donor;
    const std::string acceptor = selection.distance.acceptor.empty()
                                         ? selection.acceptor.dye
                                         : selection.distance.acceptor;
    const std::string libname_1 = selection.distance.libname_1.empty()
                                          ? selection.donor.library
                                          : selection.distance.libname_1;
    const std::string libname_2 = selection.distance.libname_2.empty()
                                          ? selection.acceptor.library
                                          : selection.distance.libname_2;
    return RotamerFRET(protein, residues, chains, donor, acceptor, libname_1,
                       libname_2, temperature, electrostatic, "lj", true, 0.5,
                       1.0, fixed_R0, r0, "", 0.05, output_prefix);
}

IMPBFF_END_NAMESPACE
