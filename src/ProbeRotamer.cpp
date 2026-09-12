#include <IMP/bff/ProbeRotamer.h>
#include <IMP/bff/ProbeDataPaths.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/Text.h>
#include <IMP/bff/internal/json.h>
#include <IMP/bff/RotamerScoring.h>
#include <IMP/bff/ProbeSampling.h>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>


// -------- from RotamerEnergy.cpp --------
/**
 * (formerly RotamerEnergy.cpp, now a section of this file)
 * \brief The interaction energy of each rotamer with its protein.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */


IMPBFF_BEGIN_NAMESPACE

std::vector<double> probe_rotamer_interaction_energies(
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
                    if (potential == PROBE_ROTAMER_POTENTIAL_GAUSS) {
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
        const std::vector<std::string> listed = probe_rotamer_drot_catalog(path);
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

std::string probe_rotamer_library_registry() { return registry().dump(); }

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

std::string probe_rotamer_library_metadata(const std::string& name) {
    return entry_for(normalize_library_name(name), name).dump();
}

std::string probe_rotamer_library_metadata_for_path(const std::string& path) {
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

std::string resolve_probe_rotamer_library_path(const std::string& name,
                                         const std::string& lib_dir) {
    if (rs_exists(name)) return name;

    const nlohmann::json meta =
            nlohmann::json::parse(probe_rotamer_library_metadata(name));
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

std::vector<std::string> infer_probe_rotamer_resnames(
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

ProbeRotamerLibrary load_probe_rotamer_library(const std::string& name,
                                    const std::string& lib_dir) {
    const std::string locator = resolve_probe_rotamer_library_path(name, lib_dir);
    // A locator addresses one library inside a family container
    // (`dyes.drot.pto::A48_C1R_cutoff10`); a plain path is a container of its
    // own. `read_probe_rotamer_drot` takes either, but everything here that reasons about
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
                                         ? probe_rotamer_library_metadata_for_path(
                                                   container)
                                         : probe_rotamer_library_metadata(name);
    nlohmann::json meta = nlohmann::json::parse(metadata, NULL, false);
    if (meta.is_discarded()) meta = nlohmann::json::object();

    const std::string dir = rs_dirname(container);
    const std::string file_name = rs_basename(container);
    ProbeRotamerLibrary lib;
    std::string sidecar_pdb;

    if (internal::ends_with(file_name, ".drot") ||
        internal::ends_with(file_name, ".drot.pto")) {
        // The shipped form (PRD-118): the template, the Z-matrix and the
        // per-conformer internal coordinates in one container, so nothing
        // beside it is read -- names, residues and elements are all in it.
        lib = read_probe_rotamer_drot(locator);
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
        lib.resnames = infer_probe_rotamer_resnames(lib.atom_names, lib.metadata);
    }
    return lib;
}

IMPBFF_END_NAMESPACE

// -------- from ProbeRotamerEnsemble.cpp --------
/**
 * (formerly ProbeRotamerEnsemble.cpp, now a section of this file)
 * \brief A rotamer library placed and screened at one labelling site.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */




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

}  // namespace

// --------------------------------------------------------------------------
// ProbeRotamerEnsemble
// --------------------------------------------------------------------------

ProbeRotamerEnsemble::ProbeRotamerEnsemble(
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

void ProbeRotamerEnsemble::get_atoms(double** out_view, int* n_out_view) const {
    internal::copy_to_view(atoms_, out_view, n_out_view);
}

void ProbeRotamerEnsemble::get_energies(double** out_view, int* n_out_view) const {
    internal::copy_to_view(energies_, out_view, n_out_view);
}

void ProbeRotamerEnsemble::get_centres(double** out_view, int* n_out_view) const {
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

void ProbeRotamerEnsemble::get_weights(double** out_view, int* n_out_view) const {
    const int n = get_n_points();
    double* buffer = internal::new_double_view(static_cast<std::size_t>(n),
                                               out_view, n_out_view);
    if (buffer == NULL) return;
    for (int i = 0; i < n; ++i) buffer[i] = points_[i * 4 + 3];
}

int ProbeRotamerEnsemble::get_n_atoms() const {
    const int n = get_n_points();
    if (n == 0) return 0;
    return static_cast<int>(atoms_.size() / (static_cast<std::size_t>(n) * 3));
}

double ProbeRotamerEnsemble::get_effective_sample_size() const {
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

// --------------------------------------------------------------------------
// Placing a library at a site
// --------------------------------------------------------------------------

ProbeRotamerEnsemble ProbeRotamerEnsemble::from_frame(
        const ProteinFrame& frame, const std::string& chain, int residue,
        const ProbeRotamerLibrary& library, const ProbeRotamerSiteOptions& options,
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

    const RotamerScoreResult score = get_rotamer_score(
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

    return ProbeRotamerEnsemble(points, ca, mu, name, params, rotamers,
                           library.atom_names, library.resnames,
                           score.energies, score.partition, library_name,
                           chain, residue);
}

ProbeRotamerEnsemble ProbeRotamerEnsemble::from_site(const std::string& structure,
                                           const std::string& chain,
                                           int residue,
                                           const std::string& library,
                                           const ProbeRotamerSiteOptions& options,
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
                      load_probe_rotamer_library(library), options, position_name);
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
