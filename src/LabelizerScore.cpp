/**
 * \file LabelizerScore.cpp
 * \brief The Labelizer label-site score: fitted tables, the seven per-residue
 *        parameters, and the weighted geometric mean.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/LabelizerScore.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <utility>
#include <IMP/bff/ProbeDataPaths.h>
#include <IMP/bff/internal/json.h>

#include <map>

IMPBFF_BEGIN_NAMESPACE

namespace {

std::string ls_table_path(const std::string& name) {
    return get_data_path("labelizer/probabilities/" + name +
                         "_P_l_after_s.json");
}

//! Distance between two atoms of a structure.
double ls_dist(const LabelizerStructure& s, int a, int b) {
    const double dx = s.xyz[3 * a + 0] - s.xyz[3 * b + 0];
    const double dy = s.xyz[3 * a + 1] - s.xyz[3 * b + 1];
    const double dz = s.xyz[3 * a + 2] - s.xyz[3 * b + 2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

//! The clamped-ramp exposure factor the tryptophan and charge terms share.
/*!
    Two linear ramps, multiplied (`tryptophan_proximity.py:106`, duplicated
    verbatim at `charge_environment.py:143`). Note that the `SLOPE_OFFSET` of 5
    in the first denominator means the factor bottoms out at 0.25 rather than
    0 at the top of its range -- that is the reference's arithmetic, not a
    transcription slip.
*/
double ls_hse_factor(double hse_up, double hse_down) {
    const double sum_hse = hse_up + hse_down;
    double sum_fact;
    if (sum_hse < 30.0) sum_fact = 1.0;
    else if (sum_hse > 45.0) sum_fact = 0.0;
    else sum_fact = 1.0 - (sum_hse - 30.0) * (1.0 / (45.0 - 30.0 + 5.0));

    double hse1_fact;
    if (hse_up < 10.0) hse1_fact = 1.0;
    else if (hse_up > 20.0) hse1_fact = 0.0;
    else hse1_fact = 1.0 - (hse_up - 10.0) * (1.0 / 10.0);

    return sum_fact * hse1_fact;
}

//! Formal charge per residue, as the reference counts it.
/*! `charge_environment.py:34`: histidine is counted as fully charged. */
double ls_formal_charge(const std::string& comp_id) {
    if (comp_id == "ARG" || comp_id == "LYS" || comp_id == "HIS") return 1.0;
    if (comp_id == "ASP" || comp_id == "GLU") return -1.0;
    return 0.0;
}

//! Residue average mass, Dalton -- the observable of `N_CR2_Mass`.
double ls_residue_mass(const std::string& c) {
    static const std::map<std::string, double> m = [] {
        std::map<std::string, double> v;
        const char* n[20] = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY",
                             "HIS","ILE","LEU","LYS","MET","PHE","PRO","SER",
                             "THR","TRP","TYR","VAL"};
        const double d[20] = {71.08,156.19,114.10,115.09,103.14,128.13,129.12,
                              57.05,137.14,113.16,113.16,128.17,131.19,147.18,
                              97.12,87.08,101.10,186.21,163.18,99.13};
        for (int i = 0; i < 20; ++i) v[n[i]] = d[i];
        return v;
    }();
    std::map<std::string, double>::const_iterator it = m.find(c);
    return it == m.end() ? -1.0 : it->second;
}

//! Heavy atoms in the side chain -- the observable of `I_CR4_N_Sidechain`.
double ls_sidechain_atoms(const std::string& c) {
    static const std::map<std::string, double> m = [] {
        std::map<std::string, double> v;
        const char* n[20] = {"GLY","ALA","SER","CYS","THR","VAL","PRO","LEU",
                             "ILE","ASN","ASP","MET","GLN","GLU","LYS","HIS",
                             "ARG","PHE","TYR","TRP"};
        const double d[20] = {0,1,2,2,3,3,3,4,4,4,4,4,5,5,5,6,7,7,8,10};
        for (int i = 0; i < 20; ++i) v[n[i]] = d[i];
        return v;
    }();
    std::map<std::string, double>::const_iterator it = m.find(c);
    return it == m.end() ? -1.0 : it->second;
}

//! Charge/polarity class -- the observable of `C_CR3_Charge`.
/*! Its four keys are `+`, `-`, `P` and `NP`. The grouping is the conventional
    one and is stated rather than inherited: the reference ships this table but
    implements no code for it (`cysteine_resemblance.py:79` raises for anything
    but `C_CR1_Name`), so there is nothing to reproduce. */
std::string ls_charge_class(const std::string& c) {
    if (c == "ARG" || c == "LYS" || c == "HIS") return "+";
    if (c == "ASP" || c == "GLU") return "-";
    if (c == "SER" || c == "THR" || c == "ASN" || c == "GLN" ||
        c == "CYS" || c == "TYR") return "P";
    return "NP";
}

}  // namespace

// ---------------------------------------------------------------------------
// Tables
// ---------------------------------------------------------------------------

const LabelizerTable& labelizer_load_table(const std::string& name) {
    static std::map<std::string, LabelizerTable> cache;
    std::map<std::string, LabelizerTable>::iterator it = cache.find(name);
    if (it != cache.end()) return it->second;

    const std::string path = ls_table_path(name);
    std::ifstream in(path.c_str());
    if (!in) {
        IMP_THROW("labelizer_load_table: no fitted table named '" << name << "' at "
                  << path, IOException);
    }
    nlohmann::json j;
    in >> j;

    LabelizerTable t;
    t.name = name;
    // The reference reads the domain off the first character of the name:
    // `C_` categorical, `N_`/`I_` numeric (labeling_parameter.py:198).
    t.categorical = !name.empty() && name[0] == 'C';
    if (t.categorical) {
        for (nlohmann::json::const_iterator k = j.begin(); k != j.end(); ++k) {
            t.by_key[k.key()] = k.value().get<double>();
        }
    } else {
        std::vector<std::pair<double, double> > rows;
        for (nlohmann::json::const_iterator k = j.begin(); k != j.end(); ++k) {
            rows.push_back(std::make_pair(std::atof(k.key().c_str()),
                                          k.value().get<double>()));
        }
        std::sort(rows.begin(), rows.end());
        for (std::size_t i = 0; i < rows.size(); ++i) {
            t.bins.push_back(rows[i].first);
            t.values.push_back(rows[i].second);
        }
    }
    cache[name] = t;
    return cache[name];
}

std::vector<std::string> labelizer_available_tables() {
    // The shipped set, named rather than globbed so a missing file is an
    // error at load rather than a silently shorter list.
    static const char* names[] = {
        "C_CR1_Name", "C_CR3_Charge", "C_CS6_Cys_In_Variety", "C_SS1_SS",
        "C_SS4_SS-1", "C_SS5_SS-2", "C_SS6_SS+1", "C_SS7_SS+2",
        "I_CR4_N_Sidechain", "I_CS1_Color", "I_CS5_Variety_Length",
        "I_SE4_HSE1_10A", "I_SE5_HSE2_10A", "I_SE6_HSE1_13A", "I_SE7_HSE2_13A",
        "I_SE8_HSE1_16A", "I_SE9_HSE2_16A", "N_CR2_Mass", "N_CS2_Score",
        "N_CS3_Lower_Score", "N_CS4_Upper_Score",
        "N_ME11_Methionin_Exclusion_Dummy", "N_SE1_RSA_Wilke",
        "N_SE10_CB_SURFACE_DIST", "N_SE11_MEAN_SURFACE_DIST", "N_SE2_RSA_Sander",
        "N_SE3_RSA_Miller", "N_SS2_Phi", "N_SS3_Psi"};
    return std::vector<std::string>(names, names + 29);
}

double labelizer_lookup(const LabelizerTable& table, double value) {
    if (table.bins.empty()) return 0.0;
    // Nearest bin centre, no interpolation (labeling_parameter.py:184).
    std::size_t best = 0;
    double best_d = std::fabs(table.bins[0] - value);
    for (std::size_t i = 1; i < table.bins.size(); ++i) {
        const double d = std::fabs(table.bins[i] - value);
        if (d < best_d) { best_d = d; best = i; }
    }
    return table.values[best];
}

double labelizer_lookup_key(const LabelizerTable& table, const std::string& key) {
    std::map<std::string, double>::const_iterator it = table.by_key.find(key);
    if (it == table.by_key.end()) {
        IMP_THROW("labelizer_lookup_key: '" << key << "' is not in table "
                  << table.name, ValueException);
    }
    return it->second;
}

char labelizer_one_letter(const std::string& comp_id) {
    static const char* const pairs[20] = {
        "ALA A", "ARG R", "ASN N", "ASP D", "CYS C", "GLN Q", "GLU E",
        "GLY G", "HIS H", "ILE I", "LEU L", "LYS K", "MET M", "PHE F",
        "PRO P", "SER S", "THR T", "TRP W", "TYR Y", "VAL V"};
    if (comp_id.size() != 3) return 'X';
    for (int i = 0; i < 20; ++i) {
        if (comp_id.compare(0, 3, pairs[i], 3) == 0) return pairs[i][4];
    }
    return 'X';
}

// ---------------------------------------------------------------------------
// The model
// ---------------------------------------------------------------------------

std::vector<LabelizerParameter> labelizer_model_paper() {
    std::vector<LabelizerParameter> m;
    m.push_back(LabelizerParameter("cs", "N_CS2_Score", 1));
    m.push_back(LabelizerParameter("se", "N_SE11_MEAN_SURFACE_DIST", 1));
    m.push_back(LabelizerParameter("tp", "", 0));
    m.push_back(LabelizerParameter("cr", "C_CR1_Name", 1));
    m.push_back(LabelizerParameter("ss", "C_SS1_SS", 1));
    m.push_back(LabelizerParameter("ce", "", 0));
    return m;
}

std::string labelizer_score_type(const std::string& tag) {
    if (tag == "cs") return "conservation";
    if (tag == "se") return "solvent_exposure";
    if (tag == "ss") return "secondary_structure";
    if (tag == "ce") return "charge_environment";
    if (tag == "tp") return "tryptophan_proximity";
    if (tag == "cr") return "cysteine_resemblance";
    if (tag == "me") return "methionine_exclusion";
    if (tag == "combined") return "combined";
    IMP_THROW("labelizer_score_type: unknown parameter tag '" << tag << "'",
              ValueException);
}

// ---------------------------------------------------------------------------
// The parameters
// ---------------------------------------------------------------------------

std::vector<LabelizerScore> labelizer_parameter_scores(
        const LabelizerStructure& s, const std::vector<LabelizerParameter>& model,
        const LabelizerOptions& options,
        const std::map<std::string, double>& conservation) {
    const std::size_t nr = s.residues.size();
    std::vector<LabelizerScore> out;
    out.reserve(nr * model.size());
    if (nr == 0) return out;

    // Everything the terms below share, computed once. The reference rebuilds
    // the surface and the neighbour search per parameter object.
    const std::string ss = labelizer_dssp(s);

    double* hse = 0; int n_hse = 0;
    labelizer_half_sphere_exposure(s, options.hse_radius, &hse, &n_hse);

    double* depth = 0; int n_depth = 0;
    labelizer_residue_depth(s, options.probe_radius, options.n_sphere_points, &depth,
                     &n_depth);

    // The exclusion term needs the exposure of every candidate residue, which
    // the reference takes as 1/depth (methionin_exclusion.py:62).
    const std::string excluded =
            options.model == LABELIZER_MODEL_PUBLISHED ? "MET" : options.exclusion_residue;

    // Relative accessibility is a whole-structure computation and there are
    // three scales; compute a scale at most once, and only if a term asks.
    std::map<int, std::vector<double> > rsa_cache;
    struct RsaFor {
        const LabelizerStructure& s;
        const LabelizerOptions& o;
        std::map<int, std::vector<double> >& cache;
        const std::vector<double>& operator()(const std::string& table) const {
            LabelizerMaxAsa scale = LABELIZER_MAXASA_WILKE;
            if (table.find("Sander") != std::string::npos)
                scale = LABELIZER_MAXASA_SANDER;
            else if (table.find("Miller") != std::string::npos)
                scale = LABELIZER_MAXASA_MILLER;
            std::map<int, std::vector<double> >::iterator it =
                    cache.find(static_cast<int>(scale));
            if (it != cache.end()) return it->second;
            double* v = 0; int n = 0;
            labelizer_relative_solvent_accessibility(s, scale, o.probe_radius,
                                              o.n_sphere_points, &v, &n);
            std::vector<double> out(v, v + n);
            std::free(v);
            return cache[static_cast<int>(scale)] = out;
        }
    } rsa_for = {s, options, rsa_cache};

    // Same treatment for the Cbeta depth: one whole-structure surface build,
    // and only if a term asks for it.
    std::vector<double> cb_cache;
    struct CbDepthFor {
        const LabelizerStructure& s;
        const LabelizerOptions& o;
        std::vector<double>& cache;
        const std::vector<double>& operator()() const {
            if (!cache.empty()) return cache;
            double* v = 0; int n = 0;
            labelizer_cbeta_depth(s, o.probe_radius, o.n_sphere_points, &v, &n);
            cache.assign(v, v + n);
            std::free(v);
            return cache;
        }
    } cb_depth_for = {s, options, cb_cache};

    for (std::size_t p = 0; p < model.size(); ++p) {
        const LabelizerParameter& par = model[p];
        const std::string score_type = labelizer_score_type(par.tag);

        for (std::size_t i = 0; i < nr; ++i) {
            const LabelizerResidue& r = s.residues[i];
            LabelizerScore row;
            row.asym_id = r.chain;
            row.seq_id = r.seq_id;
            row.comp_id = r.comp_id;
            row.score_type = score_type;
            row.status = "unavailable";
            row.value = 0.0;

            if (r.ca < 0) {
                row.status = "unresolved";
                out.push_back(row);
                continue;
            }

            if (par.tag == "cr") {
                const LabelizerTable& t = labelizer_load_table(par.table);
                if (par.table == "C_CR1_Name") {
                    const std::string key(1, labelizer_one_letter(r.comp_id));
                    if (t.by_key.count(key)) {
                        row.value = t.by_key.find(key)->second;
                        row.status = "scored";
                    }
                } else if (par.table == "C_CR3_Charge") {
                    const std::string key = ls_charge_class(r.comp_id);
                    if (t.by_key.count(key)) {
                        row.value = t.by_key.find(key)->second;
                        row.status = "scored";
                    }
                } else if (par.table == "N_CR2_Mass") {
                    const double mass = ls_residue_mass(r.comp_id);
                    if (mass > 0.0) {
                        row.value = labelizer_lookup(t, mass);
                        row.status = "scored";
                    }
                } else if (par.table == "I_CR4_N_Sidechain") {
                    const double n = ls_sidechain_atoms(r.comp_id);
                    if (n >= 0.0) {
                        row.value = labelizer_lookup(t, n);
                        row.status = "scored";
                    }
                } else {
                    IMP_THROW("labelizer_parameter_scores: table '" << par.table
                              << "' is not an implemented observable for tag "
                              << "'cr'. Implemented: C_CR1_Name, C_CR3_Charge, "
                              << "N_CR2_Mass, I_CR4_N_Sidechain.",
                              ValueException);
                }
            } else if (par.tag == "ss") {
                const LabelizerTable& t = labelizer_load_table(par.table);
                // SS4/5/6/7 score the structure of a *neighbouring* residue,
                // so the observable is the letter at an offset -- reading the
                // letter at i and calling it SS-1 is silently wrong.
                int offset = 0;
                bool known = true;
                if (par.table == "C_SS1_SS") offset = 0;
                else if (par.table == "C_SS4_SS-1") offset = -1;
                else if (par.table == "C_SS5_SS-2") offset = -2;
                else if (par.table == "C_SS6_SS+1") offset = 1;
                else if (par.table == "C_SS7_SS+2") offset = 2;
                else known = false;

                if (!known) {
                    IMP_THROW("labelizer_parameter_scores: table '" << par.table
                              << "' is not an implemented observable for tag "
                              << "'ss'. Implemented: C_SS1_SS, C_SS4_SS-1, "
                              << "C_SS5_SS-2, C_SS6_SS+1, C_SS7_SS+2. The phi "
                              << "and psi tables (N_SS2_Phi, N_SS3_Psi) are "
                              << "not: their bin centres span -153 to +333 "
                              << "degrees, which is neither the -180..180 nor "
                              << "the 0..360 convention, and the reference "
                              << "implements no code for them to check "
                              << "against (secondary_structure.py:156).",
                              ValueException);
                }
                const int j = static_cast<int>(i) + offset;
                // A neighbour across a chain break or off the end is not a
                // neighbour; the position simply has no value.
                bool ok = j >= 0 && j < static_cast<int>(nr);
                if (ok && offset != 0) {
                    const int lo = std::min(static_cast<int>(i), j);
                    const int hi = std::max(static_cast<int>(i), j);
                    for (int k = lo + 1; k <= hi && ok; ++k) {
                        if (s.residues[k].chain != s.residues[k - 1].chain ||
                            s.residues[k].seq_id != s.residues[k - 1].seq_id + 1)
                            ok = false;
                    }
                }
                if (ok) {
                    const std::string key(1, ss[j]);
                    if (t.by_key.count(key)) {
                        row.value = t.by_key.find(key)->second;
                        row.status = "scored";
                    }
                }
            } else if (par.tag == "cs") {
                if (par.table != "N_CS2_Score") {
                    // The other conservation tables score fields a ConSurf
                    // grade file carries beside the score -- the confidence
                    // bounds, the colour bin, the variety string. The importer
                    // reads the score only, so there is nothing to look them
                    // up with, and guessing would be worse than refusing. The
                    // reference refuses too (conservation_score.py:153).
                    IMP_THROW("labelizer_parameter_scores: table '" << par.table
                              << "' needs ConSurf fields labelizer_read_consurf does "
                              << "not import (confidence bounds, colour bin, "
                              << "variety). Only N_CS2_Score is implemented, "
                              << "as in the reference.", ValueException);
                }
                const std::string key = labelizer_residue_key(r.chain, r.seq_id);
                std::map<std::string, double>::const_iterator g =
                        conservation.find(key);
                if (g != conservation.end()) {
                    row.value = labelizer_lookup(labelizer_load_table(par.table), g->second);
                    row.status = "scored";
                }
            } else if (par.tag == "se") {
                const LabelizerTable& t = labelizer_load_table(par.table);
                double observable = 0.0;
                bool have = false;
                if (par.table.find("HSE") != std::string::npos) {
                    // The `I_SE*` tables read the up-count at 10, 13 or 16 A.
                    observable = hse[2 * i + 0];
                    have = true;
                } else if (par.table.find("RSA") != std::string::npos) {
                    observable = rsa_for(par.table)[i];
                    have = observable >= 0.0;
                } else if (par.table == "N_SE10_CB_SURFACE_DIST") {
                    // The Cbeta's own distance to the surface, not the
                    // residue's mean -- a different observable on the same
                    // surface, and its bins start at 0.69 where the mean's
                    // start at 1.46.
                    observable = cb_depth_for()[i];
                    have = observable >= 0.0;
                } else if (par.table == "N_SE11_MEAN_SURFACE_DIST") {
                    // The published default, clipped at 4.0
                    // (solvent_exposure.py:157).
                    observable = std::min(depth[i], 4.0);
                    have = true;
                } else {
                    IMP_THROW("labelizer_parameter_scores: table '" << par.table
                              << "' is not an implemented observable for tag "
                              << "'se'. Implemented: N_SE11_MEAN_SURFACE_DIST, "
                              << "N_SE10_CB_SURFACE_DIST, N_SE1/2/3_RSA_*, "
                              << "I_SE4..I_SE9_HSE*.", ValueException);
                }
                if (have) {
                    row.value = labelizer_lookup(t, observable);
                    row.status = "scored";
                }
            } else if (par.tag == "me") {
                if (!par.table.empty() &&
                    par.table != "N_ME11_Methionin_Exclusion_Dummy") {
                    IMP_THROW("labelizer_parameter_scores: the exclusion term is "
                              << "hard-coded (0.001/0.999) and reads no table; "
                              << "'" << par.table << "' would be ignored.",
                              ValueException);
                }
                // Within the contact distance of an exposed excluded residue,
                // measured C-alpha to C-alpha (methionin_exclusion.py:77).
                // 0.001/0.999 rather than 0/1: zero would trip the zero-veto
                // and -1 means an error, so the reference uses a thousandfold
                // penalty that keeps the geometric mean defined.
                bool too_close = false;
                for (std::size_t j = 0; j < nr && !too_close; ++j) {
                    if (j == i) continue;
                    const LabelizerResidue& q = s.residues[j];
                    if (q.comp_id != excluded || q.ca < 0) continue;
                    if (ls_dist(s, r.ca, q.ca) > options.exclusion_distance)
                        continue;
                    const double d = std::max(depth[j], 1e-10);
                    if (1.0 / d > options.exclusion_exposure) too_close = true;
                }
                row.value = too_close ? 0.001 : 0.999;
                row.status = "scored";
            } else if (par.tag == "tp") {
                static const double radii[5] = {5, 10, 20, 30, 40};
                static const double weights[5] = {100, 40, 15, 5, 2};
                double shell[5] = {0, 0, 0, 0, 0};
                for (int k = 0; k < 5; ++k) {
                    for (std::size_t j = 0; j < nr; ++j) {
                        if (j == i) continue;
                        const LabelizerResidue& q = s.residues[j];
                        if (q.comp_id != "TRP" || q.ca < 0) continue;
                        if (ls_dist(s, r.ca, q.ca) > radii[k]) continue;
                        shell[k] += ls_hse_factor(hse[2 * j + 0], hse[2 * j + 1]);
                    }
                }
                // Make the shells disjoint (tryptophan_proximity.py:76).
                for (int k = 4; k > 0; --k)
                    for (int m = k - 1; m >= 0; --m) shell[k] -= shell[m];
                double initial = 0.0;
                for (int k = 0; k < 5; ++k) initial += weights[k] * shell[k];
                double v;
                if (initial > 25.0) v = 0.0;
                else if (initial < 5.0) v = 1.0;
                else v = 1.0 - (initial - 5.0) / 25.0;
                row.value = v;
                row.status = "scored";
            } else if (par.tag == "ce") {
                static const double radii[4] = {5, 10, 15, 20};
                static const double weights[4] = {100, 40, 20, 10};
                double plus[4] = {0, 0, 0, 0}, minus[4] = {0, 0, 0, 0};
                for (int k = 0; k < 4; ++k) {
                    for (std::size_t j = 0; j < nr; ++j) {
                        if (j == i) continue;
                        const LabelizerResidue& q = s.residues[j];
                        if (q.ca < 0) continue;
                        const double c = ls_formal_charge(q.comp_id);
                        if (c == 0.0) continue;
                        if (ls_dist(s, r.ca, q.ca) > radii[k]) continue;
                        if (c > 0) plus[k] += 1.0; else minus[k] -= 1.0;
                    }
                }
                for (int k = 3; k > 0; --k)
                    for (int m = k - 1; m >= 0; --m) {
                        plus[k] -= plus[m];
                        minus[k] -= minus[m];
                    }
                // The reference computes an exposure-weighted charge and then
                // scores the *unweighted* counts (charge_environment.py:119);
                // reproduced, so the published number is the published number.
                double initial = 0.0;
                for (int k = 0; k < 4; ++k) {
                    const double net = plus[k] + minus[k];
                    const double sum = std::fabs(plus[k]) + std::fabs(minus[k]);
                    const double ratio = sum > 0.0 ? net / sum : 0.0;
                    initial += weights[k] * ratio;
                }
                initial /= 220.0;
                row.value = 1.0 - (initial + 1.0) / 2.0;
                row.status = "scored";
            } else {
                IMP_THROW("labelizer_parameter_scores: unknown tag '" << par.tag << "'",
                          ValueException);
            }
            out.push_back(row);
        }
    }
    std::free(hse);
    std::free(depth);
    return out;
}

// ---------------------------------------------------------------------------
// The combination
// ---------------------------------------------------------------------------

std::vector<LabelizerScore> labelizer_labeling_score(
        const std::vector<LabelizerScore>& parameter_scores,
        const std::vector<LabelizerParameter>& model, const LabelizerOptions& options) {
    // Group by position, keeping the order positions first appear in.
    std::vector<std::string> order;
    std::map<std::string, std::map<std::string, const LabelizerScore*> > by_pos;
    std::map<std::string, const LabelizerScore*> exemplar;
    for (std::size_t i = 0; i < parameter_scores.size(); ++i) {
        const LabelizerScore& r = parameter_scores[i];
        const std::string key = labelizer_residue_key(r.asym_id, r.seq_id);
        if (!by_pos.count(key)) { order.push_back(key); exemplar[key] = &r; }
        by_pos[key][r.score_type] = &r;
    }

    std::vector<LabelizerScore> out;
    out.reserve(order.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
        const std::string& key = order[i];
        const LabelizerScore& ex = *exemplar[key];
        LabelizerScore row;
        row.asym_id = ex.asym_id;
        row.seq_id = ex.seq_id;
        row.comp_id = ex.comp_id;
        row.score_type = "combined";
        row.status = "unavailable";
        row.value = 0.0;

        bool err = false, zero = false;
        std::vector<double> repeated;
        for (std::size_t p = 0; p < model.size(); ++p) {
            const std::string st = labelizer_score_type(model[p].tag);
            std::map<std::string, const LabelizerScore*>::const_iterator it =
                    by_pos[key].find(st);
            if (it == by_pos[key].end() || it->second->status != "scored") {
                err = true;
                continue;
            }
            const double v = it->second->value;
            // The published zero-veto looks at every term, including the ones
            // whose weight is zero (labeling_score.py:166). Weight zero does
            // not mean ignored, and that is a defect, not a convention.
            if (v == 0.0) {
                if (options.model == LABELIZER_MODEL_PUBLISHED || model[p].weight > 0)
                    zero = true;
            }
            for (int w = 0; w < model[p].weight; ++w) repeated.push_back(v);
        }

        if (err) {
            row.status = "unavailable";
        } else if (zero) {
            row.value = 0.0;
            row.status = "scored";
        } else if (repeated.empty()) {
            row.status = "unavailable";
        } else {
            double prod = 1.0;
            for (std::size_t k = 0; k < repeated.size(); ++k) prod *= repeated[k];
            row.value = std::pow(prod, 1.0 / static_cast<double>(repeated.size()));
            row.status = "scored";
        }
        out.push_back(row);
    }
    return out;
}

std::vector<LabelizerScore> labelizer_score_structure(const std::string& pdb_path,
                                        const std::vector<LabelizerParameter>& model,
                                        const LabelizerOptions& options,
                                        const std::string& conservation_path) {
    const LabelizerStructure s = labelizer_read_structure(pdb_path);
    std::map<std::string, double> conservation;
    if (!conservation_path.empty()) {
        conservation = labelizer_read_consurf(conservation_path);
    }
    std::vector<LabelizerScore> out = labelizer_parameter_scores(s, model, options, conservation);
    const std::vector<LabelizerScore> combined = labelizer_labeling_score(out, model, options);
    out.insert(out.end(), combined.begin(), combined.end());
    return out;
}

std::map<std::string, double> labelizer_combined_by_key(
        const std::vector<LabelizerScore>& scores) {
    std::map<std::string, double> out;
    for (std::size_t i = 0; i < scores.size(); ++i) {
        if (scores[i].score_type != "combined") continue;
        if (scores[i].status != "scored") continue;
        out[labelizer_residue_key(scores[i].asym_id, scores[i].seq_id)] =
                scores[i].value;
    }
    return out;
}

IMPBFF_END_NAMESPACE
