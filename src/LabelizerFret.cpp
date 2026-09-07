/**
 * \file LabelizerFret.cpp
 * \brief The Labelizer FRET pair score: which two labelling sites make the
 *        most informative FRET assay.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/LabelizerFret.h>

#include <IMP/bff/AV.h>
#include <IMP/bff/AVBuilder.h>
#include <IMP/bff/AVModel.h>
#include <IMP/bff/StatesDistance.h>
#include <IMP/bff/bff_config.h>
#include <IMP/bff/internal/OutputView.h>

#include <IMP/bff/Base.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <set>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! The van der Waals radius the reference gives every atom (`config.py`).
const double LF_COLLISION_RADIUS = 1.7;

struct LfPoint {
    double x, y, z;
    LfPoint() : x(0), y(0), z(0) {}
    LfPoint(double a, double b, double c) : x(a), y(b), z(c) {}
};

double lf_distance(const LfPoint& a, const LfPoint& b) {
    const double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

//! A dye at a site: its mean position, and its cloud when it has one.
/*! The cloud is kept, not reduced to its mean, because the pair distance the
    reference reports is a property of the two *distributions* -- see
    #IMP::bff::LlFretOptions::distance_type. */
struct LfSite {
    int residue;
    std::string key;
    LfPoint mean;
    bool placed;
    //! `States` rather than `AccessibleVolume`: the pair distance needs only
    //! the cloud, and the density grid an AV also carries is a large thing to
    //! copy per site.
    States av;
    bool has_av;
    LfSite() : residue(-1), placed(false), has_av(false) {}
};

//! The distance between two placed dyes.
/*! Two clouds give the distance the options ask for; two points can only give
    the distance between them, which is what the reference's point models do. */
double lf_pair_distance(const LfSite& a, const LfSite& b,
                        const LlFretOptions& options) {
    if (a.has_av && b.has_av) {
        return model_distance(a.av, b.av, options.distance_type,
                              options.forster_radius);
    }
    return lf_distance(a.mean, b.mean);
}

}  // namespace

// ---------------------------------------------------------------------------
// The scores
// ---------------------------------------------------------------------------

double ll_joined_label_score(const std::vector<double>& scores, LlModel model) {
    if (scores.empty()) return 0.0;
    double prod = 1.0;
    for (std::size_t i = 0; i < scores.size(); ++i) prod *= scores[i];
    if (prod < 0.0) return 0.0;
    // The published exponent is a flat 0.5 whatever the count
    // (fret_score.py:220); the geometric mean is the commented-out line below
    // it. They agree for two scores and not for four.
    const double exponent = (model == LL_MODEL_PUBLISHED)
                                    ? 0.5
                                    : 1.0 / static_cast<double>(scores.size());
    return std::pow(prod, exponent);
}

double ll_pair_score_single(double joined_label_score, double distance,
                            double forster_radius) {
    const double e = fret_efficiency(distance, forster_radius);
    return joined_label_score * (1.0 - 2.0 * std::fabs(e - 0.5));
}

double ll_pair_score_double(double joined_label_score, double distance_1,
                            double distance_2, double forster_radius) {
    const double e1 = fret_efficiency(distance_1, forster_radius);
    const double e2 = fret_efficiency(distance_2, forster_radius);
    return joined_label_score * std::fabs(e1 - e2);
}

double ll_pair_score_negative_control(double joined_label_score,
                                      double distance_1, double distance_2,
                                      double forster_radius) {
    const double e1 = fret_efficiency(distance_1, forster_radius);
    const double e2 = fret_efficiency(distance_2, forster_radius);
    const double mid = 1.0 - std::fabs(1.0 - (e1 + e2));
    const double flat = std::max(0.0, 1.0 - 20.0 * std::fabs(e1 - e2));
    return joined_label_score * mid * flat;
}

// ---------------------------------------------------------------------------
// Where the dye is
// ---------------------------------------------------------------------------

void ll_alpha_cone_mean_position(const LlStructure& s, int residue,
                                 const LlFretOptions& options,
                                 double** out_view, int* n_out_view) {
    double* cb = 0;
    int n_cb = 0;
    ll_cbeta_position(s, residue, &cb, &n_cb);
    if (n_cb != 3) {
        std::free(cb);
        internal::new_double_view(0, out_view, n_out_view);
        return;
    }
    const LfPoint c(cb[0], cb[1], cb[2]);
    std::free(cb);

    // The centroid of everything within the cone radius of the attachment.
    double sx = 0, sy = 0, sz = 0;
    std::size_t n = 0;
    const double r2 = options.alpha_cone_radius * options.alpha_cone_radius;
    for (std::size_t a = 0; a < s.size(); ++a) {
        const double dx = s.xyz[3 * a + 0] - c.x;
        const double dy = s.xyz[3 * a + 1] - c.y;
        const double dz = s.xyz[3 * a + 2] - c.z;
        if (dx * dx + dy * dy + dz * dz > r2) continue;
        sx += s.xyz[3 * a + 0]; sy += s.xyz[3 * a + 1]; sz += s.xyz[3 * a + 2];
        ++n;
    }
    if (n == 0) {
        internal::new_double_view(0, out_view, n_out_view);
        return;
    }
    const LfPoint centroid(sx / n, sy / n, sz / n);
    const double d = lf_distance(c, centroid);
    if (d < 1e-9) {
        internal::new_double_view(0, out_view, n_out_view);
        return;
    }

    double off = options.alpha_cone_offset;
    if (off < 0.0) {
        // b + 0.54*ll - 0.0225*ll^2 (fret_score.py:344): an empirical fit of
        // how far the linker reaches, not a derived quantity.
        const double r_min = std::min(options.r1, std::min(options.r2, options.r3));
        const double b = std::max(LF_COLLISION_RADIUS,
                                  2.0 * r_min - LF_COLLISION_RADIUS);
        const double ll = options.linker_length;
        off = b + 0.54 * ll - 0.0225 * ll * ll;
    }
    const double reach = (0.75 - d / options.alpha_cone_radius) *
                         (options.linker_length + off);
    const LfPoint dir((c.x - centroid.x) / d, (c.y - centroid.y) / d,
                      (c.z - centroid.z) / d);

    double* b = internal::new_double_view(3, out_view, n_out_view);
    if (b) {
        b[0] = c.x + reach * dir.x;
        b[1] = c.y + reach * dir.y;
        b[2] = c.z + reach * dir.z;
    }
}

void ll_probe_mean_position(const LlStructure& s, const std::string& pdb_path,
                          int residue, const LlFretOptions& options,
                          double** out_view, int* n_out_view) {
    if (residue < 0 || residue >= static_cast<int>(s.residues.size())) {
        internal::new_double_view(0, out_view, n_out_view);
        return;
    }
    if (options.probe_model == PROBE_MODEL_CBETA) {
        ll_cbeta_position(s, residue, out_view, n_out_view);
        return;
    }
    if (options.probe_model == PROBE_MODEL_ALPHA_CONE) {
        ll_alpha_cone_mean_position(s, residue, options, out_view, n_out_view);
        return;
    }

    // The real cloud, through the module's own builder. A site whose volume
    // comes back empty is a site no dye fits at -- the reference calls that
    // "fewer than a thousand grid points" and skips it (fret_score.py:708);
    // here an empty cloud says the same thing without a magic count.
    const LlResidue& r = s.residues[residue];
    const std::string atom = r.cb >= 0 ? "CB" : "CA";
    AccessibleVolume av = compute_av_from_structure(
            pdb_path, r.chain, r.seq_id, atom, options.linker_length,
            options.linker_width, options.r1, options.r2, options.r3,
            options.grid_resolution);
    double* mp = 0;
    int n_mp = 0;
    av.get_mean_position(&mp, &n_mp);
    if (n_mp != 3) {
        std::free(mp);
        internal::new_double_view(0, out_view, n_out_view);
        return;
    }
    double* b = internal::new_double_view(3, out_view, n_out_view);
    if (b) { b[0] = mp[0]; b[1] = mp[1]; b[2] = mp[2]; }
    std::free(mp);
}

// ---------------------------------------------------------------------------
// The pairs
// ---------------------------------------------------------------------------

namespace {

//! Every site above the label-score threshold, with its dye position placed.
/*! The placement happens once per site. The reference places the inner site
    again for every outer site, which is where its exact mode spends its time. */
std::vector<LfSite> lf_place_sites(
        const LlStructure& s, const std::string& pdb_path,
        const std::map<std::string, double>& label_scores,
        const LlFretOptions& options, ProbeModel probe_model) {
    LlFretOptions o = options;
    o.probe_model = probe_model;
    std::vector<LfSite> out;
    for (std::size_t i = 0; i < s.residues.size(); ++i) {
        const LlResidue& r = s.residues[i];
        const std::string key = ll_residue_key(r.chain, r.seq_id);
        std::map<std::string, double>::const_iterator it = label_scores.find(key);
        if (it == label_scores.end()) continue;
        if (it->second < options.label_score_threshold) continue;

        LfSite site;
        site.residue = static_cast<int>(i);
        site.key = key;

        if (probe_model == PROBE_MODEL_ACCESSIBLE_VOLUME) {
            // Keep the cloud: the pair distance is a property of the two
            // distributions, not of their two mean positions.
            const std::string atom = r.cb >= 0 ? "CB" : "CA";
            const AccessibleVolume built = compute_av_from_structure(
                    pdb_path, r.chain, r.seq_id, atom, options.linker_length,
                    options.linker_width, options.r1, options.r2, options.r3,
                    options.grid_resolution);

            double* pts = 0; int n_pts = 0;
            built.get_points(&pts, &n_pts);
            std::vector<double> cloud(pts, pts + std::max(0, n_pts));
            std::free(pts);
            double* ap = 0; int n_ap = 0;
            built.get_attachment_point(&ap, &n_ap);
            const std::vector<double> anchor(ap, ap + std::max(0, n_ap));
            std::free(ap);

            site.av = States(cloud, anchor);
            double* mp = 0;
            int n_mp = 0;
            site.av.get_mean_position(&mp, &n_mp);
            if (n_mp == 3) {
                site.mean = LfPoint(mp[0], mp[1], mp[2]);
                site.placed = true;
                site.has_av = site.av.get_n_points() > 0;
            }
            std::free(mp);
        } else {
            double* p = 0;
            int n = 0;
            ll_probe_mean_position(s, pdb_path, static_cast<int>(i), o, &p, &n);
            if (n == 3) {
                site.mean = LfPoint(p[0], p[1], p[2]);
                site.placed = true;
            }
            std::free(p);
        }
        if (site.placed) out.push_back(site);
    }
    return out;
}

bool lf_by_score(const LlPairScore& a, const LlPairScore& b) {
    return a.value > b.value;
}

//! Does this ordered pair of chains match the requested donor/acceptor pair?
/*! With neither chain set every pair is kept, which is what this module did
    before the fields existed. */
bool lf_chains_wanted(const LlFretOptions& options, const std::string& first,
                      const std::string& second) {
    if (options.donor_chain.empty() && options.acceptor_chain.empty())
        return true;
    if (!options.donor_chain.empty() && first != options.donor_chain)
        return false;
    if (!options.acceptor_chain.empty() && second != options.acceptor_chain)
        return false;
    return true;
}

//! The chain of \p state_1 that carries the same site in state 2.
std::string lf_mapped_chain(const LlFretOptions& options,
                            const std::string& chain) {
    std::map<std::string, std::string>::const_iterator it =
            options.chain_map.find(chain);
    return it == options.chain_map.end() ? chain : it->second;
}

}  // namespace

std::vector<LlPairScore> ll_pair_scores(
        const std::string& pdb_path,
        const std::map<std::string, double>& label_scores,
        const LlFretOptions& options) {
    const LlStructure s = ll_read_structure(pdb_path);
    const std::vector<LfSite> sites =
            lf_place_sites(s, pdb_path, label_scores, options, options.probe_model);

    std::vector<LlPairScore> out;
    for (std::size_t i = 0; i < sites.size(); ++i) {
        for (std::size_t j = i + 1; j < sites.size(); ++j) {
            const LlResidue& ri = s.residues[sites[i].residue];
            const LlResidue& rj = s.residues[sites[j].residue];
            // Only i<j is enumerated, so a donor/acceptor chain pair has to be
            // accepted in either order -- (A5, B7) and (B5, A7) are both
            // "donor in A, acceptor in B" labellings of distinct site pairs.
            if (!lf_chains_wanted(options, ri.chain, rj.chain) &&
                !lf_chains_wanted(options, rj.chain, ri.chain))
                continue;
            std::vector<double> ls;
            ls.push_back(label_scores.find(sites[i].key)->second);
            ls.push_back(label_scores.find(sites[j].key)->second);

            LlPairScore p;
            p.asym_id_1 = ri.chain; p.seq_id_1 = ri.seq_id;
            p.asym_id_2 = rj.chain; p.seq_id_2 = rj.seq_id;
            p.joined_label_score = ll_joined_label_score(ls, options.model);
            p.distance = lf_pair_distance(sites[i], sites[j], options);
            p.distance_2 = std::numeric_limits<double>::quiet_NaN();
            p.probe_model = options.probe_model;
            p.value = ll_pair_score_single(p.joined_label_score, p.distance,
                                           options.forster_radius);
            out.push_back(p);
        }
    }
    std::sort(out.begin(), out.end(), lf_by_score);

    // Refine the best few with the expensive dye model. Only the sites those
    // pairs actually use are rebuilt.
    if (options.n_refine > 0 && options.refine_probe_model != options.probe_model) {
        const std::size_t n =
                std::min(static_cast<std::size_t>(options.n_refine), out.size());
        std::set<std::string> wanted;
        for (std::size_t k = 0; k < n; ++k) {
            wanted.insert(ll_residue_key(out[k].asym_id_1, out[k].seq_id_1));
            wanted.insert(ll_residue_key(out[k].asym_id_2, out[k].seq_id_2));
        }
        std::map<std::string, double> subset;
        for (std::set<std::string>::const_iterator it = wanted.begin();
             it != wanted.end(); ++it) {
            subset[*it] = label_scores.find(*it)->second;
        }
        const std::vector<LfSite> refined = lf_place_sites(
                s, pdb_path, subset, options, options.refine_probe_model);
        std::map<std::string, std::size_t> by_key;
        for (std::size_t k = 0; k < refined.size(); ++k) {
            by_key[refined[k].key] = k;
        }
        LlFretOptions refine_options = options;
        refine_options.probe_model = options.refine_probe_model;
        for (std::size_t k = 0; k < n; ++k) {
            const std::string k1 = ll_residue_key(out[k].asym_id_1, out[k].seq_id_1);
            const std::string k2 = ll_residue_key(out[k].asym_id_2, out[k].seq_id_2);
            if (!by_key.count(k1) || !by_key.count(k2)) continue;
            out[k].distance = lf_pair_distance(refined[by_key[k1]],
                                               refined[by_key[k2]],
                                               refine_options);
            out[k].probe_model = options.refine_probe_model;
            out[k].value = ll_pair_score_single(out[k].joined_label_score,
                                                out[k].distance,
                                                options.forster_radius);
        }
        std::sort(out.begin(), out.end(), lf_by_score);
    }
    return out;
}

std::vector<LlPairScore> ll_pair_scores_two_states(
        const std::string& pdb_path_1, const std::string& pdb_path_2,
        const std::map<std::string, double>& label_scores_1,
        const std::map<std::string, double>& label_scores_2,
        const LlFretOptions& options) {
    const LlStructure s1 = ll_read_structure(pdb_path_1);
    const LlStructure s2 = ll_read_structure(pdb_path_2);
    const std::vector<LfSite> a =
            lf_place_sites(s1, pdb_path_1, label_scores_1, options, options.probe_model);
    const std::vector<LfSite> b =
            lf_place_sites(s2, pdb_path_2, label_scores_2, options, options.probe_model);

    // Positions are matched across the two files by CHAIN AND residue number.
    //
    // The reference matches on the residue number alone (fret_score.py:512),
    // and so did this function: `std::map<int, std::size_t>` keyed on
    // `seq_id`. On a multimer that collides -- residue 100 of chains A, B and
    // C are one key, the last one inserted wins, and an inter-protomer pair
    // gets the SAME site on both sides of the second state, so `distance_2`
    // comes back 0 and the score becomes jls*|E(d1) - 1|. Silently, for every
    // pair. A homodimer screen, where the only measurable distance is
    // i(A)-i(B), could not be done at all.
    //
    // #LlFretOptions::chain_map covers the case the old comment was reaching
    // for -- two files that name the same chain differently.
    std::map<std::string, std::size_t> by_site;
    for (std::size_t k = 0; k < b.size(); ++k) {
        const LlResidue& r = s2.residues[b[k].residue];
        by_site[ll_residue_key(r.chain, r.seq_id)] = k;
    }

    std::vector<LlPairScore> out;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const LlResidue& ri = s1.residues[a[i].residue];
        const std::string ki =
                ll_residue_key(lf_mapped_chain(options, ri.chain), ri.seq_id);
        if (!by_site.count(ki)) continue;
        for (std::size_t j = i + 1; j < a.size(); ++j) {
            const LlResidue& rj = s1.residues[a[j].residue];
            if (!lf_chains_wanted(options, ri.chain, rj.chain) &&
                !lf_chains_wanted(options, rj.chain, ri.chain))
                continue;
            const std::string kj =
                    ll_residue_key(lf_mapped_chain(options, rj.chain), rj.seq_id);
            if (!by_site.count(kj)) continue;

            const LfSite& bi = b[by_site[ki]];
            const LfSite& bj = b[by_site[kj]];
            std::vector<double> ls;
            ls.push_back(label_scores_1.find(a[i].key)->second);
            ls.push_back(label_scores_1.find(a[j].key)->second);
            ls.push_back(label_scores_2.find(bi.key)->second);
            ls.push_back(label_scores_2.find(bj.key)->second);

            LlPairScore p;
            p.asym_id_1 = ri.chain; p.seq_id_1 = ri.seq_id;
            p.asym_id_2 = rj.chain; p.seq_id_2 = rj.seq_id;
            p.joined_label_score = ll_joined_label_score(ls, options.model);
            p.distance = lf_pair_distance(a[i], a[j], options);
            p.distance_2 = lf_pair_distance(bi, bj, options);
            p.probe_model = options.probe_model;
            p.value = ll_pair_score_double(p.joined_label_score, p.distance,
                                           p.distance_2, options.forster_radius);
            out.push_back(p);
        }
    }
    std::sort(out.begin(), out.end(), lf_by_score);
    return out;
}

int ll_cbeta_difference_map(const LlStructure& s1, const LlStructure& s2,
                            double** out_view, int* n_out_view,
                            const std::string& chain) {
    std::map<int, LfPoint> p1, p2;
    for (int pass = 0; pass < 2; ++pass) {
        const LlStructure& s = pass == 0 ? s1 : s2;
        std::map<int, LfPoint>& into = pass == 0 ? p1 : p2;
        for (std::size_t i = 0; i < s.residues.size(); ++i) {
            // The matrix is indexed by residue number alone, so on a multimer
            // every chain lands on the same row and the last one wins. Taking
            // one chain is the only reading of this map that means anything.
            if (!chain.empty() && s.residues[i].chain != chain) continue;
            double* cb = 0;
            int n = 0;
            ll_cbeta_position(s, static_cast<int>(i), &cb, &n);
            if (n == 3) into[s.residues[i].seq_id] = LfPoint(cb[0], cb[1], cb[2]);
            std::free(cb);
        }
    }
    if (p1.empty() || p2.empty()) {
        internal::new_double_view(0, out_view, n_out_view);
        return 0;
    }
    const int lo = std::min(p1.begin()->first, p2.begin()->first);
    const int hi = std::max(p1.rbegin()->first, p2.rbegin()->first);
    const std::size_t n = static_cast<std::size_t>(hi - lo + 1);

    double* m = internal::new_double_view(n * n, out_view, n_out_view);
    if (!m) return lo;
    for (int i = lo; i <= hi; ++i) {
        for (int j = lo; j <= hi; ++j) {
            if (!p1.count(i) || !p1.count(j) || !p2.count(i) || !p2.count(j))
                continue;
            // A residue missing from either structure stays zero: there is no
            // change to report, and inventing one would be worse.
            m[static_cast<std::size_t>(i - lo) * n + (j - lo)] =
                    lf_distance(p2[i], p2[j]) - lf_distance(p1[i], p1[j]);
        }
    }
    return lo;
}

IMPBFF_END_NAMESPACE
