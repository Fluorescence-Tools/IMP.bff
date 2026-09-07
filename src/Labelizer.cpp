/**
 * \file Labelizer.cpp
 * \brief The native Labelizer: features, site score, pair score, the .pto store.
 *
 * Sections in the order of IMP/bff/Labelizer.h; each is marked with the file
 * it came from.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

// -------- from LabelizerFeatures.cpp --------
/**
 * (formerly LabelizerFeatures.cpp, now a section of this file)
 * \brief The structural quantities a label-site score is computed from:
 *        secondary structure, half-sphere exposure, relative solvent
 *        accessibility, residue depth, and an imported conservation grade.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/Labelizer.h>

#include <IMP/bff/AVBuilder.h>
#include <IMP/bff/SolventAccessibleSurface.h>
#include <IMP/bff/bff_config.h>
#include <IMP/bff/internal/OutputView.h>

#include <IMP/bff/Base.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <set>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

namespace {

// ---------------------------------------------------------------------------
// small vector helpers -- three doubles, no allocation
// ---------------------------------------------------------------------------

struct LfVec {
    double x, y, z;
    LfVec() : x(0), y(0), z(0) {}
    LfVec(double a, double b, double c) : x(a), y(b), z(c) {}
};

inline LfVec lf_at(const std::vector<double>& xyz, int i) {
    return LfVec(xyz[3 * i + 0], xyz[3 * i + 1], xyz[3 * i + 2]);
}
inline LfVec lf_sub(const LfVec& a, const LfVec& b) {
    return LfVec(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline LfVec lf_add(const LfVec& a, const LfVec& b) {
    return LfVec(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline LfVec lf_scale(const LfVec& a, double s) {
    return LfVec(a.x * s, a.y * s, a.z * s);
}
inline double lf_dot(const LfVec& a, const LfVec& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline LfVec lf_cross(const LfVec& a, const LfVec& b) {
    return LfVec(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                 a.x * b.y - a.y * b.x);
}
inline double lf_norm(const LfVec& a) { return std::sqrt(lf_dot(a, a)); }
inline LfVec lf_unit(const LfVec& a) {
    const double n = lf_norm(a);
    return n > 1e-12 ? lf_scale(a, 1.0 / n) : LfVec(0, 0, 0);
}
inline double lf_dist(const LfVec& a, const LfVec& b) {
    return lf_norm(lf_sub(a, b));
}

std::string lf_strip(const std::string& s) {
    std::size_t b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    std::size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

// ---------------------------------------------------------------------------
// a uniform grid, so the surface work is not quadratic
// ---------------------------------------------------------------------------

//! Bucket points on a lattice of the query cutoff, and walk the 27 neighbours.
class LfGrid {
public:
    LfGrid(const std::vector<double>& pts, double cell)
        : pts_(pts), cell_(cell > 1e-6 ? cell : 1.0) {
        const std::size_t n = pts.size() / 3;
        if (n == 0) return;
        lo_ = lf_at(pts, 0);
        LfVec hi = lo_;
        for (std::size_t i = 1; i < n; ++i) {
            const LfVec p = lf_at(pts, static_cast<int>(i));
            lo_.x = std::min(lo_.x, p.x); hi.x = std::max(hi.x, p.x);
            lo_.y = std::min(lo_.y, p.y); hi.y = std::max(hi.y, p.y);
            lo_.z = std::min(lo_.z, p.z); hi.z = std::max(hi.z, p.z);
        }
        nx_ = std::max(1, static_cast<int>((hi.x - lo_.x) / cell_) + 1);
        ny_ = std::max(1, static_cast<int>((hi.y - lo_.y) / cell_) + 1);
        nz_ = std::max(1, static_cast<int>((hi.z - lo_.z) / cell_) + 1);
        cells_.resize(static_cast<std::size_t>(nx_) * ny_ * nz_);
        for (std::size_t i = 0; i < n; ++i) {
            cells_[index(lf_at(pts, static_cast<int>(i)))].push_back(
                    static_cast<int>(i));
        }
    }

    //! Indices of every point within \p r of \p q, appended to \p out.
    void within(const LfVec& q, double r, std::vector<int>& out) const {
        out.clear();
        if (cells_.empty()) return;
        const double r2 = r * r;
        int cx, cy, cz;
        coords(q, cx, cy, cz);
        const int span = static_cast<int>(r / cell_) + 1;
        for (int i = std::max(0, cx - span); i <= std::min(nx_ - 1, cx + span); ++i)
            for (int j = std::max(0, cy - span); j <= std::min(ny_ - 1, cy + span); ++j)
                for (int k = std::max(0, cz - span); k <= std::min(nz_ - 1, cz + span); ++k) {
                    const std::vector<int>& c =
                            cells_[(static_cast<std::size_t>(i) * ny_ + j) * nz_ + k];
                    for (std::size_t m = 0; m < c.size(); ++m) {
                        if (lf_dot(lf_sub(q, lf_at(pts_, c[m])),
                                   lf_sub(q, lf_at(pts_, c[m]))) < r2) {
                            out.push_back(c[m]);
                        }
                    }
                }
    }

    //! Squared distance to the nearest point, growing the search until one is
    //! found. Returns a large number when the grid is empty.
    double nearest2(const LfVec& q) const {
        if (cells_.empty()) return 1e30;
        double best = 1e30;
        std::vector<int> hit;
        for (double r = cell_; r < 1e4; r *= 2.0) {
            within(q, r, hit);
            for (std::size_t m = 0; m < hit.size(); ++m) {
                const LfVec d = lf_sub(q, lf_at(pts_, hit[m]));
                best = std::min(best, lf_dot(d, d));
            }
            // A hit inside r is only certainly the nearest once the search
            // radius covers it; one more doubling is enough to be sure.
            if (best < r * r) {
                within(q, std::sqrt(best) + cell_, hit);
                for (std::size_t m = 0; m < hit.size(); ++m) {
                    const LfVec d = lf_sub(q, lf_at(pts_, hit[m]));
                    best = std::min(best, lf_dot(d, d));
                }
                return best;
            }
        }
        return best;
    }

private:
    void coords(const LfVec& q, int& cx, int& cy, int& cz) const {
        cx = std::min(nx_ - 1, std::max(0, static_cast<int>((q.x - lo_.x) / cell_)));
        cy = std::min(ny_ - 1, std::max(0, static_cast<int>((q.y - lo_.y) / cell_)));
        cz = std::min(nz_ - 1, std::max(0, static_cast<int>((q.z - lo_.z) / cell_)));
    }
    std::size_t index(const LfVec& q) const {
        int cx, cy, cz;
        coords(q, cx, cy, cz);
        return (static_cast<std::size_t>(cx) * ny_ + cy) * nz_ + cz;
    }

    const std::vector<double>& pts_;
    double cell_;
    LfVec lo_;
    int nx_, ny_, nz_;
    std::vector<std::vector<int> > cells_;
};

}  // namespace

// ---------------------------------------------------------------------------
// The residue view
// ---------------------------------------------------------------------------

const std::vector<std::string>& ll_standard_residues() {
    static const std::vector<std::string> names = [] {
        const char* n[20] = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                             "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                             "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
        return std::vector<std::string>(n, n + 20);
    }();
    return names;
}

std::string ll_residue_key(const std::string& chain, int seq_id) {
    std::ostringstream os;
    os << chain << seq_id;
    return os.str();
}

LlStructure ll_read_structure(const std::string& pdb_path, bool protein_only,
                              int model) {
    const std::vector<PDBAtomRecord> records = read_pdb_records(pdb_path);
    if (records.empty()) {
        IMP_THROW("ll_read_structure: no coordinates in " << pdb_path,
                  IOException);
    }

    // read_pdb_records() does not honour MODEL/ENDMDL -- it returns every
    // ATOM record in the file, so a 20-model NMR entry arrives twenty times
    // over. The models repeat the same atoms in the same order, so a model
    // boundary is the point at which the first record's identity recurs.
    // Splitting here rather than adding a second parser keeps one PDB reader
    // in the package, which is the point of PRD-117's consolidation.
    std::vector<std::size_t> starts;
    starts.push_back(0);
    for (std::size_t i = 1; i < records.size(); ++i) {
        if (records[i].chain == records[0].chain &&
            records[i].resseq == records[0].resseq &&
            records[i].atom_name == records[0].atom_name) {
            starts.push_back(i);
        }
    }
    starts.push_back(records.size());
    const int n_models = static_cast<int>(starts.size()) - 1;
    if (model < 0 || model >= n_models) {
        IMP_THROW("ll_read_structure: " << pdb_path << " has " << n_models
                  << " model(s); model " << model << " was asked for",
                  ValueException);
    }
    const std::size_t first = starts[model], last = starts[model + 1];

    std::set<std::string> standard(ll_standard_residues().begin(),
                                   ll_standard_residues().end());

    LlStructure s;
    std::string cur_chain;
    int cur_seq = 0;
    bool have_residue = false;
    for (std::size_t i = first; i < last; ++i) {
        const PDBAtomRecord& r = records[i];
        const std::string res = lf_strip(r.res_name);
        if (protein_only && standard.find(res) == standard.end()) continue;

        if (!have_residue || r.chain != cur_chain || r.resseq != cur_seq) {
            LlResidue nr;
            nr.chain = r.chain;
            nr.seq_id = r.resseq;
            nr.comp_id = res;
            s.residues.push_back(nr);
            cur_chain = r.chain;
            cur_seq = r.resseq;
            have_residue = true;
        }
        const int idx = static_cast<int>(s.vdw.size());
        s.xyz.push_back(r.x);
        s.xyz.push_back(r.y);
        s.xyz.push_back(r.z);
        s.vdw.push_back(r.vdw_radius);
        const std::string name = lf_strip(r.atom_name);
        s.atom_name.push_back(name);

        LlResidue& cr = s.residues.back();
        cr.atoms.push_back(idx);
        if (name == "N") cr.n = idx;
        else if (name == "CA") cr.ca = idx;
        else if (name == "C") cr.c = idx;
        else if (name == "O") cr.o = idx;
        else if (name == "CB") cr.cb = idx;
    }
    return s;
}

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------

namespace {

//! The Cbeta of a residue, real or virtual. Returns false when unplaceable.
bool lf_cbeta(const LlStructure& s, int ri, LfVec& out) {
    const LlResidue& r = s.residues[ri];
    if (r.cb >= 0) {
        out = lf_at(s.xyz, r.cb);
        return true;
    }
    if (r.ca < 0 || r.n < 0 || r.c < 0) return false;
    // Bio.PDB.HSExposure's glycine construction, as fret_score.py:293 uses it:
    // the N vector centred on CA, rotated -120 degrees about the CA->C axis.
    const LfVec ca = lf_at(s.xyz, r.ca);
    const LfVec n = lf_sub(lf_at(s.xyz, r.n), ca);
    const LfVec axis = lf_unit(lf_sub(lf_at(s.xyz, r.c), ca));
    const double theta = -120.0 * M_PI / 180.0;
    const double ct = std::cos(theta), st = std::sin(theta);
    // Rodrigues' rotation of n about axis by theta.
    const LfVec rot = lf_add(lf_add(lf_scale(n, ct), lf_scale(lf_cross(axis, n), st)),
                             lf_scale(axis, lf_dot(axis, n) * (1.0 - ct)));
    out = lf_add(ca, rot);
    return true;
}

}  // namespace

void ll_cbeta_position(const LlStructure& s, int residue, double** out_view,
                       int* n_out_view) {
    LfVec cb;
    if (residue < 0 || residue >= static_cast<int>(s.residues.size()) ||
        !lf_cbeta(s, residue, cb)) {
        internal::new_double_view(0, out_view, n_out_view);
        return;
    }
    double* b = internal::new_double_view(3, out_view, n_out_view);
    if (b) { b[0] = cb.x; b[1] = cb.y; b[2] = cb.z; }
}

void ll_half_sphere_exposure(const LlStructure& s, double radius,
                             double** out_view, int* n_out_view) {
    const std::size_t nr = s.residues.size();
    double* out = internal::new_double_view(2 * nr, out_view, n_out_view);
    if (!out || nr == 0) return;

    // The CA positions, so the neighbour search is over residues not atoms.
    std::vector<double> ca_xyz;
    std::vector<int> ca_res;
    ca_xyz.reserve(3 * nr);
    for (std::size_t i = 0; i < nr; ++i) {
        if (s.residues[i].ca < 0) continue;
        const LfVec p = lf_at(s.xyz, s.residues[i].ca);
        ca_xyz.push_back(p.x); ca_xyz.push_back(p.y); ca_xyz.push_back(p.z);
        ca_res.push_back(static_cast<int>(i));
    }
    LfGrid grid(ca_xyz, radius);
    std::vector<int> hit;

    for (std::size_t i = 0; i < nr; ++i) {
        const LlResidue& r = s.residues[i];
        LfVec cb;
        if (r.ca < 0 || !lf_cbeta(s, static_cast<int>(i), cb)) continue;
        const LfVec ca = lf_at(s.xyz, r.ca);
        const LfVec dir = lf_unit(lf_sub(cb, ca));
        grid.within(ca, radius, hit);
        double up = 0.0, down = 0.0;
        for (std::size_t k = 0; k < hit.size(); ++k) {
            const int rj = ca_res[hit[k]];
            if (rj == static_cast<int>(i)) continue;
            const LfVec d = lf_sub(lf_at(ca_xyz, hit[k]), ca);
            if (lf_dot(d, dir) > 0.0) up += 1.0; else down += 1.0;
        }
        out[2 * i + 0] = up;
        out[2 * i + 1] = down;
    }
}

// ---------------------------------------------------------------------------
// Secondary structure -- Kabsch and Sander (1983)
// ---------------------------------------------------------------------------

const double LL_DSSP_HBOND_ENERGY = -0.5;

std::string ll_dssp(const LlStructure& s) {
    const int nr = static_cast<int>(s.residues.size());
    std::string ss(static_cast<std::size_t>(nr), '-');
    if (nr == 0) return ss;

    // Backbone, and the amide hydrogen DSSP places rather than reads.
    std::vector<LfVec> N(nr), CA(nr), C(nr), O(nr), H(nr);
    std::vector<char> ok(nr, 0), has_h(nr, 0);
    for (int i = 0; i < nr; ++i) {
        const LlResidue& r = s.residues[i];
        if (r.n < 0 || r.ca < 0 || r.c < 0 || r.o < 0) continue;
        N[i] = lf_at(s.xyz, r.n);
        CA[i] = lf_at(s.xyz, r.ca);
        C[i] = lf_at(s.xyz, r.c);
        O[i] = lf_at(s.xyz, r.o);
        ok[i] = 1;
    }
    // Adjacent in the array is not adjacent in the chain. A residue the reader
    // dropped -- an unnatural amino acid, a modified residue -- and an
    // unresolved loop in a crystal structure both leave a gap in the author
    // numbering, and inferring a peptide bond across one produces a confident
    // wrong assignment rather than an error. `linked[i]` is the only thing
    // below that may assume i-1 and i are one bond apart.
    std::vector<char> linked(nr, 0);
    for (int i = 1; i < nr; ++i) {
        linked[i] = (s.residues[i].chain == s.residues[i - 1].chain &&
                     s.residues[i].seq_id == s.residues[i - 1].seq_id + 1)
                            ? 1
                            : 0;
    }
    // Residues i .. i+n are one unbroken stretch of chain.
    struct Run {
        const std::vector<char>& linked;
        int nr;
        bool operator()(int i, int n) const {
            if (i < 0 || i + n >= nr) return false;
            for (int k = i + 1; k <= i + n; ++k)
                if (!linked[k]) return false;
            return true;
        }
    } run = {linked, nr};

    for (int i = 1; i < nr; ++i) {
        if (!ok[i] || !ok[i - 1]) continue;
        if (!linked[i]) continue;
        if (s.residues[i].comp_id == "PRO") continue;  // no amide hydrogen
        const LfVec co = lf_unit(lf_sub(C[i - 1], O[i - 1]));
        H[i] = lf_add(N[i], co);
        has_h[i] = 1;
    }

    // hbond[i][j]: the C=O of i donates to the N-H of j.
    std::vector<std::vector<char> > hb(nr, std::vector<char>(nr, 0));
    const double q1q2f = 0.42 * 0.20 * 332.0;
    for (int i = 0; i < nr; ++i) {
        if (!ok[i]) continue;
        for (int j = 0; j < nr; ++j) {
            if (i == j || !ok[j] || !has_h[j]) continue;
            if (std::abs(i - j) < 2) continue;
            if (lf_dist(CA[i], CA[j]) > 9.0) continue;
            const double rON = lf_dist(O[i], N[j]);
            const double rCH = lf_dist(C[i], H[j]);
            const double rOH = lf_dist(O[i], H[j]);
            const double rCN = lf_dist(C[i], N[j]);
            if (rON < 1e-3 || rCH < 1e-3 || rOH < 1e-3 || rCN < 1e-3) continue;
            const double e = q1q2f * (1.0 / rON + 1.0 / rCH - 1.0 / rOH - 1.0 / rCN);
            if (e < LL_DSSP_HBOND_ENERGY) hb[i][j] = 1;
        }
    }

    // Chain-contiguous n-turns.
    std::vector<char> t3(nr, 0), t4(nr, 0), t5(nr, 0);
    for (int i = 0; i < nr; ++i) {
        for (int n = 3; n <= 5; ++n) {
            if (!run(i, n)) continue;
            const int j = i + n;
            if (!hb[i][j]) continue;
            if (n == 3) t3[i] = 1; else if (n == 4) t4[i] = 1; else t5[i] = 1;
        }
    }

    // Bridges, then ladders. A ladder of one bridge is B, longer is E.
    //
    // Ladders are joined through **beta-bulges**, which is not a refinement:
    // without it a bulged sheet fragments into isolated bridges, and on
    // maltose-binding protein that cost fourteen strand residues against the
    // reference (E -> B and E -> S). A bulge is a step of one residue on one
    // strand against up to five on the other, so two bridges of the same type
    // belong to one ladder when they advance together, one of them by exactly
    // one residue.
    struct Bridge { int i, j; bool parallel; };
    std::vector<Bridge> bridges;
    for (int i = 1; i + 1 < nr; ++i) {
        for (int j = i + 3; j + 1 < nr; ++j) {
            const bool anti = (hb[i][j] && hb[j][i]) ||
                              (hb[i - 1][j + 1] && hb[j - 1][i + 1]);
            const bool para = (hb[i - 1][j] && hb[j][i + 1]) ||
                              (hb[j - 1][i] && hb[i][j + 1]);
            if (!anti && !para) continue;
            Bridge br;
            br.i = i; br.j = j; br.parallel = para && !anti ? true : !anti ? true : false;
            br.parallel = para && !anti;
            bridges.push_back(br);
        }
    }

    // Union-find over bridges: same component means same ladder.
    std::vector<int> parent(bridges.size());
    for (std::size_t k = 0; k < parent.size(); ++k) parent[k] = static_cast<int>(k);
    struct Find {
        std::vector<int>& p;
        int operator()(int x) const {
            while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
            return x;
        }
    } find = {parent};

    for (std::size_t a = 0; a < bridges.size(); ++a) {
        for (std::size_t b = a + 1; b < bridges.size(); ++b) {
            if (bridges[a].parallel != bridges[b].parallel) continue;
            const int di = bridges[b].i - bridges[a].i;
            // Antiparallel strands run against each other, so j decreases as
            // i increases; parallel strands advance together.
            const int dj = bridges[a].parallel ? bridges[b].j - bridges[a].j
                                               : bridges[a].j - bridges[b].j;
            if (di < 1 || dj < 1 || di > 5 || dj > 5) continue;
            if (di != 1 && dj != 1) continue;   // a bulge steps 1 on one side
            const int ra = find(static_cast<int>(a)), rb = find(static_cast<int>(b));
            if (ra != rb) parent[ra] = rb;
        }
    }

    std::map<int, int> ladder_size;
    for (std::size_t k = 0; k < bridges.size(); ++k) {
        ladder_size[find(static_cast<int>(k))] += 1;
    }

    std::vector<char> in_bridge(nr, 0), in_ladder(nr, 0);
    for (std::size_t k = 0; k < bridges.size(); ++k) {
        const bool extended = ladder_size[find(static_cast<int>(k))] > 1;
        const int ends[2] = {bridges[k].i, bridges[k].j};
        for (int e = 0; e < 2; ++e) {
            in_bridge[ends[e]] = 1;
            if (extended) in_ladder[ends[e]] = 1;
        }
    }

    // Bends: the curvature of the CA trace over +/- two residues.
    std::vector<char> bend(nr, 0);
    for (int i = 2; i + 2 < nr; ++i) {
        if (!ok[i - 2] || !ok[i] || !ok[i + 2]) continue;
        if (!run(i - 2, 4)) continue;
        const LfVec a = lf_unit(lf_sub(CA[i], CA[i - 2]));
        const LfVec b = lf_unit(lf_sub(CA[i + 2], CA[i]));
        const double c = std::max(-1.0, std::min(1.0, lf_dot(a, b)));
        if (std::acos(c) * 180.0 / M_PI > 70.0) bend[i] = 1;
    }

    // Assignment, in DSSP's order of precedence: H B E G I T S.
    std::vector<char> out(nr, '-');
    // The pi-helix goes first, and that is a version difference, not a
    // preference: DSSP 2 and earlier let a 4-helix win wherever a 4-turn and
    // a 5-turn overlap, and DSSP 3.0 reversed it (Touw 2015). The reference
    // CSVs were produced against the later behaviour -- assigning H first
    // put alpha where the reference has pi on exactly one 4-residue stretch
    // of 1DDB and nowhere else -- so this reproduces what the published
    // numbers were computed with.
    // A 3-10 or pi helix is three or five residues; a shorter leftover is not
    // one. Assigning the remainder after an alpha helix has taken the rest
    // produced one- and two-residue `G` stubs at helix C-termini where the
    // reference has `T` -- six of them on maltose-binding protein.
    struct SpanFree {
        const std::vector<char>& out;
        int nr;
        bool operator()(int i, int n) const {
            if (i < 0 || i + n > nr) return false;
            for (int k = i; k < i + n; ++k)
                if (out[k] != '-') return false;
            return true;
        }
    } span_free = {out, nr};

    for (int i = 1; i < nr; ++i)
        if (t5[i - 1] && t5[i] && run(i, 4) && span_free(i, 5))
            for (int k = i; k < i + 5; ++k) out[k] = 'I';
    for (int i = 1; i < nr; ++i)
        if (t4[i - 1] && t4[i] && run(i, 3))
            for (int k = i; k < std::min(nr, i + 4); ++k)
                if (out[k] == '-') out[k] = 'H';
    for (int i = 0; i < nr; ++i)
        if (in_bridge[i] && out[i] == '-') out[i] = in_ladder[i] ? 'E' : 'B';
    for (int i = 1; i < nr; ++i)
        if (t3[i - 1] && t3[i] && run(i, 2) && span_free(i, 3))
            for (int k = i; k < i + 3; ++k) out[k] = 'G';
    for (int i = 0; i < nr; ++i) {
        if (out[i] != '-') continue;
        // A residue is in a turn when it lies strictly *between* the two
        // residues an n-turn hydrogen-bonds: DSSP marks k+1..k+n-1, not the
        // donor k or the acceptor k+n themselves. Including k marked the
        // turn's own first residue and cost nine residues of agreement
        // against the reference on 1DDB.
        bool turn = false;
        for (int n = 3; n <= 5 && !turn; ++n)
            for (int k = std::max(0, i - n + 1); k < i && !turn; ++k) {
                if (k + n <= i) continue;
                if ((n == 3 && t3[k]) || (n == 4 && t4[k]) || (n == 5 && t5[k]))
                    turn = true;
            }
        if (turn) out[i] = 'T';
    }
    for (int i = 0; i < nr; ++i)
        if (out[i] == '-' && bend[i]) out[i] = 'S';

    for (int i = 0; i < nr; ++i) ss[static_cast<std::size_t>(i)] = out[i];
    return ss;
}

// ---------------------------------------------------------------------------
// Solvent exposure
// ---------------------------------------------------------------------------

std::map<std::string, double> ll_max_asa(LlMaxAsa scale) {
    const char* res[20] = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                           "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                           "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
    // Tien 2013 theoretical; Kabsch and Sander 1983; Miller 1987.
    static const double wilke[20] = {129, 274, 195, 193, 167, 225, 223, 104,
                                     224, 197, 201, 236, 224, 240, 159, 155,
                                     172, 285, 263, 174};
    static const double sander[20] = {106, 248, 157, 163, 135, 198, 194, 84,
                                      184, 169, 164, 205, 188, 197, 136, 130,
                                      142, 227, 222, 142};
    static const double miller[20] = {113, 241, 158, 151, 140, 189, 183, 85,
                                      194, 182, 180, 211, 204, 218, 143, 122,
                                      146, 259, 229, 160};
    const double* v = scale == LL_MAXASA_SANDER ? sander
                    : scale == LL_MAXASA_MILLER ? miller : wilke;
    std::map<std::string, double> m;
    for (int i = 0; i < 20; ++i) m[res[i]] = v[i];
    return m;
}

void ll_residue_sasa(const LlStructure& s, double probe_radius,
                     int n_sphere_points, double** out_view, int* n_out_view) {
    const std::size_t nr = s.residues.size();
    double* out = internal::new_double_view(nr, out_view, n_out_view);
    if (!out || nr == 0) return;

    std::vector<int> all(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) all[i] = static_cast<int>(i);
    const std::vector<double> pts = sphere_points(n_sphere_points);
    const std::vector<double> atom_asa =
            solvent_accessible_surface_area_per_atom(s.xyz, s.vdw, all, pts,
                                                     probe_radius);
    for (std::size_t i = 0; i < nr; ++i) {
        double total = 0.0;
        const std::vector<int>& a = s.residues[i].atoms;
        for (std::size_t k = 0; k < a.size(); ++k) total += atom_asa[a[k]];
        out[i] = total;
    }
}

void ll_relative_solvent_accessibility(const LlStructure& s, LlMaxAsa scale,
                                       double probe_radius, int n_sphere_points,
                                       double** out_view, int* n_out_view) {
    const std::size_t nr = s.residues.size();
    double* out = internal::new_double_view(nr, out_view, n_out_view);
    if (!out || nr == 0) return;

    double* area = 0;
    int n_area = 0;
    ll_residue_sasa(s, probe_radius, n_sphere_points, &area, &n_area);
    const std::map<std::string, double> maxasa = ll_max_asa(scale);
    for (std::size_t i = 0; i < nr; ++i) {
        std::map<std::string, double>::const_iterator it =
                maxasa.find(s.residues[i].comp_id);
        // A residue type the scale does not cover has no relative value.
        // Negative says "not computed"; zero would say "fully buried".
        out[i] = (it == maxasa.end() || it->second <= 0.0)
                         ? -1.0
                         : area[i] / it->second;
    }
    std::free(area);
}

namespace {

//! The solvent-excluded surface as a flat point list.
/*! Rolls the probe centre over every atom and keeps the positions no atom
    occludes, then pulls each back by one probe radius onto the contact
    surface. Shared by the two depth measures so the expensive part happens
    once per call rather than once per observable. */
std::vector<double> lf_build_surface(const LlStructure& s, double probe_radius,
                                     int n_sphere_points) {
    const std::vector<double> unit = sphere_points(n_sphere_points);
    const std::size_t n_unit = unit.size() / 3;
    double r_max = 0.0;
    for (std::size_t a = 0; a < s.size(); ++a) r_max = std::max(r_max, s.vdw[a]);

    LfGrid atoms(s.xyz, r_max + probe_radius);
    std::vector<int> hit;
    std::vector<double> surface;
    surface.reserve(3 * s.size() * n_unit / 8);

    for (std::size_t a = 0; a < s.size(); ++a) {
        const LfVec c = lf_at(s.xyz, static_cast<int>(a));
        const double ra = s.vdw[a] + probe_radius;
        atoms.within(c, ra + r_max + probe_radius, hit);
        for (std::size_t j = 0; j < n_unit; ++j) {
            const LfVec dir(unit[3 * j + 0], unit[3 * j + 1], unit[3 * j + 2]);
            const LfVec p = lf_add(c, lf_scale(dir, ra));
            bool free_point = true;
            for (std::size_t k = 0; k < hit.size(); ++k) {
                const std::size_t b = static_cast<std::size_t>(hit[k]);
                if (b == a) continue;
                const double rb = s.vdw[b] + probe_radius;
                const LfVec d = lf_sub(p, lf_at(s.xyz, static_cast<int>(b)));
                if (lf_dot(d, d) < rb * rb) { free_point = false; break; }
            }
            if (!free_point) continue;
            const LfVec q = lf_add(c, lf_scale(dir, s.vdw[a]));
            surface.push_back(q.x); surface.push_back(q.y); surface.push_back(q.z);
        }
    }
    return surface;
}

}  // namespace

void ll_cbeta_depth(const LlStructure& s, double probe_radius,
                    int n_sphere_points, double** out_view, int* n_out_view) {
    const std::size_t nr = s.residues.size();
    double* out = internal::new_double_view(nr, out_view, n_out_view);
    if (!out || nr == 0 || s.size() == 0) return;

    const std::vector<double> surface =
            lf_build_surface(s, probe_radius, n_sphere_points);
    if (surface.empty()) {
        for (std::size_t i = 0; i < nr; ++i) out[i] = 0.0;
        return;
    }
    LfGrid surf(surface, 3.0);
    for (std::size_t i = 0; i < nr; ++i) {
        LfVec cb;
        if (!lf_cbeta(s, static_cast<int>(i), cb)) { out[i] = -1.0; continue; }
        out[i] = std::sqrt(surf.nearest2(cb));
    }
}

void ll_residue_depth(const LlStructure& s, double probe_radius,
                      int n_sphere_points, double** out_view, int* n_out_view) {
    const std::size_t nr = s.residues.size();
    double* out = internal::new_double_view(nr, out_view, n_out_view);
    if (!out || nr == 0 || s.size() == 0) return;

    const std::vector<double> surface =
            lf_build_surface(s, probe_radius, n_sphere_points);
    if (surface.empty()) {
        for (std::size_t i = 0; i < nr; ++i) out[i] = 0.0;
        return;
    }
    LfGrid surf(surface, 3.0);
    std::vector<double> atom_depth(s.size(), 0.0);
    for (std::size_t a = 0; a < s.size(); ++a) {
        atom_depth[a] = std::sqrt(surf.nearest2(lf_at(s.xyz, static_cast<int>(a))));
    }
    for (std::size_t i = 0; i < nr; ++i) {
        const std::vector<int>& a = s.residues[i].atoms;
        double total = 0.0;
        for (std::size_t k = 0; k < a.size(); ++k) total += atom_depth[a[k]];
        out[i] = a.empty() ? 0.0 : total / static_cast<double>(a.size());
    }
}

void ll_global_charge(const LlStructure& s, double** out_view,
                      int* n_out_view) {
    double* out = internal::new_double_view(3, out_view, n_out_view);
    if (!out) return;
    double net = 0.0, positive = 0.0, negative = 0.0;
    for (std::size_t i = 0; i < s.residues.size(); ++i) {
        const std::string& c = s.residues[i].comp_id;
        // The `ce` term's table, so the two cannot disagree about what is
        // charged. Histidine counts as +1 here; see the header.
        double q = 0.0;
        if (c == "ARG" || c == "LYS" || c == "HIS") q = 1.0;
        else if (c == "ASP" || c == "GLU") q = -1.0;
        if (q == 0.0) continue;
        net += q;
        if (q > 0.0) positive += q; else negative += q;
    }
    out[0] = net;
    out[1] = positive;
    out[2] = negative;
}

// ---------------------------------------------------------------------------
// Conservation, imported
// ---------------------------------------------------------------------------

std::map<std::string, double> ll_read_consurf(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) {
        IMP_THROW("ll_read_consurf: cannot read " << path, IOException);
    }
    std::map<std::string, double> grades;

    const bool is_pdb = path.size() > 4 &&
                        path.compare(path.size() - 4, 4, ".pdb") == 0;
    std::string line;
    if (is_pdb) {
        // The normalised grade rides in the B-factor column, and the first
        // atom of a residue carries it (conservation_score.py:107).
        std::set<std::string> seen;
        while (std::getline(in, line)) {
            if (line.compare(0, 4, "ATOM") != 0 || line.size() < 66) continue;
            const std::string chain = lf_strip(line.substr(21, 1));
            const int seq = std::atoi(line.substr(22, 4).c_str());
            const std::string key = ll_residue_key(chain, seq);
            if (!seen.insert(key).second) continue;
            grades[key] = std::atof(line.substr(60, 6).c_str());
        }
        return grades;
    }

    // A `.grades` table: fifteen header lines, then whitespace-separated
    // columns whose third is `ALA123:A` and whose fourth is the score
    // (conservation_score.py:129).
    int skipped = 0;
    while (std::getline(in, line)) {
        if (skipped < 15) { ++skipped; continue; }
        std::istringstream ls(line);
        std::vector<std::string> f;
        std::string tok;
        while (ls >> tok) f.push_back(tok);
        if (f.size() < 4) continue;
        const std::string& pos = f[2];
        const std::size_t colon = pos.find(':');
        if (colon == std::string::npos || colon < 4) continue;
        const std::string seq = pos.substr(3, colon - 3);
        const std::string chain = pos.substr(colon + 1);
        grades[ll_residue_key(chain, std::atoi(seq.c_str()))] =
                std::atof(f[3].c_str());
    }
    return grades;
}

IMPBFF_END_NAMESPACE

// -------- from LabelizerScore.cpp --------
/**
 * (formerly LabelizerScore.cpp, now a section of this file)
 * \brief The Labelizer label-site score: fitted tables, the seven per-residue
 *        parameters, and the weighted geometric mean.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/DataPaths.h>
#include <IMP/bff/internal/json.h>


#include <map>

IMPBFF_BEGIN_NAMESPACE

namespace {

std::string ls_table_path(const std::string& name) {
    return get_data_path("labelizer/probabilities/" + name +
                         "_P_l_after_s.json");
}

//! Distance between two atoms of a structure.
double ls_dist(const LlStructure& s, int a, int b) {
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

const LlTable& ll_load_table(const std::string& name) {
    static std::map<std::string, LlTable> cache;
    std::map<std::string, LlTable>::iterator it = cache.find(name);
    if (it != cache.end()) return it->second;

    const std::string path = ls_table_path(name);
    std::ifstream in(path.c_str());
    if (!in) {
        IMP_THROW("ll_load_table: no fitted table named '" << name << "' at "
                  << path, IOException);
    }
    nlohmann::json j;
    in >> j;

    LlTable t;
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

std::vector<std::string> ll_available_tables() {
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

double ll_lookup(const LlTable& table, double value) {
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

double ll_lookup_key(const LlTable& table, const std::string& key) {
    std::map<std::string, double>::const_iterator it = table.by_key.find(key);
    if (it == table.by_key.end()) {
        IMP_THROW("ll_lookup_key: '" << key << "' is not in table "
                  << table.name, ValueException);
    }
    return it->second;
}

char ll_one_letter(const std::string& comp_id) {
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

std::vector<LlParameter> ll_model_paper() {
    std::vector<LlParameter> m;
    m.push_back(LlParameter("cs", "N_CS2_Score", 1));
    m.push_back(LlParameter("se", "N_SE11_MEAN_SURFACE_DIST", 1));
    m.push_back(LlParameter("tp", "", 0));
    m.push_back(LlParameter("cr", "C_CR1_Name", 1));
    m.push_back(LlParameter("ss", "C_SS1_SS", 1));
    m.push_back(LlParameter("ce", "", 0));
    return m;
}

std::string ll_score_type(const std::string& tag) {
    if (tag == "cs") return "conservation";
    if (tag == "se") return "solvent_exposure";
    if (tag == "ss") return "secondary_structure";
    if (tag == "ce") return "charge_environment";
    if (tag == "tp") return "tryptophan_proximity";
    if (tag == "cr") return "cysteine_resemblance";
    if (tag == "me") return "methionine_exclusion";
    if (tag == "combined") return "combined";
    IMP_THROW("ll_score_type: unknown parameter tag '" << tag << "'",
              ValueException);
}

// ---------------------------------------------------------------------------
// The parameters
// ---------------------------------------------------------------------------

std::vector<LlScore> ll_parameter_scores(
        const LlStructure& s, const std::vector<LlParameter>& model,
        const LlOptions& options,
        const std::map<std::string, double>& conservation) {
    const std::size_t nr = s.residues.size();
    std::vector<LlScore> out;
    out.reserve(nr * model.size());
    if (nr == 0) return out;

    // Everything the terms below share, computed once. The reference rebuilds
    // the surface and the neighbour search per parameter object.
    const std::string ss = ll_dssp(s);

    double* hse = 0; int n_hse = 0;
    ll_half_sphere_exposure(s, options.hse_radius, &hse, &n_hse);

    double* depth = 0; int n_depth = 0;
    ll_residue_depth(s, options.probe_radius, options.n_sphere_points, &depth,
                     &n_depth);

    // The exclusion term needs the exposure of every candidate residue, which
    // the reference takes as 1/depth (methionin_exclusion.py:62).
    const std::string excluded =
            options.model == LL_MODEL_PUBLISHED ? "MET" : options.exclusion_residue;

    // Relative accessibility is a whole-structure computation and there are
    // three scales; compute a scale at most once, and only if a term asks.
    std::map<int, std::vector<double> > rsa_cache;
    struct RsaFor {
        const LlStructure& s;
        const LlOptions& o;
        std::map<int, std::vector<double> >& cache;
        const std::vector<double>& operator()(const std::string& table) const {
            LlMaxAsa scale = LL_MAXASA_WILKE;
            if (table.find("Sander") != std::string::npos)
                scale = LL_MAXASA_SANDER;
            else if (table.find("Miller") != std::string::npos)
                scale = LL_MAXASA_MILLER;
            std::map<int, std::vector<double> >::iterator it =
                    cache.find(static_cast<int>(scale));
            if (it != cache.end()) return it->second;
            double* v = 0; int n = 0;
            ll_relative_solvent_accessibility(s, scale, o.probe_radius,
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
        const LlStructure& s;
        const LlOptions& o;
        std::vector<double>& cache;
        const std::vector<double>& operator()() const {
            if (!cache.empty()) return cache;
            double* v = 0; int n = 0;
            ll_cbeta_depth(s, o.probe_radius, o.n_sphere_points, &v, &n);
            cache.assign(v, v + n);
            std::free(v);
            return cache;
        }
    } cb_depth_for = {s, options, cb_cache};

    for (std::size_t p = 0; p < model.size(); ++p) {
        const LlParameter& par = model[p];
        const std::string score_type = ll_score_type(par.tag);

        for (std::size_t i = 0; i < nr; ++i) {
            const LlResidue& r = s.residues[i];
            LlScore row;
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
                const LlTable& t = ll_load_table(par.table);
                if (par.table == "C_CR1_Name") {
                    const std::string key(1, ll_one_letter(r.comp_id));
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
                        row.value = ll_lookup(t, mass);
                        row.status = "scored";
                    }
                } else if (par.table == "I_CR4_N_Sidechain") {
                    const double n = ls_sidechain_atoms(r.comp_id);
                    if (n >= 0.0) {
                        row.value = ll_lookup(t, n);
                        row.status = "scored";
                    }
                } else {
                    IMP_THROW("ll_parameter_scores: table '" << par.table
                              << "' is not an implemented observable for tag "
                              << "'cr'. Implemented: C_CR1_Name, C_CR3_Charge, "
                              << "N_CR2_Mass, I_CR4_N_Sidechain.",
                              ValueException);
                }
            } else if (par.tag == "ss") {
                const LlTable& t = ll_load_table(par.table);
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
                    IMP_THROW("ll_parameter_scores: table '" << par.table
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
                    IMP_THROW("ll_parameter_scores: table '" << par.table
                              << "' needs ConSurf fields ll_read_consurf does "
                              << "not import (confidence bounds, colour bin, "
                              << "variety). Only N_CS2_Score is implemented, "
                              << "as in the reference.", ValueException);
                }
                const std::string key = ll_residue_key(r.chain, r.seq_id);
                std::map<std::string, double>::const_iterator g =
                        conservation.find(key);
                if (g != conservation.end()) {
                    row.value = ll_lookup(ll_load_table(par.table), g->second);
                    row.status = "scored";
                }
            } else if (par.tag == "se") {
                const LlTable& t = ll_load_table(par.table);
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
                    IMP_THROW("ll_parameter_scores: table '" << par.table
                              << "' is not an implemented observable for tag "
                              << "'se'. Implemented: N_SE11_MEAN_SURFACE_DIST, "
                              << "N_SE10_CB_SURFACE_DIST, N_SE1/2/3_RSA_*, "
                              << "I_SE4..I_SE9_HSE*.", ValueException);
                }
                if (have) {
                    row.value = ll_lookup(t, observable);
                    row.status = "scored";
                }
            } else if (par.tag == "me") {
                if (!par.table.empty() &&
                    par.table != "N_ME11_Methionin_Exclusion_Dummy") {
                    IMP_THROW("ll_parameter_scores: the exclusion term is "
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
                    const LlResidue& q = s.residues[j];
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
                        const LlResidue& q = s.residues[j];
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
                        const LlResidue& q = s.residues[j];
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
                IMP_THROW("ll_parameter_scores: unknown tag '" << par.tag << "'",
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

std::vector<LlScore> ll_labeling_score(
        const std::vector<LlScore>& parameter_scores,
        const std::vector<LlParameter>& model, const LlOptions& options) {
    // Group by position, keeping the order positions first appear in.
    std::vector<std::string> order;
    std::map<std::string, std::map<std::string, const LlScore*> > by_pos;
    std::map<std::string, const LlScore*> exemplar;
    for (std::size_t i = 0; i < parameter_scores.size(); ++i) {
        const LlScore& r = parameter_scores[i];
        const std::string key = ll_residue_key(r.asym_id, r.seq_id);
        if (!by_pos.count(key)) { order.push_back(key); exemplar[key] = &r; }
        by_pos[key][r.score_type] = &r;
    }

    std::vector<LlScore> out;
    out.reserve(order.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
        const std::string& key = order[i];
        const LlScore& ex = *exemplar[key];
        LlScore row;
        row.asym_id = ex.asym_id;
        row.seq_id = ex.seq_id;
        row.comp_id = ex.comp_id;
        row.score_type = "combined";
        row.status = "unavailable";
        row.value = 0.0;

        bool err = false, zero = false;
        std::vector<double> repeated;
        for (std::size_t p = 0; p < model.size(); ++p) {
            const std::string st = ll_score_type(model[p].tag);
            std::map<std::string, const LlScore*>::const_iterator it =
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
                if (options.model == LL_MODEL_PUBLISHED || model[p].weight > 0)
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

std::vector<LlScore> ll_score_structure(const std::string& pdb_path,
                                        const std::vector<LlParameter>& model,
                                        const LlOptions& options,
                                        const std::string& conservation_path) {
    const LlStructure s = ll_read_structure(pdb_path);
    std::map<std::string, double> conservation;
    if (!conservation_path.empty()) {
        conservation = ll_read_consurf(conservation_path);
    }
    std::vector<LlScore> out = ll_parameter_scores(s, model, options, conservation);
    const std::vector<LlScore> combined = ll_labeling_score(out, model, options);
    out.insert(out.end(), combined.begin(), combined.end());
    return out;
}

std::map<std::string, double> ll_combined_by_key(
        const std::vector<LlScore>& scores) {
    std::map<std::string, double> out;
    for (std::size_t i = 0; i < scores.size(); ++i) {
        if (scores[i].score_type != "combined") continue;
        if (scores[i].status != "scored") continue;
        out[ll_residue_key(scores[i].asym_id, scores[i].seq_id)] =
                scores[i].value;
    }
    return out;
}

IMPBFF_END_NAMESPACE

// -------- from LabelizerFret.cpp --------
/**
 * (formerly LabelizerFret.cpp, now a section of this file)
 * \brief The Labelizer FRET pair score: which two labelling sites make the
 *        most informative FRET assay.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/AV.h>
#include <IMP/bff/AVModel.h>
#include <IMP/bff/StatesDistance.h>


#include <limits>

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

// -------- from LabelizerIO.cpp --------
/**
 * (formerly LabelizerIO.cpp, now a section of this file)
 * \brief A scored structure as one `.mmfdb.pto` container.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */


#include <IMP/bff/Pto.h>



IMPBFF_BEGIN_NAMESPACE

const char* const LL_PTO_README = "README";
const char* const LL_PTO_STRUCTURE = "structure.pdb";
const char* const LL_PTO_SCORES = "label_scores.json";
const char* const LL_PTO_PAIRS = "label_pairs.json";
const char* const LL_PTO_MODEL = "label_model.json";

namespace {

//! The README that goes in first, so the file explains itself.
/*! Not documentation *about* the format: the format telling a reader what it
    is, in the file, in language that survives this library not building. */
const char* const LL_README_TEXT =
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

std::string ll_read_file(const std::string& path) {
    std::ifstream in(path.c_str(), std::ios::binary);
    if (!in) {
        IMP_THROW("ll_write_pto: cannot read " << path, IOException);
    }
    std::ostringstream os;
    os << in.rdbuf();
    return os.str();
}

void ll_add(PtoWriter& w, const std::string& name, const std::string& kind,
            const std::string& encoding, const std::string& bytes) {
    w.add(name, kind, encoding, bytes.data(), bytes.size());
}

//! The tags of one object, flattened into the JSON its payload carries.
/*! PtoWriter frames bytes and carries a kind and an encoding; it has no tag
    element. Rather than add one -- which would change the container format,
    which this profile does not do -- the profile tags of an object ride in a
    `_tags` member of the object's own JSON. A payload that is not JSON keeps
    its tags in the model object instead. */
nlohmann::json ll_tags_json(const std::vector<MfdbTag>& tags) {
    nlohmann::json j = nlohmann::json::object();
    for (std::size_t i = 0; i < tags.size(); ++i) j[tags[i].item] = tags[i].value;
    return j;
}

nlohmann::json ll_score_row(const LlScore& s) {
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

std::vector<MfdbColumn> ll_score_columns() {
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

std::vector<MfdbColumn> ll_pair_columns() {
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

const char* ll_dye_model_name(ProbeModel m) {
    switch (m) {
        case PROBE_MODEL_CBETA: return "cbeta";
        case PROBE_MODEL_ALPHA_CONE: return "alpha_cone";
        case PROBE_MODEL_ACCESSIBLE_VOLUME: return "accessible_volume";
    }
    return "unknown";
}

std::string ll_object_text(const std::string& path, const std::string& name) {
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

void ll_write_pto(const std::string& path, const std::string& pdb_path,
                  const std::vector<LlScore>& scores,
                  const std::vector<LlPairScore>& pairs,
                  const std::string& settings_json) {
    const std::string structure = ll_read_file(pdb_path);
    const std::string structure_sum = mfdb_checksum(structure);

    // The scores, tidy: one row per (position, score_type).
    nlohmann::json score_rows = nlohmann::json::array();
    for (std::size_t i = 0; i < scores.size(); ++i) {
        score_rows.push_back(ll_score_row(scores[i]));
    }
    nlohmann::json score_doc;
    score_doc["rows"] = score_rows;

    nlohmann::json pair_rows = nlohmann::json::array();
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        const LlPairScore& p = pairs[i];
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
        r["probe_model"] = ll_dye_model_name(p.probe_model);
        pair_rows.push_back(r);
    }
    nlohmann::json pair_doc;
    pair_doc["rows"] = pair_rows;

    // The model object carries the settings, the container-level conformance
    // tags, the operation, and the provenance edges -- everything that is
    // about the file rather than about a row.
    nlohmann::json model_doc;
    model_doc["_container"] = ll_tags_json(mfdb_container_tags());
    model_doc["_operation"] = ll_tags_json(mfdb_operation_tags(
            "analysis", "labelizer", settings_json, "IMP.bff",
            get_module_version()));
    model_doc["settings"] = nlohmann::json::parse(settings_json);

    nlohmann::json edges = nlohmann::json::array();
    edges.push_back(ll_tags_json(mfdb_edge_tags(
            LL_PTO_STRUCTURE, LL_PTO_SCORES, "derived_from")));
    if (!pairs.empty()) {
        // The pair table is coarser than the score table -- one row per two
        // positions -- so the edge names the columns that join them rather
        // than leaving it to position.
        edges.push_back(ll_tags_json(mfdb_edge_tags(
                LL_PTO_SCORES, LL_PTO_PAIRS, "maps_rows_of",
                "_mmfdb_label_score.seq_id", "seq_id_1")));
    }
    model_doc["_edges"] = edges;

    const std::string score_bytes = score_doc.dump(1);
    const std::string pair_bytes = pair_doc.dump(1);

    model_doc["_artifacts"] = nlohmann::json::object();
    model_doc["_artifacts"][LL_PTO_STRUCTURE] = ll_tags_json(mfdb_artifact_tags(
            LL_PTO_STRUCTURE, "processed_data", "text", "", -1, structure_sum));
    model_doc["_artifacts"][LL_PTO_SCORES] = ll_tags_json(mfdb_artifact_tags(
            LL_PTO_SCORES, "parameter_table", "json", "label_site",
            static_cast<long>(scores.size()), mfdb_checksum(score_bytes),
            ll_score_columns()));
    if (!pairs.empty()) {
        model_doc["_artifacts"][LL_PTO_PAIRS] = ll_tags_json(mfdb_artifact_tags(
                LL_PTO_PAIRS, "analysis_result", "json", "pair",
                static_cast<long>(pairs.size()), mfdb_checksum(pair_bytes),
                ll_pair_columns()));
    }
    const std::string model_bytes = model_doc.dump(1);

    PtoWriter w(path);
    ll_add(w, LL_PTO_README, "readme", "text", LL_README_TEXT);
    ll_add(w, LL_PTO_STRUCTURE, "label.structure", "text", structure);
    ll_add(w, LL_PTO_SCORES, "label.scores", "json", score_bytes);
    if (!pairs.empty()) {
        ll_add(w, LL_PTO_PAIRS, "label.pairs", "json", pair_bytes);
    }
    ll_add(w, LL_PTO_MODEL, "label.model", "json", model_bytes);
    w.close();
}

// ---------------------------------------------------------------------------
// Reading
// ---------------------------------------------------------------------------

std::vector<LlScore> ll_read_pto_scores(const std::string& path) {
    const std::string text = ll_object_text(path, LL_PTO_SCORES);
    if (text.empty()) {
        IMP_THROW("ll_read_pto_scores: " << path << " carries no "
                  << LL_PTO_SCORES, IOException);
    }
    const nlohmann::json j = nlohmann::json::parse(text);
    std::vector<LlScore> out;
    const nlohmann::json& rows = j.at("rows");
    for (nlohmann::json::const_iterator it = rows.begin(); it != rows.end();
         ++it) {
        LlScore s;
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

std::vector<LlPairScore> ll_read_pto_pairs(const std::string& path) {
    const std::string text = ll_object_text(path, LL_PTO_PAIRS);
    std::vector<LlPairScore> out;
    if (text.empty()) return out;
    const nlohmann::json j = nlohmann::json::parse(text);
    const nlohmann::json& rows = j.at("rows");
    for (nlohmann::json::const_iterator it = rows.begin(); it != rows.end();
         ++it) {
        LlPairScore p;
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

std::string ll_extract_pto_structure(const std::string& path,
                                     const std::string& out_pdb_path) {
    PtoReader r(path);
    const int i = r.find(LL_PTO_STRUCTURE);
    if (i < 0) {
        IMP_THROW("ll_extract_pto_structure: " << path
                  << " carries no structure", IOException);
    }
    const std::vector<unsigned char> d = r.data(r.objects()[i]);
    const std::string bytes(d.begin(), d.end());
    const std::string got = mfdb_checksum(bytes);

    // What the container says it should be. Verification is explicit: opening
    // a file never hashes a payload, extracting one always does.
    const std::string model = ll_object_text(path, LL_PTO_MODEL);
    if (!model.empty()) {
        const nlohmann::json j = nlohmann::json::parse(model);
        if (j.count("_artifacts") &&
            j["_artifacts"].count(LL_PTO_STRUCTURE)) {
            const nlohmann::json& a = j["_artifacts"][LL_PTO_STRUCTURE];
            if (a.count("_mmfdb_artifact.checksum")) {
                const std::string want =
                        a["_mmfdb_artifact.checksum"].get<std::string>();
                if (want != got) {
                    IMP_THROW("ll_extract_pto_structure: " << path
                              << " does not hold the bytes it says it does: "
                              << "recorded " << want << ", found " << got,
                              IOException);
                }
            }
        }
    }
    std::ofstream out(out_pdb_path.c_str(), std::ios::binary);
    if (!out) {
        IMP_THROW("ll_extract_pto_structure: cannot write " << out_pdb_path,
                  IOException);
    }
    out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    return got;
}

std::string ll_read_pto_settings(const std::string& path) {
    const std::string model = ll_object_text(path, LL_PTO_MODEL);
    if (model.empty()) return "";
    const nlohmann::json j = nlohmann::json::parse(model);
    return j.count("settings") ? j["settings"].dump() : std::string();
}

std::string ll_settings_json(const std::vector<LlParameter>& model,
                             const LlOptions& options,
                             const LlFretOptions& fret_options,
                             const std::string& conservation_path) {
    nlohmann::json j;
    nlohmann::json terms = nlohmann::json::array();
    for (std::size_t i = 0; i < model.size(); ++i) {
        nlohmann::json t;
        t["tag"] = model[i].tag;
        t["score_type"] = ll_score_type(model[i].tag);
        t["table"] = model[i].table;
        t["weight"] = model[i].weight;
        terms.push_back(t);
    }
    j["model"] = terms;
    j["arithmetic"] =
            options.model == LL_MODEL_PUBLISHED ? "published" : "corrected";
    j["probe_radius"] = options.probe_radius;
    j["n_sphere_points"] = options.n_sphere_points;
    j["exclusion_distance"] = options.exclusion_distance;
    j["exclusion_exposure"] = options.exclusion_exposure;
    j["exclusion_residue"] = options.exclusion_residue;
    j["hse_radius"] = options.hse_radius;
    j["conservation_path"] = conservation_path;

    nlohmann::json f;
    f["probe_model"] = ll_dye_model_name(fret_options.probe_model);
    f["refine_probe_model"] = ll_dye_model_name(fret_options.refine_probe_model);
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
