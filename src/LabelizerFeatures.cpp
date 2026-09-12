/**
 * \file LabelizerFeatures.cpp
 * \brief The structural quantities a label-site score is computed from:
 *        secondary structure, half-sphere exposure, relative solvent
 *        accessibility, residue depth, and an imported conservation grade.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/LabelizerFeatures.h>
#include <IMP/bff/ProbeAccessibleVolumeBuilder.h>
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

const std::vector<std::string>& labelizer_standard_residues() {
    static const std::vector<std::string> names = [] {
        const char* n[20] = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                             "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                             "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
        return std::vector<std::string>(n, n + 20);
    }();
    return names;
}

std::string labelizer_residue_key(const std::string& chain, int seq_id) {
    std::ostringstream os;
    os << chain << seq_id;
    return os.str();
}

LabelizerStructure labelizer_read_structure(const std::string& pdb_path, bool protein_only,
                              int model) {
    const std::vector<PDBAtomRecord> records = read_pdb_records(pdb_path);
    if (records.empty()) {
        IMP_THROW("labelizer_read_structure: no coordinates in " << pdb_path,
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
        IMP_THROW("labelizer_read_structure: " << pdb_path << " has " << n_models
                  << " model(s); model " << model << " was asked for",
                  ValueException);
    }
    const std::size_t first = starts[model], last = starts[model + 1];

    std::set<std::string> standard(labelizer_standard_residues().begin(),
                                   labelizer_standard_residues().end());

    LabelizerStructure s;
    std::string cur_chain;
    int cur_seq = 0;
    bool have_residue = false;
    for (std::size_t i = first; i < last; ++i) {
        const PDBAtomRecord& r = records[i];
        const std::string res = lf_strip(r.res_name);
        if (protein_only && standard.find(res) == standard.end()) continue;

        if (!have_residue || r.chain != cur_chain || r.resseq != cur_seq) {
            LabelizerResidue nr;
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

        LabelizerResidue& cr = s.residues.back();
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
bool lf_cbeta(const LabelizerStructure& s, int ri, LfVec& out) {
    const LabelizerResidue& r = s.residues[ri];
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

void labelizer_cbeta_position(const LabelizerStructure& s, int residue, double** out_view,
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

void labelizer_half_sphere_exposure(const LabelizerStructure& s, double radius,
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
        const LabelizerResidue& r = s.residues[i];
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

const double LABELIZER_DSSP_HBOND_ENERGY = -0.5;

std::string labelizer_dssp(const LabelizerStructure& s) {
    const int nr = static_cast<int>(s.residues.size());
    std::string ss(static_cast<std::size_t>(nr), '-');
    if (nr == 0) return ss;

    // Backbone, and the amide hydrogen DSSP places rather than reads.
    std::vector<LfVec> N(nr), CA(nr), C(nr), O(nr), H(nr);
    std::vector<char> ok(nr, 0), has_h(nr, 0);
    for (int i = 0; i < nr; ++i) {
        const LabelizerResidue& r = s.residues[i];
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
            if (e < LABELIZER_DSSP_HBOND_ENERGY) hb[i][j] = 1;
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

std::map<std::string, double> labelizer_max_asa(LabelizerMaxAsa scale) {
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
    const double* v = scale == LABELIZER_MAXASA_SANDER ? sander
                    : scale == LABELIZER_MAXASA_MILLER ? miller : wilke;
    std::map<std::string, double> m;
    for (int i = 0; i < 20; ++i) m[res[i]] = v[i];
    return m;
}

void labelizer_residue_sasa(const LabelizerStructure& s, double probe_radius,
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

void labelizer_relative_solvent_accessibility(const LabelizerStructure& s, LabelizerMaxAsa scale,
                                       double probe_radius, int n_sphere_points,
                                       double** out_view, int* n_out_view) {
    const std::size_t nr = s.residues.size();
    double* out = internal::new_double_view(nr, out_view, n_out_view);
    if (!out || nr == 0) return;

    double* area = 0;
    int n_area = 0;
    labelizer_residue_sasa(s, probe_radius, n_sphere_points, &area, &n_area);
    const std::map<std::string, double> maxasa = labelizer_max_asa(scale);
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
std::vector<double> lf_build_surface(const LabelizerStructure& s, double probe_radius,
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

void labelizer_cbeta_depth(const LabelizerStructure& s, double probe_radius,
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

void labelizer_residue_depth(const LabelizerStructure& s, double probe_radius,
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

void labelizer_global_charge(const LabelizerStructure& s, double** out_view,
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

std::map<std::string, double> labelizer_read_consurf(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) {
        IMP_THROW("labelizer_read_consurf: cannot read " << path, IOException);
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
            const std::string key = labelizer_residue_key(chain, seq);
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
        grades[labelizer_residue_key(chain, std::atoi(seq.c_str()))] =
                std::atof(f[3].c_str());
    }
    return grades;
}

IMPBFF_END_NAMESPACE
