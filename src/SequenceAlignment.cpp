/**
 * \file SequenceAlignment.cpp
 * \brief Local sequence alignment with affine gaps.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/SequenceAlignment.h>

#include <IMP/bff/ProbeDataPaths.h>
#include <IMP/bff/internal/Text.h>
#include <IMP/bff/internal/json.h>

#include <IMP/bff/IMPCompatibility.h>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <set>
#include <sstream>


#include <algorithm>
#include <limits>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! The filled Gotoh matrices, plus where the best local score sits.
struct Tables {
    std::vector<double> main;
    std::vector<signed char> trace;   // 0 stop, 1 diagonal, 2 up, 3 left
    double best;
    int best_i, best_j;
};

Tables fill(const std::string& q, const std::string& t,
            double match, double mismatch, double gap_open, double gap_extend) {
    const int n = (int) q.size(), m = (int) t.size();
    const int w = m + 1;
    const double neg = -std::numeric_limits<double>::infinity();

    Tables tb;
    tb.main.assign((size_t) (n + 1) * w, 0.0);
    tb.trace.assign((size_t) (n + 1) * w, 0);
    tb.best = 0.0; tb.best_i = 0; tb.best_j = 0;

    // Only the previous row of each gap matrix is ever read, so they are two
    // rows rather than two (n+1)x(m+1) planes.
    std::vector<double> gap_q_prev((size_t) w, neg), gap_q_cur((size_t) w, neg);
    std::vector<double> gap_t_cur((size_t) w, neg);

    for (int i = 1; i <= n; ++i) {
        const char qi = q[i - 1];
        gap_q_cur.assign((size_t) w, neg);
        gap_t_cur.assign((size_t) w, neg);
        for (int j = 1; j <= m; ++j) {
            const double up   = std::max(tb.main[(size_t) (i - 1) * w + j] + gap_open,
                                         gap_q_prev[j] + gap_extend);
            const double left = std::max(tb.main[(size_t) i * w + (j - 1)] + gap_open,
                                         gap_t_cur[j - 1] + gap_extend);
            gap_q_cur[j] = up;
            gap_t_cur[j] = left;

            const double diag = tb.main[(size_t) (i - 1) * w + (j - 1)]
                              + (qi == t[j - 1] ? match : mismatch);
            double cell = 0.0;
            cell = std::max(cell, diag);
            cell = std::max(cell, up);
            cell = std::max(cell, left);
            tb.main[(size_t) i * w + j] = cell;

            // precedence: stop, then diagonal, then a gap in the template,
            // then a gap in the query
            signed char step;
            if (cell == 0.0)      step = 0;
            else if (cell == diag) step = 1;
            else if (cell == up)   step = 2;
            else                   step = 3;
            tb.trace[(size_t) i * w + j] = step;

            if (cell > tb.best) { tb.best = cell; tb.best_i = i; tb.best_j = j; }
        }
        gap_q_prev.swap(gap_q_cur);
    }
    return tb;
}
}

double smith_waterman_score(const std::string& query, const std::string& templ,
                            double match, double mismatch,
                            double gap_open, double gap_extend) {
    if (query.empty() || templ.empty()) return 0.0;
    return fill(query, templ, match, mismatch, gap_open, gap_extend).best;
}

std::vector<AlignedBlock> smith_waterman(const std::string& query,
                                         const std::string& templ,
                                         double match, double mismatch,
                                         double gap_open, double gap_extend) {
    std::vector<AlignedBlock> out;
    if (query.empty() || templ.empty()) return out;

    const int w = (int) templ.size() + 1;
    const Tables tb = fill(query, templ, match, mismatch, gap_open, gap_extend);
    if (tb.best <= 0.0) return out;

    // Walk back to the first zero, collecting runs of diagonal steps.
    int i = tb.best_i, j = tb.best_j;
    int run_q_end = i, run_t_end = j;
    while (i > 0 && j > 0 && tb.trace[(size_t) i * w + j] != 0) {
        const signed char step = tb.trace[(size_t) i * w + j];
        if (step == 1) { --i; --j; continue; }
        if (run_q_end != i || run_t_end != j)
            out.push_back(AlignedBlock(i, run_q_end, j, run_t_end));
        if (step == 2) --i; else --j;
        run_q_end = i; run_t_end = j;
    }
    if (run_q_end != i || run_t_end != j)
        out.push_back(AlignedBlock(i, run_q_end, j, run_t_end));

    std::reverse(out.begin(), out.end());
    return out;
}

namespace {

//! The best local alignment window of `query` in `templ`, one-based.
/*! Returns `(start, end, identity)`; identity is matched residues over the
    aligned length, which is what a caller thresholds on. */
void align_best_window(const std::string& query, const std::string& templ,
                       int& start, int& end, double& identity) {
    start = 0;
    end = 0;
    identity = 0.0;
    const std::vector<AlignedBlock> blocks = smith_waterman(query, templ);
    if (blocks.empty()) return;

    int matched = 0, aligned = 0;
    for (std::size_t i = 0; i < blocks.size(); ++i) {
        const AlignedBlock& b = blocks[i];
        const int length = b.query_end - b.query_start;
        aligned += length;
        for (int k = 0; k < length; ++k) {
            const std::size_t q = b.query_start + k, t = b.template_start + k;
            if (q < query.size() && t < templ.size() && query[q] == templ[t]) {
                ++matched;
            }
        }
    }
    start = blocks.front().query_start + 1;
    end = blocks.back().query_end;
    identity = aligned > 0 ? static_cast<double>(matched) / aligned : 0.0;
}

}  // namespace

std::string fp_library_json() {
    const std::string path = get_probe_data_dir() + "/fp_library.json";
    std::ifstream in(path.c_str());
    if (!in) IMP_THROW("cannot read the FP library " << path, IOException);
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

std::vector<SequenceSegment> get_fp_domains(const std::string& sequence,
                                             double min_identity) {
    nlohmann::json library = nlohmann::json::parse(fp_library_json(), NULL,
                                                   false);
    std::vector<SequenceSegment> hits;
    if (library.is_discarded() || !library.is_object()) return hits;

    for (nlohmann::json::const_iterator it = library.begin();
         it != library.end(); ++it) {
        if (!it.value().is_object() || !it.value().contains("sequence")) {
            continue;
        }
        const std::string templ = it.value()["sequence"].get<std::string>();

        // A motif narrows the search to a window; without one the whole
        // sequence is aligned, which is the same answer for more work.
        std::size_t seed = std::string::npos;
        std::size_t seed_length = 0;
        if (it.value().contains("motifs") && it.value()["motifs"].is_array()) {
            const nlohmann::json& motifs = it.value()["motifs"];
            for (nlohmann::json::const_iterator m = motifs.begin();
                 m != motifs.end(); ++m) {
                if (!m->is_string()) continue;
                const std::string motif = m->get<std::string>();
                const std::size_t at = sequence.find(motif);
                if (at != std::string::npos && motif.size() > seed_length) {
                    seed = at;
                    seed_length = motif.size();
                }
            }
        }

        int offset = 0;
        std::string window = sequence;
        if (seed != std::string::npos) {
            offset = static_cast<int>(seed > 50 ? seed - 50 : 0);
            const std::size_t stop = std::min(sequence.size(), seed + 300);
            window = sequence.substr(offset, stop - offset);
        }

        int start = 0, end = 0;
        double identity = 0.0;
        align_best_window(window, templ, start, end, identity);
        if (identity < min_identity || start <= 0) continue;

        SequenceSegment hit;
        hit.kind = "fp";
        hit.name = it.key();
        hit.start = offset + start;
        hit.end = offset + end;
        hit.identity = identity;
        hits.push_back(hit);
    }

    // Best identity first, then the longer domain, then the earlier one.
    std::stable_sort(hits.begin(), hits.end(),
                     [](const SequenceSegment& a, const SequenceSegment& b) {
                         if (a.identity != b.identity) {
                             return a.identity > b.identity;
                         }
                         if ((a.end - a.start) != (b.end - b.start)) {
                             return (a.end - a.start) < (b.end - b.start);
                         }
                         return a.start < b.start;
                     });

    // A residue belongs to one domain: a hit overlapping one already taken is
    // the same protein found twice, or a worse match for the same stretch.
    std::set<int> covered;
    std::vector<SequenceSegment> out;
    for (std::size_t i = 0; i < hits.size(); ++i) {
        bool overlaps = false;
        for (int r = hits[i].start; r <= hits[i].end && !overlaps; ++r) {
            overlaps = covered.count(r) > 0;
        }
        if (overlaps) continue;
        for (int r = hits[i].start; r <= hits[i].end; ++r) covered.insert(r);
        out.push_back(hits[i]);
    }
    return out;
}

std::map<int, double> parse_plddt_from_pdb(const std::string& pdb_path,
                                           const std::string& chain_id) {
    std::ifstream in(pdb_path.c_str());
    if (!in) IMP_THROW("cannot read " << pdb_path, IOException);
    std::map<int, double> sum;
    std::map<int, int> count;
    std::string line;
    while (std::getline(in, line)) {
        if (line.compare(0, 4, "ATOM") != 0) continue;
        if (line.size() < 66) continue;
        const std::string chain = internal::trimmed(line.substr(21, 1));
        if (!chain_id.empty() && !chain.empty() && chain != chain_id) continue;
        const std::string residue = internal::trimmed(line.substr(22, 4));
        const std::string b = internal::trimmed(line.substr(60, 6));
        char* end = NULL;
        const long index = std::strtol(residue.c_str(), &end, 10);
        if (end == residue.c_str() || *end != '\0') continue;
        const double value = std::strtod(b.c_str(), &end);
        if (*end != '\0') continue;
        sum[static_cast<int>(index)] += value;
        count[static_cast<int>(index)] += 1;
    }
    std::map<int, double> out;
    for (std::map<int, double>::const_iterator it = sum.begin();
         it != sum.end(); ++it) {
        out[it->first] = it->second / count[it->first];
    }
    return out;
}

std::vector<SequenceSegment> segments_from_plddt(
        int sequence_length, const std::map<int, double>& plddt,
        const std::vector<SequenceSegment>& fp_domains, double rigid_threshold,
        int min_rigid_length) {
    std::vector<SequenceSegment> out;
    if (sequence_length <= 0) return out;

    std::vector<char> label(sequence_length, 'L');
    for (int i = 0; i < sequence_length; ++i) {
        std::map<int, double>::const_iterator score = plddt.find(i + 1);
        if (score != plddt.end() && score->second >= rigid_threshold) {
            label[i] = 'R';
        }
    }
    // A detected fluorescent protein is a body whatever its pLDDT says.
    for (std::size_t d = 0; d < fp_domains.size(); ++d) {
        for (int i = fp_domains[d].start - 1;
             i < fp_domains[d].end && i < sequence_length; ++i) {
            if (i >= 0) label[i] = 'R';
        }
    }
    // A handful of confident residues between two linkers is not a body that
    // can be moved as one.
    for (int i = 0; i < sequence_length;) {
        if (label[i] != 'R') {
            ++i;
            continue;
        }
        int j = i;
        while (j < sequence_length && label[j] == 'R') ++j;
        if (j - i < min_rigid_length) {
            for (int k = i; k < j; ++k) label[k] = 'L';
        }
        i = j;
    }

    for (int i = 0; i < sequence_length;) {
        const char current = label[i];
        int j = i;
        while (j < sequence_length && label[j] == current) ++j;

        SequenceSegment segment;
        segment.kind = current == 'R' ? "core" : "linker";
        segment.name = segment.kind;
        segment.start = i + 1;
        segment.end = j;
        // A segment that is mostly one fluorescent protein *is* that protein.
        for (std::size_t d = 0; d < fp_domains.size(); ++d) {
            const int from = std::max(i + 1, fp_domains[d].start);
            const int to = std::min(j, fp_domains[d].end);
            if (from <= to && (to - from + 1) > 0.5 * (j - i)) {
                segment.kind = "fp";
                segment.name = fp_domains[d].name;
                segment.identity = fp_domains[d].identity;
                break;
            }
        }
        out.push_back(segment);
        i = j;
    }
    return out;
}

IMPBFF_END_NAMESPACE
