/**
 * \file SequenceAlignment.cpp
 * \brief Local sequence alignment with affine gaps.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/SequenceAlignment.h>

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

            // the same precedence the Python had: stop, then diagonal, then a
            // gap in the template, then a gap in the query
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

IMPBFF_END_NAMESPACE
