/**
 * \file IMP/bff/SequenceAlignment.h
 * \brief Local sequence alignment with affine gaps.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_SEQUENCEALIGNMENT_H
#define IMPBFF_SEQUENCEALIGNMENT_H

#include <IMP/bff/bff_config.h>

#include <string>
#include <utility>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! One aligned block: half-open `[start, end)` in each sequence.
struct IMPBFFEXPORT AlignedBlock {
    int query_start;
    int query_end;
    int template_start;
    int template_end;

    AlignedBlock() : query_start(0), query_end(0),
                     template_start(0), template_end(0) {}
    AlignedBlock(int qs, int qe, int ts, int te)
        : query_start(qs), query_end(qe), template_start(ts), template_end(te) {}
};

//! Smith-Waterman local alignment with affine gaps, Gotoh's formulation.
/** Three matrices: one for a residue pair and one per sequence for a gap
    continuing in it, which is what makes an affine penalty exact rather than a
    per-position approximation. Replaces Biopython's `PairwiseAligner` in local
    mode with the same scoring, so imp.bff needs nothing beyond what IMP brings.

    Returns the aligned blocks in ascending order, in the spelling Biopython's
    `aligned` uses. */
IMPBFFEXPORT std::vector<AlignedBlock> smith_waterman(
        const std::string& query, const std::string& templ,
        double match = 2.0, double mismatch = -1.0,
        double gap_open = -5.0, double gap_extend = -1.0);

//! The best local-alignment score of `query` against `templ`.
IMPBFFEXPORT double smith_waterman_score(
        const std::string& query, const std::string& templ,
        double match = 2.0, double mismatch = -1.0,
        double gap_open = -5.0, double gap_extend = -1.0);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_SEQUENCEALIGNMENT_H
