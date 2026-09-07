/**
 * \file IMP/bff/SequenceAlignment.h
 * \brief Local sequence alignment with affine gaps.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_SEQUENCEALIGNMENT_H
#define IMPBFF_SEQUENCEALIGNMENT_H

#include <IMP/bff/bff_config.h>

#include <IMP/bff/Base.h>

#include <map>
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

    IMP_SHOWABLE_INLINE(AlignedBlock,
                        out << "AlignedBlock([" << query_start << ", "
                            << query_end << ") -> [" << template_start << ", "
                            << template_end << "))");
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

//! One stretch of a sequence, and what it is.
struct IMPBFFEXPORT SequenceSegment {
    //! `fp`, `core` or `linker`.
    std::string kind;
    //! The fluorescent protein's name for an `fp` segment; the kind otherwise.
    std::string name;
    //! One-based, inclusive.
    int start, end;
    //! Sequence identity against the template, for an `fp` segment.
    double identity;

    SequenceSegment() : start(0), end(0), identity(0) {}

    IMP_SHOWABLE_INLINE(SequenceSegment,
                        out << "SequenceSegment(" << kind << " " << name << " "
                            << start << "-" << end << ")");
};
IMP_VALUES(SequenceSegment, SequenceSegments);

//! The bundled fluorescent-protein library (`data/cgprobe/fp_library.json`).
/*! A copy of fpsim's registry, as JSON text. **fpsim is the home of this
    algorithm**; the copy exists so `imp_bff dye label-fp` works without it,
    and the drift between them is tracked in chisurf PRD-96. Do not extend it
    here. */
IMPBFFEXPORT std::string fp_library_json();

//! Which fluorescent proteins a fusion sequence contains, and where.
/*!
    Each library entry is looked for twice over: first by its motifs, which
    narrows the search to a window around the match and makes the alignment
    cheap, and otherwise by aligning the whole sequence. Hits are taken best
    identity first and a hit overlapping one already taken is dropped -- a
    residue belongs to one domain.

    \param[in] sequence the fusion's one-letter sequence
    \param[in] min_identity below this a hit is not a domain
    \return one segment per domain, `kind` `fp`, in the order taken
*/
IMPBFFEXPORT std::vector<SequenceSegment> find_fp_domains(
        const std::string& sequence, double min_identity = 0.35);

//! Per-residue pLDDT from an AlphaFold PDB's B-factor column.
/*!
    \param[in] pdb_path the structure
    \param[in] chain_id the chain; empty reads every chain
    \return one score per residue number, averaged over its atoms
*/
IMPBFFEXPORT std::map<int, double> parse_plddt_from_pdb(
        const std::string& pdb_path, const std::string& chain_id = "A");

//! Cut a sequence into rigid cores, linkers and fluorescent proteins.
/*!
    A residue is rigid when its pLDDT reaches \p rigid_threshold, and a
    detected fluorescent protein is rigid whatever its pLDDT says. A rigid run
    shorter than \p min_rigid_length becomes linker: a handful of confident
    residues between two linkers is not a body that can be moved as one.

    \param[in] sequence_length how many residues
    \param[in] plddt per-residue scores (#parse_plddt_from_pdb)
    \param[in] fp_domains what #find_fp_domains found
    \param[in] rigid_threshold the pLDDT a rigid residue reaches
    \param[in] min_rigid_length the shortest run that stays rigid
    \return the segments, in sequence order, covering every residue
*/
IMPBFFEXPORT std::vector<SequenceSegment> segments_from_plddt(
        int sequence_length, const std::map<int, double>& plddt,
        const std::vector<SequenceSegment>& fp_domains,
        double rigid_threshold = 70.0, int min_rigid_length = 12);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_SEQUENCEALIGNMENT_H
