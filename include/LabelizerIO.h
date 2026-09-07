/**
 *  \file IMP/bff/LabelizerIO.h
 *  \brief A scored structure as one `.mmfdb.pto` container: the structure it
 *         was computed from, the scores, and how they were produced.
 *
 * The reference writes six CSVs, a JSON heat map, four PDBs whose B-factor
 * column carries a dimensionless score, and a zip of the lot
 * (`labelizer.py:423`). Those files cannot say what a row is, cannot say what
 * a column is in, and cannot say what produced them — and one of them, the
 * conservation PDB, is written to the same path it was read from, which is how
 * the shipped 1DDB example came to feed its own output back in as its input.
 *
 * This writes **one file**. It holds the structure verbatim, the scores as
 * tables whose columns are named by `_mmfdb_label_score` items, and the
 * complete settings; every controlled value is a term from the MMFDB
 * dictionary, checked by `Pto.h`'s profile section before it is written.
 *
 * The container is #IMP::bff::PtoWriter's and the vocabulary is
 * `Pto.h`'s. Nothing here invents either.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_LABELIZERIO_H
#define IMPBFF_LABELIZERIO_H

#include <IMP/bff/bff_config.h>

#include <IMP/bff/LabelizerFret.h>
#include <IMP/bff/LabelizerScore.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! The object names a label container uses, so a reader can ask for one.
/*! Every kind this package writes is under the `label.` prefix. `drot.` and
    `rot.bbdep.` belong to the rotamer containers; the namespaces are disjoint
    by agreement so one walker can read both. */
IMPBFFEXPORT extern const char* const LL_PTO_README;
//! The input structure, byte for byte as it was read.
IMPBFFEXPORT extern const char* const LL_PTO_STRUCTURE;
//! The per-position scores, one row per (position, score_type).
IMPBFFEXPORT extern const char* const LL_PTO_SCORES;
//! The pair scores, one row per pair.
IMPBFFEXPORT extern const char* const LL_PTO_PAIRS;
//! The complete settings the run used.
IMPBFFEXPORT extern const char* const LL_PTO_MODEL;

//! Write a scored structure as one PTO.MFDB container.
/*!
    The structure goes in **first and verbatim** and is never rewritten, so it
    can be recovered byte for byte and checked against its recorded SHA-256.
    A README object precedes it explaining, in the file, how to walk the
    framing and get the structure back out — a container outlives the software
    that wrote it, and the person who needs it most is the one for whom this
    library will not build.

    Absence is preserved: a position whose score was not computed is written
    with its `status` and **no value**, never a sentinel.

    \param[in] path the container to write; `.mmfdb.pto` by convention, though
               conformance is stated by the tags inside and not by the name
    \param[in] pdb_path the structure that was scored, embedded verbatim
    \param[in] scores rows from #IMP::bff::ll_score_structure
    \param[in] pairs rows from #IMP::bff::ll_pair_scores, or empty
    \param[in] settings_json the complete settings, as JSON text; a partial
               record is worse than none because it looks reproducible
    \throw IOException when the structure cannot be read or the container
           cannot be written
    \throw ValueException when a value is not a dictionary term
*/
IMPBFFEXPORT void ll_write_pto(const std::string& path,
                               const std::string& pdb_path,
                               const std::vector<LlScore>& scores,
                               const std::vector<LlPairScore>& pairs,
                               const std::string& settings_json);

//! Read the per-position scores back out of a container.
/*! \param[in] path a container written by #ll_write_pto
    \return the rows, in the order they were written
    \throw IOException when the file is not a PTO document or holds no scores */
IMPBFFEXPORT std::vector<LlScore> ll_read_pto_scores(const std::string& path);

//! Read the pair scores back out of a container.
/*! \return the rows, or empty when the container carries none */
IMPBFFEXPORT std::vector<LlPairScore> ll_read_pto_pairs(
        const std::string& path);

//! Recover the embedded structure, verifying it against its checksum.
/*!
    \param[in] path the container
    \param[in] out_pdb_path where to write the recovered structure
    \return the SHA-256 that was recorded and matched
    \throw IOException when the container has no structure, or when the bytes
           do not hash to what the container says they should
*/
IMPBFFEXPORT std::string ll_extract_pto_structure(
        const std::string& path, const std::string& out_pdb_path);

//! The settings JSON the run recorded.
IMPBFFEXPORT std::string ll_read_pto_settings(const std::string& path);

//! The settings of a scoring run, as the JSON the container stores.
/*! Every field of both options structs, so the run is reproducible from the
    file alone. */
IMPBFFEXPORT std::string ll_settings_json(
        const std::vector<LlParameter>& model, const LlOptions& options,
        const LlFretOptions& fret_options, const std::string& conservation_path);

IMPBFF_END_NAMESPACE

#endif /* IMPBFF_LABELIZERIO_H */
