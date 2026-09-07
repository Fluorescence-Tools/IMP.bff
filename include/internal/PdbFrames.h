/**
 *  \file IMP/bff/internal/PdbFrames.h
 *  \brief A PDB file as #IMP::bff::ProteinFrame values, read the way IMP reads it.
 *
 * The core's replacement for `IMP::atom::read_pdb` / `read_multimodel_pdb`
 * followed by the connection layer's `protein_frame_from_hierarchy` (PRD-137
 * step 5c). It reproduces IMP's rules, and test/io/test_pdb_frames_match_imp.py
 * checks it against them on every shipped PDB:
 *
 * - a `MODEL` record starts a new frame at the next accepted atom; `ENDMDL`
 *   does nothing; atoms before the first `MODEL` are the first frame;
 * - `NonWater` accepts an `ATOM`/`HETATM` line whose alternate-location
 *   indicator is blank or `A` and whose residue name (columns 18-20, raw) is
 *   neither `HOH` nor `DOD`; `All` accepts every line;
 * - the atom type is the four-character name field trimmed for `ATOM` (`UNK`
 *   when empty) and `"HET:"` plus the raw field for `HETATM`; the atom name is
 *   the second whitespace-separated token of `"Atom <type> of residue <n>"`,
 *   which is what the hierarchy path reads back out of the particle name;
 * - the residue name is columns 18-20 trimmed (`UNK` when empty), the chain is
 *   column 22 as a one-character string, the residue index is columns 23-26.
 *
 * Internal: the public doors are #IMP::bff::load_protein_frames and
 * #IMP::bff::load_structure.
 */
#ifndef IMPBFF_INTERNAL_PDBFRAMES_H
#define IMPBFF_INTERNAL_PDBFRAMES_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/HierarchyFrame.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_INTERNAL_NAMESPACE

enum PdbSelect { PDB_NON_WATER, PDB_ALL };

//! The frames of a PDB file. `max_frames` < 0 means all; `first_model_only` stops at the second MODEL.
IMPBFFEXPORT std::vector<ProteinFrame> read_pdb_frames(
        const std::string& path, PdbSelect select, bool first_model_only,
        int max_frames);

IMPBFF_END_INTERNAL_NAMESPACE

#endif  // IMPBFF_INTERNAL_PDBFRAMES_H
