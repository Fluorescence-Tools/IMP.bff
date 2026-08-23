/**
 *  \file IMP/bff/RotamerSite.h
 *  \brief Placing a rotamer library on a residue: the backbone frame, the
 *         atom selectors, and the library registry.
 *
 * The transforms and the registry resolution were Python in
 * `representation/rotamer.py`. They are kernels and file logic, not glue:
 * the backbone frame is nine multiplications a caller should not re-type,
 * the selector matching has the same ambiguity rules as the scoring one, and
 * the registry is a JSON lookup with a cutoff-suffix grammar.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_ROTAMERSITE_H
#define IMPBFF_ROTAMERSITE_H

#include <IMP/bff/bff_config.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

// --------------------------------------------------------------------------
// The backbone frame
// --------------------------------------------------------------------------

//! CA, N, C coordinates of `(chain, residue)` in a protein frame.
/*!
    The first atom of each name that matches: same chain when both name one
    (case-blind), and the residue number when the frame carries one. A missing
    backbone atom raises -- a frame that cannot place a library is not a frame
    a caller can guess around.

    \param[in] coords flat, three per atom
    \param[in] atom_names,chain_ids one per atom; chain_ids may be empty
    \param[in] residue_indices one per atom, -1 where unknown; may be empty
    \param[in] chain the chain asked for; empty matches any
    \param[in] residue the residue number asked for
    \param[out] out_view,n_out_view nine values: CA, then N, then C
    \throw ValueException when CA, N or C is missing
*/
IMPBFFEXPORT void resolve_backbone_site(
        const std::vector<double>& coords,
        const std::vector<std::string>& atom_names,
        const std::vector<std::string>& chain_ids,
        const std::vector<int>& residue_indices,
        const std::string& chain, int residue,
        double** out_view, int* n_out_view);

//! The site frame's rotation: rows are x (CA->N), y (in the N-CA-C plane),
//! z = x cross y.
/*! \param[out] out_view,n_out_view nine values, row-major */
IMPBFFEXPORT void backbone_rotation(
        const std::vector<double>& ca, const std::vector<double>& n,
        const std::vector<double>& c,
        double** out_view, int* n_out_view);

//! Library coordinates into the backbone frame at CA.
/*!
    \param[in] coords flat, `n_rotamers * n_atoms * 3`, the library's own frame
    \param[in] ca,n,c the site's backbone atoms, three values each
    \param[out] out_view,n_out_view flat, the input's shape, in the protein's
               frame: `rotated + ca`
*/
IMPBFFEXPORT void transform_library_to_site(
        const std::vector<double>& coords,
        const std::vector<double>& ca, const std::vector<double>& n,
        const std::vector<double>& c,
        double** out_view, int* n_out_view);

// --------------------------------------------------------------------------
// Atom selectors
// --------------------------------------------------------------------------

//! Indices of the atoms named by FRETpredict selectors.
/*!
    A selector is `NAME` or `NAME and resname RES`. The resname clause is
    **honoured** when \p resnames is given: atom names repeat between the dye
    residue and its linker -- `C13` is in both `A48` and `C1R` -- and a
    selector that ignores the residue takes whichever comes first in the atom
    ordering, which is the dye today only by luck. Without \p resnames the
    clause cannot be checked and the first name match is returned.

    \param[in] atom_names one per atom, uppercased internally
    \param[in] selectors one entry per wanted atom
    \param[in] resnames one per atom, or empty to skip the clause check
    \return one index per selector entry, in order
    \throw ValueException naming the selector when one matches no atom
*/
IMPBFFEXPORT std::vector<int> selector_atom_indices(
        const std::vector<std::string>& atom_names,
        const std::vector<std::string>& selectors,
        const std::vector<std::string>& resnames = std::vector<std::string>());

// `selector_resnames(selectors)` lives in Scoring.h beside the rest of the
// selector mini-language.

// --------------------------------------------------------------------------
// The bundled library registry
// --------------------------------------------------------------------------

//! The registry key of a library name: the name without its cutoff suffix.
/*! `'AlexaFluor 488 C1R cutoff30'` -> `'AlexaFluor 488 C1R'`. */
IMPBFFEXPORT std::string normalize_library_name(const std::string& name);

//! The cutoff of a library name, or -1 when the name carries none.
IMPBFFEXPORT int library_name_cutoff(const std::string& name);

//! One registry entry of `libraries.json`, plus the resolved spelling.
/*! The returned JSON object is the registry's own entry with `name` (the
    key), `library_name` (the name as asked) and `cutoff` added -- the shape
    the Python `rotamer_library_metadata` returned. The registry is read once
    and cached. An unknown library raises. */
IMPBFFEXPORT std::string rotamer_library_metadata(const std::string& name);

//! The library file stem the metadata and cutoff select.
/*! The registry's `filename` carries FRETpredict's default cutoff; a name
    that asks for another cutoff replaces it. */
IMPBFFEXPORT std::string library_filename(const std::string& metadata_json,
                                          int cutoff);

//! Resolve a library name (or an explicit path) to a library file.
/*!
    The FRETpredict library files (module data, `data/rotamer_library`) are
    canonical: `<stem>.pdb` + `<stem>_cutoff<N>.bcif` (+ weights) per cutoff.
    They are tried first so the *requested cutoff* is the one loaded -- the
    RMF templates hold only the cutoff-30 clustering, so resolving every name
    to `<stem>.rmf3` silently returned the wrong library for cutoff10/20.

    \param[in] name a registry name, or a path that exists
    \param[in] lib_dir an extra directory to search after the canonical one;
               empty searches the module's template directory only
    \return the existing file path
    \throw IOException when nothing matches, or when only the cutoff-30 RMF
           template exists for another cutoff (that mismatch once answered
           quietly)
*/
IMPBFFEXPORT std::string resolve_rotamer_library_path(
        const std::string& name, const std::string& lib_dir = "");

IMPBFF_END_NAMESPACE

#endif //IMPBFF_ROTAMERSITE_H
