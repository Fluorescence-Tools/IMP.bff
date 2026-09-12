/**
 *  \file IMP/bff/StripMask.h
 *  \brief The fps `strip_mask` dialect: which atoms a selection removes.
 *
 * fps.json positions are calibrated for a structure whose attachment residue
 * does not wall in its own probe, so the side chain goes -- minus the
 * attachment atom, and the backbone stays -- before anything is measured. That
 * convention is what every `allowed_sphere_radius` in a shipped fps.json is
 * calibrated against, and it is spelled as a selection so that a declared mask
 * and a consumer's default go through one parser.
 *
 * The grammar is #IMP::bff::SelectionExpression, whole: `chain A and resid 132
 * and not name CA+CB+C+N+O`, `chain A and resid 115 and not name CA CB C N O`,
 * `(resid 132 and not name CA+CB+C+N+O) or resname HOH SOL WAT`, and anything
 * else the language accepts. A mask this cannot read is **refused loudly**,
 * because the alternative is computing against obstacles the document said to
 * remove.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_STRIPMASK_H
#define IMPBFF_STRIPMASK_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/SelectionExpression.h>


#include <IMP/bff/IMPCompatibility.h>

#include <set>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! The backbone a labelling site keeps when its side chain is stripped.
IMPBFFEXPORT std::vector<std::string> backbone_atom_names();

//! The obstacles a volume sees: `(x, y, z, r)` per atom, flat.
/*! The array form of a strip, for the kernels that take obstacles as numbers
    rather than as a hierarchy -- #IMP::bff::get_av and everything built on
    it. Four doubles per atom, in hierarchy order, **every** atom present: the
    ones \p mask names carry a **radius of zero**, which is what a strip is
    everywhere in this module. A rasteriser asks `distance < radius`, so an
    atom of no size blocks nothing, and the array keeps row-for-atom
    correspondence with the structure it came from.

    \param[in] hierarchy the structure
    \param[in] mask a selection expression; empty keeps every atom
    \param[in] keep an atom to keep even when the mask names it, as
               `chain/resi/name` -- the attachment atom, whose removal would
               let a volume grow through its own anchor. Empty keeps nothing
               extra.
    \throw ValueException when the mask cannot be read */

//! What a mask removes from a structure, atom by atom.
/*! The diagnostic a caller needs to answer "did my mask do what I meant":
    how many atoms it names, which residues they belong to, and -- the case
    worth catching -- whether it names none at all. */
struct IMPBFFEXPORT StripReport {
    //! The mask, as given.
    std::string mask;
    //! How many atoms it selects.
    int n_selected;
    //! How many atoms the structure has.
    int n_atoms;
    //! `chain/resi/name` of each selected atom, in structure order.
    std::vector<std::string> atoms;
    //! `chain/resi` of each residue it touches, without repeats.
    std::vector<std::string> residues;

    StripReport() : n_selected(0), n_atoms(0) {}

    IMP_SHOWABLE_INLINE(StripReport,
                        out << "StripReport(" << n_selected << "/" << n_atoms
                            << " atoms, " << residues.size() << " residues)");
};
IMP_VALUES(StripReport, StripReports);

//! What \p mask removes from \p hierarchy.
/*! \throw ValueException when the mask cannot be read. A mask that reads and
           selects nothing is **not** an error -- it is reported as
           `n_selected == 0`, which is the answer the caller asked for. */

//! Parse an fps `strip_mask`.
/*! \throw ValueException on a syntax error or an unanswerable keyword; see
           #IMP::bff::SelectionExpression. */
IMPBFFEXPORT SelectionExpression parse_strip_mask(const std::string& mask);

//! Mask that strips residue `(chain, resseq)` except \p keep_atom_names.
/*! An empty \p chain matches any chain. Spelled in the fps dialect so a
    declared mask and a consumer default round-trip through the same parser. */
IMPBFFEXPORT std::string site_strip_mask(
        const std::string& chain, int resseq,
        const std::vector<std::string>& keep_atom_names);

//! The AV default strip for an attachment site, as an fps `strip_mask`.
/*! The attachment residue's side chain minus the attachment atom; the backbone
    stays. E.g. `chain A and resid 36 and not name N+CA+C+O+CB`. */
IMPBFFEXPORT std::string default_strip_mask(const std::string& chain,
                                            int resseq,
                                            const std::string& atom_name);

//! The PDB lines with the atoms selected by \p mask removed.
/*!
    \param[in] lines the file's lines, verbatim
    \param[in] mask an fps `strip_mask`
    \param[in] keep_chain,keep_resseq,keep_atom_name the attachment atom, kept
               whatever the mask selects; an empty chain matches any, and an
               empty atom name protects nothing
*/
IMPBFFEXPORT std::vector<std::string> strip_pdb_lines(
        const std::vector<std::string>& lines, const std::string& mask,
        const std::string& keep_chain = "", int keep_resseq = 0,
        const std::string& keep_atom_name = "");

//! A copy of \p pdb_path with the atoms the mask selects removed.
/*!
    An empty \p strip_mask means the default strip. The attachment atom is
    always kept, whatever the mask selects, because the attachment is resolved
    by `(chain, residue, atom name)` from this file.

    **Cached on (path, mtime, size, site, mask).** One temporary file per
    distinct combination per process; a stripped copy is read once per AV and
    an fps.json builds dozens.

    \return the stripped copy's path, or \p pdb_path itself on an I/O failure —
            a cloud computed against the unstripped structure beats no cloud
*/
IMPBFFEXPORT std::string stripped_pdb_for(const std::string& pdb_path,
                                          const std::string& chain, int resseq,
                                          const std::string& atom_name,
                                          const std::string& strip_mask = "");

IMPBFF_END_NAMESPACE

#endif //IMPBFF_STRIPMASK_H
