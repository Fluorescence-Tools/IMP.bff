/**
 *  \file IMP/bff/StripMask.h
 *  \brief The fps `strip_mask` dialect: which atoms a selection removes.
 *
 * fps.json positions are calibrated for a structure whose attachment residue
 * does not wall in its own dye, so the side chain goes — minus the attachment
 * atom, and the backbone stays — before anything is measured. That convention
 * is what every `allowed_sphere_radius` in a shipped fps.json is calibrated
 * against, and it is spelled as a PyMOL-style selection so a declared mask and
 * a consumer's default round-trip through one parser.
 *
 * The grammar is the dialect fps documents actually carry, and nothing more:
 * an `and`-chain of `chain <id>`, `resid <n>` (`resi` accepted), and one name
 * term, `name A+B+...` or `not name A+B+...`. `+` is the list separator —
 * PyMOL's own; a space-separated list is a parse error there and here.
 * Anything else (`or`, parentheses, other keywords) is **refused loudly**: a
 * mask this cannot read must not be silently ignored or approximated, because
 * the alternative is computing against obstacles the document said to remove.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_STRIPMASK_H
#define IMPBFF_STRIPMASK_H

#include <IMP/bff/bff_config.h>

#include <IMP/value_macros.h>
#include <IMP/showable_macros.h>

#include <set>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! The backbone a labelling site keeps when its side chain is stripped.
IMPBFFEXPORT std::vector<std::string> backbone_atom_names();

//! A parsed `strip_mask`: which (chain, resseq, atom name) it selects.
class IMPBFFEXPORT StripSelection {
    std::string chain_;               //!< empty = any chain
    std::vector<int> resids_;         //!< empty = any residue
    std::set<std::string> names_;     //!< empty and !has_names_ = no name term
    bool has_names_;
    bool negate_;                     //!< the name term was `not name ...`

public:
    StripSelection() : has_names_(false), negate_(false) {}
    StripSelection(const std::string& chain, const std::vector<int>& resids,
                   const std::vector<std::string>& names, bool has_names,
                   bool negate);

    std::string get_chain() const { return chain_; }
    std::vector<int> get_resids() const { return resids_; }
    std::vector<std::string> get_names() const;
    bool get_has_names() const { return has_names_; }
    bool get_negate() const { return negate_; }

    //! True when an atom is selected by the mask.
    bool matches(const std::string& chain, int resseq,
                 const std::string& name) const;

    IMP_SHOWABLE_INLINE(StripSelection,
                        out << "StripSelection(chain=" << chain_ << ")");
};
IMP_VALUES(StripSelection, StripSelections);

//! Parse an fps `strip_mask`.
/*! \throw ValueException for anything outside the dialect */
IMPBFFEXPORT StripSelection parse_strip_mask(const std::string& mask);

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
