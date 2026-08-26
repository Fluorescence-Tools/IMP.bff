/**
 *  \file IMP/bff/ProbeAttachment.h
 *  \brief Putting a probe on a residue: the backbone frame, the site, the move.
 *
 * A labelling site is a residue, and what a label is attached *to* is that
 * residue's backbone frame -- the origin at CA, x toward N, z out of the
 * peptide plane. Everything here follows from that one construction: a probe
 * drawn in the canonical frame is placed by transforming it into the site's,
 * and two structures are superposed by aligning one site's frame onto
 * another's.
 *
 * \note That frame existed **twice**: as `backbone_rotation` in
 * `RotamerSite.h`, which the rotamer libraries are placed with, and again in
 * the Python labelling layer, which the explicit dyes were placed with. The
 * same three lines of cross products, in two languages, for the two halves of
 * the same package. There is one now, and `backbone_rotation` is a view of it.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_PROBEATTACHMENT_H
#define IMPBFF_PROBEATTACHMENT_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/DyeLibrary.h>
#include <IMP/bff/StripMask.h>

#include <IMP/algebra/Transformation3D.h>
#include <IMP/atom/Atom.h>
#include <IMP/atom/Hierarchy.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! The backbone frame at a residue, from its three backbone atoms.
/*!
    Origin at CA, x along CA->N, z perpendicular to the peptide plane, y
    completing the right-handed set. Returned as the transformation *into* the
    frame, which is what places a probe drawn in the canonical frame onto this
    residue.

    \param[in] ca,n,c the residue's backbone positions
*/
IMPBFFEXPORT IMP::algebra::Transformation3D backbone_frame(
        const IMP::algebra::Vector3D& ca, const IMP::algebra::Vector3D& n,
        const IMP::algebra::Vector3D& c);

//! The same, looked up in a structure.
/*! \param[in] hierarchy the structure
    \param[in] chain_id the chain; empty matches any
    \param[in] residue the residue number
    \throw ValueException naming the atom when N, CA or C is missing -- a
           residue that cannot define a frame cannot carry a label, and
           guessing which atom stood in for it is how a probe ends up somewhere
           plausible and wrong */
IMPBFFEXPORT IMP::algebra::Transformation3D backbone_frame(
        IMP::atom::Hierarchy hierarchy, const std::string& chain_id,
        int residue);

//! The backbone atoms of a labelling site: CA, N, C, in that order.
/*! \throw ValueException naming the atom when one is missing */
IMPBFFEXPORT IMP::ParticlesTemp resolve_probe_site(
        IMP::atom::Hierarchy hierarchy, const std::string& chain_id,
        int residue);

//! The atoms of \p hierarchy a strip mask selects.
/*! \param[in] mask an fps `strip_mask`, e.g.
               `chain A and resid 36 and not name N+CA+C+O` */
IMPBFFEXPORT IMP::ParticlesTemp select_atoms(IMP::atom::Hierarchy hierarchy,
                                             const std::string& mask);

//! Which rows of a parallel atom table a mask *keeps*.
/*!
    For a caller holding coordinates as an array rather than a hierarchy: the
    three identifying columns come in, a keep flag per row goes out, and the
    caller slices its own array with it. The slicing is numpy's business; which
    rows to slice is the mask's.

    \param[in] chains,resseqs,names one per atom, parallel
    \param[in] mask an fps `strip_mask`
    \throw ValueException when the three columns disagree in length
*/
IMPBFFEXPORT std::vector<int> strip_keep_mask(
        const std::vector<std::string>& chains,
        const std::vector<int>& resseqs,
        const std::vector<std::string>& names, const std::string& mask);

//! Remove the atoms a mask selects, in place; returns how many went.
/*! The atoms are detached from their residue and destroyed. A caller that
    wants to keep the original clones it first (`IMP::atom::create_clone`),
    which is what the Python did behind an `inplace=False` default -- and a
    default that silently copies a whole structure is a default that hides
    what a call costs. */
IMPBFFEXPORT int strip_hierarchy(IMP::atom::Hierarchy hierarchy,
                                 const std::string& mask);

//! The keep-set at a labelling site for an explicit probe: the backbone.
/*! `N`, `CA`, `C`, `O` and a terminal `OXT`. **CB is stripped**, unlike the
    accessible-volume convention: an explicit linker is built off CA and
    replaces the whole side chain, so a CB left behind clashes with the
    linker's first atom -- as true of a spin label's linker as of a dye's.
    Each consumer owns its default; the engine takes the mask (PRD-106). */
IMPBFFEXPORT std::vector<std::string> site_keep_atom_names();

//! Strip the side chain at a labelling site, in place; returns how many went.
IMPBFFEXPORT int strip_sidechain_at_site(
        IMP::atom::Hierarchy hierarchy, const std::string& chain_id,
        int residue,
        const std::vector<std::string>& keep_atom_names =
                std::vector<std::string>());

//! Move a probe drawn in the canonical frame onto a site's frame.
/*! \param[in] probe the molecule, transformed in place -- a dye, a spin
               label, anything drawn about a backbone
    \param[in] frame the site's frame (#backbone_frame) */
IMPBFFEXPORT void place_probe(IMP::atom::Hierarchy probe,
                              const IMP::algebra::Transformation3D& frame);

//! The same, from the site's three backbone positions.
IMPBFFEXPORT void place_probe_from_coords(IMP::atom::Hierarchy probe,
                                          const IMP::algebra::Vector3D& ca,
                                          const IMP::algebra::Vector3D& n,
                                          const IMP::algebra::Vector3D& c);

//! Superpose one structure on another by a residue's backbone frame.
/*!
    \param[in] source the structure to move, transformed in place
    \param[in] chain_id,residue which of *its* residues to align
    \param[in] ca,n,c where that residue should end up
*/
IMPBFFEXPORT void align_hierarchies(IMP::atom::Hierarchy source,
                                    const std::string& chain_id, int residue,
                                    const IMP::algebra::Vector3D& ca,
                                    const IMP::algebra::Vector3D& n,
                                    const IMP::algebra::Vector3D& c);

//! The `_flr_sample_probe_details.fluorophore_type` enumeration, verbatim.
/*! `donor`, `acceptor`, `unspecified` -- a **FRET role**, which only a
    fluorophore has. */
IMPBFFEXPORT std::vector<std::string> fluorophore_types();

//! The flrCIF item each #IMP::bff::ProbePosition field is written as, as JSON.
/*! `asym_id`, `seq_id`, `comp_id` and `atom_id` are the dictionary's names
    for what the fps dialect calls `chain_identifier`, `residue_seq_number`,
    the residue name and `atom_name`. Keeping the mapping here means a writer
    and a reader cannot disagree about which item a field is. */
IMPBFFEXPORT std::string probe_position_flrcif_items();

//! The flrCIF items a #IMP::bff::Quencher maps to, as JSON.
/*! The identifiers exist in the dictionary; the photophysics does not. A null
    value marks a field that is bff-native: checked across all ten `.dic`
    files in the stack, nothing matches "quench". */
IMPBFFEXPORT std::string quencher_flrcif_items();

//! Where a probe is attached: flrCIF's `_flr_poly_probe_position`.
/*! The name follows the dictionary. flrCIF says **probe** throughout
    (`_flr_probe_list`, `_flr_poly_probe_position`,
    `_flr_sample_probe_details`) and never "label" for this; the five places
    the dictionary does write `label` mean mmCIF's `label_asym_id` /
    `label_seq_id` naming convention -- which is precisely what the fields
    below *are*, so calling the type `Label` would have named it after its own
    field prefix. */
/*!
    The *position* half of a labelling site and nothing else. How the probe's
    accessible volume is computed -- linker length, radii, grid resolution --
    is a representation parameter and belongs to whichever representation is
    used; conflating the two is why an `AV` ended up doubling as "the dye"
    (PRD-113).

    **A label is not a dye.** A spin label at a cysteine, a quencher
    tryptophan and a fluorophore are, structurally, the same object: a
    residue, an attachment atom, a linker and a rotamer library. Nothing here
    distinguishes them and nothing here should -- the difference is in what is
    *measured* from them, which is spectra for FRET (`dye`) and a dipole for
    DEER, and lives with the measurement rather than with the site.
*/
struct IMPBFFEXPORT ProbePosition {
    //! The chain, the residue number, and the attachment atom.
    std::string asym_id;
    int seq_id;
    std::string atom_id;
    //! The residue name, when known.
    std::string comp_id;
    //! What is attached, by name: `AlexaFluor488`, `MTSSL`, `TRP`.
    std::string probe;
    //! The photophysics, when the probe is a fluorophore this package knows --
    //! the value, not its name: a label carrying only a name cannot derive an
    //! R0 or a correlation time without going back to the library for what it
    //! already had. Empty `name` means "not a fluorophore, or not known", and
    //! a site with an unknown probe is still a usable site.
    Dye dye;
    //! The FRET role: `donor`, `acceptor` or `unspecified`. A probe that is
    //! not a fluorophore has none, and asking for one is an error rather than
    //! a value nobody reads.
    std::string fluorophore_type;
    //! The author's name for the position, when the file gives one.
    std::string auth_name;
    //! The residue was mutated, or chemically modified, to carry the label.
    bool mutation_flag, modification_flag;
    std::string description;

    ProbePosition()
        : seq_id(0), atom_id("CB"), fluorophore_type("unspecified"),
          mutation_flag(false), modification_flag(false) {}
    /*! \param[in] asym_id,seq_id,atom_id where the probe is attached
        \param[in] probe what is attached, by name
        \param[in] fluorophore_type the FRET role, when there is one
        \throw ValueException for a role outside the enumeration */
    ProbePosition(const std::string& asym_id, int seq_id,
          const std::string& atom_id = "CB", const std::string& probe = "",
          const std::string& fluorophore_type = "unspecified");

    //! Give the label the photophysics of a fluorophore.
    /*! Sets `dye` and `probe` together, so the two cannot disagree about what
        is attached. */
    void set_dye(const Dye& dye);

    //! Two labels are the same when they name the same position of the same
    //! probe: what a round trip through a file has to preserve.
    bool operator==(const ProbePosition& other) const;
    bool operator!=(const ProbePosition& other) const { return !(*this == other); }

    //! `(asym_id, seq_id, atom_id)` -- what identifies the position.
    std::string get_key() const;
    //! The position half of an fps-dialect dict, as JSON.
    std::string get_source_info() const;

    IMP_SHOWABLE_INLINE(ProbePosition,
                        out << "ProbePosition(" << asym_id << seq_id << "." << atom_id
                            << ", " << fluorophore_type << ")");
};
IMP_VALUES(ProbePosition, ProbePositions);

//! A label from the fps-dialect position dict the AV builder takes.
/*! Reads the *position* fields and ignores the AV parameters in the same
    dict, which is the split #IMP::bff::ProbePosition exists to make. */
IMPBFFEXPORT ProbePosition probe_position_from_source_info(const std::string& source_info_json,
                                          const Dye& dye = Dye());

//! One probe and the site it goes on.
/*! \note The hierarchy and the transformation are held privately and handed
    out **by value**. They are IMP value types, and a public member of one
    would have SWIG return a non-const reference to it, which IMP's wrapper
    layer refuses -- for the good reason that a caller mutating it would be
    editing this object's insides through a borrowed handle. */
class IMPBFFEXPORT ProbeAttachment {
    IMP::atom::Hierarchy probe_;
    std::string chain_;
    int residue_;
    IMP::algebra::Transformation3D frame_;
    int n_stripped_;

public:
    ProbeAttachment() : residue_(0), n_stripped_(0) {}
    //! \param[in] probe the molecule to attach -- a dye, a spin label, any
    //! residue-attached group
    ProbeAttachment(IMP::atom::Hierarchy probe, const std::string& chain,
                    int residue)
        : probe_(probe), chain_(chain), residue_(residue), n_stripped_(0) {}

    IMP::atom::Hierarchy get_probe() const { return probe_; }
    std::string get_chain() const { return chain_; }
    int get_residue() const { return residue_; }
    //! Where the site's backbone frame turned out to be; #attach_probes fills
    //! it, so a caller can see where each probe was put.
    IMP::algebra::Transformation3D get_frame() const { return frame_; }
    void set_frame(const IMP::algebra::Transformation3D& f) { frame_ = f; }
    //! How many side-chain atoms were removed for it.
    int get_n_stripped() const { return n_stripped_; }
    void set_n_stripped(int n) { n_stripped_ = n; }

    IMP_SHOWABLE_INLINE(ProbeAttachment,
                        out << "ProbeAttachment(" << chain_ << residue_ << ")");
};
IMP_VALUES(ProbeAttachment, ProbeAttachments);

//! Attach probes to a structure at the sites named, in place.
/*!
    \param[in] protein the labelled structure
    \param[in] attachments the probes and their sites
    \param[in] strip_site_sidechain remove each site's side chain first
    \return the attachments, each with the frame its probe was placed on and
            how many atoms were stripped for it
    \throw ValueException when a site has no backbone frame
*/
IMPBFFEXPORT std::vector<ProbeAttachment> attach_probes(
        IMP::atom::Hierarchy protein,
        const std::vector<ProbeAttachment>& attachments,
        bool strip_site_sidechain = false);

//! The backbone-dependent C-beta at a labelling residue.
/*!
    #IMP::rotamer::RotamerCalculator queries a Dunbrack library for the
    side-chain geometry that the residue's actual phi/psi imply, and this asks
    it for one atom: the C-beta an attachment is anchored at. For a
    crystallographic structure the difference from the deposited C-beta is
    small; for an NMR or modelled one, where C-beta may be absent or badly
    placed, it is a physically consistent anchor rather than a coordinate.

    This was a `%pythoncode` def behind a lazy door, because `IMP.rotamer` was
    not one of this module's modules. It is one now.

    \param[in] hierarchy the structure
    \param[in] chain_id,resnum which residue
    \param[in] rl_path a Dunbrack rotamer library file; empty queries an empty
               library, which answers with the structure's own C-beta
    \param[in] prob_threshold cumulative probability cut passed to the
               calculator
    \return the position
    \throw ValueException when the residue is not in the structure, or has no
           C-beta at all (glycine, or an incomplete model) -- there is no
           anchor to report and no coordinate that would stand for one
*/
IMPBFFEXPORT IMP::algebra::Vector3D get_anchor_cb_position(
        IMP::atom::Hierarchy hierarchy, std::string chain_id, int resnum,
        std::string rl_path = "", double prob_threshold = 0.01);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_PROBEATTACHMENT_H
