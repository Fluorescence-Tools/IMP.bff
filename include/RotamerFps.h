/**
 *  \file IMP/bff/RotamerFps.h
 *  \brief fps.json and rotamer ensembles: reading a labelled pair out of a
 *         file, and writing predicted distances back into one.
 *
 * The fps.json format names the same thing several ways -- a position's chain
 * is `chain_identifier`, `chain` or `segid`, a distance's donor is
 * `position1_name`, `donor_position`, `donor_position_name` or
 * `dye1_position` -- because it has been written by several programs. The
 * aliases are the format's, so reading them is this module's job and not
 * every caller's.
 *
 * Entries cross the language boundary as JSON text, which is what the rest of
 * the fps layer does (#IMP::bff::read_fps_json, #IMP::bff::write_fps_json):
 * an entry carries whatever keys its writer put there, and a typed struct
 * would either lose them or have to grow a field per program. What *is* typed
 * is the part this module reasons about -- which chain, which residue, which
 * library, which dye.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_ROTAMERFPS_H
#define IMPBFF_ROTAMERFPS_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/RotamerEnsemble.h>

#include <map>
#include <string>
#include <utility>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! A labelling position of an fps.json file.
struct IMPBFFEXPORT RotamerPosition {
    std::string name;
    //! Empty when the file names none, which matches any chain.
    std::string chain;
    int residue;
    //! The anchor atom; `CA` unless the file says otherwise.
    std::string atom_name;
    //! The dye, and the rotamer library; either may be empty.
    std::string dye, library;
    //! `donor`, `acceptor`, ... as the file spells it; may be empty.
    std::string role;

    RotamerPosition() : residue(0), atom_name("CA") {}

    IMP_SHOWABLE_INLINE(RotamerPosition,
                        out << "RotamerPosition(" << name << ", " << chain
                            << residue << ")");
};
IMP_VALUES(RotamerPosition, RotamerPositions);

//! A FRET distance between two labelling positions of an fps.json file.
struct IMPBFFEXPORT RotamerDistance {
    std::string name;
    std::string donor_position, acceptor_position;
    //! The dyes and their libraries, taken from the distance entry when it
    //! names them and from the two positions when it does not.
    std::string donor, acceptor, libname_1, libname_2;

    RotamerDistance() {}

    IMP_SHOWABLE_INLINE(RotamerDistance,
                        out << "RotamerDistance(" << name << ": "
                            << donor_position << " -> " << acceptor_position
                            << ")");
};
IMP_VALUES(RotamerDistance, RotamerDistances);

//! One distance of an fps.json file, with both its positions and the document.
struct IMPBFFEXPORT RotamerFpsSelection {
    RotamerPosition donor, acceptor;
    RotamerDistance distance;
    //! The whole document, as JSON text: every position, every distance, and
    //! the score sets and extra sections merged.
    std::string positions, distances, extra;

    RotamerFpsSelection() : positions("{}"), distances("{}"), extra("{}") {}

    IMP_SHOWABLE_INLINE(RotamerFpsSelection,
                        out << "RotamerFpsSelection(" << distance.name << ")");
};
IMP_VALUES(RotamerFpsSelection, RotamerFpsSelections);

//! Parse one fps.json position entry.
/*! \param[in] name the entry's key
    \param[in] payload_json the entry */
IMPBFFEXPORT RotamerPosition rotamer_position_from_payload(
        const std::string& name, const std::string& payload_json);

//! Parse one fps.json distance entry, filling from its positions.
/*! A distance that names no dye or library takes the ones its positions name,
    which is how an fps file written position-first still describes a pair.

    \param[in] name the entry's key
    \param[in] payload_json the entry
    \param[in] positions_json every position of the document */
IMPBFFEXPORT RotamerDistance rotamer_distance_from_payload(
        const std::string& name, const std::string& payload_json,
        const std::string& positions_json);

//! Read one FRET distance and its two positions out of an fps.json file.
/*!
    \param[in] path the fps.json file
    \param[in] distance_name which distance; empty takes the first
    \throw ValueException when the file has no distances, when \p
           distance_name is not one of them, or when a distance names a
           position the file does not have
*/
IMPBFFEXPORT RotamerFpsSelection read_rotamer_fps(
        const std::string& path, const std::string& distance_name = "");

//! The fps.json entry of a rotamer position (`simulation_type` `R1`).
/*!
    \param[in] chain,residue,library what the position is
    \param[in] atom_name the anchor atom
    \param[in] dye the dye's name; empty leaves the key out
    \param[in] temperature K; NaN leaves the key out
    \param[in] electrostatic 1 true, 0 false, -1 leaves the key out
    \param[in] potential `lj` or `gauss`; empty leaves the key out
    \return the entry, as JSON text

    A key that is left out is not the same as a key written empty: the reader
    falls back on its own default for the first and takes the second at its
    word.
*/
IMPBFFEXPORT std::string rotamer_position_payload(
        const std::string& chain, int residue, const std::string& library,
        const std::string& atom_name = "CA", const std::string& dye = "",
        double temperature = std::numeric_limits<double>::quiet_NaN(),
        int electrostatic = -1, const std::string& potential = "");

//! The `R1` entries describing existing ensembles, keyed as they are.
/*! Each ensemble's own `params` -- the record of how it was screened -- fills
    the temperature, the potential and the electrostatic flag, so a file
    written from ensembles says how to reproduce them. */
IMPBFFEXPORT std::string rotamer_positions_payload(
        const std::map<std::string, RotamerEnsemble>& ensembles,
        const std::string& atom_name = "CA");

//! Predicted fps.json distance entries between rotamer ensembles.
/*!
    \param[in] ensembles the screened ensembles, keyed by position name
    \param[in] pairs which of them to measure, donor first
    \param[in] forster_radius \f$R_0\f$ at \f$\kappa^2 = 2/3\f$, A
    \param[in] distance_type `RDAMean` (\f$\langle R_{DA}\rangle\f$),
               `RDAMeanE` (the FRET-averaged one) or `Rmp` (between mean
               positions) -- from the full pair matrix, with no sampling
    \param[in] error the error bar to write; negative uses \p error_fraction
    \param[in] error_fraction of the distance, when no \p error is given
    \param[in] kappa2 `isotropic` puts \f$\kappa^2 = 2/3\f$ on every pair,
               which is the fps.json and #IMP::bff::ProbeNetworkRestraint
               convention and is what makes the number comparable with an
               AV's; `dipoles` uses the ensembles' own per-pair
               \f$\kappa^2\f$, which is the orientation-resolved answer
    \return the entries, keyed `<donor>_<acceptor>`, as JSON text
    \throw ValueException for an unknown \p distance_type or \p kappa2, or a
           pair naming an ensemble that is not there
*/
IMPBFFEXPORT std::string distances_from_ensembles(
        const std::map<std::string, RotamerEnsemble>& ensembles,
        const std::vector<std::pair<std::string, std::string> >& pairs,
        double forster_radius,
        const std::string& distance_type = "RDAMeanE", double error = -1.0,
        double error_fraction = 0.05,
        const std::string& kappa2 = "isotropic");

//! One ensemble per fps.json position that names a rotamer library.
/*!
    A position's `rotamer_library` (or `library_name` / `libname` / `library`)
    field selects the library; \p library_map, keyed by position name,
    overrides it or supplies one for a position that names none -- which is
    how an AV-only fps file is screened as rotamers. Positions with neither
    are skipped rather than guessed at.

    The structure is read once and each library at most once, however many
    positions share them.

    \param[in] fps_json the fps.json path
    \param[in] structure the PDB the positions refer to
    \param[in] library_map position name -> library name
    \param[in] options the screening parameters, for every position alike
    \param[in] frame_index which model of a multi-MODEL PDB
*/
IMPBFFEXPORT std::map<std::string, RotamerEnsemble> rotamer_ensembles_from_fps(
        const std::string& fps_json, const std::string& structure,
        const std::map<std::string, std::string>& library_map =
                std::map<std::string, std::string>(),
        const RotamerSiteOptions& options = RotamerSiteOptions(),
        int frame_index = 0);

//! Write (or merge into) an fps.json with rotamer (`R1`) positions.
/*!
    \param[in] path the file to write
    \param[in] positions_json,distances_json the entries, as JSON objects
    \param[in] merge_into an existing fps.json whose positions and distances
               are kept; entries of the same name are replaced. Its score sets
               and extra sections are carried over untouched.
    \param[in] validate check against the fps.json schema before writing
    \throw ValueException when validation fails
*/
IMPBFFEXPORT void write_rotamer_fps(const std::string& path,
                                    const std::string& positions_json,
                                    const std::string& distances_json = "{}",
                                    const std::string& merge_into = "",
                                    bool validate = true);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_ROTAMERFPS_H
