/**
 *  \file IMP/bff/FPSIO.h
 *  \brief fps.json, and the legacy C# FPS `.txt` files it replaced.
 *
 * The one reader and writer of the format #IMP::bff::fps_json_schema defines.
 * fps.json had three readers and no definition until PRD-97; files are now
 * checked against the schema rather than against another parser.
 *
 * The legacy C# `.txt` formats are **read only**. A format nobody can still read
 * is data that has been lost, and a decade of measurements live in those files —
 * so they are readable, and deliberately not writable: nothing should produce
 * another one.
 *
 * Everything here speaks JSON **text**, not typed records. A position carries an
 * open set of keys whose types differ per key and whose unknown members are
 * preserved on the way through; that is what a JSON object is, and turning it
 * into a struct here would either drop the unknown keys or reinvent the
 * object.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_FPSIO_H
#define IMPBFF_FPSIO_H

#include <IMP/bff/bff_config.h>

#include <IMP/bff/Base.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! The four top-level sections of an fps.json document, as JSON text.
struct IMPBFFEXPORT FPSDocument {
    //! `{name: {param: value}}`; may be `{}`.
    std::string positions, distances;
    //! The `"χ²"` section: `{name: {"distances": [...], ...}}`.
    std::string score_sets;
    //! Every other top-level key, e.g. `Evaluators`.
    std::string extra;
    //! Molecule names in order of appearance; only the legacy `.txt` door fills it.
    std::vector<std::string> molecules;

    FPSDocument()
        : positions("{}"), distances("{}"), score_sets("{}"), extra("{}") {}

    IMP_SHOWABLE_INLINE(FPSDocument, out << "FPSDocument()");
};
IMP_VALUES(FPSDocument, FPSDocuments);

//! Read a C# FPS labelling-positions `.txt` file.
/*!
    Each line is `name molecule dye type ...`, the tail depending on the type:
    `AV1` carries linker length, width, one radius and an **atom serial**; `AV3`
    three radii and an atom serial; `XYZ` a fixed coordinate triple. Lines that
    are blank, commented, or too short for their type are skipped — the format
    has no version and no schema, so a line that does not parse is not an error,
    it is a line from a dialect this does not know.

    The atom serial is resolved to chain, residue number and atom name by
    scanning the molecule's PDB. Without one the serial is carried in
    `residue_seq_number` as a proxy, which is what the caller gets when no PDB
    for that molecule was found.

    \param[in] path the positions file
    \param[in] pdb_paths PDBs to resolve atom serials against, matched to a
               molecule by basename. Empty scans the positions file's own
               directory for `*.pdb`.
*/
IMPBFFEXPORT FPSDocument read_old_lps_txt(
        const std::string& path,
        const std::vector<std::string>& pdb_paths = std::vector<std::string>());

//! Read a C# FPS experimental-distances `.txt` file. \return the JSON text
/*! A lone token on the first line is the file-wide `distance_type`; anything
    else is the first distance and the file is re-read from the top. */
IMPBFFEXPORT std::string read_old_distances_txt(const std::string& path);

//! The distance names of a legacy C# FPS distances `.txt`, **in file order**.
/*!
    The order matters and is otherwise lost. FPS's project file stores its
    distance selection as a `Boolean[]` parallel to its `DistanceList`, i.e.
    to this file's line order (`ProjectData.SelectedDistances`), while
    #read_old_distances_txt hands back a JSON *object* whose keys come out
    sorted -- so resolving a positional selection through that reader silently
    pairs each flag with the wrong distance. This is the order to index with.

    Names are built exactly as #read_old_distances_txt builds its keys,
    `position1_position2`, and the same lines are skipped, so the two readings
    of one file cannot disagree about what a distance is called.

    \throw IOException when the file cannot be read
*/
IMPBFFEXPORT std::vector<std::string> read_old_distances_order(
        const std::string& path);

//! Load an fps.json file, or a legacy C# `.txt` pair.
/*!
    A path that does not end in `.json` is read as a legacy positions file, and
    a `Distances.txt` beside it is read for the distances.

    \param[in] validate check against the schema and throw listing the
           violations when the file does not conform
    \throw ValueException when \p validate is set and the file does not conform
    \throw IOException when the file cannot be read
*/
IMPBFFEXPORT FPSDocument read_fps_json(
        const std::string& path,
        const std::vector<std::string>& pdb_paths = std::vector<std::string>(),
        bool validate = false);

//! Write an fps.json file.
/*! \param[in] score_sets_json the `"χ²"` section; `{}` omits the key entirely
    \param[in] validate check the assembled payload before writing
    \throw ValueException when \p validate is set and the payload does not
           conform — refusing to write a non-conforming file */
IMPBFFEXPORT void write_fps_json(const std::string& path,
                                 const std::string& positions_json,
                                 const std::string& distances_json,
                                 const std::string& score_sets_json = "{}",
                                 const std::string& extra_json = "{}",
                                 bool validate = false);

//! Keep only the positions the C++ AV scorer understands, and their distances.
/*!
    Rotamer-ensemble positions (`simulation_type == "R1"`) are Python-only:
    #IMP::bff::ProbeNetworkRestraint never reads `simulation_type` and would score
    one as an AV1 with its AV parameters. A distance survives only if **both**
    its ends do.
*/
IMPBFFEXPORT FPSDocument fps_positions_for_docking(
        const std::string& positions_json,
        const std::string& distances_json = "{}");

//! The `Evaluators` array of an fps.json file, as JSON text.
/*! `[]` for a missing or unparseable file — an evaluator list is optional, and
    a caller asking for one is asking whether there is one. */
IMPBFFEXPORT std::string read_evaluators_json(const std::string& path);

//! Replace the `Evaluators` key of an fps.json file, keeping the rest.
IMPBFFEXPORT void write_evaluators_json(const std::string& path,
                                        const std::string& evaluators_json);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_FPSIO_H
