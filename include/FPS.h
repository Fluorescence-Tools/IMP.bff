#ifndef IMPBFF_FPS_H
#define IMPBFF_FPS_H

/**
 *  \file IMP/bff/FPS.h
 *  \brief fps.json: what a file may say, and reading and writing one.
 *
 * Two former headers: the **schema** (formerly `FPSSchema.h`) -- the field
 * tables, the three dialects, the validator -- and the **doors** (formerly
 * `FPSIO.h`) -- fps.json in and out, and the legacy C# FPS `.txt` readers.
 *
 * `FPSProject.h` and `FPSExport.h` are not here: both include `Docking.h`,
 * which is IMP-side (PRD-137's connection layer), and this header is what
 * the IMP-free readers depend on.
 *
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

// -------- from FPSSchema.h --------
/**
 *  (formerly IMP/bff/FPSSchema.h, now a section of this file)
 *  \brief What an fps.json file may say, and whether one says it.
 *
 * fps.json describes labelling positions and distance measurements on a
 * structure. It had three readers and no definition until PRD-97; the field
 * tables here **are** the definition, and `data/fps_json_schema.json` is derived
 * from them by #IMP::bff::fps_json_schema (a test regenerates it and fails on
 * drift).
 *
 * Three dialects share the format and are recorded per field:
 *
 * - `network` — ChiSurf Positions/Distances files;
 * - `flat` — imp.bff template files;
 * - `csfps` — what converting a legacy C# FPS `.txt` produces.
 *
 * Validation is self-contained: no JSON-Schema library is involved, and files
 * are checked against these tables rather than against another parser. Unknown
 * keys are a **warning**, never an error — readers and writers preserve them, so
 * an old file with an extra key still validates.
 *
 * References: the flrCIF extension dictionary (mmCIF `ihm_flr` extension),
 * categories `_flr_FPS_AV_parameter`, `_flr_FPS_global_parameter`,
 * `_flr_FPS_mean_probe_position` and `_flr_fret_distance_restraint`.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/bff_config.h>

#include <IMP/bff/Base.h>

#include <map>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! The JSON types an fps.json field may hold.
enum FPSFieldType {
    FPS_STRING = 0,
    FPS_NUMBER = 1,
    FPS_INTEGER = 2,
    FPS_BOOLEAN = 3,
    FPS_ARRAY = 4,
    //! A nested object: the project's parameter blocks, and nothing else yet.
    FPS_OBJECT = 5
};

//! The schema version `fps_json_schema()` stamps.
IMPBFFEXPORT std::string fps_schema_version();

//! One field of an fps.json object.
struct IMPBFFEXPORT FPSField {
    std::string name;
    FPSFieldType type;
    //! The value a reader may assume when the key is absent, as JSON text.
    /*! Empty means *no* default: the key is then required for the dialects in
        \p dialects. `"null"` is not used — an absent default is the empty
        string, so that a field whose default really is JSON null stays
        expressible. */
    std::string default_json;
    //! The canonical flrCIF item name, or empty where flrCIF has no such item.
    std::string flrcif;
    //! The local definition, for fields flrCIF does not define.
    std::string authored;
    //! Which historical dialect(s) carried the key.
    std::vector<std::string> dialects;
    //! The permitted values, when the field is an enumeration.
    std::vector<std::string> enum_values;
    bool required;

    FPSField() : type(FPS_STRING), required(false) {}

    IMP_SHOWABLE_INLINE(FPSField, out << "FPSField(" << name << ")");
};
IMP_VALUES(FPSField, FPSFields);

//! Errors and warnings from checking a payload against the tables.
/*! An empty \p errors means the file conforms. A \p warning flags an unknown
    field or a legal-but-suspicious value — an AV3 position with a zero third
    radius, say, which behaves like an AV1. */
struct IMPBFFEXPORT FPSValidation {
    std::vector<std::string> errors, warnings;

    bool get_is_valid() const { return errors.empty(); }

    IMP_SHOWABLE_INLINE(FPSValidation,
                        out << "FPSValidation(" << errors.size() << " errors, "
                            << warnings.size() << " warnings)");
};
IMP_VALUES(FPSValidation, FPSValidations);

//! The `simulation_type` values a *file* may carry.
/*! AV1: one dye radius. AV3: three. XYZ: a fixed mean position with no volume
    simulation, carried by legacy C# FPS conversions. R1: a rotamer ensemble —
    a FRETpredict-style library placed in the residue's backbone frame and
    Boltzmann-screened against the structure. */
IMPBFFEXPORT std::vector<std::string> fps_simulation_types();

//! The `simulation_type` values the **C++ AV scorer** understands.
/*! R1 is missing, and that is the point: `IMP::bff::ProbeNetworkRestraint` never
    reads `simulation_type` and would score an R1 position as an AV1 with its AV
    parameters. Filter with fps_positions_for_docking() first. */
IMPBFFEXPORT std::vector<std::string> fps_av_simulation_types();

//! `distance_type` values and their flrCIF spellings.
/*! `RDAMean` = \f$\langle R_{DA}\rangle\f$, `RDAMeanE` =
    \f$\langle R_{DA}\rangle_E\f$, `Rmp` = \f$R_{mp}\f$. */
IMPBFFEXPORT std::map<std::string, std::string> fps_distance_types();

IMPBFFEXPORT std::vector<FPSField> fps_position_fields();
IMPBFFEXPORT std::vector<FPSField> fps_distance_fields();
IMPBFFEXPORT std::vector<FPSField> fps_score_set_fields();

//! The fields an fps.json carries at the top level, beside its three sections.
/*! `strip_mask` here is the document's own: what the *structure* carries that
    no volume should be blocked by -- waters, ions, a crystallisation
    additive. A position's `strip_mask` is about the *labelling site*, and the
    two are a union. */
IMPBFFEXPORT std::vector<FPSField> fps_document_fields();

//! The fields of the `Project` section: what binds one run together (PRD-121 G1).
/*! Three of them (`parameters`, `conversion`, `av`) are nested objects whose
    own tables are #fps_mode_parameter_fields, #fps_conversion_fields and
    #fps_av_global_fields. #IMP::bff::FPSProject is the typed form. */
IMPBFFEXPORT std::vector<FPSField> fps_project_fields();

//! One mode's block of `Project.parameters` -- FPS's `FPSParameters`.
IMPBFFEXPORT std::vector<FPSField> fps_mode_parameter_fields();

//! `Project.conversion` -- FPS's `ConversionParameters`.
IMPBFFEXPORT std::vector<FPSField> fps_conversion_fields();

//! `Project.av` -- FPS's `AVGlobalParameters`.
IMPBFFEXPORT std::vector<FPSField> fps_av_global_fields();

//! Check one `Project` section, given as JSON text.
IMPBFFEXPORT FPSValidation validate_project(const std::string& project_json,
                                            const std::string& name =
                                                    "Project");

//! The strip mask a position is computed with: the document's, then its own.
/*! \param[in] document_mask the file's top-level `strip_mask`, or empty
    \param[in] position_mask the position's own `strip_mask`, or empty
    \return an expression selecting the union of the two, or empty when both
            are. Two masks are combined with `or` rather than concatenated,
            so each keeps its own precedence. */
IMPBFFEXPORT std::string combined_strip_mask(const std::string& document_mask,
                                             const std::string& position_mask);

//! The JSON-Schema document for the network-dialect fps.json, as JSON text.
/*! `data/fps_json_schema.json` is this function's output. */
IMPBFFEXPORT std::string fps_json_schema();

//! Check one position object, given as JSON text.
IMPBFFEXPORT FPSValidation validate_position(const std::string& position_json,
                                             const std::string& name =
                                                     "position");

//! Check one distance object, given as JSON text.
/*! \param[in] position_names names a `position1_name`/`position2_name` may
           refer to; empty skips the cross-reference check, which is what a
           single-object check wants */
IMPBFFEXPORT FPSValidation validate_distance(
        const std::string& distance_json,
        const std::string& name = "distance",
        const std::vector<std::string>& position_names =
                std::vector<std::string>());

//! Check a whole parsed payload, given as JSON text.
IMPBFFEXPORT FPSValidation fps_schema_validate(const std::string& payload_json);

IMPBFF_END_NAMESPACE

// -------- from FPSIO.h --------
/**
 *  (formerly IMP/bff/FPSIO.h, now a section of this file)
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


#endif  // IMPBFF_FPS_H
