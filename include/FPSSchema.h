/**
 *  \file IMP/bff/FPSSchema.h
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
#ifndef IMPBFF_FPSSCHEMA_H
#define IMPBFF_FPSSCHEMA_H

#include <IMP/bff/bff_config.h>

#include <IMP/showable_macros.h>
#include <IMP/value_macros.h>

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
    FPS_ARRAY = 4
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
/*! R1 is missing, and that is the point: `IMP::bff::AVNetworkRestraint` never
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

#endif //IMPBFF_FPSSCHEMA_H
