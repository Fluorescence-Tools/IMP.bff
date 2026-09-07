/**
 *  \file IMP/bff/PtoProfile.h
 *  \brief PTO.MFDB — the vocabulary layer over the PTO container: which terms
 *         a tag may carry, what a column means, and where a result came from.
 *
 * #IMP::bff::PtoWriter frames bytes and deliberately knows nothing about what
 * they are. This says what a conforming file puts *beside* them: an artifact
 * kind, a row grain, a column's units, and a provenance edge back to whatever
 * the result was computed from. It is the PTO.MFDB profile of PTO 1.0, and it
 * is not specific to any one payload — the labelizer writes label-score tables
 * through it, and a rotamer container could carry its licence and citation
 * through the same tags.
 *
 * Three decisions are worth stating, because each rules out a tempting
 * alternative.
 *
 * **Tags are additive, never a decoding dependency.** A reader that wants the
 * bytes needs #IMP::bff::PtoReader and nothing here. That property is what lets
 * a container written by this package be read by a walker in another language
 * with no shared code, and it is worth more than any validation this file could
 * enforce at read time.
 *
 * **The controlled vocabularies are compiled in, not read from a dictionary at
 * run time.** They come from `mmfdb_flr_ext.dic`, which lives in a sibling
 * repository that is not a dependency of this one; parsing it at run time would
 * make writing a file depend on a checkout being present. This is the pattern
 * the package already uses for flrCIF names in
 * #IMP::bff::probe_flrcif_items and `FPS.h`, and it comes with the same
 * obligation: `test/label/test_pto_profile.py` re-reads the dictionary and
 * fails when these lists drift from it. A term is *checked*, never *recalled*.
 *
 * **Absence is a value.** `_mmfdb_label_score` says so in as many words: a
 * position with no score carries a `status` saying why and no number at all.
 * The alternative — the reference's `-1` for excluded and `0` for no
 * contribution, in the same column as real scores — cannot be read back.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_PTOPROFILE_H
#define IMPBFF_PTOPROFILE_H

#include <IMP/bff/bff_config.h>

#include <IMP/bff/Base.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! The profile this package writes, and the version of it.
IMPBFFEXPORT extern const char* const MFDB_PROFILE;
//! The profile version written into `_mmfdb_container.profile_version`.
IMPBFFEXPORT extern const char* const MFDB_PROFILE_VERSION;
//! The lowest profile version a reader must understand to read the file.
IMPBFFEXPORT extern const char* const MFDB_PROFILE_READ_VERSION;
//! The `_dictionary.version` of `mmfdb_flr_ext.dic` these terms came from.
IMPBFFEXPORT extern const char* const MFDB_DICTIONARY_VERSION;

// ---------------------------------------------------------------------------
// One tag
// ---------------------------------------------------------------------------

//! A tag: an mmCIF item name and its value.
/*! The item name *is* the tag name, so an object's tags read as the dictionary
    row they are and no second naming convention exists to be translated. */
struct IMPBFFEXPORT MfdbTag {
    //! An mmCIF item name, e.g. `"_mmfdb_artifact.artifact_kind"`.
    std::string item;
    //! Its value.
    std::string value;

    MfdbTag() {}
    MfdbTag(const std::string& i, const std::string& v) : item(i), value(v) {}
    IMP_SHOWABLE_INLINE(MfdbTag, out << item << " = " << value);
};
IMP_VALUES(MfdbTag, MfdbTags);

// ---------------------------------------------------------------------------
// The controlled vocabularies
// ---------------------------------------------------------------------------

//! `_mmfdb_artifact.artifact_kind`.
IMPBFFEXPORT const std::vector<std::string>& mfdb_artifact_kinds();
//! `_mmfdb_artifact.data_format`.
IMPBFFEXPORT const std::vector<std::string>& mfdb_data_formats();
//! `_mmfdb_artifact.row_grain` — what one row of a tabular artifact *is*.
IMPBFFEXPORT const std::vector<std::string>& mfdb_row_grains();
//! `_mmfdb_edge.relationship_type`.
IMPBFFEXPORT const std::vector<std::string>& mfdb_relationship_types();
//! `_mmfdb_column.units`.
IMPBFFEXPORT const std::vector<std::string>& mfdb_units();
//! `_mmfdb_operation.operation_type` — the subset this package writes.
/*! The dictionary's enumeration is much larger and mostly photon-level; a
    reader must expect terms that are not here. `clustering` covers grouping
    conformations — a rotamer library built from a trajectory, MSM microstates,
    #IMP::bff::cluster_frames_leader — and is distinct from the dictionary's
    `ndxplorer_clustering`, which is one tool's burst-selection step. */
IMPBFFEXPORT const std::vector<std::string>& mfdb_operation_types();
//! `_mmfdb_optical_property.property_name` — the canonical dye scalars.
/*! The stored data is not confined to these: the reference spectra database
    holds 19882 rows under 76 names written by importers that never agreed.
    This is what a *writer* uses; a reader must expect others. */
IMPBFFEXPORT const std::vector<std::string>& mfdb_optical_property_names();
//! `_mmfdb_optical_property.unit`.
IMPBFFEXPORT const std::vector<std::string>& mfdb_optical_property_units();
//! `_mmfdb_spectrum.spectrum_type` — dyes *and* the optics they are measured
//! through, since a filter's transmission is the same shape as an emission.
IMPBFFEXPORT const std::vector<std::string>& mfdb_spectrum_types();
//! `_mmfdb_spectrum.wavelength_unit`. `nm`, not `nanometres` — see the item.
IMPBFFEXPORT const std::vector<std::string>& mfdb_spectrum_wavelength_units();
//! `_mmfdb_spectrum.intensity_unit`.
IMPBFFEXPORT const std::vector<std::string>& mfdb_spectrum_intensity_units();
//! `_mmfdb_probe.probe_type` -- what a probe *is*, as opposed to
//! `_flr_probe_list.probe_origin`, which is how it got onto the molecule.
IMPBFFEXPORT const std::vector<std::string>& mfdb_probe_types();
//! `_mmfdb_label_score.score_type`.
IMPBFFEXPORT const std::vector<std::string>& mfdb_label_score_types();
//! `_mmfdb_label_score.status`.
IMPBFFEXPORT const std::vector<std::string>& mfdb_label_score_statuses();
//! `_mmfdb_label_score.definition`.
IMPBFFEXPORT const std::vector<std::string>& mfdb_label_score_definitions();

//! Is \p value in the enumeration \p item declares?
/*! \param[in] item an mmCIF item name this profile knows an enumeration for
    \param[in] value the candidate
    \return false also when \p item has no enumeration here — an unknown item
            cannot vouch for a value */
IMPBFFEXPORT bool mfdb_is_term(const std::string& item,
                               const std::string& value);

//! #mfdb_is_term, or throw.
/*! \throw ValueException naming the item, the value and what was allowed. A
    writer that cannot find a word does not coin one; it adds one to the
    dictionary first, and this is what makes that unavoidable. */
IMPBFFEXPORT void mfdb_check_term(const std::string& item,
                                  const std::string& value);

// ---------------------------------------------------------------------------
// Columns
// ---------------------------------------------------------------------------

//! What one column of a tabular artifact holds.
/*!
    Carried with the column rather than with the file, so a caller reading two
    columns out of a large table still learns what they are. Units are a
    `_mmfdb_column.units` term: **no unit means the unit is unknown**, and
    `dimensionless` is a positive claim for a ratio that genuinely has none.
*/
struct IMPBFFEXPORT MfdbColumn {
    //! `_mmfdb_column.name`.
    std::string name;
    //! `_mmfdb_column.units`, a #mfdb_units term, or empty for unknown.
    std::string units;
    //! `_mmfdb_column.item` — the mmCIF item this column corresponds to.
    std::string item;
    //! `_mmfdb_column.description`, free text.
    std::string description;

    MfdbColumn() {}
    MfdbColumn(const std::string& n, const std::string& u,
               const std::string& i, const std::string& d = "")
        : name(n), units(u), item(i), description(d) {}
    IMP_SHOWABLE_INLINE(MfdbColumn, out << "MfdbColumn(" << name << ")");
};
IMP_VALUES(MfdbColumn, MfdbColumns);

//! The columns as the JSON a `columns` tag carries.
/*! \throw ValueException when a units value is not a #mfdb_units term */
IMPBFFEXPORT std::string mfdb_columns_json(
        const std::vector<MfdbColumn>& columns);

// ---------------------------------------------------------------------------
// The tag sets
// ---------------------------------------------------------------------------

//! The container-level tags that declare conformance.
/*! A reader decides whether a file conforms by reading these, never by its
    name: a renamed file is still conformant and a `.mmfdb.pto` without them is
    not. */
IMPBFFEXPORT std::vector<MfdbTag> mfdb_container_tags();

//! The tags of one artifact.
/*!
    \param[in] artifact_id the identity that survives export, unlike the
               container-local `FileUID`
    \param[in] artifact_kind a #mfdb_artifact_kinds term
    \param[in] data_format a #mfdb_data_formats term
    \param[in] row_grain a #mfdb_row_grains term, or empty for a non-table
    \param[in] row_count rows, or negative to omit
    \param[in] checksum the SHA-256 of the payload, lowercase hexadecimal
    \param[in] columns the column descriptors, or empty
    \throw ValueException when any controlled value is not a term
*/
IMPBFFEXPORT std::vector<MfdbTag> mfdb_artifact_tags(
        const std::string& artifact_id, const std::string& artifact_kind,
        const std::string& data_format, const std::string& row_grain,
        long row_count, const std::string& checksum,
        const std::vector<MfdbColumn>& columns = std::vector<MfdbColumn>());

//! The tags of the operation that produced the derived artifacts.
/*!
    `settings_json` must be **complete**. A partial settings record is worse
    than none, because it looks reproducible and is not; `settings_hash` is the
    identity of a run, so re-running with the same settings replaces an
    artifact in place and changing one produces a new artifact.

    \param[in] operation_type a #mfdb_operation_types term
    \param[in] algorithm which implementation, where the type says only what
    \param[in] settings_json the complete settings, as JSON text
    \param[in] software_package,software_version who computed it
    \throw ValueException when `operation_type` is not a term
*/
IMPBFFEXPORT std::vector<MfdbTag> mfdb_operation_tags(
        const std::string& operation_type, const std::string& algorithm,
        const std::string& settings_json, const std::string& software_package,
        const std::string& software_version);

//! One provenance edge: what a derived artifact came from, and how.
/*!
    Three facts, three tags, because asking a file how a result relates to what
    it came from must not return an integer: the parent, the relation, and —
    when the two tables count different things — the columns that join them.

    \param[in] source_node_id the parent artifact
    \param[in] target_node_id the derived artifact
    \param[in] relationship_type a #mfdb_relationship_types term
    \param[in] source_row_column,target_row_column the join, or empty
    \throw ValueException when `relationship_type` is not a term
*/
IMPBFFEXPORT std::vector<MfdbTag> mfdb_edge_tags(
        const std::string& source_node_id, const std::string& target_node_id,
        const std::string& relationship_type,
        const std::string& source_row_column = "",
        const std::string& target_row_column = "");

//! Who made the data an artifact holds, and on what terms.
/*!
    Attribution is **per artifact, not per document**, and that is not a
    hypothesis: a family container can hold libraries under different licences
    at once — the shipped dye and spin-label rotamers are GPL-3.0-only while the
    converted Dunbrack side-chain table is "free for academic use" — and a
    document-level licence would be wrong the moment two of them share a file.
*/
struct IMPBFFEXPORT MfdbAttribution {
    //! Upstream authors of the *data*, as they name themselves. Free text.
    std::string author;
    //! A DOI where one exists, otherwise enough reference to find it.
    std::string citation;
    //! An SPDX identifier — `"GPL-3.0-only"`, `"MIT"` — where one exists.
    std::string license;
    //! The conditions in words, for terms SPDX has no identifier for.
    std::string terms;
    //! Where the terms actually live.
    std::string terms_url;
    //! What this was made from, outside this database: a file, a release.
    std::string source;
    //! The intermediary it reached us through, when the licences differ along
    //! the chain — Dunbrack-2010 arrives via FASPR, which is MIT while the
    //! library is not.
    std::string redistributed_via;

    IMP_SHOWABLE_INLINE(MfdbAttribution,
                        out << "MfdbAttribution(" << (license.empty() ? terms
                                                                      : license)
                            << ")");
};
IMP_VALUES(MfdbAttribution, MfdbAttributions);

//! The attribution of one artifact, as tags.
/*!
    **At least one of `license` and `terms` is required, and neither is
    preferred.** A great deal of scientific data is released on conditions SPDX
    has no identifier for, and coining a `LicenseRef` for "free for academic
    use" would state something the upstream did not; refusing the artifact
    because its terms are unusual is how a notice comes to be dropped instead
    of carried. Empty fields are omitted rather than written blank.

    \param[in] attribution who, from where, on what terms
    \throw ValueException when neither a license nor terms is given
*/
IMPBFFEXPORT std::vector<MfdbTag> mfdb_attribution_tags(
        const MfdbAttribution& attribution);

//! The SHA-256 of a payload, lowercase hexadecimal.
/*! `_mmfdb_artifact.checksum`, with `checksum_algorithm` `"sha256"`.
    Opening a file never hashes a payload; verification is explicit. */
IMPBFFEXPORT std::string mfdb_checksum(const std::string& bytes);

IMPBFF_END_NAMESPACE

#endif /* IMPBFF_PTOPROFILE_H */
