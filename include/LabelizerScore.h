/**
 *  \file IMP/bff/LabelizerScore.h
 *  \brief The Labelizer label-site score: the fitted likelihood tables, the
 *         seven per-residue parameters, and the weighted geometric mean that
 *         combines them.
 *
 * A 1:1 port of the model published as Gebhardt *et al.*, *Nat. Commun.* **16**,
 * 3305 (2025), whose reference implementation is the Python package
 * `labelizer`. The structural inputs are in
 * #IMP::bff::ll_read_structure and its neighbours in `LabelizerFeatures.h`;
 * this header is the model on top of them.
 *
 * Two things about it are easy to get wrong and are therefore stated here.
 *
 * **A parameter score is a likelihood ratio, not a probability.** Each table
 * holds \f$P(\mathrm{labelable} \mid s) / P(\mathrm{labelable})\f$ binned over
 * the observable, so values run from 0 to about 3.9 and the combined score is
 * unbounded above. A cysteine on a surface loop exceeds 2 easily, and the
 * `LS_THRESHOLD` of 0.5 the pair layer filters on is a threshold on that
 * unbounded quantity.
 *
 * **The published implementation has defects, and they are reproduced by
 * default.** #LlModel selects; see its documentation for what differs. The
 * default is #LL_MODEL_PUBLISHED, so a caller who asks for nothing gets the
 * numbers in the paper.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_LABELIZERSCORE_H
#define IMPBFF_LABELIZERSCORE_H

#include <IMP/bff/bff_config.h>

#include <IMP/bff/LabelizerFeatures.h>

#include <IMP/bff/Base.h>

#include <map>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

// ---------------------------------------------------------------------------
// Which arithmetic
// ---------------------------------------------------------------------------

//! Reproduce the published implementation, or its intended arithmetic.
/*!
    Three defects in the reference change a number, and all three are
    reproduced under #LL_MODEL_PUBLISHED because the paper's values were
    computed with them.

    1. The joined label score of a pair takes `prod() ** 0.5` where the
       geometric mean is `prod() ** (1/N)` (`fret_score.py:220`, with the
       correct line commented out directly beneath it). For a two-residue,
       one-conformation pair the two agree; for the four values of a
       two-conformation pair the published form is the geometric mean
       *squared*, so the single- and double-conformation scores are not on the
       same scale.
    2. A parameter whose weight is zero is still evaluated, and if its value is
       exactly `0.0` it zeroes the whole labeling score
       (`labeling_score.py:166`). Weight zero does not mean "ignored".
    3. The excluded amino acid of the exclusion term is configurable, written
       to `bad_aa` and read from `self.aa`, so it is always methionine
       (`labelizer.py:233`).

    Two further reference defects are not selectable because they make the code
    raise rather than compute, and are ported as fixed: the tryptophan and
    charge terms call their base `set_up` without a table name and fail on load
    (`tryptophan_proximity.py:40`, `charge_environment.py:63`), and the
    accessible-volume guard tests its import the wrong way round
    (`fret_score.py:121`). Both are noted where they occur.
*/
enum LlModel {
    //! The arithmetic the paper's numbers were computed with. The default.
    LL_MODEL_PUBLISHED = 0,
    //! The arithmetic the reference's own documentation describes.
    LL_MODEL_CORRECTED = 1
};

// ---------------------------------------------------------------------------
// The fitted tables
// ---------------------------------------------------------------------------

//! One fitted likelihood-ratio table, as the reference ships it.
/*!
    The reference encodes the domain in the first character of the table name:
    `C_` is categorical and its keys are strings, `N_` and `I_` are numeric and
    their keys are bin centres (`labeling_parameter.py:198`).

    A numeric lookup is **nearest bin centre, with no interpolation**
    (`labeling_parameter.py:184`), which is reproduced exactly: interpolating
    would be smoother and would not be this model.
*/
struct IMPBFFEXPORT LlTable {
    //! The table name, e.g. `"C_SS1_SS"`.
    std::string name;
    //! True when the keys are strings rather than bin centres.
    bool categorical;
    //! Categorical entries.
    std::map<std::string, double> by_key;
    //! Numeric bin centres, ascending, and their values.
    std::vector<double> bins;
    std::vector<double> values;

    LlTable() : categorical(false) {}
    IMP_SHOWABLE_INLINE(LlTable, out << "LlTable(" << name << ")");
};
IMP_VALUES(LlTable, LlTables);

//! Load a fitted table by name from the shipped data directory.
/*!
    Cached on the name: a scoring run reads each table once however many
    residues it scores.

    \param[in] name e.g. `"C_CR1_Name"`, without the `_P_l_after_s.json` suffix
    \return the table
    \throw IOException when no such table ships
*/
IMPBFFEXPORT const LlTable& ll_load_table(const std::string& name);

//! Every fitted table that ships, by name.
IMPBFFEXPORT std::vector<std::string> ll_available_tables();

//! Look a numeric observable up in a table: the value of the nearest bin.
IMPBFFEXPORT double ll_lookup(const LlTable& table, double value);

//! Look a categorical observable up in a table.
/*! \throw ValueException when the key is not in the table */
IMPBFFEXPORT double ll_lookup_key(const LlTable& table, const std::string& key);

//! The one-letter code of a three-letter residue name, or `'X'`.
IMPBFFEXPORT char ll_one_letter(const std::string& comp_id);

// ---------------------------------------------------------------------------
// The model
// ---------------------------------------------------------------------------

//! One term of the label score: which parameter, which table, what weight.
struct IMPBFFEXPORT LlParameter {
    //! The reference's two-letter tag: `cs se ss ce tp cr me`.
    std::string tag;
    //! The fitted table this term reads, empty when the term is hard-coded.
    std::string table;
    //! The weight, an **integer repeat count** in the geometric mean.
    int weight;

    LlParameter() : weight(0) {}
    LlParameter(const std::string& t, const std::string& tb, int w)
        : tag(t), table(tb), weight(w) {}
    IMP_SHOWABLE_INLINE(LlParameter,
                        out << "LlParameter(" << tag << ", w=" << weight << ")");
};
IMP_VALUES(LlParameter, LlParameters);

//! The published model: conservation, solvent exposure, cysteine resemblance
//! and secondary structure at weight 1; tryptophan and charge at weight 0.
/*! `label_score_model_paper`, `labelizer.py:40`. The tryptophan and charge
    terms are switched off in the paper, not merely down-weighted. */
IMPBFFEXPORT std::vector<LlParameter> ll_model_paper();

//! What the scoring run may be told, beyond the model itself.
struct IMPBFFEXPORT LlOptions {
    //! Published or corrected arithmetic.
    LlModel model;
    //! Rolling-probe radius for the surface terms, Angstrom.
    double probe_radius;
    //! Unit-sphere samples per atom for the surface terms.
    int n_sphere_points;
    //! Exclusion term: the contact distance between C-alpha atoms, Angstrom.
    double exclusion_distance;
    //! Exclusion term: the exposure above which the excluded residue counts.
    double exclusion_exposure;
    //! Exclusion term: which residue is excluded. Ignored under
    //! #LL_MODEL_PUBLISHED, which is always methionine — see #LlModel.
    std::string exclusion_residue;
    //! Half-sphere radius the tryptophan and charge terms weight exposure by.
    double hse_radius;

    LlOptions()
        : model(LL_MODEL_PUBLISHED), probe_radius(1.4), n_sphere_points(590),
          exclusion_distance(6.0), exclusion_exposure(0.4),
          exclusion_residue("MET"), hse_radius(13.0) {}
};

// ---------------------------------------------------------------------------
// The result
// ---------------------------------------------------------------------------

//! One score at one position — a row of `_mmfdb_label_score`.
/*!
    The field names are the MMFDB dictionary's, so a table of these *is* the
    artifact it will be written as, with no second naming convention in
    between.

    A position that was not scored carries **no value**: `status` says why, and
    `value` must not be read. The reference writes `-1` for "excluded" and `0`
    for "no contribution" into the same column, which is why a consumer of its
    CSVs cannot tell a deliberately excluded residue from a genuinely zero
    score, nor either from a real value that happens to be zero.
*/
struct IMPBFFEXPORT LlScore {
    //! `_mmfdb_label_score.asym_id` — the chain, as the file names it.
    std::string asym_id;
    //! `_mmfdb_label_score.seq_id` — the residue number.
    int seq_id;
    //! `_mmfdb_label_score.comp_id` — the residue name.
    std::string comp_id;
    //! `_mmfdb_label_score.score_type` — a dictionary term, e.g.
    //! `"solvent_exposure"`, `"combined"`.
    std::string score_type;
    //! `_mmfdb_label_score.value`. Meaningful only when `status == "scored"`.
    double value;
    //! `_mmfdb_label_score.status` — `scored`, `excluded`, `unresolved` or
    //! `unavailable`.
    std::string status;

    LlScore() : seq_id(0), value(0.0), status("unavailable") {}
    IMP_SHOWABLE_INLINE(LlScore, out << "LlScore(" << asym_id << seq_id << " "
                                     << score_type << "=" << value << ")");
};
IMP_VALUES(LlScore, LlScores);

//! The dictionary `score_type` a two-letter reference tag names.
/*! `cs -> conservation`, `se -> solvent_exposure`, `ss -> secondary_structure`,
    `ce -> charge_environment`, `tp -> tryptophan_proximity`,
    `cr -> cysteine_resemblance`, `me -> methionine_exclusion`. The mapping is
    the one `mmfdb_flr_ext.dic` writes down in the description of
    `_mmfdb_label_score.score_type`; the tags themselves are that program's
    private naming and are not dictionary values. */
IMPBFFEXPORT std::string ll_score_type(const std::string& tag);

// ---------------------------------------------------------------------------
// The parameters
// ---------------------------------------------------------------------------

//! Every per-residue parameter score of one structure, tidy.
/*!
    One row per (position, parameter) for each parameter in \p model, whatever
    its weight — the reference evaluates weight-zero terms too, and under
    #LL_MODEL_PUBLISHED that is load-bearing.

    \param[in] s the structure
    \param[in] model the terms to evaluate; see #ll_model_paper
    \param[in] options probe radii, exclusion settings, and the arithmetic
    \param[in] conservation grades by residue key, from #ll_read_consurf; an
               empty map leaves every conservation row `unavailable`
    \return the parameter scores
*/
IMPBFFEXPORT std::vector<LlScore> ll_parameter_scores(
        const LlStructure& s, const std::vector<LlParameter>& model,
        const LlOptions& options,
        const std::map<std::string, double>& conservation);

//! The combined label score per position, from parameter scores.
/*!
    The weighted geometric mean, with the weight acting as an **integer repeat
    count** rather than an exponent (`labeling_score.py:187`): a term of weight
    3 is multiplied in three times and the root is over the total count. For
    integer weights this equals the exponent form; the reference offers no
    other, and `harmonic` raises there.

    Two statuses propagate rather than becoming numbers. A position where any
    term is `unavailable` is `unavailable`. A position where any term is
    exactly zero scores zero — the published zero-veto, which is how the
    exclusion term vetoes a site.

    \param[in] parameter_scores the output of #ll_parameter_scores
    \param[in] model the same model, for the weights
    \param[in] options for #LlModel
    \return one `combined` row per position
*/
IMPBFFEXPORT std::vector<LlScore> ll_labeling_score(
        const std::vector<LlScore>& parameter_scores,
        const std::vector<LlParameter>& model, const LlOptions& options);

//! Score a structure end to end: features, parameters, and the combination.
/*!
    The convenience door, and what the program and the examples call.

    \param[in] pdb_path the structure
    \param[in] model the terms; see #ll_model_paper
    \param[in] options the settings
    \param[in] conservation_path a ConSurf `.grades` or B-factor PDB, or empty
    \return the parameter rows followed by the `combined` rows
    \throw IOException when a file cannot be read
*/
IMPBFFEXPORT std::vector<LlScore> ll_score_structure(
        const std::string& pdb_path, const std::vector<LlParameter>& model,
        const LlOptions& options, const std::string& conservation_path = "");

//! The combined score per residue key, for the pair layer.
/*! \param[in] scores rows from #ll_score_structure or #ll_labeling_score
    \return `"<chain><seq_id>" -> value`, only for positions with `status`
            `scored` */
IMPBFFEXPORT std::map<std::string, double> ll_combined_by_key(
        const std::vector<LlScore>& scores);

IMPBFF_END_NAMESPACE

#endif /* IMPBFF_LABELIZERSCORE_H */
