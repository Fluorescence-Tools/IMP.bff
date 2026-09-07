#ifndef IMPBFF_LABELIZER_H
#define IMPBFF_LABELIZER_H

/**
 *  \file IMP/bff/Labelizer.h
 *  \brief The native Labelizer: label-site scores and the FRET pair score built on them.
 *
 * Four former headers, in the order the score is computed:
 *
 * 1. **Features** (formerly `LabelizerFeatures.h`) -- the structural
 *    quantities a label-site score is computed from.
 * 2. **Score** (formerly `LabelizerScore.h`) -- the fitted likelihood tables
 *    and the weighted geometric mean over the features.
 * 3. **FRET pair** (formerly `LabelizerFret.h`) -- which two labelling sites
 *    make the best pair, over the site scores and the distances between them.
 * 4. **IO** (formerly `LabelizerIO.h`) -- a scored structure as one
 *    `.mmfdb.pto` container.
 *
 * The `Ll` prefix and the `ll_` free functions mark what came from the
 * Labelizer package, so a reader can tell the ported model from the module's
 * own physics without consulting a PRD.
 *
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

// -------- from LabelizerFeatures.h --------
/**
 *  (formerly IMP/bff/LabelizerFeatures.h, now a section of this file)
 *  \brief The structural quantities a label-site score is computed from:
 *         secondary structure, half-sphere exposure, relative solvent
 *         accessibility, residue depth, and an imported conservation grade.
 *
 * These were the parts of the Labelizer package (Gebhardt *et al.*,
 * *Nat. Commun.* **16**, 3305, 2025) that did not run: `ss` shelled out to
 * DSSP and refused to start on macOS at all
 * (`secondary_structure.py:126`, `RuntimeError("Unknown platform")`), while
 * `se` and `me` shelled out to MSMS, whose bundled binaries are 32-bit
 * ppc/i386 Mach-O. A score model whose default solvent-exposure term cannot be
 * evaluated is not a model anyone can check, so the two external programs are
 * replaced here by native kernels and the difference is *measured* rather than
 * asserted — see `okf/validation/labelizer_ab.md`.
 *
 * Nothing here re-implements a kernel the module already has: the accessible
 * surface comes from #IMP::bff::solvent_accessible_surface_area, the structure
 * from #IMP::bff::read_pdb_records.
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

// ---------------------------------------------------------------------------
// The residue view
// ---------------------------------------------------------------------------

//! One residue, and where its atoms are in the flat arrays of #LlStructure.
/*! The backbone indices are resolved once, at read time, because every kernel
    below wants them and each would otherwise re-scan the atom list. A missing
    atom is `-1` and stays `-1`: a residue without a CA cannot be scored, and
    substituting a neighbouring atom yields a plausible, wrong answer. */
struct IMPBFFEXPORT LlResidue {
    //! Author chain identifier, as the PDB writes it.
    std::string chain;
    //! Author residue number.
    int seq_id;
    //! Three-letter residue name, e.g. `"TRP"`.
    std::string comp_id;
    //! Indices into #LlStructure::xyz (as triples) and the parallel arrays.
    std::vector<int> atoms;
    //! Backbone atom indices, or -1 when the atom is absent.
    int n, ca, c, o, cb;

    LlResidue() : seq_id(0), n(-1), ca(-1), c(-1), o(-1), cb(-1) {}
    IMP_SHOWABLE_INLINE(LlResidue,
                        out << "LlResidue(" << chain << seq_id << ":"
                            << comp_id << ")");
};
IMP_VALUES(LlResidue, LlResidues);

//! A structure as the scoring layer sees it: flat coordinates, grouped residues.
struct IMPBFFEXPORT LlStructure {
    //! `3N` coordinates, xyzxyz…
    std::vector<double> xyz;
    //! `N` van der Waals radii.
    std::vector<double> vdw;
    //! `N` atom names, stripped of padding.
    std::vector<std::string> atom_name;
    //! Residues in file order.
    std::vector<LlResidue> residues;

    //! Number of atoms.
    std::size_t size() const { return vdw.size(); }
};

//! The twenty standard amino acids, three-letter, uppercase.
IMPBFFEXPORT const std::vector<std::string>& ll_standard_residues();

//! Read a structure and group it into residues.
/*!
    Reads through #IMP::bff::read_pdb_records, so the parse and its cache are
    shared with the AV builder.

    \param[in] pdb_path the structure
    \param[in] protein_only keep only the residues named by
               #ll_standard_residues. This is how the reference's
               `remove_hetatoms` (`pdbhelper.py:16`, which detaches every
               residue whose hetflag is not blank) is reproduced: the record
               type is not carried on a #IMP::bff::PDBAtomRecord, and for a
               protein the two tests differ only on modified residues such as
               MSE, which the reference also drops.
    \param[in] model for a multi-model file, which model to take. The reference
               takes `structure[0]` unconditionally (`pdbhelper.py:35`).
    \return the grouped structure
    \throw IOException when the file has no coordinates
*/
IMPBFFEXPORT LlStructure ll_read_structure(const std::string& pdb_path,
                                           bool protein_only = true,
                                           int model = 0);

// ---------------------------------------------------------------------------
// Geometry the reference needed and biopython supplied
// ---------------------------------------------------------------------------

//! The Cbeta position of a residue, real or reconstructed for glycine.
/*!
    A glycine has no Cbeta, and the reference builds a virtual one rather than
    skipping the residue: the N vector is centred on CA, rotated by -120
    degrees about the CA->C axis, and added back to CA
    (`fret_score.py:293`, itself lifted from `Bio/PDB/HSExposure.py`).

    \param[in] s the structure
    \param[in] residue index into `s.residues`
    \param[out] out_view,n_out_view three coordinates, or zero-length when the
                residue lacks the backbone atoms to place one
*/
IMPBFFEXPORT void ll_cbeta_position(const LlStructure& s, int residue,
                                    double** out_view, int* n_out_view);

//! Half-sphere exposure about Cbeta, the `(up, down)` neighbour counts.
/*!
    For each residue the CA-CB direction splits a sphere of the given radius in
    two; `up` counts the CA atoms of other residues on the side the side chain
    points to, `down` the rest. It is a cheap burial measure that needs no
    surface, which is why the reference offers it as an alternative to MSMS
    (`solvent_exposure.py:141`, `Bio.PDB.HSExposure.HSExposureCB`).

    \param[in] s the structure
    \param[in] radius sphere radius in Angstrom; the reference uses 13.0 for
               the tryptophan and charge terms and 10/13/16 for the `I_SE*`
               solvent-exposure tables
    \param[out] out_view,n_out_view `2 * n_residues` values, up then down for
                each residue in `s.residues` order
*/
IMPBFFEXPORT void ll_half_sphere_exposure(const LlStructure& s, double radius,
                                          double** out_view, int* n_out_view);

// ---------------------------------------------------------------------------
// Secondary structure
// ---------------------------------------------------------------------------

//! DSSP secondary structure, one 8-state code per residue.
/*!
    Kabsch and Sander (1983), implemented natively. The electrostatic hydrogen
    bond energy between the C=O of residue `i` and the N-H of residue `j` is

        E = 332 * 0.42 * 0.20 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN)   kcal/mol

    with a bond declared below #LL_DSSP_HBOND_ENERGY. The amide hydrogen is not
    in the file and is placed the way DSSP places it: on N, along the direction
    opposite the preceding residue's C=O.

    The eight states are `H` (alpha helix), `B` (isolated bridge), `E` (strand),
    `G` (3-10 helix), `I` (pi helix), `T` (turn), `S` (bend) and `-` (none),
    assigned in that order of precedence.

    \param[in] s the structure
    \return one character per residue, in `s.residues` order
*/
IMPBFFEXPORT std::string ll_dssp(const LlStructure& s);

//! The hydrogen-bond energy below which a bond is declared, kcal/mol.
IMPBFFEXPORT extern const double LL_DSSP_HBOND_ENERGY;

// ---------------------------------------------------------------------------
// Solvent exposure
// ---------------------------------------------------------------------------

//! Which maximum-accessibility scale a relative accessibility is taken against.
enum LlMaxAsa {
    //! Tien *et al.* (2013) theoretical maxima. The reference's `N_SE1`.
    LL_MAXASA_WILKE = 0,
    //! Kabsch and Sander (1983). The reference's `N_SE2`.
    LL_MAXASA_SANDER = 1,
    //! Miller *et al.* (1987). The reference's `N_SE3`.
    LL_MAXASA_MILLER = 2
};

//! The maximum solvent-accessible area per residue type, Angstrom squared.
IMPBFFEXPORT std::map<std::string, double> ll_max_asa(LlMaxAsa scale);

//! Solvent-accessible surface area per residue, Angstrom squared.
/*!
    A Shrake-Rupley sum over the residue's atoms, through the module's existing
    #IMP::bff::solvent_accessible_surface_area.

    \param[in] s the structure
    \param[in] probe_radius rolling-probe radius; 1.4 is the DSSP convention
    \param[in] n_sphere_points sampling density per atom
    \param[out] out_view,n_out_view one area per residue
*/
IMPBFFEXPORT void ll_residue_sasa(const LlStructure& s,
                                  double probe_radius, int n_sphere_points,
                                  double** out_view, int* n_out_view);

//! Relative solvent accessibility per residue, a fraction.
/*!
    #ll_residue_sasa divided by the residue's maximum. A residue type absent
    from the scale yields a negative value, which the caller must treat as
    "not computed" rather than as zero exposure.

    \param[in] s the structure
    \param[in] scale which maximum-accessibility table
    \param[in] probe_radius,n_sphere_points passed to #ll_residue_sasa
    \param[out] out_view,n_out_view one fraction per residue
*/
IMPBFFEXPORT void ll_relative_solvent_accessibility(
        const LlStructure& s, LlMaxAsa scale, double probe_radius,
        int n_sphere_points, double** out_view, int* n_out_view);

//! Mean distance from a residue's atoms to the solvent-excluded surface.
/*!
    The reference's default solvent-exposure term is MSMS residue depth
    (`N_SE11_MEAN_SURFACE_DIST`, `solvent_exposure.py:157`), and MSMS does not
    run on any machine this package is developed on. This computes the same
    quantity from a numerically-built solvent-excluded surface: the probe
    centre is rolled over the atoms on the golden-spiral point set
    #IMP::bff::sphere_points, the points that no atom occludes are kept, and
    each is pulled back onto the surface by one probe radius.

    It is *not* MSMS and is not claimed to be. What it is measured against is
    the reference's own published output for 1DDB; the agreement is recorded in
    `okf/validation/labelizer_ab.md`.

    \param[in] s the structure
    \param[in] probe_radius rolling-probe radius, Angstrom
    \param[in] n_sphere_points points per atom used to build the surface
    \param[out] out_view,n_out_view one depth per residue, Angstrom
*/
IMPBFFEXPORT void ll_residue_depth(const LlStructure& s, double probe_radius,
                                   int n_sphere_points, double** out_view,
                                   int* n_out_view);

//! Distance from each residue's Cbeta to the solvent-excluded surface.
/*!
    The observable of the `N_SE10_CB_SURFACE_DIST` table, and a different
    quantity from #ll_residue_depth: that averages over every atom of the
    residue, this measures the one atom a label is attached at. The bins differ
    accordingly -- the Cbeta table starts at 0.69 Angstrom where the mean-atom
    table starts at 1.46.

    Glycine has no Cbeta and gets the reconstructed one from
    #ll_cbeta_position, so every residue has a value.

    \param[in] s the structure
    \param[in] probe_radius,n_sphere_points as #ll_residue_depth
    \param[out] out_view,n_out_view one distance per residue, Angstrom;
                negative where no Cbeta could be placed
*/
IMPBFFEXPORT void ll_cbeta_depth(const LlStructure& s, double probe_radius,
                                 int n_sphere_points, double** out_view,
                                 int* n_out_view);

//! The whole protein's formal charge: net, positive and negative.
/*!
    A charged dye is not indifferent to the charge of what it is attached to,
    and the reference computes this (`charge_environment.py:160`) — though
    nothing there consumes it, so it is offered rather than used.

    Formal charges only, from the same table the `ce` term counts with:
    ARG, LYS and HIS `+1`, ASP and GLU `-1`, everything else zero. Histidine
    counted as fully charged is the reference's choice and is wrong at
    physiological pH, where it is mostly neutral; it is reproduced because the
    charge environment term depends on it and the two must agree.

    The reference also offers a variant reading a per-atom partial charge out
    of the PDB **occupancy** column, for a PQR file loaded as a PDB. That is
    not ported: `IMP::bff::PDBAtomRecord` does not carry occupancy, and reading
    a charge out of a field that means "fraction of the site this atom
    occupies" is the same category of abuse as the score-in-the-B-factor this
    package's container exists to end. A caller wanting partial charges should
    read a PQR as a PQR.

    \param[in] s the structure
    \param[out] out_view,n_out_view three values: net, positive, negative.
                The negative total is reported as a negative number, so the
                three sum as `net = positive + negative`.
*/
IMPBFFEXPORT void ll_global_charge(const LlStructure& s, double** out_view,
                                   int* n_out_view);

// ---------------------------------------------------------------------------
// Conservation, imported
// ---------------------------------------------------------------------------

//! A conservation grade per residue, keyed `"<chain><seq_id>"`.
/*!
    Nothing here computes conservation: an alignment and a rate estimate are a
    different program, and the reference imports them too
    (`conservation_score.py`, whose ConSurf call lives in the web backend, not
    in the package). Both shapes the reference reads are supported.

    \param[in] path a ConSurf `.grades` table, or a PDB carrying the normalised
               grade in the B-factor column
    \return grade by residue key; an empty map when the file holds none
    \throw IOException when the file cannot be read
*/
IMPBFFEXPORT std::map<std::string, double> ll_read_consurf(
        const std::string& path);

//! The residue key the score tables and the reference CSVs are indexed by.
/*! `"<chain><seq_id>"`, e.g. `"A123"` — the reference's own convention
    (`labeling_parameter.py:129`), kept so its output files can be compared
    row for row. */
IMPBFFEXPORT std::string ll_residue_key(const std::string& chain, int seq_id);

IMPBFF_END_NAMESPACE

 /* IMPBFF_LABELIZERFEATURES_H */

// -------- from LabelizerScore.h --------
/**
 *  (formerly IMP/bff/LabelizerScore.h, now a section of this file)
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

 /* IMPBFF_LABELIZERSCORE_H */

// -------- from LabelizerFret.h --------
/**
 *  (formerly IMP/bff/LabelizerFret.h, now a section of this file)
 *  \brief The Labelizer FRET pair score: which two labelling sites make the
 *         most informative FRET assay.
 *
 * The per-residue score in `LabelizerScore.h` says whether a dye can be put
 * somewhere. This says whether a *pair* of places is worth measuring: a pair is
 * good when both sites are labelable and the distance between the dyes sits
 * where FRET is sensitive to it — near \f$R_0\f$ for one conformation, or far
 * apart between two conformations.
 *
 * Three ways of placing the dye are offered, and they are a cost ladder rather
 * than alternatives. #PROBE_MODEL_CBETA puts it on the Cbeta and costs nothing;
 * #PROBE_MODEL_ALPHA_CONE is the reference's analytic approximation of where the dye
 * mean ends up; #PROBE_MODEL_ACCESSIBLE_VOLUME builds the real cloud. The reference
 * screens every pair with a cheap method and rebuilds only the best few hundred
 * exactly (`fret_score.py:537`), and #ll_pair_scores does the same.
 *
 * **On the reference's `weights` switch, which does not arise here.** Its
 * deployed backend calls `calc_fret_score(..., weights=False)`
 * (`backendapp/tasks.py:167`), running `removeWeights` to set every occupied
 * voxel of the LabelLib grid to 1 — a parameter its own documentation calls
 * "I DON'T KNOW" (`labelizer.py:279`). Measured here: #IMP::bff::compute_av
 * already returns a uniformly weighted cloud (2327 points at a real site, one
 * distinct weight, exactly 1), so weighted and binarised are the same cloud
 * and there is nothing to switch. An option for it was written, measured to
 * change no distance by more than 0, and removed.
 *
 * **The accessible volume is the module's own.** The reference calls LabelLib,
 * which is banned here (owner rule, 2026-08-11); the cloud comes from
 * #IMP::bff::compute_av and the distance from #IMP::bff::model_distance, so
 * there is one AV implementation in this package and this is not it.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */



#include <limits>

IMPBFF_BEGIN_NAMESPACE

//! How the dye position a pair distance is measured between is obtained.
enum ProbeModel {
    //! The Cbeta itself. No dye model at all; the reference's `SIMPLE`.
    PROBE_MODEL_CBETA = 0,
    //! The reference's `GEBHARDT` analytic mean position; see
    //! #ll_alpha_cone_mean_position.
    PROBE_MODEL_ALPHA_CONE = 1,
    //! A real accessible volume, through #IMP::bff::compute_av. The
    //! reference's `KALININ`, with its LabelLib call replaced.
    PROBE_MODEL_ACCESSIBLE_VOLUME = 2
};

//! The dye and pair settings a FRET pair score needs.
struct IMPBFFEXPORT LlFretOptions {
    //! Which dye model the screen uses.
    ProbeModel probe_model;
    //! Which dye model the refinement uses; the reference refines with AVs.
    ProbeModel refine_probe_model;
    //! How many top-scoring pairs are recomputed with \p refine_probe_model.
    //! Zero refines nothing.
    int n_refine;
    //! Forster radius, Angstrom. Derive it with
    //! #IMP::bff::forster_radius rather than guessing.
    double forster_radius;
    //! Sites whose combined label score is below this are not paired at all
    //! (`config.py`, `LS_THRESHOLD`).
    double label_score_threshold;
    //! Linker length and width, Angstrom.
    double linker_length, linker_width;
    //! The three dye radii, Angstrom.
    double r1, r2, r3;
    //! AV grid spacing, Angstrom.
    double grid_resolution;
    //! Which pair distance the accessible-volume model reports: `"Rmp"` (the
    //! distance between the two mean positions), `"RDAMean"`
    //! (\f$\langle R_{DA}\rangle\f$) or `"RDAMeanE"` (\f$R_E\f$, the
    //! FRET-averaged distance).
    /*!
        **`RDAMeanE` is the default because it is what the reference reports.**
        Its `effDistance` inverts the mean *efficiency*
        (`label_lib_functions.py:138`), and that is not the distance between
        the mean positions. The ordering is
        \f$R_{mp} \le R_E \le \langle R_{DA}\rangle\f$: averaging a distance
        over two distributions cannot fall below the distance between their
        means (Jensen), and the \f$1/R^6\f$ weighting of the efficiency
        average pulls \f$R_E\f$ back towards the near side. Measured on real
        sites here the spread is **2.2 to 3.0 Angstrom**, and since the pair
        score peaks sharply at \f$R = R_0\f$ that is a large error in the
        quantity being ranked.

        Ignored by #PROBE_MODEL_CBETA and #PROBE_MODEL_ALPHA_CONE, which place a *point*
        and so have only one distance to report — as the reference's `SIMPLE`
        and `GEBHARDT` also do.
    */
    std::string distance_type;
    //! Neighbour radius the alpha cone measures burial over
    //! (`config.py`, `AV_CALC_RADIUS`).
    double alpha_cone_radius;
    //! The alpha cone's offset; negative derives it as the reference does.
    double alpha_cone_offset;
    //! Published or corrected arithmetic; see #LlModel.
    LlModel model;
    //! Restrict pairs to one donor chain and one acceptor chain.
    /*!
        Empty (the default) pairs every labelable site with every other,
        whatever chain they are in, which is what this module did before these
        fields existed. Set both and only pairs with the first site in
        \p donor_chain and the second in \p acceptor_chain survive.

        This is what makes a **homodimer** screen possible. A single-cysteine
        mutation puts the same residue on both protomers, so the only
        intermolecular distance that exists is i(A)-i(B); with the chains
        unset those pairs are drowned in intra-protomer ones, and on a
        structure with more than one copy of the assembly they are mixed with
        lattice contacts that are not biology at all. The reference Python
        carried exactly this argument (`fret_score.py:479`, filter at `:728`)
        and it was not ported.
    */
    std::string donor_chain, acceptor_chain;
    //! How the chains of the second conformation map onto the first.
    /*!
        Empty (the default) means the two files use the same chain names, which
        is the usual case. Give it when they do not: `{"A", "E"}, {"B", "F"}`
        reads state 1's chain A as state 2's chain E.
    */
    std::map<std::string, std::string> chain_map;

    LlFretOptions()
        : probe_model(PROBE_MODEL_ALPHA_CONE),
          refine_probe_model(PROBE_MODEL_ACCESSIBLE_VOLUME), n_refine(300),
          forster_radius(52.0), label_score_threshold(0.5),
          linker_length(20.0), linker_width(4.5), r1(9.0), r2(4.5), r3(1.5),
          grid_resolution(1.5), distance_type("RDAMeanE"),
          alpha_cone_radius(20.0), alpha_cone_offset(-1.0),
          model(LL_MODEL_PUBLISHED) {}
};

//! One scored pair of labelling sites.
struct IMPBFFEXPORT LlPairScore {
    //! The two positions, as `_flr_poly_probe_position` names a position.
    std::string asym_id_1, asym_id_2;
    int seq_id_1, seq_id_2;
    //! The FRET pair score.
    double value;
    //! The dye-dye distance, Angstrom. With two conformations, the first.
    double distance;
    //! The second conformation's distance; NaN when there is only one.
    double distance_2;
    //! The joined label score of the two sites.
    double joined_label_score;
    //! Which dye model produced `distance`, so a refined row is
    //! distinguishable from a screened one.
    ProbeModel probe_model;

    LlPairScore()
        // `distance_2` defaults to NaN, not 0: zero is a perfectly plausible
        // distance and would be read as an efficiency of 1, so "there is no
        // second conformation" must not look like "the two dyes coincide".
        : seq_id_1(0), seq_id_2(0), value(0.0), distance(0.0),
          distance_2(std::numeric_limits<double>::quiet_NaN()),
          joined_label_score(0.0), probe_model(PROBE_MODEL_CBETA) {}
    IMP_SHOWABLE_INLINE(LlPairScore,
                        out << "LlPairScore(" << asym_id_1 << seq_id_1 << "_"
                            << asym_id_2 << seq_id_2 << "=" << value << ")");
};
IMP_VALUES(LlPairScore, LlPairScores);

// ---------------------------------------------------------------------------
// The scores
// ---------------------------------------------------------------------------

//! The combined label score of the sites of a pair.
/*!
    Under #LL_MODEL_PUBLISHED this is `prod(scores) ** 0.5`, which is the
    reference's arithmetic (`fret_score.py:220`) and is the geometric mean only
    when there are two scores. With the four scores of a two-conformation pair
    it is the geometric mean *squared*, so the one- and two-conformation pair
    scores are not on a common scale. #LL_MODEL_CORRECTED uses
    `prod(scores) ** (1/N)`, which is the line the reference has commented out
    immediately below the one it runs.

    \param[in] scores the label scores of the sites, two or four of them
    \param[in] model which arithmetic
*/
IMPBFFEXPORT double ll_joined_label_score(const std::vector<double>& scores,
                                          LlModel model = LL_MODEL_PUBLISHED);

//! The pair score of one conformation: best when the efficiency is one half.
/*! \f$jls \cdot (1 - 2|E - 0.5|)\f$ (`fret_score.py:263`). A pair is most
    informative where the transfer efficiency responds most steeply to the
    distance, which is at \f$R \approx R_0\f$. */
IMPBFFEXPORT double ll_pair_score_single(double joined_label_score,
                                         double distance,
                                         double forster_radius);

//! The pair score of two conformations: best when the efficiency changes most.
/*! \f$jls \cdot |E(d_1) - E(d_2)|\f$ (`fret_score.py:225`). */
IMPBFFEXPORT double ll_pair_score_double(double joined_label_score,
                                         double distance_1, double distance_2,
                                         double forster_radius);

//! The negative-control pair score: a pair that should *not* change.
/*!
    \f$jls \cdot (1 - |1 - (E_1 + E_2)|) \cdot \max(0, 1 - 20|E_1 - E_2|)\f$
    (`fret_score.py:243`) — rewarding a mid-range mean efficiency and punishing
    any change in it. The reference computes this and then discards it, the
    line that would have stored it being commented out at `fret_score.py:528`;
    it is exposed here because a control pair is a real experimental need and
    the formula was already written.
*/
IMPBFFEXPORT double ll_pair_score_negative_control(double joined_label_score,
                                                   double distance_1,
                                                   double distance_2,
                                                   double forster_radius);

// ---------------------------------------------------------------------------
// Where the dye is
// ---------------------------------------------------------------------------

//! The reference's analytic estimate of the dye's mean position.
/*!
    An "alpha cone": the dye mean is placed along the outward direction — from
    the centroid of the atoms near the site towards the Cbeta — at a distance
    that shrinks as the site becomes buried (`fret_score.py:326`):

    \f[
        \mathrm{mean} = c_\beta + \left(\tfrac{3}{4} - \frac{d}{R}\right)
                        (\ell + \mathrm{off})\,\hat{n}
    \f]

    where \f$d\f$ is the distance from the Cbeta to the centroid of the atoms
    within \f$R\f$ of it, \f$\ell\f$ the linker length, and the offset, when
    not given, is \f$b + 0.54\ell - 0.0225\ell^2\f$ with
    \f$b = \max(1.7,\, 2\min(r_i) - 1.7)\f$ — an empirical fit of how far the
    linker actually reaches.

    It costs one neighbour search where an accessible volume costs a grid
    search, which is what makes screening every pair of a protein affordable.

    \param[in] s the structure
    \param[in] residue index into `s.residues`
    \param[in] options the linker and cone settings
    \param[out] out_view,n_out_view three coordinates, or zero-length when the
                residue has no placeable Cbeta
*/
IMPBFFEXPORT void ll_alpha_cone_mean_position(const LlStructure& s, int residue,
                                              const LlFretOptions& options,
                                              double** out_view,
                                              int* n_out_view);

//! The dye mean position at one site, by whichever model is asked for.
/*!
    \param[in] s the structure
    \param[in] pdb_path the same structure on disk; the accessible volume
               builder is a structure-level front door and reads it
    \param[in] residue index into `s.residues`
    \param[in] options which model, and its parameters
    \param[out] out_view,n_out_view three coordinates, or zero-length when the
                site cannot carry a dye — a buried site whose volume comes back
                empty returns nothing rather than a fabricated position
*/
IMPBFFEXPORT void ll_probe_mean_position(const LlStructure& s,
                                       const std::string& pdb_path, int residue,
                                       const LlFretOptions& options,
                                       double** out_view, int* n_out_view);

// ---------------------------------------------------------------------------
// The pairs
// ---------------------------------------------------------------------------

//! Score every pair of labelable sites of one structure.
/*!
    Sites whose combined label score is below
    #LlFretOptions::label_score_threshold are dropped before pairing, so the
    cost is quadratic in the *labelable* sites rather than in the residues.

    Pairs are screened with #LlFretOptions::probe_model and the top
    #LlFretOptions::n_refine are then recomputed with
    #LlFretOptions::refine_probe_model. **Each site's dye position is computed
    once and reused**; the reference rebuilds the inner site's accessible
    volume inside the inner loop (`fret_score.py:653`), which makes its exact
    mode quadratic in AV builds rather than linear. That is a pure cost
    difference and changes no number, so it is not behind #LlModel.

    \param[in] pdb_path the structure
    \param[in] label_scores combined scores by residue key, from
               #IMP::bff::ll_combined_by_key
    \param[in] options the dye, the pair filter and the arithmetic
    \return the pairs, highest score first
*/
IMPBFFEXPORT std::vector<LlPairScore> ll_pair_scores(
        const std::string& pdb_path,
        const std::map<std::string, double>& label_scores,
        const LlFretOptions& options);

//! Score every pair across two conformations of the same molecule.
/*!
    The pair score is the *change* in transfer efficiency, so this is the mode
    that answers "which pair would report on this conformational change".
    Positions are matched by residue number between the two structures, which
    is the reference's convention (`fret_score.py:512`) and means the two files
    must share a numbering.

    \param[in] pdb_path_1,pdb_path_2 the two conformations
    \param[in] label_scores_1,label_scores_2 combined scores for each
    \param[in] options the dye, the pair filter and the arithmetic
    \return the pairs, highest score first
*/
IMPBFFEXPORT std::vector<LlPairScore> ll_pair_scores_two_states(
        const std::string& pdb_path_1, const std::string& pdb_path_2,
        const std::map<std::string, double>& label_scores_1,
        const std::map<std::string, double>& label_scores_2,
        const LlFretOptions& options);

//! The change in Cbeta-Cbeta distance between two conformations.
/*!
    A residue-by-residue map of how far apart two positions move, which is what
    a conformational change looks like before any dye is involved
    (`fret_score.py:920`).

    The matrix is indexed by residue number, so it describes **one chain**.
    Leaving \p chain empty on a multimer folds every chain onto the same index
    and the map is meaningless -- give the chain.

    \param[in] s1,s2 the two conformations
    \param[in] chain which chain to map; empty takes every residue, which is
               only correct for a single-chain structure
    \param[out] out_view,n_out_view a square matrix, row-major, over the union
                of the residue numbers present; a position missing from either
                structure is zero
    \return the residue number the matrix starts at, so a row index can be
            turned back into a position
*/
IMPBFFEXPORT int ll_cbeta_difference_map(const LlStructure& s1,
                                         const LlStructure& s2,
                                         double** out_view, int* n_out_view,
                                         const std::string& chain = "");

IMPBFF_END_NAMESPACE

 /* IMPBFF_LABELIZERFRET_H */

// -------- from LabelizerIO.h --------
/**
 *  (formerly IMP/bff/LabelizerIO.h, now a section of this file)
 *  \brief A scored structure as one `.mmfdb.pto` container: the structure it
 *         was computed from, the scores, and how they were produced.
 *
 * The reference writes six CSVs, a JSON heat map, four PDBs whose B-factor
 * column carries a dimensionless score, and a zip of the lot
 * (`labelizer.py:423`). Those files cannot say what a row is, cannot say what
 * a column is in, and cannot say what produced them — and one of them, the
 * conservation PDB, is written to the same path it was read from, which is how
 * the shipped 1DDB example came to feed its own output back in as its input.
 *
 * This writes **one file**. It holds the structure verbatim, the scores as
 * tables whose columns are named by `_mmfdb_label_score` items, and the
 * complete settings; every controlled value is a term from the MMFDB
 * dictionary, checked by `Pto.h`'s profile section before it is written.
 *
 * The container is #IMP::bff::PtoWriter's and the vocabulary is
 * `Pto.h`'s. Nothing here invents either.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */



IMPBFF_BEGIN_NAMESPACE

//! The object names a label container uses, so a reader can ask for one.
/*! Every kind this package writes is under the `label.` prefix. `drot.` and
    `rot.bbdep.` belong to the rotamer containers; the namespaces are disjoint
    by agreement so one walker can read both. */
IMPBFFEXPORT extern const char* const LL_PTO_README;
//! The input structure, byte for byte as it was read.
IMPBFFEXPORT extern const char* const LL_PTO_STRUCTURE;
//! The per-position scores, one row per (position, score_type).
IMPBFFEXPORT extern const char* const LL_PTO_SCORES;
//! The pair scores, one row per pair.
IMPBFFEXPORT extern const char* const LL_PTO_PAIRS;
//! The complete settings the run used.
IMPBFFEXPORT extern const char* const LL_PTO_MODEL;

//! Write a scored structure as one PTO.MFDB container.
/*!
    The structure goes in **first and verbatim** and is never rewritten, so it
    can be recovered byte for byte and checked against its recorded SHA-256.
    A README object precedes it explaining, in the file, how to walk the
    framing and get the structure back out — a container outlives the software
    that wrote it, and the person who needs it most is the one for whom this
    library will not build.

    Absence is preserved: a position whose score was not computed is written
    with its `status` and **no value**, never a sentinel.

    \param[in] path the container to write; `.mmfdb.pto` by convention, though
               conformance is stated by the tags inside and not by the name
    \param[in] pdb_path the structure that was scored, embedded verbatim
    \param[in] scores rows from #IMP::bff::ll_score_structure
    \param[in] pairs rows from #IMP::bff::ll_pair_scores, or empty
    \param[in] settings_json the complete settings, as JSON text; a partial
               record is worse than none because it looks reproducible
    \throw IOException when the structure cannot be read or the container
           cannot be written
    \throw ValueException when a value is not a dictionary term
*/
IMPBFFEXPORT void ll_write_pto(const std::string& path,
                               const std::string& pdb_path,
                               const std::vector<LlScore>& scores,
                               const std::vector<LlPairScore>& pairs,
                               const std::string& settings_json);

//! Read the per-position scores back out of a container.
/*! \param[in] path a container written by #ll_write_pto
    \return the rows, in the order they were written
    \throw IOException when the file is not a PTO document or holds no scores */
IMPBFFEXPORT std::vector<LlScore> ll_read_pto_scores(const std::string& path);

//! Read the pair scores back out of a container.
/*! \return the rows, or empty when the container carries none */
IMPBFFEXPORT std::vector<LlPairScore> ll_read_pto_pairs(
        const std::string& path);

//! Recover the embedded structure, verifying it against its checksum.
/*!
    \param[in] path the container
    \param[in] out_pdb_path where to write the recovered structure
    \return the SHA-256 that was recorded and matched
    \throw IOException when the container has no structure, or when the bytes
           do not hash to what the container says they should
*/
IMPBFFEXPORT std::string ll_extract_pto_structure(
        const std::string& path, const std::string& out_pdb_path);

//! The settings JSON the run recorded.
IMPBFFEXPORT std::string ll_read_pto_settings(const std::string& path);

//! The settings of a scoring run, as the JSON the container stores.
/*! Every field of both options structs, so the run is reproducible from the
    file alone. */
IMPBFFEXPORT std::string ll_settings_json(
        const std::vector<LlParameter>& model, const LlOptions& options,
        const LlFretOptions& fret_options, const std::string& conservation_path);

IMPBFF_END_NAMESPACE

 /* IMPBFF_LABELIZERIO_H */

#endif  // IMPBFF_LABELIZER_H
