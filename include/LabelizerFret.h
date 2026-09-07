/**
 *  \file IMP/bff/LabelizerFret.h
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
#ifndef IMPBFF_LABELIZERFRET_H
#define IMPBFF_LABELIZERFRET_H

#include <IMP/bff/bff_config.h>

#include <IMP/bff/LabelizerFeatures.h>
#include <IMP/bff/LabelizerScore.h>

#include <IMP/bff/Base.h>

#include <limits>
#include <map>
#include <string>
#include <vector>

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

#endif /* IMPBFF_LABELIZERFRET_H */
