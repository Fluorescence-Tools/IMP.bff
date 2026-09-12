/**
 *  \file IMP/bff/ProbePairSelection.h
 *  \brief Greedy probe selection — which probe pair (or, for a homo-oligomer, which
 *         labelling site) to measure next.
 *
 * Experiment planning, after Olga: given an ensemble of candidate structures and
 * a set of measurable pairs, repeatedly add the pair that most reduces the
 * expected RMSD between the true structure and the one the measurements would
 * pick out. Greedy because each step takes the locally best pair and never
 * revisits it — which is what makes it tractable, and what makes the scoring
 * step the whole cost. When the experiment is a homo-oligomer with a
 * statistical labelling mix, the unit changes: sites are mutated once and all
 * their cross-protomer combinations are measured, so
 * #IMP::bff::select_probe_positions greedies over sites and scores each step
 * by the pair set the sites imply.
 *
 * The cost is in the scoring step. Every remaining candidate pair has to be
 * scored against every ordered pair of ensemble frames, and the natural
 * vectorised form builds a `(candidates, frames, frames)` array — at a thousand
 * frames and a hundred candidates, 800 MB that exists only to be reduced away.
 * The kernel here never materialises it: the reduction runs inside the loop, so
 * the memory is `O(frames)` per candidate and the arithmetic is the same.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_PROBEPAIRSELECTION_H
#define IMPBFF_PROBEPAIRSELECTION_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Chi-squared right-tail probability \f$Q(\nu/2, \chi^2/2)\f$, elementwise.
/*!
    The regularized upper incomplete gamma. Because \f$\nu/2\f$ is an integer or
    a half-integer, a few-term closed form serves for all but the largest
    \f$\nu\f$, where a continued fraction takes over.

    \param[in] chisq values, any length
    \param[in] ndof degrees of freedom
    \return the right-tail probability, elementwise
*/
IMPBFFEXPORT std::vector<double> chi2_right_tail(
        const std::vector<double>& chisq, int ndof);

//! Expected mean RMSD after adding each candidate pair.
/*!
    For candidate \f$k\f$ the accumulated \f$\chi^2\f$ gains
    \f$(e_{ki} - e_{kj})^2 / \sigma^2\f$, the right-tail probability of the
    result weights the pairwise RMSDs, and the column means are averaged.

    \param[in] rmsds pairwise RMSD between frames, flat `n * n`
    \param[in] chi2 accumulated chi-squared so far, flat `n * n`
    \param[in] e_add per-frame predicted measurement of each candidate, flat `m * n`
    \param[in] inv_err_sq \f$1/\sigma^2\f$ of the new measurement
    \param[in] ndof degrees of freedom after the addition
    \param[in] diag_weight Olga's diagonal correction on the denominator
    \param[in] n_frames,n_candidates shapes
    \return one expected RMSD per candidate
*/
IMPBFFEXPORT std::vector<double> expected_rmsd_after_adding(
        const std::vector<double>& rmsds,
        const std::vector<double>& chi2,
        const std::vector<double>& e_add,
        double inv_err_sq, int ndof, double diag_weight,
        int n_frames, int n_candidates);

//! Expected mean RMSD given the chi-squared accumulated so far.
/*!
    Olga's `rmsdMeanMean(...)`: the right-tail probability of `chi2` weights the
    pairwise RMSDs, and the column means are averaged.

    \param[in] rmsds pairwise RMSD between frames, flat `n * n`
    \param[in] chi2 accumulated chi-squared, flat `n * n`
    \param[in] ndof degrees of freedom
    \param[in] diag_weight Olga's diagonal correction on the denominator
    \param[in] n_frames shape
    \return the expected mean RMSD
*/
IMPBFFEXPORT double expected_rmsd(
        const std::vector<double>& rmsds, const std::vector<double>& chi2,
        int ndof, double diag_weight, int n_frames);

//! Olga-style greedy informative pair selection.
/*!
    Repeatedly adds the candidate pair that leaves the smallest expected mean
    RMSD, then reports the precision decay over the selection.

    `predicted_measurements` is expected to be finite: Olga's GUI does its NaN filtering before
    calling the selector, and so must a caller here.

    The result is published as two managed numpy views: the selected pair
    indices in selection order, and the expected mean RMSD after each of them
    (same length). Shapes are part of the contract, stated once here rather
    than by a Python wrapper that split a struct into two arrays.

    \param[in] predicted_measurements,n_frames,n_pairs predicted measurement per frame and pair (FRET, EPR or PRE)
    \param[in] rmsds,n_rmsd_rows,n_rmsd_cols pairwise RMSD between frames
    \param[in] measurement_error expected absolute measurement error, in the prediction units
    \param[in] max_pairs how many to select; capped at `n_pairs`
    \param[in] unique_only a candidate may be selected at most once
    \param[in] diag_weight Olga's diagonal correction on the denominator
    \param[out] out_pairs,n_out_pairs selected pair indices (int view)
    \param[out] out_decay,n_out_decay expected mean RMSD after each (double view)
*/
IMPBFFEXPORT void select_probe_pairs(
        double* predicted_measurements, int n_frames, int n_pairs,
        double* rmsds, int n_rmsd_rows, int n_rmsd_cols,
        double measurement_error, int max_pairs,
        bool unique_only = true, double diag_weight = 0.99,
        int** out_pairs = 0, int* n_out_pairs = 0,
        double** out_decay = 0, int* n_out_decay = 0);

//! Olga-style greedy informative *labelling site* selection.
/*!
    For a monomer, a FRET pair is a double mutant and pairs are the unit worth
    greedy over. For a homo-oligomer it is not: the labelling mix is
    statistical, each site is mutated once, and every cross-protomer
    combination of the chosen sites is measurable — in a dimer, sites {1, 2}
    measure 1:1, 2:2, 1:2 and 2:1. A pair-wise greedy would credit each of
    those as its own double mutant and buy, per mutation, roughly half the
    experiment that is actually there.

    So this selector greedies over **sites**. Each step scores a candidate
    site by the expected mean RMSD left by *all* pairs the enlarged site set
    implies — a pair is implied when both of its sites are chosen — not by the
    single pair the monomeric selector would add. The arithmetic per pair is
    Olga's unchanged: the implied pairs' \f$\chi^2\f$ gains accumulate together,
    and the right-tail probability of the sum weights the pairwise RMSDs.

    Which pairs exist is the caller's statement about the biology, passed as
    `pair_sites`: row \f$p\f$ holds the two site indices pair \f$p\f$ connects
    (a homodimer's inter-protomer 1:2 pair is `(0, 1)`, its 1:1 pair `(0, 0)`).
    The selector is agnostic to the oligomer order — a trimer or tetramer
    simply offers more rows per site — and never re-selects a site: a residue
    cannot be mutated twice, so there is no `unique_only` to turn off.

    The result is published as two managed numpy views: the selected site
    indices in selection order, and the expected mean RMSD after each of them
    (same length). The number of *measurements* a step implies is
    \f$|\{p : \mathrm{both\ sites\ chosen}\}|\f$, which the caller can count
    from `pair_sites` directly.

    \param[in] predicted_measurements,n_frames,n_pairs predicted measurement per frame and pair (FRET, EPR or PRE)
    \param[in] rmsds,n_rmsd_rows,n_rmsd_cols pairwise RMSD between frames
    \param[in] pair_sites,n_site_rows,n_site_cols the two site indices per
               pair, `n_pairs * 2`, values in `[0, n_sites)`
    \param[in] measurement_error expected absolute measurement error, in the prediction units
    \param[in] max_sites how many sites to select; capped at the site count
    \param[in] diag_weight Olga's diagonal correction on the denominator
    \param[out] out_sites,n_out_sites selected site indices (int view)
    \param[out] out_decay,n_out_decay expected mean RMSD after each (double
                view)
*/
IMPBFFEXPORT void select_probe_positions(
        double* predicted_measurements, int n_frames, int n_pairs,
        double* rmsds, int n_rmsd_rows, int n_rmsd_cols,
        int* pair_sites, int n_site_rows, int n_site_cols,
        double measurement_error, int max_sites,
        double diag_weight = 0.99,
        int** out_sites = 0, int* n_out_sites = 0,
        double** out_decay = 0, int* n_out_decay = 0);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_PROBEPAIRSELECTION_H
