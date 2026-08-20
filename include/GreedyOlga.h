/**
 *  \file IMP/bff/GreedyOlga.h
 *  \brief Greedy Olga — which FRET pair to measure next.
 *
 * Experiment planning, after Olga: given an ensemble of candidate structures and
 * a set of measurable pairs, repeatedly add the pair that most reduces the
 * expected RMSD between the true structure and the one the measurements would
 * pick out. Greedy because each step takes the locally best pair and never
 * revisits it — which is what makes it tractable, and what makes the scoring
 * step the whole cost.
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
#ifndef IMPBFF_GREEDYOLGA_H
#define IMPBFF_GREEDYOLGA_H

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
    \param[in] e_add per-frame efficiency of each candidate, flat `m * n`
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

//! What a greedy selection leaves behind: the pairs, and the precision decay.
struct IMPBFFEXPORT GreedyPairSelection {
    //! Indices into the pair axis of `effs`, in selection order.
    std::vector<int> pairs;
    //! Expected mean RMSD after each step, same length as `pairs`.
    std::vector<double> decay;
};

//! Olga-style greedy informative pair selection.
/*!
    Repeatedly adds the candidate pair that leaves the smallest expected mean
    RMSD, then reports the precision decay over the selection.

    `effs` is expected to be finite: Olga's GUI does its NaN filtering before
    calling the selector, and so must a caller here.

    \param[in] effs,n_frames,n_pairs FRET efficiency per frame and pair
    \param[in] rmsds,n_rmsd_rows,n_rmsd_cols pairwise RMSD between frames
    \param[in] err expected absolute error in FRET efficiency
    \param[in] max_pairs how many to select; capped at `n_pairs`
    \param[in] unique_only a candidate may be selected at most once
    \param[in] diag_weight Olga's diagonal correction on the denominator
    \return the selected pairs and the expected mean RMSD after each of them
*/
IMPBFFEXPORT GreedyPairSelection select_informative_pairs(
        double* effs, int n_frames, int n_pairs,
        double* rmsds, int n_rmsd_rows, int n_rmsd_cols,
        double err, int max_pairs,
        bool unique_only = true, double diag_weight = 0.99);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_GREEDYOLGA_H
