/**
 *  \file IMP/bff/LifetimeSpectrum.h
 *  \brief The experiment-neutral output: (amplitude, rate constant) pairs.
 *
 * A forward model's answer is a set of deactivation rate constants and how much
 * population carries each. Everything an instrument adds -- convolution with an
 * IRF, pileup, counting statistics, a TAC's bin width, a background -- is
 * applied to that, outside this package. Emitting a binned curve instead bakes
 * a bin width and a time range into the model's output and makes the
 * instrument's job ambiguous.
 *
 * A spectrum is exact whenever each species keeps its rate for the whole
 * excited-state lifetime -- one entry per accessible-volume point, per rotamer,
 * per conformer. When the dye moves fast enough to average over rates the decay
 * is not a sum of exponentials at all, and a spectrum is then the static
 * approximation to it; see #DiffusionSolver for the case where that fails.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_LIFETIMESPECTRUM_H
#define IMPBFF_LIFETIMESPECTRUM_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Evaluate \f$F(t) = \sum_i a_i e^{-k_i t}\f$ -- **unconvolved**.
/*!
    No instrument response, no pileup, no counting noise: this is the model's
    curve, and folding an IRF into it belongs to whatever owns the instrument.

    \param[in] amplitudes per species
    \param[in] rate_constants per species, 1/ns
    \param[in] time the abscissa, ns
    \return \p time.size() values
*/
IMPBFFEXPORT void lifetime_spectrum_decay(
        const std::vector<double>& amplitudes,
        const std::vector<double>& rate_constants,
        const std::vector<double>& time,
        double** out_view, int* n_out_view);

//! Reduce many species to few, preserving the first two moments exactly.
/*!
    An accessible volume gives one species per point -- tens of thousands -- and
    almost none of them are distinguishable. Species are binned on a uniform
    grid in the rate constant; each bin returns the summed amplitude and the
    **amplitude-weighted mean rate** of its members.

    That choice makes \f$\sum_i a_i\f$ and \f$\sum_i a_i k_i\f$ exact, so the
    total population and the initial slope \f$F'(0)\f$ are preserved whatever
    \p n_bins is. Higher moments are not: the curvature is under-stated by an
    amount that falls as the square of the bin width, so the reduction is a
    genuine approximation and the caller chooses how much of one.

    \param[in] amplitudes per species
    \param[in] rate_constants per species, 1/ns
    \param[in] n_bins bins to distribute the rate range over
    \param[out] out_rates the amplitude-weighted mean rate of each non-empty bin
    \return the summed amplitude of each non-empty bin; empty bins are dropped,
            so the result may be shorter than \p n_bins
*/
IMPBFFEXPORT std::vector<double> lifetime_spectrum_coarse_grain(
        const std::vector<double>& amplitudes,
        const std::vector<double>& rate_constants,
        int n_bins,
        std::vector<double>& out_rates);

//! Weighted histogram of the outer product of two distributions.
/*!
    The distribution of \f$a_i b_j\f$ with weight \f$w_i v_j\f$ -- what you
    need when one quantity multiplies another and you want the result's
    histogram. Convolving a distance distribution with a \f$\kappa^2\f$
    distance-ratio distribution is exactly this: every \f$R_{DA}\f$ times every
    ratio, weighted by both.

    Done without forming the outer product. That matters at the sizes this is
    reached with: the array version allocates \f$n \times m\f$ twice, once for
    the products and once for the weights, purely to hand them to a histogram.

    Bin edges are uniform from \p lo to \p hi. Values outside are dropped; the
    top edge is closed, so a value exactly at \p hi lands in the last bin rather
    than nowhere.

    \param[in] a,weights_a the first distribution and its weights
    \param[in] b,weights_b the second
    \param[in] n_bins bins between \p lo and \p hi
    \param[in] lo,hi histogram range
    \return \p n_bins accumulated weights
*/
IMPBFFEXPORT std::vector<double> outer_product_histogram(
        const std::vector<double>& a, const std::vector<double>& weights_a,
        const std::vector<double>& b, const std::vector<double>& weights_b,
        int n_bins, double lo, double hi);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_LIFETIMESPECTRUM_H
