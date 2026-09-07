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
#include <IMP/bff/InteractionTerms.h>

#include <IMP/bff/Base.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE


//! A set of deactivation channels and the population that carries each.
/*!
    The package's output contract, as an object rather than two loose arrays.
    It is a **value**: copyable, comparable by content, and immutable once
    built -- amplitudes and rate constants are validated in the constructor and
    never change afterwards, so nothing downstream has to re-check them.

    Amplitudes are **not** normalised. An unnormalised spectrum still carries
    the relative answer, and normalising on construction would silently discard
    a quantum yield. Call normalized() when that is what is wanted.

    Rate constants, not lifetimes, in 1/ns. Rates are what the physics adds
    (\f$k = k_{rad} + k_{PET} + k_{FRET} + \dots\f$), a zero rate is a
    legitimate non-decaying population, and an infinite lifetime is an awkward
    number to carry around.

    \p exact records whether each species genuinely holds its rate for the whole
    excited-state lifetime -- one species per accessible-volume point, per
    rotamer, per conformer. When the dye reorganises fast enough to average over
    rates during the lifetime the decay is not a sum of exponentials at all, and
    no set of amplitudes reproduces it. Which of the two a number came from is
    part of the answer, so it is recorded rather than assumed.
*/
class IMPBFFEXPORT LifetimeSpectrum {
    std::vector<double> amplitudes_;
    std::vector<double> rate_constants_;
    bool exact_;

public:
    //! \throws IMP::ValueException on a length mismatch or a negative rate
    /*! One constructor, not two: an overloaded one cannot carry keyword
        arguments through SWIG, and the 28 lines of Python that used to
        reconstruct them by hand are what this default-argument form
        replaces. Default-constructed is the empty spectrum the value-vector
        template needs. */
    LifetimeSpectrum(const std::vector<double>& amplitudes =
                             std::vector<double>(),
                     const std::vector<double>& rate_constants =
                             std::vector<double>(),
                     bool exact = true);

    //! Population fraction per species. A numpy view over the object's buffer.
    void get_amplitudes(double** out_view, int* n_out_view) const;
    //! Total deactivation rate per species, 1/ns.
    void get_rate_constants(double** out_view, int* n_out_view) const;
    //! \f$1/k\f$ in ns; `inf` for a non-decaying species.
    void get_lifetimes(double** out_view, int* n_out_view) const;

    bool get_exact() const { return exact_; }
    int get_n_species() const { return static_cast<int>(amplitudes_.size()); }
    double get_total_amplitude() const;

    //! \f$\langle\tau\rangle_x = \sum a_i \tau_i / \sum a_i\f$ -- over *molecules*.
    /*! This is the one that enters a FRET efficiency, because efficiency is
        defined per molecule. Zero for an empty population; infinite if any
        species does not decay. */
    double get_species_averaged_lifetime() const;

    //! \f$\langle\tau\rangle_f = \sum a_i \tau_i^2 / \sum a_i \tau_i\f$ -- over *photons*.
    /*! What a long-lived species dominates. Equal to the species average only
        for a single exponential; the gap between them measures the sample's
        heterogeneity, and reporting one under the other's name is a common way
        to be wrong by tens of percent. */
    double get_intensity_averaged_lifetime() const;

    //! The same spectrum with amplitudes summing to 1 (unchanged if empty).
    LifetimeSpectrum normalized() const;

    //! \f$F(t) = \sum_i a_i e^{-k_i t}\f$ on the caller's axis. **Unconvolved.**
    void decay(const std::vector<double>& time,
               double** out_view, int* n_out_view) const;

    //! Reduce many species to few, preserving \f$\sum a\f$ and \f$\sum a k\f$.
    /*! The result is marked inexact unless nothing was actually merged. */
    LifetimeSpectrum coarse_grain(int n_bins = 128) const;

    IMP_SHOWABLE_INLINE(LifetimeSpectrum,
                        out << "LifetimeSpectrum(" << get_n_species()
                            << " species, <tau>_x = "
                            << get_species_averaged_lifetime() << " ns, "
                            << (exact_ ? "exact" : "approximate") << ")");
};
IMP_VALUES(LifetimeSpectrum, LifetimeSpectrums);

//! \f$E = 1 - \langle\tau\rangle_{x,DA} / \langle\tau\rangle_{x,D}\f$.
/*!
    The **species** average, not the intensity average. Efficiency is defined
    per molecule -- the fraction of excitations that transfer -- and the
    intensity average weights each molecule by how many photons it emitted,
    which is exactly the thing transfer suppresses. Using it gives an efficiency
    that is systematically too low, by tens of percent on a heterogeneous sample.

    \param[in] donor_only the donor without an acceptor
    \param[in] donor_acceptor the same donor with one
    \return efficiency; 0 if the donor-only lifetime is zero or not finite
*/
//! One species per rate constant.
/*!
    The static limit taken literally: if a population of weight \f$w_i\f$ decays
    at \f$k_i\f$ and keeps that rate, then \f$F(t) = \sum_i w_i e^{-k_i t}\f$ is
    the decay, with no fitting and no approximation. Whether that premise holds
    is the caller's to know, and `exact` records the answer.

    \param[in] rates total deactivation rate per species, 1/ns
    \param[in] weights population per species; uniform when empty
    \param[in] exact see LifetimeSpectrum
*/
IMPBFFEXPORT LifetimeSpectrum lifetime_spectrum_from_rates(
        const std::vector<double>& rates,
        const std::vector<double>& weights = std::vector<double>(),
        bool exact = true);

//! The spectrum of a state ensemble under a set of interaction terms.
/*!
    One species per state -- per accessible-volume point, per rotamer, per
    conformer -- carrying the summed rate of every channel acting on it. This is
    the reduction that connects a representation and the photophysics to an
    experiment-neutral answer, and it is representation-agnostic for the same
    reason the terms are: it consumes states.

    \warning This is the **static** limit. It is exact when each state holds its
    rate for the whole excited-state lifetime, and wrong when the dye
    reorganises fast enough to average over rates -- then the decay is not a sum
    of exponentials at all and the population has to be propagated instead
    (#IMP::bff::GridDiffusionSolver, or the Brownian walk). Pass `exact = false`
    when using it outside that limit, so the spectrum says what it is.

    \param[in] terms the interaction terms to sum
    \param[in] first the dye whose states are being rated; its weights are the
               populations unless `weights` overrides them
    \param[in] second the second participant, for the 2-body terms
    \param[in] weights population per state; the first participant's own
               weights when empty, which is what an accessible volume's
               occupancy already is
    \param[in] exact see LifetimeSpectrum
*/
IMPBFFEXPORT LifetimeSpectrum lifetime_spectrum_from_states(
        const InteractionTerms& terms, const States& first,
        const States& second = IMP::bff::States(),
        const std::vector<double>& weights = std::vector<double>(),
        bool exact = true);

IMPBFFEXPORT double fret_efficiency_from_lifetimes(
        const LifetimeSpectrum& donor_only,
        const LifetimeSpectrum& donor_acceptor);

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

//! \f$R_{app} = R_{DA} \cdot (R_{app}/R_{DA})\f$, as a multiplicative convolution.
/*!
    The front door over #outer_product_histogram a FRET caller actually uses:
    bin \p r_da against \p r_ratio, but report only the bins that carry weight
    (the threshold is \f$10^{-10}\f$ of the peak) as `(centres, hist)` on one
    floor with the range:
    \f$[\min R_{DA} \cdot \min r, \max R_{DA} \cdot \max r]\f$. A range of zero
    width is widened by half a bin on each side, as numpy's
    `np.histogram`-with-a-point-range contract does. Empty input returns two
    empty views.

    \param[in] r_da,amp the donor distance distribution and its amplitudes
    \param[in] r_ratio,weights_ratio the \f$R_{app}/R_{DA}\f$ distribution
    \param[in] n_bins bins on the output axis
    \param[out] out_centres,n_centres bin centres that carry weight
    \param[out] out_hist,n_hist the corresponding weights
*/
IMPBFFEXPORT void convolve_distance_with_k2_ratio(
        const std::vector<double>& r_da, const std::vector<double>& amp,
        const std::vector<double>& r_ratio,
        const std::vector<double>& weights_ratio, int n_bins,
        double** out_centres, int* n_centres, double** out_hist, int* n_hist);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_LIFETIMESPECTRUM_H
