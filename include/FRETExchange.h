/**
 *  \file IMP/bff/FRETExchange.h
 *  \brief FRET when the labels exchange between states during the lifetime.
 *
 * The usual averages assume one of two limits: either each state holds its rate
 * for the whole excited-state lifetime (*static*), or it reorganises fast enough
 * to average over them (*dynamic*). Between the two, the answer depends on the
 * competition between fluorescence decay and conformational transition, and
 * neither limit is right.
 *
 * #IMP::bff::fret_efficiency_exact_kinetic solves that competition exactly, by
 * integrating the master equation of the labelled population; the three limits
 * are in #IMP::bff::fret_efficiency_regimes for comparison.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_FRETEXCHANGE_H
#define IMPBFF_FRETEXCHANGE_H

#include <IMP/bff/bff_config.h>

#include <IMP/bff/IMPCompatibility.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! The three averaging limits of a FRET pair, side by side.
struct IMPBFFEXPORT FRETRegimes {
    //! Each state keeps its rate for the whole lifetime: average the efficiency.
    double static_efficiency;
    //! Fast enough to average \f$\kappa^2\f$, then convert once.
    double dynamic;
    //! Fast enough to average the *rate*, then convert once.
    double dynamic_plus;
    double kappa2_avg;

    FRETRegimes()
        : static_efficiency(0), dynamic(0), dynamic_plus(0),
          kappa2_avg(2.0 / 3.0) {}

    IMP_SHOWABLE_INLINE(FRETRegimes,
                        out << "FRETRegimes(static=" << static_efficiency
                            << ", dynamic=" << dynamic << ")");
};
IMP_VALUES(FRETRegimes, FRETRegimesList);

//! The exact FRET efficiency of an exchanging ensemble.
/*!
    Solves the master equation of the labelled population rather than assuming a
    limit: with \f$P\f$ the transition probability matrix over \p dt, the rate
    matrix is \f$M = (P - I)/dt\f$ and the integrated populations
    \f$G_j = \int_0^\infty p_j\,dt\f$ solve
    \f$(k_{rad} I + K_{FRET} - M)^T G = w\f$. The efficiency is
    \f$\sum_i k_{FRET,i} G_i\f$ — the fraction of excitations that leave by
    transfer.

    **The transpose is load-bearing.** `P[i][j]` is the probability of
    \f$i \to j\f$, so populations evolve as a *row* vector and \f$G\f$ solves the
    transposed system. Solving \f$A G = w\f$ instead is right only for a
    symmetric \f$M\f$, and for a real transition matrix it gives the wrong
    fast-exchange limit.

    \param[in] p_matrix transition probabilities over \p dt, flat `n * n`,
               rows summing to 1
    \param[in] fret_rates per state, 1/ns
    \param[in] tau0 donor lifetime without acceptor, ns
    \param[in] dt the time step `p_matrix` was measured over, ns
    \param[in] weights the stationary distribution; empty takes the eigenvector
               of \f$P^T\f$ at eigenvalue 1
*/
IMPBFFEXPORT double fret_efficiency_exact_kinetic(
        const std::vector<double>& p_matrix,
        const std::vector<double>& fret_rates, double tau0, double dt = 1.0,
        const std::vector<double>& weights = std::vector<double>());

//! The exact efficiency for two exchanging ensembles.
/*!
    The joint kinetics is the Kronecker product of the two transition matrices,
    which is \f$(n_d n_a)^2\f$ — expensive for large libraries and the reason
    this is not the default.

    \param[in] dist_matrix,kappa2_matrix flat `nd * na`
    \param[in] p_d,p_a transition matrices, flat `nd * nd` and `na * na`
    \param[in] weights_d,weights_a stationary distributions
    \param[in] forster_radius \f$R_0\f$ at \f$\kappa^2 = 2/3\f$, A
    \param[in] tau0 donor lifetime, ns
    \param[in] dt the transition time step, ns
*/
IMPBFFEXPORT double fret_efficiency_exact_kinetic_pair(
        const std::vector<double>& dist_matrix,
        const std::vector<double>& kappa2_matrix,
        const std::vector<double>& p_d, const std::vector<double>& p_a,
        const std::vector<double>& weights_d,
        const std::vector<double>& weights_a, double forster_radius = 52.0,
        double tau0 = 4.0, double dt = 0.1);

//! The three limits, for a pair of weighted state sets.
/*!
    \param[in] dist_matrix,kappa2_matrix flat `n_donor * n_acceptor`
    \param[in] weights_d,weights_a per-state weights
    \param[in] forster_radius \f$R_0\f$ at \f$\kappa^2 = 2/3\f$, A
*/
IMPBFFEXPORT FRETRegimes fret_efficiency_regimes(
        const std::vector<double>& dist_matrix,
        const std::vector<double>& kappa2_matrix,
        const std::vector<double>& weights_d,
        const std::vector<double>& weights_a, double forster_radius = 52.0);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_FRETEXCHANGE_H
