/**
 *  \file IMP/bff/PolymerChain.h
 *  \brief End-to-end distance distributions of ideal and worm-like chains.
 *
 * The linker between an attachment point and a dye is a short polymer, and its
 * end-to-end distribution is what an accessible volume approximates
 * geometrically. These give it analytically instead.
 *
 * Ported from Python by PRD-113: numba is a prototyping tool in this package,
 * not a runtime dependency, so every numerical kernel is C++.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_POLYMERCHAIN_H
#define IMPBFF_POLYMERCHAIN_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! RMS end-to-end distance of an ideal (Gaussian) chain, \f$b\sqrt{N}\f$.
/*!
    \param[in] segment_length Kuhn segment length in Angstrom
    \param[in] number_of_segments number of Kuhn segments
*/
IMPBFFEXPORT double gaussian_chain_ree(
        double segment_length, int number_of_segments);

//! Radial distribution of an ideal chain.
/*!
    \f$P(r) = 4\pi r^2 (3/2\pi\langle r^2\rangle)^{3/2}\exp(-3r^2/2\langle r^2\rangle)\f$

    Not normalised on the given axis unless that axis covers the mass.

    \param[in] distances the r axis
    \param[in] segment_length Kuhn segment length
    \param[in] number_of_segments number of Kuhn segments
*/
IMPBFFEXPORT void gaussian_chain(
        const std::vector<double>& distances,
        double segment_length,
        int number_of_segments,
        double** out_view, int* n_out_view
);

//! Radial distribution of a worm-like chain.
/*!
    The multi-piece analytical solution of Becker, Rosa and Everaers
    (Eur Phys J E 32:53-69, 2010). \f$\kappa\f$ is the dimensionless
    persistence-length ratio, and the expression branches at
    \f$\kappa = 0.125\f$.

    Values at or beyond the contour length stay zero: the chain cannot be longer
    than itself, and the closed form diverges there.

    \param[in] distances the r axis
    \param[in] kappa dimensionless persistence length
    \param[in] chain_length contour length; 0 takes the largest r on the axis
    \param[in] normalize divide by the sum
    \param[in] distance multiply by \f$r^2\f$, giving a distance distribution
        rather than a density in space
    \param[out] out_view,n_out_view the distribution as a managed view
*/
IMPBFFEXPORT void worm_like_chain(
        const std::vector<double>& distances,
        double kappa,
        double chain_length = 0.0,
        bool normalize = true,
        bool distance = true,
        double** out_view = 0, int* n_out_view = 0
);

//! Worm-like chain broadened by the dye linkers at each end.
/*!
    Convolves :func:`worm_like_chain` with a Gaussian of width \p sigma, which
    stands for the finite reach and flexibility of the linkers -- the chain
    distribution is between the *attachment points*, and what is measured is
    between the *dyes*.

    \param[in] distances the r axis
    \param[in] kappa dimensionless persistence length
    \param[in] chain_length total length; 0 takes the largest r on the axis
    \param[in] sigma linker broadening in Angstrom
    \param[in] normalize divide by the sum
    \param[out] out_view,n_out_view the broadened ring as a managed view
*/
IMPBFFEXPORT void worm_like_chain_linker(
        const std::vector<double>& distances,
        double kappa,
        double chain_length = 0.0,
        double sigma = 6.0,
        bool normalize = true,
        double** out_view = 0, int* n_out_view = 0
);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_POLYMERCHAIN_H
