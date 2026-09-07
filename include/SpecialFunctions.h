/**
 *  \file IMP/bff/SpecialFunctions.h
 *  \brief Special functions the distance distributions need.
 *
 * Only what the polymer distributions call for. In particular the modified
 * Bessel function \f$I_0\f$, which the worm-like chain needs and which is
 * **deliberately the polynomial approximation** rather than the exact
 * function -- see i0().
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_SPECIALFUNCTIONS_H
#define IMPBFF_SPECIALFUNCTIONS_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Modified Bessel function \f$I_0(x)\f$, for any real \p x.
/*!
    The Abramowitz & Stegun polynomial approximation (Numerical Recipes'
    `bessi0`), correct to about 1e-7 -- **not** the exact function, and that is
    the point rather than a shortcut.

    `std::cyl_bessel_i(0, x)` and `scipy.special.i0` are correct to machine
    precision, so substituting either would move every fitted worm-like-chain
    distribution by more than the optimiser's tolerance and quietly invalidate
    published fits. ChiSurf kept the polynomial for that reason; this is the
    same nine coefficients in the same order, so the two agree exactly rather
    than approximately.

    \f$I_0\f$ is **even**, which is the property the worm-like chain depends
    on: its argument is negative, and \f$I_0(-x) = I_0(x)\f$ grows where
    \f$\exp(-x)\f$ decays. Transcribing this as `exp` is not a small error --
    see the note on worm_like_chain().

    \param[in] x the argument
*/
IMPBFFEXPORT double i0(double x);

//! i0() over an axis.
/*!
    \param[in] x the arguments
    \param[out] out_view,n_out_view \f$I_0\f$ at each point, as a managed view
*/
IMPBFFEXPORT void i0_array(
        const std::vector<double>& x,
        double** out_view = 0, int* n_out_view = 0
);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_SPECIALFUNCTIONS_H
