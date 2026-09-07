/**
 *  \file IMP/bff/PolymerChain.h
 *  \brief End-to-end distance distributions of ideal and worm-like chains.
 *
 * The linker between an attachment point and a dye is a short polymer, and its
 * end-to-end distribution is what an accessible volume approximates
 * geometrically. These give it analytically instead.
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

//! End-to-end distance distribution of an Ising two-state Gaussian chain.
/*!
    The tractable Ising-worm-like-chain form for a partially structured chain:
    every residue is **structured** (S) or **unstructured** (U) under a
    nearest-neighbour Ising Hamiltonian -- \p coupling is the cooperativity
    \f$J\f$, \p field the bias \f$h\f$ towards S -- and each residue
    contributes a Gaussian bond of mean-square extension \f$b_S^2\f$ or
    \f$b_U^2\f$.

    Because the bonds are Gaussian the characteristic function factorises per
    residue, so the Boltzmann-weighted end-to-end \f$\varphi(k)\f$ is an exact
    2x2 **transfer-matrix product**,

    \f$\varphi(k) = \mathbf{1}^\top[\prod_i M_i(k)]\mathbf{1} / \varphi(0)\f$,
    \f$M(k)_{\sigma\sigma'} = W(\sigma,\sigma')e^{-k^2b_{\sigma'}^2/6}\f$,

    and \f$P(R)\f$ follows from the isotropic inverse transform
    \f$P(R) = (2R/\pi)\int_0^\infty k\sin(kR)\varphi(k)\,dk\f$, taken on a
    \p n_k -point trapezoidal grid. With \f$b_S = b_U\f$ it reduces to
    gaussian_chain().

    The product is formed **as written**, without rescaling between residues.
    That is deliberate rather than careless: \f$\varphi\f$ is used only as the
    ratio \f$\varphi(k)/\varphi(0)\f$, so a per-residue rescaling would cancel
    only if it were the same at every \p k -- and a per-\p k rescaling silently
    changes the answer. The cost is that a long chain with a large \p coupling
    can overflow to a non-finite \f$\varphi\f$; those entries end up zero, as
    they do in the numpy original this reproduces.

    \param[in] distances the r axis
    \param[in] number_of_residues residues (bonds) between the dyes
    \param[in] b_structured RMS bond contribution of a structured residue
    \param[in] b_unstructured RMS bond contribution of an unstructured residue
    \param[in] coupling Ising nearest-neighbour coupling \f$J\f$
    \param[in] field Ising field \f$h\f$; positive biases towards structured
    \param[in] n_k number of k points in the inverse transform
    \param[out] out_view,n_out_view the normalised distribution as a view
*/
IMPBFFEXPORT void ising_chain(
        const std::vector<double>& distances,
        int number_of_residues,
        double b_structured,
        double b_unstructured,
        double coupling = 1.5,
        double field = 0.0,
        int n_k = 2000,
        double** out_view = 0, int* n_out_view = 0
);

//! Radial distribution of a self-avoiding walk with Flory exponent \p nu.
/*!
    The des Cloizeaux form (Zheng et al., JACS 2018), the standard model for
    disordered and unfolded chains in single-molecule FRET:

    \f$P(r) \propto r^{2+\theta}\exp[-(r/r_0)^{\delta}]\f$, with
    \f$\theta = (\gamma-1)/\nu\f$ and \f$\delta = 1/(1-\nu)\f$, and the
    scale \f$r_0\f$ fixed so that \f$\sqrt{\langle r^2\rangle}\f$ is
    \p r_rms. With \p nu = 0.5 and \p gamma_exp = 1 it reduces to
    gaussian_chain().

    \param[in] distances the r axis, positive
    \param[in] r_rms target root-mean-square inter-dye distance
    \param[in] nu Flory exponent; ~0.588 expanded, 0.5 theta, < 0.4 collapsed
    \param[in] gamma_exp SAW susceptibility exponent; 1.0 is ideal
    \param[out] out_view,n_out_view the normalised density as a managed view
*/
IMPBFFEXPORT void saw_nu(
        const std::vector<double>& distances,
        double r_rms,
        double nu = 0.588,
        double gamma_exp = 1.1615,
        double** out_view = 0, int* n_out_view = 0
);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_POLYMERCHAIN_H
