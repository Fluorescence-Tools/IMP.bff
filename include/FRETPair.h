/**
 *  \file IMP/bff/FRETPair.h
 *  \brief Distances, orientation factors and efficiencies over all pairs of two
 *         labelled ensembles.
 *
 * Two rotamer libraries, two accessible-volume clouds, or one of each: the
 * observable is a sum over every donor state against every acceptor state. That
 * is an \f$N_1 \times N_2\f$ problem, and the vectorised form pays for it twice
 * — once in an \f$(N_1, N_2, 3)\f$ array of separation vectors, and again in the
 * five full-size temporaries the efficiency expression walks through.
 *
 * These compute the same matrices in one pass each.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_FRETPAIR_H
#define IMPBFF_FRETPAIR_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Separation and \f$\kappa^2\f$ for every donor/acceptor state pair.
/*!
    With transition dipoles, \f$\kappa^2\f$ is
    \f$(\hat\mu_D\cdot\hat\mu_A - 3(\hat\mu_D\cdot\hat r)(\hat\mu_A\cdot\hat r))^2\f$.
    Without them it is the isotropic 2/3 everywhere, which is what an accessible
    volume can say: a point cloud carries no orientation.

    The separation vectors are formed one at a time and discarded. The Python
    this replaces built the whole \f$(N_1, N_2, 3)\f$ array and returned it in
    its result dictionary, where **nothing ever read it** — it was the largest
    allocation in the call and it existed for no consumer.

    \param[in] points1,n_points1 donor state centres, flat `n1 * 3`
    \param[in] points2,n_points2 acceptor state centres, flat `n2 * 3`
    \param[in] mu1,n_mu1 donor transition dipoles, flat `n1 * 3`, or length 0
    \param[in] mu2,n_mu2 acceptor dipoles, or length 0
    \param[in] n1,n2 state counts
    \param[out] out_view,n_out_view two concatenated `n1 * n2` blocks: the
                separations, then \f$\kappa^2\f$. Handed back as a **numpy view
                over the kernel's own buffer**, allocated with `malloc` because
                numpy releases it with `free`. The parameter is named `out_view`
                rather than `output` on purpose: `output` is claimed by both the
                managed and the unmanaged typemap, and binding the unmanaged one
                would leak the whole array on every call. On allocation failure
                the view comes back empty rather than as a null pointer. Returning a `std::vector` instead makes SWIG build one
                Python float per element and numpy walk them back: **66 ns each**,
                which on a 400x350 pair matrix is 18 ms of a 20 ms call against
                about 1 ms of arithmetic.

    Dipoles need not be normalised; they are normalised here, and a zero-length
    one contributes zero to every dot product.
*/
IMPBFFEXPORT void fret_pair_matrices(
        double* points1, int n_points1,
        double* points2, int n_points2,
        double* mu1, int n_mu1,
        double* mu2, int n_mu2,
        int n1, int n2,
        double** out_view, int* n_out_view);

//! Per-pair transfer efficiency and rate ratio.
/*!
    \f$E = 1/(1 + \tfrac{2}{3}(R/R_0)^6/\kappa^2)\f$ and
    \f$k_{FRET}/k_{rad} = \tfrac{3}{2}\kappa^2 (R_0/R)^6\f$.

    The degenerate cases are kept exactly as the Python had them, because they
    are reachable: coincident states make \f$(R/R_0)^6\f$ zero, and an
    orthogonal dipole pair makes \f$\kappa^2\f$ zero, so the efficiency can come
    out `nan` or `+inf`. Both mean complete transfer and both become 1.

    \param[in] r,n_r separations, flat `n1 * n2`
    \param[in] kappa2,n_kappa2 orientation factors, same shape
    \param[in] forster_radius \f$R_0\f$ at \f$\kappa^2 = 2/3\f$
    \return two concatenated blocks: the efficiencies, then the rate ratios.
            The rate ratios are **not** sanitised — an infinite rate at zero
            separation is true, and the caller decides what to do with it.
*/
IMPBFFEXPORT std::vector<double> fret_pair_efficiency_matrices(
        double* r, int n_r,
        double* kappa2, int n_kappa2,
        double forster_radius);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_FRETPAIR_H
