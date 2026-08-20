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

#include <IMP/value_macros.h>
#include <IMP/showable_macros.h>

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
                Python float per element and numpy walk them back: **35-40 ns each**,
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
    \param[out] out_view,n_out_view two concatenated blocks: the efficiencies,
                then the rate ratios. The rate ratios are **not** sanitised — an infinite rate at zero
            separation is true, and the caller decides what to do with it.
*/
IMPBFFEXPORT void fret_pair_efficiency_matrices(
        double* r, int n_r,
        double* kappa2, int n_kappa2,
        double forster_radius,
        double** out_view, int* n_out_view);

//! The pair geometry of two labelled ensembles: every (i, j) at once.
/*!
    `R` and `kappa2` are `n1 * n2`, row-major; `weight` is the normalised outer
    product of the two weight vectors. `kappa2_avg` is the weight-averaged
    orientation factor.

    The separation vectors used to come back too, as `(n1, n2, 3)`. **Nothing
    ever read them** -- they were the largest allocation in the call, built for
    no consumer -- so each one is formed and discarded.
*/
struct IMPBFFEXPORT FRETPairGeometry {
    std::vector<double> R, kappa2, weight;
    double kappa2_avg;
    int n1, n2;

    FRETPairGeometry() : kappa2_avg(2.0 / 3.0), n1(0), n2(0) {}

    void get_R(double** out_view, int* n_out_view) const;
    void get_kappa2(double** out_view, int* n_out_view) const;
    void get_weight(double** out_view, int* n_out_view) const;

    IMP_SHOWABLE_INLINE(FRETPairGeometry,
                        out << "FRETPairGeometry(" << n1 << " x " << n2 << ")");
};
IMP_VALUES(FRETPairGeometry, FRETPairGeometries);

//! What a pair geometry transfers, in the three averaging limits.
/*!
    `static_` averages the per-pair efficiency; `dynamic1` uses the ensemble
    \f$\langle\kappa^2\rangle\f$ and then averages; `dynamic2` averages the
    *rate* and converts once. They differ by which quantity the ensemble is
    fast compared with, and reporting one under another's name is a common way
    to be wrong by tens of percent.
*/
struct IMPBFFEXPORT FRETPairEfficiencies {
    double static_efficiency, dynamic1, dynamic2, kappa2_avg;
    //! \f$R_0\f$ the efficiencies were computed at, A.
    double forster_radius;
    std::vector<double> E, rate_ratio, k_fret;
    //! The geometry these came from, carried so the value is self-contained --
    //! every consumer reads `R` and `kappa2` beside `E`, and holding the
    //! geometry separately is how the two came to be passed around in one dict.
    std::vector<double> R, kappa2, weight;
    int n1, n2;

    FRETPairEfficiencies()
        : static_efficiency(0), dynamic1(0), dynamic2(0),
          kappa2_avg(2.0 / 3.0), forster_radius(0), n1(0), n2(0) {}

    void get_E(double** out_view, int* n_out_view) const;
    void get_rate_ratio(double** out_view, int* n_out_view) const;
    //! Per-pair FRET rate, 1/ns; empty unless a lifetime was given.
    void get_k_fret(double** out_view, int* n_out_view) const;
    void get_R(double** out_view, int* n_out_view) const;
    void get_kappa2(double** out_view, int* n_out_view) const;
    void get_weight(double** out_view, int* n_out_view) const;

    IMP_SHOWABLE_INLINE(FRETPairEfficiencies,
                        out << "FRETPairEfficiencies(E=" << static_efficiency
                            << ")");
};
IMP_VALUES(FRETPairEfficiencies, FRETPairEfficienciesList);

//! Distances, \f$\kappa^2\f$ and pair weights over all (i, j).
/*!
    \param[in] points1,points2 flat `(n, 3)` centres, A
    \param[in] weights1,weights2 one per point, normalised or not
    \param[in] mu1,mu2 flat `(n, 3)` transition dipoles; empty means the
               isotropic 2/3 everywhere, which is what an AV cloud has
*/
IMPBFFEXPORT FRETPairGeometry fret_pair_geometry(
        const std::vector<double>& points1, const std::vector<double>& weights1,
        const std::vector<double>& points2, const std::vector<double>& weights2,
        const std::vector<double>& mu1 = std::vector<double>(),
        const std::vector<double>& mu2 = std::vector<double>());

//! The efficiencies of a pair geometry.
/*!
    \param[in] geometry from #IMP::bff::fret_pair_geometry
    \param[in] forster_radius \f$R_0\f$ for \f$\kappa^2 = 2/3\f$, same units
    \param[in] tau0 donor lifetime, ns; negative leaves `k_fret` empty
*/
IMPBFFEXPORT FRETPairEfficiencies fret_pair_efficiencies(
        const FRETPairGeometry& geometry, double forster_radius,
        double tau0 = -1.0);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_FRETPAIR_H
