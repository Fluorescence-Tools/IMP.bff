#ifndef IMPBFF_FRET_H
#define IMPBFF_FRET_H

/**
 *  \file IMP/bff/FRET.h
 *  \brief FRET between two labels: every pair of two clouds, and the rate along one walk.
 *
 * Two former headers: the **pair** (formerly `FRETPair.h`) -- distances,
 * orientation factors and efficiencies over all pairs of two point clouds --
 * and the **rate trace** (formerly `FRETRateTrace.h`) -- the FRET rate along
 * a dye trajectory, from the orientation factor it carries.
 *
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

// -------- from FRETPair.h --------
/**
 *  (formerly IMP/bff/FRETPair.h, now a section of this file)
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
#include <IMP/bff/bff_config.h>

#include <IMP/bff/Base.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Separation and \f$\kappa^2\f$ for every donor/acceptor state pair.
/*!
    With transition dipoles, \f$\kappa^2\f$ is
    \f$(\hat\mu_D\cdot\hat\mu_A - 3(\hat\mu_D\cdot\hat r)(\hat\mu_A\cdot\hat r))^2\f$.
    Without them it is the isotropic 2/3 everywhere, which is what an accessible
    volume can say: a point cloud carries no orientation.

    The separation vectors are formed one at a time and discarded, rather than
    materialised as an \f$(N_1, N_2, 3)\f$ array: that array is the largest
    allocation in the call and no consumer of the result reads it.

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

    The degenerate cases are handled explicitly, because they are reachable:
    coincident states make \f$(R/R_0)^6\f$ zero, and an
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

//! `fret_pair_geometry` followed by `fret_pair_efficiencies`, in one call.
/*! Convenience for a caller who wants the efficiencies and does not need the
    intermediate geometry. */
IMPBFFEXPORT FRETPairEfficiencies fret_pair_distribution(
        const std::vector<double>& points1, const std::vector<double>& weights1,
        const std::vector<double>& points2, const std::vector<double>& weights2,
        double forster_radius,
        const std::vector<double>& mu1 = std::vector<double>(),
        const std::vector<double>& mu2 = std::vector<double>(),
        double tau0 = -1.0);

IMPBFF_END_NAMESPACE

 //IMPBFF_FRETPAIR_H

// -------- from FRETRateTrace.h --------
/**
 *  (formerly IMP/bff/FRETRateTrace.h, now a section of this file)
 *  \brief FRET rate along a dye trajectory.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/OrientationFactor.h>


IMPBFF_BEGIN_NAMESPACE

//! FRET rate per frame, donor trajectory against a static acceptor cloud.
/*!
    **The arithmetic mean of rates**, averaged over the acceptor cloud -- the
    fast-exchange limit, in which the acceptor re-randomises within the donor's
    excited-state lifetime. ``fret_map`` accumulates the mean transfer *time*
    instead and inverts it, which is the static limit. Different physics, not
    two spellings of one thing.

    \param[in] trajectory donor positions, flat, three per frame
    \param[in] acceptor_points acceptor cloud, flat, three per point
    \param[in] R0 Forster radius
    \param[in] tau0 unquenched donor lifetime
    \param[in] r_min2 squared closest approach; pairs nearer than this are
        clamped, because \f$1/r^6\f$ diverges and two dyes cannot interpenetrate
    \param[in] kappa2_scale multiplier on \f$R_0^6\f$ from the orientation factor
    \param[out] out_view,n_out_view the rate at each frame, as a numpy view over
                the kernel's buffer -- a trajectory runs to millions of frames
                and a returned `std::vector` costs ~35-40 ns each to hand back
                (internal/OutputView.h carries the measured table).
                See internal/OutputView.h.
*/
IMPBFFEXPORT void fret_rate_trace_kernel(
        const std::vector<double>& trajectory,
        const std::vector<double>& acceptor_points,
        double R0,
        double tau0,
        double r_min2,
        double kappa2_scale,
        double** out_view, int* n_out_view
);

//! FRET rate per frame from two trajectories, paired frame by frame.
/*!
    Both dyes resolved in time. The two walks are independent, so pairing frame
    *i* with frame *i* samples the joint distribution -- provided both were
    produced on the same time base.

    \param[in] donor,acceptor positions, flat, three per frame, same length
    \param[in] R0,tau0,r_min2,kappa2_scale as above
    \param[out] out_view,n_out_view the rate at each frame, as above
*/
IMPBFFEXPORT void fret_rate_pair_trace_kernel(
        const std::vector<double>& donor,
        const std::vector<double>& acceptor,
        double R0,
        double tau0,
        double r_min2,
        double kappa2_scale,
        double** out_view, int* n_out_view
);

//! The largest physically meaningful orientation factor: collinear dipoles.
extern IMPBFFEXPORT const double MAX_KAPPA2;

//! Turn an orientation factor into a multiplier on \f$R_0^6\f$.
/*!
    A published Forster radius is quoted *at* the isotropic 2/3, so the transfer
    rate carries the ratio \f$\kappa^2/(2/3)\f$ -- the \f$1.5\kappa^2\f$ of
    the literature -- and passing 2/3 leaves the rate as it was.

    **The diffusing-sphere model does not compute \f$\kappa^2\f$, and cannot**:
    its dye is a structureless point in a volume, with no transition dipole to
    orient. This makes the *assumption* explicit and adjustable instead of
    silently baked into \f$R_0\f$. kappa2_from_dipoles() does compute it, but it
    takes the two dipole vectors, which is exactly what a point has not got.

    There is **no sentinel for "unspecified"**, and deliberately: the isotropic
    default is the value #IMP::bff::kappa2_isotropic() names, so a caller with
    nothing better to say passes that. NaN is an *error* here rather than a
    stand-in for it -- a NaN reaching this from a failed orientation calculation
    would otherwise be silently answered with the isotropic average.

    \throw ValueException for NaN, or outside [0, MAX_KAPPA2]
*/
IMPBFFEXPORT double kappa2_scale(double kappa2);

// The two trace functions below take \p kappa2 without a default: it precedes
// the output view, and C++ will not have a defaulted parameter before an
// undefaulted one. Pass kappa2_isotropic() for "no orientation information".

//! Per-frame FRET rate (1/ns) along a donor trajectory, acceptor averaged.
/*!
    \f$k(t) = (1/\tau_0)\,\langle (R_0/|r_D(t) - r_A|)^6\rangle_A\,
    (\kappa^2/(2/3))\f$

    This resolves the **donor** in time -- the point of simulating its diffusion
    at all -- while treating the acceptor as sampling its own volume quickly
    compared with the donor's excited-state lifetime. That is the *fast-acceptor
    limit*, and it is not valid for an immobilised acceptor; use
    fret_rate_pair_trace() then.

    \param[in] trajectory donor positions, flat, three per frame
    \param[in] acceptor_points acceptor cloud, flat, three per point
    \param[in] r_min closest approach of the two dye centres, A. Both clouds are
               centre positions and can overlap, where \f$(R_0/r)^6\f$ would
               diverge; two dye spheres cannot interpenetrate.
    \param[in] max_acceptor_points cap on the points averaged over, taken by an
               even stride. The average converges long before the full cloud,
               and the cost is `n_frames * n_acceptor`.
    \throw ValueException for an empty cloud, a non-positive \p tau0, or a
           \p kappa2 kappa2_scale() refuses
*/
IMPBFFEXPORT void fret_rate_trace(
        const std::vector<double>& trajectory,
        const std::vector<double>& acceptor_points,
        double R0,
        double tau0,
        double r_min,
        int max_acceptor_points,
        double kappa2,
        double** out_view, int* n_out_view);

//! Per-frame FRET rate (1/ns) from *two* trajectories, both resolved in time.
/*!
    The general case. fret_rate_trace() averages the acceptor over its volume
    instead, which is the fast-acceptor limit. The two do not agree, and the
    reason is **occupancy**, not the averaging order: a cloud weights every
    accessible voxel the same, while a trajectory is weighted by where the dye
    actually spends its time, and a sticky dye dwells near the surface. Which
    way that moves the efficiency depends on the geometry.

    **Unequal lengths are refused, not truncated.** A trajectory is a
    *concatenation of one walk per excitation*, so an arbitrary cut splits a
    walk and pairs one excitation's tail against another's head. Unequal lengths
    are the normal case -- the two dyes sit at different sites and are quenched
    differently -- and the lever is the number of *walks*, not `t_max` or
    `t_step`.

    \throw ValueException for an empty trajectory, a non-positive \p tau0, or
           a frame-count mismatch
*/
IMPBFFEXPORT void fret_rate_pair_trace(
        const std::vector<double>& donor_trajectory,
        const std::vector<double>& acceptor_trajectory,
        double R0,
        double tau0,
        double r_min,
        double kappa2,
        double** out_view, int* n_out_view);

IMPBFF_END_NAMESPACE

 //IMPBFF_FRETRATETRACE_H

#endif  // IMPBFF_FRET_H
