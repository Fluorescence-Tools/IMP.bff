/**
 *  \file IMP/bff/AVDistance.h
 *  \brief Distances and reductions over accessible-volume point clouds.
 *
 * These take **point arrays** rather than the ``AV`` decorators the
 * `av_distance` family in AV.h takes, so they serve a rotamer library, a
 * coarse-grained ensemble or an MD trajectory as readily as an accessible
 * volume -- anything that supplies states.
 *
 * Ported from Python by PRD-113: numba is a prototyping tool in this package,
 * not a runtime dependency, so every numerical kernel is C++.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_AVDISTANCE_H
#define IMPBFF_AVDISTANCE_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! What split_contact_volume() found in a voxel.
enum AVVoxel {
    AV_VOXEL_EMPTY = 0,       //!< no density: outside the accessible volume
    AV_VOXEL_CONTACT = 1,     //!< inside a slow centre's sphere
    AV_VOXEL_FREE = 2         //!< accessible, and away from every slow centre
};

//! Weighted mean position of a point cloud.
/*!
    \param[in] points flat, four per point: x, y, z, weight
    \return three coordinates; the origin for an empty or zero-weight cloud
*/
IMPBFFEXPORT std::vector<double> points_weighted_mean(
        const std::vector<double>& points);

//! Random distance and weight-product samples between two clouds.
/*!
    Draws one point from each cloud per sample. **The two clouds are
    independent**, so any pairing samples the joint distribution.

    The clouds arrive straight from numpy, one `(n, 4)` array per cloud via
    IN_ARRAY2. The result is published as an `(n_samples, 2)` managed numpy
    view: column 0 is the distance, column 1 the weight product. Both shapes
    are part of the contract -- they used to be applied by Python `.ravel()` /
    `.reshape()` over the flat buffers, and the caller's picture of the data is
    not something a second side should state a second time.

    \param[in] p1,p2 flat clouds, four per point. Via IN_ARRAY2, carrying the
                (n, 4) shape; the values are read as `4 * n` in row order.
    \param[in] n_samples samples to draw
    \param[in] seed for reproducibility
    \param[out] output,n_output1,n_output2 `(n_samples, 2)` view. The buffer
                is malloc'd and numpy adopts it.
*/
IMPBFFEXPORT void random_distances(
        double* p1, int n_p1, int n_p1c,
        double* p2, int n_p2, int n_p2c,
        int n_samples,
        int seed,
        double** output, int* n_output1, int* n_output2
);

//! Weighted mean inter-point distance \f$\langle R_{DA}\rangle\f$.
IMPBFFEXPORT double average_distance(
        const std::vector<double>& p1,
        const std::vector<double>& p2,
        int n_samples = 50000,
        int seed = 0
);

//! FRET-averaged distance \f$\langle R_{DA}\rangle_E\f$.
/*!
    Averages the *efficiency* over the sampled pairs and converts back, which is
    not the same as averaging the distance -- \f$1/r^6\f$ weights close pairs far
    more heavily, so this is always the shorter of the two.

    \return 0 when the mean efficiency saturates at 1, infinity when it reaches 0
*/
IMPBFFEXPORT double mean_fret_distance(
        const std::vector<double>& p1,
        const std::vector<double>& p2,
        double forster_radius = 52.0,
        int n_samples = 50000,
        int seed = 0
);

//! The four distance statistics of a weighted distance sample.
/*!
    One reduction rather than four passes, and one place where the conventions
    live: \f$\langle R_{DA}\rangle\f$, \f$R_E\f$, \f$\langle E\rangle\f$ and
    the width of the distance distribution were computed in three modules that
    disagreed at the limits.

    Note what is **not** here: \f$R_{mp}\f$, the distance between the clouds'
    mean positions. It is not a function of the distribution of pair distances
    -- \f$|\langle a\rangle - \langle b\rangle|\f$ cannot be recovered from
    \f$|a - b|\f$ -- so it takes the clouds, not a sample of them. A slot in
    this tuple returned the mean distance under that name until 2026-07-28;
    on 148l E15/E90 the two were 8 % apart.

    \param[in] distances sampled pair distances
    \param[in] weights per-sample weight
    \param[in] forster_radius \f$R_0\f$
    \return four values: mean distance, \f$R_E\f$, mean efficiency, standard
            deviation. \f$R_E\f$ is 0 when the mean efficiency saturates at 1
            and infinity when it reaches 0; all four are 0 for zero total weight.
*/
IMPBFFEXPORT std::vector<double> distance_sample_statistics(
        const std::vector<double>& distances,
        const std::vector<double>& weights,
        double forster_radius = 52.0
);

//! Asymmetric chi-squared contribution of one distance restraint.
/*!
    \f$\chi^2 = (d_m - d_e)^2 / \sigma^2\f$, with \f$\sigma\f$ taken from
    whichever side of the experimental value the model falls on. A non-positive
    error contributes nothing rather than dividing by zero.
*/
IMPBFFEXPORT double chi2_score(double model_distance,
                               double experimental_distance,
                               double error_neg, double error_pos);

//! Single-pair FRET efficiency \f$1/(1 + (r/R_0)^6)\f$.
IMPBFFEXPORT double fret_efficiency(double distance,
                                    double forster_radius = 52.0);

//! The distance a single-pair FRET efficiency implies.
IMPBFFEXPORT double distance_from_fret_efficiency(
        double efficiency, double forster_radius = 52.0);

//! Occupied voxels of a density as a weighted point cloud.
/*!
    Published as an `(n, 4)` managed numpy view: x, y, z, density. The shape is
    part of the contract, stated once here rather than by a Python `.reshape()`.

    \param[in] density `(nx, ny, nz)` numpy grid, via IN_ARRAY3
    \param[in] nx,ny,nz grid shape
    \param[in] dg voxel edge
    \param[in] r0 grid anchor
    \param[in] threshold keep voxels strictly above this
    \param[out] output,n_output1,n_output2 `(n, 4)` view; the buffer is
                malloc'd and numpy adopts it.
*/
IMPBFFEXPORT void density_to_points(
        double* density, int nx, int ny, int nz,
        double dg,
        const std::vector<double>& r0,
        double threshold = 0.0,
        double** output = 0, int* n_output1 = 0, int* n_output2 = 0
);
//! Label each voxel of an accessible volume contact, free, or empty.
/*!
    A voxel is *contact* when it lies inside any slow centre's sphere. One
    label array rather than a pair of masks: the two are complementary within
    the occupied volume, so returning both invites them to disagree.

    \p r0 is **the grid anchor**: the position of voxel \f$(ng-1)/2\f$ on each
    axis, not the corner and not the attachment atom. The inverse map uses
    `floor` and the *integer* offset, both deliberately: the float corner
    differs on every even \p ng -- the normal case -- putting the stamped
    spheres half a voxel from the density they mask, and truncation would round
    a negative offset up while every other index map rounds down.

    \param[in] density flat ng^3
    \param[in] ng voxels per axis
    \param[in] dg voxel edge
    \param[in] rad per-centre slow radius
    \param[in] rs centre coordinates, flat
    \param[in] r0 grid anchor
    \param[out] output_i,dim1,dim2,dim3 an `(ng, ng, ng)` #AVVoxel view; the
                buffer is malloc'd and numpy adopts it.
*/
IMPBFFEXPORT void split_contact_volume(
        const std::vector<double>& density,
        int ng,
        double dg,
        const std::vector<double>& rad,
        const std::vector<double>& rs,
        const std::vector<double>& r0,
        int** output_i, int* dim1, int* dim2, int* dim3
);

//! The contact and free regions of an accessible volume, as two uint8 masks.
/*!
    A voxel is *contact* when it lies inside any slow centre's sphere. The two
    masks are complementary within the occupied volume; returning both invites
    them to disagree, which is why this is one call publishing the pair.

    `uint8`, not `float64`: at `ng = 92` a float64 pair costs 12.5 MB per
    labelling site against 1.6 MB. A residue scan builds one per site.

    \param[in] density `(ng, ng, ng)`, straight from numpy via IN_ARRAY3
    \param[in] radius a scalar slow radius, broadcast over every centre, or one
                value per centre
    \param[out] contact,free each `(ng, ng, ng)` uint8 views; both buffers are
                malloc'd and numpy adopts them
*/
IMPBFFEXPORT void split_contact_volume_masks(
        double* density, int ng, int ng2, int ng3,
        double dg,
        const std::vector<double>& radius,
        double* rs, int n_rs, int n_rsc,
        const std::vector<double>& r0,
        unsigned char** contact, int* contact_dim1, int* contact_dim2,
        int* contact_dim3,
        unsigned char** free, int* free_dim1, int* free_dim2, int* free_dim3
);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_AVDISTANCE_H
