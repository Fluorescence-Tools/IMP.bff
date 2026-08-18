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

    \param[in] p1,p2 flat clouds, four per point
    \param[in] n_samples samples to draw
    \param[in] seed for reproducibility
    \return flat, two per sample: distance, weight product
*/
IMPBFFEXPORT std::vector<double> random_distances(
        const std::vector<double>& p1,
        const std::vector<double>& p2,
        int n_samples,
        int seed = 0
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

//! Occupied voxels of a density as a weighted point cloud.
/*!
    \param[in] density flat, nx*ny*nz
    \param[in] nx,ny,nz grid shape
    \param[in] dg voxel edge
    \param[in] r0 grid anchor
    \param[in] threshold keep voxels strictly above this
    \return flat, four per kept point
*/
IMPBFFEXPORT std::vector<double> density_to_points(
        const std::vector<double>& density,
        int nx, int ny, int nz,
        double dg,
        const std::vector<double>& r0,
        double threshold = 0.0
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
    \return flat ng^3 of #AVVoxel
*/
IMPBFFEXPORT std::vector<int> split_contact_volume(
        const std::vector<double>& density,
        int ng,
        double dg,
        const std::vector<double>& rad,
        const std::vector<double>& rs,
        const std::vector<double>& r0
);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_AVDISTANCE_H
