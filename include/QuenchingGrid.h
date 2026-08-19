/**
 *  \file IMP/bff/QuenchingGrid.h
 *  \brief Stamping spheres of influence onto an accessible-volume grid.
 *
 * A quencher slows a dye near it and quenches it on contact. Both are stamped
 * the same way -- a sphere of a given radius around each centre, combined into
 * whatever is already there -- and differ only in how they combine: stickiness
 * *multiplies* (identity 1) and rates *add* (identity 0), because rates of
 * parallel channels add.
 *
 * Ported from Python by PRD-113: numba is a prototyping tool in this package,
 * not a runtime dependency, so every numerical kernel is C++.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_QUENCHINGGRID_H
#define IMPBFF_QUENCHINGGRID_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! How a stamped sphere combines with what is already on the grid.
enum GridCombine {
    GRID_COMBINE_MULTIPLY = 0,  //!< stickiness: factors multiply, identity 1
    GRID_COMBINE_ADD = 1        //!< rates: parallel channels add, identity 0
};

//! Voxel indices and integer radii of a set of sphere centres.
/*!
    Uses `floor`, **not** truncation. `int()` truncates toward zero, so a centre
    on the negative side of the origin would round *up* while every other map
    here rounds down -- the walk's own occupancy test and the trajectory sampler
    both floor. The two disagreed by one voxel per axis for every centre with a
    negative offset, which is half the grid.

    \param[in] rs centre coordinates, flat, three per centre
    \param[in] r0 grid origin (the attachment point), three doubles
    \param[in] dg voxel edge
    \param[in] ng voxels per axis
    \param[in] radius per-centre radius in Angstrom
    \param[out] ix0,iy0,iz0 voxel indices of each centre
    \param[out] radius_idx radius of each centre in voxels
*/
IMPBFFEXPORT void center_grid_indices(
        const std::vector<double>& rs,
        const std::vector<double>& r0,
        double dg,
        int ng,
        const std::vector<double>& radius,
        std::vector<int>& ix0,
        std::vector<int>& iy0,
        std::vector<int>& iz0,
        std::vector<int>& radius_idx
);

//! Stamp per-centre values into spheres on the accessible volume.
/*!
    Visits only the voxels inside each centre's bounding box rather than testing
    every voxel against every centre, and keeps the outer loop over x-slabs so
    the work is free of write races.

    Sphere membership is strict, \f$d^2 < r^2\f$, so a voxel exactly on the
    boundary is outside.

    Voxels where \p density is zero keep the identity: the dye cannot be there,
    so nothing is stamped.

    \param[in] density accessible-volume density, flat, ng^3
    \param[in] ng voxels per axis
    \param[in] radius per-centre radius
    \param[in] rs centre coordinates, flat
    \param[in] r0 grid origin
    \param[in] dg voxel edge
    \param[in] values per-centre value -- a factor, or a rate
    \param[in] combine GRID_COMBINE_MULTIPLY or GRID_COMBINE_ADD
*/
IMPBFFEXPORT void stamp_spheres(
        const std::vector<double>& density,
        int ng,
        const std::vector<double>& radius,
        const std::vector<double>& rs,
        const std::vector<double>& r0,
        double dg,
        const std::vector<double>& values,
        int combine
,
        double** out_view, int* n_out_view);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_QUENCHINGGRID_H
