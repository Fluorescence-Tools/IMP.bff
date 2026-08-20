/**
 * \file IMP/bff/RotamerStatistics.h
 * \brief Weighted averaging over a rotamer ensemble.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_ROTAMERSTATISTICS_H
#define IMPBFF_ROTAMERSTATISTICS_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Per-frame weights from a pair of partition functions.
/** `z` is `(n_frames, 2)` flattened -- the donor's and acceptor's partition
    function in each frame. The weight of a frame is the product, normalised
    over frames. A single `(2,)` row means one frame and weight 1.

    When every product is zero the weights are uniform rather than undefined:
    a frame in which neither dye has any accessible conformer says nothing
    about the others. */
IMPBFFEXPORT void rotamer_frame_weights(
        double* z_values, int n_frames, int n_pair,
        double** out_view, int* n_out_view);

//! Weighted mean, standard deviation and standard error of `values`.
/** Non-finite values and their weights are dropped first -- a frame where the
    dye could not be placed contributes nothing rather than poisoning the mean.
    The weights are renormalised **after** that removal, so the surviving
    frames still sum to one.

    Returns `(mean, sd, se)`; all three are NaN when nothing is finite. The
    standard error divides by the *count* of surviving frames, not by the
    effective count -- which is what the Python did. */
IMPBFFEXPORT std::vector<double> weighted_average_sd_se(
        double* values, int n_values,
        double* weights, int n_weights);

//! The effective number of contributing frames, `exp` of the entropy.
/** Zero weights are dropped, then `exp(-sum w log(w / uniform))`. Equals the
    frame count for uniform weights and falls toward 1 as the weight
    concentrates on one frame. Zero when nothing contributes. */
IMPBFFEXPORT double effective_frame_fraction(double* weights, int n_weights);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_ROTAMERSTATISTICS_H
