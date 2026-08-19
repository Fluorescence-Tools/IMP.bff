/**
 *  \file IMP/bff/TrajectoryIO.h
 *  \brief Reading rotamer-library trajectories from BinaryCIF.
 *
 * BinaryCIF is this package's trajectory format as of 2026-08-19, replacing
 * DCD and XTC. It is smaller than either -- 1.27 bytes per coordinate against
 * DCD's 4.31 and XTC's 1.60, on a 0.1 A grid -- and it is decoded by the C
 * implementation of `ihm` that IMP already vendors, so reading it costs no new
 * dependency and no build change.
 *
 * The measurements, the precision that buys the size, and the two encoding
 * traps are in `okf/validation/bcif_for_trajectories.md`. The encoder is
 * `scripts/trajectory_to_bcif.py`, which exists because python-ihm's writer
 * implements neither FixedPoint nor IntegerPacking.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_TRAJECTORYIO_H
#define IMPBFF_TRAJECTORYIO_H

#include <IMP/bff/bff_config.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Read coordinates from a BinaryCIF trajectory.
/*!
    The file stores one row per (atom, frame) with `x`, `y`, `z` columns, laid
    out atom-major -- every frame of atom 0, then atom 1 -- so that the delta
    encoding runs along an atom's own series.

    \param[in] path the `.bcif` file
    \param[in] n_atoms atoms per frame; the frame count follows from the row
               count, and a row count that is not a multiple of this is an error
               rather than a silent truncation
    \param[in] category the CIF category to read, `_rotamer_coord` by default
    \param[out] out_view,n_out_view flat `(n_frames, n_atoms, 3)` in Angstrom
    \throws IMP::IOException if the file cannot be read or does not parse
    \throws IMP::ValueException if the row count is not a multiple of \p n_atoms
*/
IMPBFFEXPORT void read_bcif_trajectory(
        const std::string& path, int n_atoms,
        const std::string& category, double** out_view, int* n_out_view);

//! How many rows a BinaryCIF trajectory holds, without decoding it all.
/*!
    Cheap enough to call first: it still parses, but it discards the
    coordinates instead of accumulating them.
*/
IMPBFFEXPORT int bcif_trajectory_rows(const std::string& path,
                                      const std::string& category);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_TRAJECTORYIO_H
