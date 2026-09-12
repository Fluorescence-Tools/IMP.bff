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

#include <IMP/bff/IMPCompatibility.h>

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

//! What a DCD header says, without loading any coordinates.
struct IMPBFFEXPORT DCDHeader {
    int n_frames, n_atoms, first_step, step_stride, charmm_version;
    double time_step;
    bool has_unit_cell;
    //! `"<"` for little-endian, `">"` for big — the struct prefix Python used.
    std::string endianness;
    //! Byte offset at which frame data begins.
    int offset;

    DCDHeader()
        : n_frames(0), n_atoms(0), first_step(0), step_stride(0),
          charmm_version(0), time_step(0.0), has_unit_cell(false),
          endianness("<"), offset(0) {}

    IMP_SHOWABLE_INLINE(DCDHeader, out << "DCDHeader(" << n_frames
                                       << " frames, " << n_atoms << " atoms)");
};
IMP_VALUES(DCDHeader, DCDHeaders);

//! Read a DCD header without loading coordinates.
/*!
    DCD carries no magic number for endianness; the convention is to read the
    leading Fortran record length and see which byte order makes it the expected
    84.

    Only what the rotamer libraries need is implemented — fixed atom counts, no
    velocity blocks, no four-dimensional trajectories. Anything outside that
    throws rather than guessing, because a trajectory silently read as the wrong
    shape is worse than one that refuses.

    \throw ValueException when the file is not a DCD this reader handles
    \throw IOException when it cannot be read
*/
IMPBFFEXPORT DCDHeader read_dcd_header(const std::string& path);

//! Read coordinates from a DCD trajectory.
/*! \param[in] max_frames stop after this many; negative reads all of them
    \param[out] out_view,n_out_view `n_frames * n_atoms * 3`, in the file's own
                units (Angstrom for the bundled libraries) */
IMPBFFEXPORT void read_dcd(const std::string& path, int max_frames = -1,
                           double** out_view = NULL, int* n_out_view = NULL);

//! Coordinates from a trajectory, whatever format it is in.
/*!
    **BinaryCIF is the format this package stores.** The rotamer libraries were
    re-encoded on 2026-08-19: 44.78 MB of DCD and XTC became 17.99 MB of
    `.bcif`, verified exact on a 0.1 A grid, and read through the C parser IMP
    already vendors. `.dcd` still reads, because the format is not gone from the
    world — a user's own library may be one. It is simply not what is shipped.

    \param[in] n_atoms required for BinaryCIF, which stores one row per
               (atom, frame) and cannot infer the split; callers have it from
               the companion PDB. Ignored for DCD, which carries its own.
    \throw ValueException for any other suffix, or for a BinaryCIF with no
           \p n_atoms
*/
IMPBFFEXPORT void read_trajectory(const std::string& path, int n_atoms,
                                  int max_frames, double** out_view,
                                  int* n_out_view);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_TRAJECTORYIO_H
