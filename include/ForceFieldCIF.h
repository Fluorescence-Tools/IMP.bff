/**
 * \file IMP/bff/ForceFieldCIF.h
 * \brief A coarse-grained force-field system in and out of mmCIF.
 *
 * Reading is IMP's parser -- the vendored `ihm_format.h` behind
 * `IMP::atom::read_mmcif`. Writing is this module's own, through the one
 * writer in `internal/Cif.h`, because IMP can write no CIF from C++: the
 * vendored ihm library is read-only, `IMP::atom` has readers only, and
 * `IMP.mmcif` is Python with no spelling for these `_ff_*` tables.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_FORCEFIELDCIF_H
#define IMPBFF_FORCEFIELDCIF_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/ProbeForceField.h>

#include <string>
#include <utility>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Read a `ProbeForceFieldSystem` from an mmCIF file.
/** Eighteen `_ff_*` categories, through the `ihm` C reader IMP already vendors
    and `libimp_atom` already exports -- the same parser the BinaryCIF
    trajectory reader uses, in text mode.

    Sites may be referred to either by full id (`site_id_1`) or by a compact
    integer (`n1`), and group membership additionally by half-open ranges. The
    compact form is resolved after the file is read, because a term may refer
    to a site number whose `_ff_site` row has not been seen yet; sites without
    a declared `site_no` are numbered by first appearance into the gaps. */
IMPBFFEXPORT ProbeForceFieldSystem read_forcefield_cif(const std::string& path);


// --------------------------------------------------------------------------
// String utilities (site-id splitting, range compression)
// --------------------------------------------------------------------------

//! Split a site id like "CX4/S1" into (prefix, serial). Returns ("", -1) on failure.
IMPBFFEXPORT std::pair<std::string, int> split_site_id(const std::string& site_id);

//! Expand a group range "A1".."A3" into ["A1","A2","A3"].
IMPBFFEXPORT std::vector<std::string> expand_site_range(
        const std::string& start_id, const std::string& end_id);

//! Compress a list of integers into sorted contiguous (start, end) pairs.
IMPBFFEXPORT std::vector<std::pair<int, int>> compress_int_ranges(
        const std::vector<int>& nos);

// --------------------------------------------------------------------------
// Force-field system writer
// --------------------------------------------------------------------------

//! Write a ProbeForceFieldSystem to mmCIF.
/*!
    Mirrors the reader in ForceFieldCIF.h. Writes the eighteen _ff_* categories
    plus _atom_site (from the component MOL2/PDB files) and _flr_probe_list.
    Site numbers are assigned and group members are written as integer ranges.
*/
IMPBFFEXPORT void write_probe_forcefield_cif(const std::string& path,
                                            const ProbeForceFieldSystem& system);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_FORCEFIELDCIF_H
