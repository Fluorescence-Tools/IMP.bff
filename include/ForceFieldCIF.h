/**
 * \file IMP/bff/ForceFieldCIF.h
 * \brief Reading a coarse-grained force-field system from mmCIF.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_FORCEFIELDCIF_H
#define IMPBFF_FORCEFIELDCIF_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/DyeForceField.h>

#include <string>

IMPBFF_BEGIN_NAMESPACE

//! Read a `DyeForceFieldSystem` from an mmCIF file.
/** Eighteen `_ff_*` categories, through the `ihm` C reader IMP already vendors
    and `libimp_atom` already exports -- the same parser the BinaryCIF
    trajectory reader uses, in text mode.

    Sites may be referred to either by full id (`site_id_1`) or by a compact
    integer (`n1`), and group membership additionally by half-open ranges. The
    compact form is resolved after the file is read, because a term may refer
    to a site number whose `_ff_site` row has not been seen yet; sites without
    a declared `site_no` are numbered by first appearance into the gaps. */
IMPBFFEXPORT DyeForceFieldSystem read_forcefield_cif(const std::string& path);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_FORCEFIELDCIF_H
