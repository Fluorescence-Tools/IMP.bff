/**
 * \file EmBridge.cpp
 * \brief The module's lattice as an IMP::em::DensityMap.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/EmBridge.h>
#include <IMP/em/DensityHeader.h>

#include <algorithm>

IMPBFF_BEGIN_NAMESPACE

IMP::em::DensityMap* create_density_map(const DensityGrid* grid) {
    IMP_USAGE_CHECK(grid, "create_density_map: no grid");
    const GridHeader* gh = grid->get_header();
    IMP_NEW(IMP::em::DensityMap, dm, ());
    // set_void_map allocates and zeroes; the spacing and origin go on after,
    // in the order IMP's own readers use, so the tops come out consistent.
    dm->set_void_map(gh->get_nx(), gh->get_ny(), gh->get_nz());
    dm->update_voxel_size(gh->get_spacing());
    dm->set_origin(gh->get_xorigin(), gh->get_yorigin(), gh->get_zorigin());
    dm->get_header_writable()->set_resolution(gh->get_resolution());
    const double* data = grid->get_data();
    std::copy(data, data + grid->get_number_of_voxels(), dm->get_data());
    return dm.release();
}

IMPBFF_END_NAMESPACE
