/**
 *  \file IMP/bff/internal/AVLatticeState.h
 *  \brief Per-AV bookkeeping for the space-fixed (lattice) evaluation path.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_INTERNAL_AV_LATTICE_STATE_H
#define IMPBFF_INTERNAL_AV_LATTICE_STATE_H

#include <IMP/bff/bff_config.h>
#include <IMP/Pointer.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/algebra/VectorD.h>
#include <IMP/bff/AVOccupancyMap.h>

#include <vector>

IMPBFF_BEGIN_INTERNAL_NAMESPACE

//! State the lattice path of AV::resample() carries between evaluations.
/** Lives beside the path map on the AV handle (shared between copies of the
    handle) -- it is *not* particle state, so a fresh handle starts cold. */
struct AVLatticeState {
    // The window: lattice index of voxel 0 and the (cubic) edge in voxels
    bool have_window = false;
    int k0[3] = {0, 0, 0};
    int n = 0;

    // What the current tiles were computed from
    bool have_result = false;
    IMP::algebra::Vector3D last_source;
    IMP::algebra::VectorD<9> last_parameter;
    unsigned long generation1 = 0;      // occupancy generation, pass 1
    unsigned long generation2 = 0;      // occupancy generation, pass 2
    IMP::algebra::Vector3D last_mean;

    // Bumped whenever the tiles change; consumers cache against it
    unsigned long result_generation = 0;

    // The (x, y, z, density) cloud of the current tiles, computed once per
    // result and shared by the mean position and the quadrature.
    std::vector<IMP::algebra::Vector4D> cloud;
    unsigned long cloud_generation = 0;
    bool cloud_valid = false;

    // Occupancy sources: private windows, one per pass.
    IMP::Pointer<AVOccupancyMap> private1;
    IMP::Pointer<AVOccupancyMap> private2;

    // Diagnostics
    long n_skip = 0, n_local = 0, n_full = 0, n_roll = 0;
};

IMPBFF_END_INTERNAL_NAMESPACE

#endif /* IMPBFF_INTERNAL_AV_LATTICE_STATE_H */
