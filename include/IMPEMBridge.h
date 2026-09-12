/**
 *  \file IMP/bff/IMPEMBridge.h
 *  \brief The module's lattice as an IMP::em::DensityMap.
 *
 * The one place the module still speaks `IMP::em` (PRD-137: the core's
 * lattice, DensityGrid, is its own, and writes MRC itself). This copies a
 * grid into an em map for whatever IMP's em module does with one -- the
 * other map formats, its readers, its fitting. Connection layer.
 *
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#ifndef IMPBFF_IMPEMBRIDGE_H
#define IMPBFF_IMPEMBRIDGE_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/DensityGrid.h>
#include <IMP/em/DensityMap.h>

IMPBFF_BEGIN_NAMESPACE

//! A new IMP::em::DensityMap holding a copy of \p grid's values, header and origin.
/*! Was PathMap::create_density_map(); that spelling stays in Python. */
IMPBFFEXPORT IMP::em::DensityMap* create_density_map(const DensityGrid* grid);

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_IMPEMBRIDGE_H
