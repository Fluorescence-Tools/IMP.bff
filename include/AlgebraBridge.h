/**
 *  \file IMP/bff/AlgebraBridge.h
 *  \brief The core's flat numbers as IMP::algebra values.
 *
 * The core plans, samples and stores in plain numbers; IMP holds rigid
 * bodies as Transformation3D. The conversions live here, in the connection
 * layer, so the core's headers name none of IMP's geometry beyond VectorD
 * (PRD-137 step 6). Formerly in RRT.h.
 *
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#ifndef IMPBFF_ALGEBRABRIDGE_H
#define IMPBFF_ALGEBRABRIDGE_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/RRT.h>
#include <IMP/algebra/Transformation3D.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! A rigid-body configuration as a transformation, and back.
/*! `(tx, ty, tz, rx, ry, rz)` with the angles in the fixed-xyz convention --
    the six numbers #grow_rigid_body_rrt plans in. */
IMPBFFEXPORT std::vector<double> transformation_to_configuration(
        const IMP::algebra::Transformation3D& transformation);

IMPBFFEXPORT IMP::algebra::Transformation3D configuration_to_transformation(
        const std::vector<double>& configuration);

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_ALGEBRABRIDGE_H
