/**
 * \file AlgebraBridge.cpp
 * \brief The core's flat numbers as IMP::algebra values (formerly in RRT.cpp).
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/AlgebraBridge.h>
#include <IMP/algebra/Rotation3D.h>

IMPBFF_BEGIN_NAMESPACE

std::vector<double> transformation_to_configuration(
        const IMP::algebra::Transformation3D& transformation) {
    const IMP::algebra::Vector3D t = transformation.get_translation();
    const IMP::algebra::FixedXYZ e =
            IMP::algebra::get_fixed_xyz_from_rotation(
                    transformation.get_rotation());
    std::vector<double> out;
    out.push_back(t[0]);
    out.push_back(t[1]);
    out.push_back(t[2]);
    out.push_back(e.get_x());
    out.push_back(e.get_y());
    out.push_back(e.get_z());
    return out;
}

IMP::algebra::Transformation3D configuration_to_transformation(
        const std::vector<double>& c) {
    if (c.size() < 6) return IMP::algebra::Transformation3D();
    return IMP::algebra::Transformation3D(
            IMP::algebra::get_rotation_from_fixed_xyz(c[3], c[4], c[5]),
            IMP::algebra::Vector3D(c[0], c[1], c[2]));
}

IMPBFF_END_NAMESPACE
