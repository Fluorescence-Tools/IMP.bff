/**
 *  \file IMP/bff/internal/Normalize.h
 *  \brief Divide a discretised density by its sum.
 *
 * Shared rather than duplicated per translation unit: IMP concatenates a
 * module's sources, so a helper in an anonymous namespace in two `.cpp` files
 * is a redefinition, not two private copies.
 *
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#ifndef IMPBFF_INTERNAL_NORMALIZE_H
#define IMPBFF_INTERNAL_NORMALIZE_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_INTERNAL_NAMESPACE

//! Divide by the sum in place, leaving an all-zero vector alone.
inline void normalize_sum(std::vector<double>& y) {
    double total = 0.0;
    for (double v : y) total += v;
    if (total > 0.0) {
        for (double& v : y) v /= total;
    }
}

IMPBFF_END_INTERNAL_NAMESPACE

#endif //IMPBFF_INTERNAL_NORMALIZE_H
