/**
 *  \file IMP/bff/internal/GridShape.h
 *  \brief Recovering a cubic grid's edge length from a flat buffer.
 *
 * Grids cross the SWIG boundary flat, because a numpy view over a kernel's
 * buffer is one-dimensional and the shape belongs to the caller's picture of
 * the data. An object that stores one therefore has to recover `ng` from the
 * length, and it has to *check* rather than assume: a non-cubic length is a
 * caller error and rounding the cube root silently would accept it.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#ifndef IMPBFF_INTERNAL_GRIDSHAPE_H
#define IMPBFF_INTERNAL_GRIDSHAPE_H

#include <IMP/bff/bff_config.h>

#include <cmath>
#include <cstddef>

IMPBFF_BEGIN_INTERNAL_NAMESPACE

//! The edge length of a cubic grid with \p n voxels, or 0 if it is not cubic.
inline int cube_side(std::size_t n) {
    if (n == 0) return 0;
    const int s = static_cast<int>(std::lround(std::cbrt(static_cast<double>(n))));
    return (static_cast<std::size_t>(s) * s * s == n) ? s : 0;
}

IMPBFF_END_INTERNAL_NAMESPACE

#endif //IMPBFF_INTERNAL_GRIDSHAPE_H
