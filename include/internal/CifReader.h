/**
 *  \file IMP/bff/internal/CifReader.h
 *  \brief Reading values out of IMP's mmCIF parser.
 *
 * Parsing is `ihm_format.h`, the C reader IMP vendors and runs its own
 * `IMP::atom::read_mmcif` with. This header is only the accessors on top of a
 * parsed keyword, so that "present, and neither `.` nor `?`" is written once
 * rather than beside every category handler.
 *
 * Internal. The public surface takes paths and returns values.
 */
#ifndef IMPBFF_INTERNAL_CIFREADER_H
#define IMPBFF_INTERNAL_CIFREADER_H

#include <IMP/bff/bff_config.h>

#include "ihm_format.h"

#include <IMP/bff/internal/Text.h>

#include <string>

IMPBFF_BEGIN_INTERNAL_NAMESPACE

//! A keyword carries a value only if it is in the file and neither `.` nor `?`.
inline bool has(ihm_keyword* k) {
    return k && k->in_file && !k->omitted && !k->unknown;
}

inline std::string txt(ihm_keyword* k) {
    return has(k) && k->data.str ? std::string(k->data.str) : std::string();
}

inline double dbl(ihm_keyword* k, double fallback) {
    return has(k) ? k->data.fval : fallback;
}

//! A missing number is NaN, for a column where every value is meaningful.
inline double dbl(ihm_keyword* k) {
    return has(k) ? k->data.fval : nan_value();
}

inline int integer(ihm_keyword* k, int fallback) {
    return has(k) ? k->data.ival : fallback;
}

inline bool flag(ihm_keyword* k, bool fallback) {
    return has(k) ? k->data.bval : fallback;
}

IMPBFF_END_INTERNAL_NAMESPACE

#endif  // IMPBFF_INTERNAL_CIFREADER_H
