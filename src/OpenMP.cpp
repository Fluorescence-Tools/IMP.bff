/**
 * \file OpenMP.cpp
 * \brief Whether this build actually threads.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/OpenMP.h>

#ifdef _OPENMP
#include <omp.h>
#endif

IMPBFF_BEGIN_NAMESPACE

bool built_with_openmp() {
#ifdef _OPENMP
    return true;
#else
    return false;
#endif
}

int openmp_thread_count() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

IMPBFF_END_NAMESPACE
