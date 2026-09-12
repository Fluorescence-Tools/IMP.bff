/**
 *  \file IMP/bff/OpenMP.h
 *  \brief Whether this build actually threads, which is not the same question
 *         as whether the code asks it to.
 *
 * The kernels in this module carry `#pragma omp parallel for` over their outer
 * loops. A pragma is a *comment* to a compiler that was not given `-fopenmp`,
 * so it costs nothing, warns about nothing, and does nothing — and the code
 * reads as parallel either way.
 *
 * That is worth being able to ask at runtime, because it changes what a
 * measurement means. Every speed figure recorded for this module was taken with
 * OpenMP **off** (IMP's CMake leaves `OpenMP_CXX_FLAGS` empty in the local
 * arm64 build), so they are single-threaded numbers against single-threaded
 * numpy, and the threading is a further factor still on the table.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_OPENMP_H
#define IMPBFF_OPENMP_H

#include <IMP/bff/bff_config.h>

IMPBFF_BEGIN_NAMESPACE

//! True when this module was compiled with OpenMP enabled.
/*! If false, every `#pragma omp` in the module is inert and the kernels run on
    one thread whatever the machine has. */
IMPBFFEXPORT bool built_with_openmp();

//! Threads the kernels will actually use: `omp_get_max_threads()`, or 1.
IMPBFFEXPORT int openmp_thread_count();

IMPBFF_END_NAMESPACE

#endif //IMPBFF_OPENMP_H
