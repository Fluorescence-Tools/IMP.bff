/**
 *  \file IMP/bff/internal/OutputView.h
 *  \brief Handing a large array back to numpy without copying it.
 *
 * The cost of crossing the SWIG boundary, per element, measured on this box on
 * 2026-08-19 (`quenched_decay` and `lifetime_spectrum_decay`, 10k–400k
 * elements, work held constant while the array size varied):
 *
 *     OUT  returned std::vector -> tuple            8–16 ns   (SWIG)
 *          ...then np.asarray over that tuple      ~27 ns     (numpy)
 *          ------------------------------------------------
 *          what the caller actually pays          35–40 ns
 *          ARGOUTVIEWM_ARRAY1                        free
 *          a SWIG proxy walked per __getitem__     ~340 ns
 *
 *     IN   ndarray -> const std::vector<double>&   32–37 ns   (default typemap)
 *          list    -> const std::vector<double>&    8–13 ns
 *          ndarray -> const std::vector<double>&    4.4 ns    (bulk-copy
 *                                                              typemap, in
 *                                                              IMP_bff.types.i)
 *          ndarray -> (double*, int) IN_ARRAY1        free
 *
 * An earlier note in this file put the return path at 66 ns/element. That was
 * too high; the decomposition above is what two kernels agree on. The ordering
 * it implies has not changed — for a large result the return path still costs
 * far more than the arithmetic, and a 400×350 pair matrix was ~20 ms of
 * marshalling against about 1 ms of work.
 *
 * numpy's `ARGOUTVIEWM_ARRAY1` typemap takes a `(pointer, length)` pair and
 * wraps the buffer in an ndarray that **owns** it. No copy, no conversion.
 *
 * Free, and dangerous. Three things make it so, and none of them announces
 * itself:
 *
 * * **The allocator must match.** numpy releases the buffer with `free`, so it
 *   has to come from `malloc`/`calloc`. `new[]` here is undefined behaviour
 *   that will not show up on this side of the boundary.
 * * **The parameter name decides which typemap binds.** `output` is claimed by
 *   *both* `ARGOUTVIEW` (numpy does not own it — the buffer leaks) and
 *   `ARGOUTVIEWM` in `IMP_bff.types.i`; whichever is declared last wins.
 *   Use `out_view`, which only the managed typemap claims.
 * * **Allocation can fail.** A null pointer becomes an ndarray over address
 *   zero, and the caller faults on first touch, a long way from the cause.
 *
 * These helpers are the one place that gets it right.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#ifndef IMPBFF_INTERNAL_OUTPUTVIEW_H
#define IMPBFF_INTERNAL_OUTPUTVIEW_H

#include <IMP/bff/bff_config.h>

#include <cstdlib>
#include <cstring>
#include <vector>

IMPBFF_BEGIN_INTERNAL_NAMESPACE

//! Allocate a zeroed buffer for numpy to adopt, and publish it.
/*!
    \param[in] n elements wanted
    \param[out] out_view,n_out_view the pair numpy's typemap consumes
    \return the buffer, or `nullptr` if allocation failed — in which case an
            **empty** view has already been published, so the caller can simply
            return and Python sees a zero-length array rather than a fault.
*/
inline double* new_double_view(std::size_t n, double** out_view, int* n_out_view) {
    double* buffer = static_cast<double*>(std::calloc(n ? n : 1, sizeof(double)));
    if (buffer == nullptr) {
        *out_view = static_cast<double*>(std::calloc(1, sizeof(double)));
        *n_out_view = 0;
        return nullptr;
    }
    *out_view = buffer;
    *n_out_view = static_cast<int>(n);
    return buffer;
}

//! The same, for an integer result.
inline int* new_int_view(std::size_t n, int** out_view, int* n_out_view) {
    int* buffer = static_cast<int*>(std::calloc(n ? n : 1, sizeof(int)));
    if (buffer == nullptr) {
        *out_view = static_cast<int*>(std::calloc(1, sizeof(int)));
        *n_out_view = 0;
        return nullptr;
    }
    *out_view = buffer;
    *n_out_view = static_cast<int>(n);
    return buffer;
}

//! Publish a copy of a vector the kernel had to build anyway.
/*! For kernels whose algorithm needs its own buffer — a double-buffered sweep
    cannot know in advance which of the two holds the answer — copying once at
    the end is a `memcpy`, microseconds against milliseconds of marshalling.
    Prefer filling the view directly where the shape of the loop allows it. */
inline void copy_to_view(const std::vector<double>& v,
                         double** out_view, int* n_out_view) {
    double* buffer = new_double_view(v.size(), out_view, n_out_view);
    if (buffer != nullptr && !v.empty()) {
        std::memcpy(buffer, v.data(), v.size() * sizeof(double));
    }
}

IMPBFF_END_INTERNAL_NAMESPACE

#endif //IMPBFF_INTERNAL_OUTPUTVIEW_H
