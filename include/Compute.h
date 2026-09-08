/**
 *  \file IMP/bff/Compute.h
 *  \brief Where a kernel runs: the CPU, or something loaded at run time.
 *
 * Some of this library's work is the kind a GPU is good at -- above all the
 * explicit propagation in `DiffusionSolver.h`, which is thousands of steps of
 * a stencil over an `ng^3` grid and costs seconds at the resolutions that
 * matter. None of it is work this library should *depend* on a GPU for: it
 * ships through conda-forge as part of IMP, where the runtime dependency list
 * is a public contract, and it has to keep running on machines that have no
 * GPU at all.
 *
 * So the accelerator is neither a build option nor a dependency. It is a
 * plugin: a shared library that is looked for when the first expensive call
 * happens, used if it is there and if it agrees to take the job, and
 * otherwise not mentioned again. #IMP::bff::get_compute_backend_name says
 * which one answered, so that nobody measures a CPU and reports a GPU.
 *
 * The plugin boundary is C, not C++: a backend may be built by a different
 * compiler than the one that built this library, and `std::vector` across
 * that line is a promise neither side can keep. Buffers are the caller's,
 * sized by the contract below.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_COMPUTE_H
#define IMPBFF_COMPUTE_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/Base.h>

#include <string>

IMPBFF_BEGIN_NAMESPACE

extern "C" {

//! Propagate a density on the grid; the whole loop, not one step.
/*!
    The whole loop deliberately: this library once ported the kernel to C++
    and left the loop in Python, and the suite became forty times slower
    because a solve is 1 000-10 000 steps and each crossed the binding. A
    backend that took one step would repeat that across a wider gap.

    \param[in] cur,d,decay,bounds the four input grids, each `ng*ng*ng`
    \param[in] ng voxels per axis
    \param[in] flux_form `FLUX_SMOLUCHOWSKI` or `FLUX_ITO`
    \param[in] n_steps steps to take
    \param[in] n_out report the surviving population every this many steps
    \param[out] out_fluorescence the reported populations, `n_fluorescence` of
                them, allocated by the caller
    \param[in] n_fluorescence `n_steps / n_out + 1`
    \param[out] out_density the final density, `ng*ng*ng`, allocated by the
                caller
    \return 0 when it did the work; anything else declines the job and the
            caller runs it on the CPU -- a backend with no device, without the
            precision the grid needs, or simply too small a problem to be
            worth the transfer says so this way, per call.
*/
typedef int (*ImpBffPropagateFn)(const double* cur, const double* d,
                                 const double* decay, const double* bounds,
                                 int ng, int flux_form, int n_steps, int n_out,
                                 double* out_fluorescence, int n_fluorescence,
                                 double* out_density);

//! What a plugin offers. Version first, so a mismatch is caught, not crashed.
struct ImpBffComputeBackend {
    //! #IMPBFF_COMPUTE_BACKEND_ABI as the plugin was compiled against it.
    int abi;
    //! For #IMP::bff::get_compute_backend_name, e.g. `wgpu:Metal:Apple M1 Pro`.
    const char* name;
    //! May be null: then this backend has nothing to say about propagation.
    ImpBffPropagateFn propagate;
};

}  // extern "C"

//! Bumped whenever #ImpBffComputeBackend changes shape.
#define IMPBFF_COMPUTE_BACKEND_ABI 1

//! The symbol a plugin exports: `const ImpBffComputeBackend* (*)(void)`.
#define IMPBFF_COMPUTE_BACKEND_SYMBOL "imp_bff_compute_backend"

//! Load a plugin and use it. The CPU stays the default until this succeeds.
/*!
    \param[in] library the plugin's path; the caller found it, because finding
               it is a question about the machine (an environment variable, a
               Python package's data directory) rather than about numerics
    \param[in] argument passed to the plugin's entry point, or empty -- the
               wgpu backend takes the path of the wgpu-native library here
    \return true when the plugin loaded, reported a matching ABI, and was
            installed. On false the CPU path is unchanged and
            #get_compute_backend_error says why.
*/
IMPBFFEXPORT bool load_compute_backend(const std::string& library,
                                       const std::string& argument = "");

//! Go back to the CPU.
IMPBFFEXPORT void reset_compute_backend();

//! `cpu`, or whatever the loaded backend calls itself.
IMPBFFEXPORT std::string get_compute_backend_name();

//! Why the last #load_compute_backend failed; empty when none has.
IMPBFFEXPORT std::string get_compute_backend_error();

//! The backend in force, or null. For the kernels, not for callers.
IMPBFFEXPORT const ImpBffComputeBackend* get_compute_backend();

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_COMPUTE_H
