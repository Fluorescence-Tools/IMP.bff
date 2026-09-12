/**
 * \file ComputeBackend.cpp
 * \brief Loading an accelerator, and doing without one.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/ComputeBackend.h>

#include <sstream>

#if IMPBFF_WITH_GPU
#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif
#endif

IMPBFF_BEGIN_NAMESPACE

namespace {

//! The backend in force, its handle, and why the last attempt failed.
/*! Function-local statics: the order in which translation units initialise is
    not ours to arrange, and a caller may ask before main() in principle. */
struct State {
    const ImpBffComputeBackend* backend;
    void* handle;
    std::string error;
    State() : backend(0), handle(0) {}
};

State& state() {
    static State s;
    return s;
}

#if IMPBFF_WITH_GPU
void close_handle(void*& handle) {
    if (!handle) return;
#ifdef _WIN32
    FreeLibrary(reinterpret_cast<HMODULE>(handle));
#else
    dlclose(handle);
#endif
    handle = 0;
}

std::string last_error() {
#ifdef _WIN32
    const DWORD code = GetLastError();
    std::ostringstream oss;
    oss << "error " << code;
    return oss.str();
#else
    const char* e = dlerror();
    return e ? std::string(e) : std::string("unknown error");
#endif
}
#endif  // IMPBFF_WITH_GPU

}  // namespace

#if !IMPBFF_WITH_GPU

// Compiled without the door: no loader, no dynamic-library call anywhere in
// this translation unit, and a refusal that says why rather than one that
// looks like a missing file.
bool load_compute_backend(const std::string& library, const std::string& argument) {
    (void)library;
    (void)argument;
    state().error =
            "this build of IMP.bff was compiled with IMPBFF_WITH_GPU=0, so it "
            "loads no compute backend; the kernels run on the CPU";
    return false;
}

#else

bool load_compute_backend(const std::string& library, const std::string& argument) {
    State& s = state();
    s.error.clear();
    if (library.empty()) {
        s.error = "no library named";
        return false;
    }

#ifdef _WIN32
    void* handle = reinterpret_cast<void*>(LoadLibraryA(library.c_str()));
#else
    // LOCAL, not GLOBAL: a backend brings its own copy of a graphics stack,
    // and putting those symbols in the global namespace is how two of them
    // start answering for each other.
    dlerror();
    void* handle = dlopen(library.c_str(), RTLD_NOW | RTLD_LOCAL);
#endif
    if (!handle) {
        s.error = "cannot load " + library + ": " + last_error();
        return false;
    }

    typedef const ImpBffComputeBackend* (*EntryFn)(const char*);
#ifdef _WIN32
    EntryFn entry = reinterpret_cast<EntryFn>(GetProcAddress(
            reinterpret_cast<HMODULE>(handle), IMPBFF_COMPUTE_BACKEND_SYMBOL));
#else
    EntryFn entry = reinterpret_cast<EntryFn>(
            dlsym(handle, IMPBFF_COMPUTE_BACKEND_SYMBOL));
#endif
    if (!entry) {
        s.error = library + " has no " + IMPBFF_COMPUTE_BACKEND_SYMBOL;
        close_handle(handle);
        return false;
    }

    const ImpBffComputeBackend* backend = entry(argument.c_str());
    if (!backend) {
        // The plugin loaded and then declined -- no device, no driver, a
        // wgpu-native it could not open. That is a normal outcome, not a bug.
        s.error = library + " found no usable device";
        close_handle(handle);
        return false;
    }
    if (backend->abi != IMPBFF_COMPUTE_BACKEND_ABI) {
        std::ostringstream oss;
        oss << library << " speaks backend ABI " << backend->abi
            << ", this build speaks " << IMPBFF_COMPUTE_BACKEND_ABI;
        s.error = oss.str();
        close_handle(handle);
        return false;
    }

    close_handle(s.handle);
    s.handle = handle;
    s.backend = backend;
    return true;
}

#endif  // IMPBFF_WITH_GPU

void reset_compute_backend() {
    State& s = state();
    s.backend = 0;
    // The handle stays open: the backend's name and function pointers may
    // still be held by something mid-call, and closing under it is worse than
    // holding a few megabytes.
}

std::string get_compute_backend_name() {
    const State& s = state();
    return (s.backend && s.backend->name) ? std::string(s.backend->name)
                                          : std::string("cpu");
}

std::string get_compute_backend_error() { return state().error; }

const ImpBffComputeBackend* get_compute_backend() { return state().backend; }

bool built_with_gpu_support() { return IMPBFF_WITH_GPU != 0; }

IMPBFF_END_NAMESPACE
