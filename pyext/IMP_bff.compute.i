/*
 * Where a kernel runs, and how the accelerator is found.
 *
 * The C++ side (Compute.h) knows how to load a plugin and nothing about where
 * one might be. Finding it is a question about the machine -- an environment
 * variable, a Python package's data directory -- so it is answered here,
 * where `import wgpu` is a natural thing to write.
 */

%include "IMP/bff/Compute.h"

%pythoncode %{
import os as _gpu_os

_GPU_TOLD = False


def _wgpu_library():
    """The wgpu-native shared library wgpu-py installs, or None.

    No operating system ships WebGPU: Windows has D3D12, macOS has Metal,
    Linux has Vulkan, and those are what WebGPU sits on. The PyPI package
    `wgpu` does ship it, for every platform this library targets, and puts it
    in a predictable place -- so if it is installed, the accelerator costs
    nothing to reach. Note the naming: Windows drops the `lib` prefix.
    """
    try:
        import wgpu
    except ImportError:
        return None
    resources = _gpu_os.path.join(_gpu_os.path.dirname(wgpu.__file__), "resources")
    for name in ("libwgpu_native-release.dylib", "libwgpu_native-release.so",
                 "wgpu_native-release.dll"):
        path = _gpu_os.path.join(resources, name)
        if _gpu_os.path.isfile(path):
            return path
    return None


def _gpu_plugin():
    """This module's own backend, which sits beside the extension."""
    here = _gpu_os.path.dirname(_gpu_os.path.abspath(__file__))
    for name in ("libimp_bff_wgpu.dylib", "libimp_bff_wgpu.so",
                 "imp_bff_wgpu.dll"):
        path = _gpu_os.path.join(here, name)
        if _gpu_os.path.isfile(path):
            return path
    return None


def enable_gpu(plugin=None, wgpu_library=None, quiet=False):
    """Look for an accelerator and use it if there is one.

    Called once at import. Returns the backend's name -- `cpu` when nothing
    was found, which is not a failure and not an error, just the answer.

    `IMP_BFF_GPU=off` stops it; `IMP_BFF_GPU_LIBRARY` names the plugin
    outright. Otherwise the plugin is looked for beside this module and the
    wgpu-native library in a wgpu-py installation.
    """
    global _GPU_TOLD
    if _gpu_os.environ.get("IMP_BFF_GPU", "").lower() in ("off", "0", "no"):
        return get_compute_backend_name()

    plugin = plugin or _gpu_os.environ.get("IMP_BFF_GPU_LIBRARY") or _gpu_plugin()
    if not plugin:
        return get_compute_backend_name()

    wgpu_library = wgpu_library or _gpu_os.environ.get("IMP_BFF_WGPU_LIBRARY") \
        or _wgpu_library() or ""
    if not wgpu_library and not _GPU_TOLD and not quiet:
        # Said once, and only when there is a plugin that could have used it:
        # a silent fallback is how somebody comes to believe they measured a
        # GPU when they measured a CPU.
        _GPU_TOLD = True
        import warnings
        warnings.warn(
            "IMP.bff found its GPU backend but no wgpu library to run it on; "
            "the kernels stay on the CPU. `pip install wgpu` provides one. "
            "Set IMP_BFF_GPU=off to stop looking.", RuntimeWarning, stacklevel=2)
        return get_compute_backend_name()

    if not load_compute_backend(plugin, wgpu_library) and not quiet:
        if not _GPU_TOLD:
            _GPU_TOLD = True
            import warnings
            warnings.warn("IMP.bff could not use its GPU backend (%s); the "
                          "kernels stay on the CPU."
                          % (get_compute_backend_error() or "no reason given"),
                          RuntimeWarning, stacklevel=2)
    return get_compute_backend_name()


try:
    enable_gpu()
except Exception:  # a backend must never stop the module from importing
    pass
%}
