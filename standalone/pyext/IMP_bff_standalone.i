/*
 * The standalone build's SWIG entry: IMP.bff without IMP (PRD-137 step 6c).
 *
 * IMP's module tooling writes a module's entry around swig.i-in -- the kernel
 * interface, every dependency's, the exception translation, numpy. This file
 * is that entry for the core alone: the preamble, the IMP_SWIG_* equivalents
 * (IMP_bff_standalone.macros.i), and then core.i verbatim -- the same file
 * the IMP build wraps before layer.i. The generated shadow module is
 * installed as IMP/bff/__init__.py next to _IMP_bff.
 */
%module(directors="1", moduleimport="from . import _IMP_bff") IMP_bff
%feature("autodoc", 1);
%feature("kwargs", 1);

%{
#define SWIG_FILE_WITH_INIT
#include <IMP/bff/bff_config.h>
#include <IMP/bff/Base.h>
#include "bff_core_headers.h"
%}

/* The SWIG side of the module's vocabulary: the export/namespace macros and
   Base.h's standalone branch (IMP_VALUES, IMP_SHOWABLE_INLINE, ...). IMP's
   kernel interface provides these to a module's SWIG; here the headers do. */
/* numpy.i (pulled in by core.i) needs the C API initialised once per module;
   IMP's kernel does this in IMP_kernel.import_numpy.i. */
%{
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
%}
%init %{
    import_array();
%}

/* Stream insertion is a C++ convenience (IMP_SHOWABLE_INLINE writes it as a
   friend); Python has show()/__str__, never __lshift__. */
%ignore operator<<;
%ignore IMP::operator<<;
%ignore IMP::bff::operator<<;

%include <std_string.i>
%implicitconv;  /* a dict is a ProbeMap, a list an FFComponentSpecList -- as under IMP's tooling */

/* The IMP_SWIG_* equivalents and, first of all, the exception translation:
   %exception covers only what is wrapped after it, and bff_config.h
   (get_data_path, which throws) comes next. */
%include "IMP_bff_standalone.macros.i"

%include <IMP/bff/bff_config.h>
/* Base.h hides its C++ exception classes from SWIG; the module makes Python
   ones (IMP_bff_standalone.macros.i), as IMP's kernel does. */
%ignore IMP::Pointer;
%ignore IMP::Object::ref;
%ignore IMP::Object::unref;
%ignore IMP::Object::release_ref;
%include <IMP/Object.h>
%include <IMP/bff/Base.h>

%pythoncode %{
# The standalone core: no IMP.atom, no IMP.Model. What is here is what
# fits fluorescence data, samples labels and reads structures on its own.
IMPBFF_STANDALONE = True
%}

%include "IMP_bff.core.i"

%pythoncode %{
def get_module_version():
    return _IMP_bff.get_module_version()

def get_module_name():
    return "IMP::bff"

def get_build():
    """Which IMP.bff this is: "core" (no IMP, this build) or "imp" (the IMP
    module, which adds the connection layer)."""
    return "core"

# ---- where the data is (PRD-137 step 6d) ----------------------------------
# A wheel ships the small, always-needed part of data/ beside this file and
# lists the rest (rotamer_library/, cgprobe/: 62 MB) in data/registry.json
# with sha256 sums; those are fetched on first use into a cache and looked
# up from there. The C++ side reads IMP_BFF_DATA as a PATH-like list, so the
# search order is: whatever the caller put in IMP_BFF_DATA, the wheel's own
# data/, the fetch cache. A conda or plain install has all of data/ in the
# directory the build was told about and nothing here applies.
import os as _os

DATA_URL = _os.environ.get("IMP_BFF_DATA_URL", "https://www.peulen.xyz/downloads/imp-bff-data/")
_PACKAGE_DATA = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "data")
_REGISTRY_FILE = _os.path.join(_PACKAGE_DATA, "registry.json")


def get_data_cache_dir():
    """Where fetched data files live: IMP_BFF_CACHE, else pooch's per-user
    cache directory for "imp.bff"."""
    d = _os.environ.get("IMP_BFF_CACHE")
    if d:
        return d
    try:
        import pooch
        return str(pooch.os_cache("imp.bff"))
    except ImportError:
        return _os.path.join(_os.path.expanduser("~"), ".cache", "imp.bff")


def get_data_registry():
    """The files fetched on demand: relative path -> "sha256:<hex>"; empty
    when this install carries all of its data."""
    try:
        with open(_REGISTRY_FILE) as fh:
            import json
            return json.load(fh)
    except OSError:
        return {}


def fetch_data(file_names=None, progressbar=False):
    """Fetch data files listed in the registry into the cache (all of them
    when `file_names` is None) and return their paths. Needs `pooch`."""
    registry = get_data_registry()
    if not registry:
        return []
    names = list(registry) if file_names is None else list(file_names)
    unknown = [n for n in names if n not in registry]
    if unknown:
        raise IOException("not in the data registry: " + ", ".join(unknown))
    import pooch
    pup = pooch.create(path=get_data_cache_dir(), base_url=DATA_URL, registry=registry)
    return [pup.fetch(n, progressbar=progressbar) for n in names]


if _os.path.isdir(_PACKAGE_DATA):
    _own = _os.environ.get("IMP_BFF_DATA", "")
    _os.environ["IMP_BFF_DATA"] = _os.pathsep.join(
        ([_own] if _own else []) + [_PACKAGE_DATA, get_data_cache_dir()])

_get_data_path_installed = get_data_path


def get_data_path(file_name):
    """The path of a data file: as installed, or fetched into the cache when
    it is one of the registry's (the first time; `fetch_data()` takes them
    all at once, e.g. before going offline)."""
    try:
        return _get_data_path_installed(file_name)
    except IOException:
        if file_name not in get_data_registry():
            raise
    fetch_data([file_name])
    return _get_data_path_installed(file_name)


def _fetch_data_main(argv=None):
    """`imp_bff_fetch_data`: fetch every registry file, or the ones named."""
    import argparse
    ap = argparse.ArgumentParser(description=_fetch_data_main.__doc__)
    ap.add_argument("names", nargs="*", help="registry entries; all when none")
    ap.add_argument("--list", action="store_true", help="print the registry and exit")
    ap.add_argument("--cache", help="fetch into this directory (default: %s)" % get_data_cache_dir())
    a = ap.parse_args(argv)
    if a.cache:
        _os.environ["IMP_BFF_CACHE"] = a.cache
    if a.list:
        for k, v in sorted(get_data_registry().items()):
            print(k, v)
        return 0
    for p in fetch_data(a.names or None, progressbar=True):
        print(p)
    return 0


def _warn_if_this_replaced_the_imp_module():
    # The pip core and the conda imp.bff package both own site-packages/IMP/bff/.
    # When an IMP kernel is importable next to this core, the pip install has
    # overwritten the IMP module's files: the connection layer's names are gone
    # until the conda package is reinstalled.
    import importlib.util
    try:
        found = importlib.util.find_spec("IMP.atom") is not None
    except Exception:
        found = False
    if found:
        import warnings
        warnings.warn(
            "IMP.bff is the IMP-free core (package 'bff'), but IMP itself is "
            "installed here: the connection layer (get_av_from_structure, the "
            "restraints, ...) is unavailable. `conda install imp.bff` restores "
            "the IMP module build.", RuntimeWarning, stacklevel=2)


_warn_if_this_replaced_the_imp_module()
%}
