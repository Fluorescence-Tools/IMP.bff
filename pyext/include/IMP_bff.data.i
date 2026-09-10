/*
 * Where the large data lives (PRD-137 step 6d): not in git, not in a wheel.
 *
 * data/ in the repository holds the small, always-needed files and
 * data/registry.json -- the path and sha256 of every file of
 * rotamer_library/ and cgprobe/ (62 MB), which are served from DATA_URL
 * (utility/data_registry.py --fetch restores them into a checkout). A wheel
 * ships the small set beside the module; a conda package ships everything
 * (its build fetches first). At run time get_data_path() looks where the
 * build installed, then in the fetch cache, and fetches a registry file the
 * first time it is asked for. The same code serves the standalone module and
 * the IMP module build: it wraps the extension's get_data_path, which both
 * shadow modules call.
 */
%pythoncode %{
import os as _os

DATA_URL = _os.environ.get("IMP_BFF_DATA_URL", "https://www.peulen.xyz/downloads/imp.bff/")


def _package_data_dir():
    """The wheel's own IMP/bff/data, or None."""
    d = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "data")
    return d if _os.path.isdir(d) else None


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


_installed_get_data_path = _IMP_bff.get_data_path


def get_data_registry():
    """The files fetched on demand: relative path -> "sha256:<hex>"."""
    import json
    candidates = []
    pkg = _package_data_dir()
    if pkg:
        candidates.append(_os.path.join(pkg, "registry.json"))
    try:
        candidates.append(_installed_get_data_path("registry.json"))
    except Exception:
        pass
    for p in candidates:
        try:
            with open(p) as fh:
                return json.load(fh)
        except OSError:
            continue
    return {}


def fetch_data(file_names=None, progressbar=False):
    """Fetch registry files into the cache (all of them when `file_names` is
    None) and return their paths. Needs `pooch`."""
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


def _fetching_get_data_path(file_name):
    try:
        return _installed_get_data_path(file_name)
    except IOException:
        cached = _os.path.join(get_data_cache_dir(), file_name)
        if _os.path.isfile(cached):
            return cached
        if file_name not in get_data_registry():
            raise
    fetch_data([file_name])
    return _os.path.join(get_data_cache_dir(), file_name)


_fetching_get_data_path.__doc__ = """The path of a data file: as installed, else
from the fetch cache, fetched first when it is one of the registry's
(`fetch_data()` takes them all at once, e.g. before going offline)."""
_IMP_bff.get_data_path = _fetching_get_data_path

# the standalone module reads IMP_BFF_DATA as a PATH-like list: the caller's
# directories, the wheel's own data, the fetch cache
if _package_data_dir():
    _own = _os.environ.get("IMP_BFF_DATA", "")
    _os.environ["IMP_BFF_DATA"] = _os.pathsep.join(
        ([_own] if _own else []) + [_package_data_dir(), get_data_cache_dir()])


# IMP's own data (top.lib, par.lib, the element table) where this build links
# IMP and carries it: the C++ side reads IMP_DATA, so a wheel points it at
# what it ships unless the caller has already chosen (PRD-139).
_IMP_DATA = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "imp_data")
if _os.path.isdir(_IMP_DATA) and not _os.environ.get("IMP_DATA"):
    _os.environ["IMP_DATA"] = _IMP_DATA


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
%}
