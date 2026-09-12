#!/usr/bin/env python3
"""data/registry.json: the data files a wheel does not carry, with sha256 sums.

The wheel (package `bff`, PRD-137) ships the small part of data/ and fetches
`rotamer_library/` and `cgprobe/` on first use from IMP.bff.DATA_URL, checked
against this registry. Regenerate it after changing either directory:

    python utility/data_registry.py            # writes data/registry.json
    python utility/data_registry.py --verify   # checks the tree against it

The two directories are not in git: they live on the download host
(DATA_URL, https://www.peulen.xyz/downloads/imp.bff/, uploaded by hand with
paths kept) and

    python utility/data_registry.py --fetch            # restores data/
    python utility/data_registry.py --fetch --cache D  # via a cache dir (CI)

brings them into a checkout, checked against the sums. Needs pooch.
test/test_data_registry.py keeps the registry and a fetched tree in step.

`data/academic/` is a third, opt-in set: files whose licence lets a user
download them for academic use but does not let this project redistribute
them (today: FASPR's `dun2010bbdep.bin`, the Dunbrack 2010 side-chain
library). They are named in the registry like everything else -- the sums
have to live somewhere -- but nothing fetches them unless BFF_ACADEMIC is
set in the environment, and a build that packages data/ never carries them
(the wheel and sdist exclude the directory outright). The flag is the
user's assertion of the terms:

    BFF_ACADEMIC=1 python utility/data_registry.py --fetch --only academic
"""

import argparse
import hashlib
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FETCHED = ("rotamer_library", "cgprobe")
#: Opt-in, licence-gated: fetched only with BFF_ACADEMIC set.
ACADEMIC = ("academic",)
REGISTRY = os.path.join(DATA, "registry.json")
DATA_URL = os.environ.get("IMP_BFF_DATA_URL", "https://www.peulen.xyz/downloads/imp.bff/")


def academic_on():
    """BFF_ACADEMIC set to anything but an explicit negative."""
    return os.environ.get("BFF_ACADEMIC", "").strip().lower() not in (
        "", "0", "false", "no")


def is_academic(name):
    return name.startswith("academic/")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return "sha256:" + h.hexdigest()


def tree():
    """The fetched tree as the checkout happens to hold it.

    `academic/` is walked only when present: a checkout that has not asked
    for the opt-in set must still regenerate a registry that keeps its sums.
    """
    out = {}
    for top in FETCHED + ACADEMIC:
        if top in ACADEMIC and not os.path.isdir(os.path.join(DATA, top)):
            continue
        for dirpath, _, files in os.walk(os.path.join(DATA, top)):
            for f in sorted(files):
                if f.startswith("."):
                    continue
                p = os.path.join(dirpath, f)
                out[os.path.relpath(p, DATA).replace(os.sep, "/")] = sha256(p)
    return out


def active(registry, only=None):
    """The registry entries a fetch takes: the two shipped sets, plus the
    academic set only under BFF_ACADEMIC, cut to `only`'s directories when
    the caller named any."""
    tops = set(FETCHED)
    if academic_on():
        tops |= set(ACADEMIC)
    if only:
        tops &= set(only)
    return {k: v for k, v in registry.items()
            if k.split("/", 1)[0] in tops}


def fetch(target, cache=None, quiet=False, only=None):
    """Bring registry files into `target` (a data/ directory), through
    `cache` when given (pooch keeps its copy there; the file is then copied
    into place), checked against the registry's sums.

    `only` restricts the fetch to files under the given top-level data/
    directories (e.g. `["cgprobe"]`) -- CI fetches just what its lane
    reads, not the 229-file rotamer_library too. `academic/` arrives only
    when BFF_ACADEMIC is set, with or without `only`."""
    import shutil
    import pooch
    with open(REGISTRY) as fh:
        registry = active(json.load(fh), only)
    pup = pooch.create(path=cache or target, base_url=DATA_URL,
                       registry=registry)
    n = 0
    for name in sorted(registry):
        src = pup.fetch(name, progressbar=not quiet)
        dst = os.path.join(target, name)
        if os.path.abspath(src) != os.path.abspath(dst):
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copyfile(src, dst)
        n += 1
    return n


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--verify", action="store_true", help="check the tree against the registry")
    ap.add_argument("--fetch", nargs="?", const=DATA, metavar="DIR",
                    help="fetch every registry file into DIR (default: data/)")
    ap.add_argument("--cache", metavar="DIR", help="with --fetch: pooch's download directory")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--only", nargs="+", metavar="DIR",
                    help="with --fetch: top-level data/ directories to take "
                         "(default: all the fetched ones, academic only under "
                         "BFF_ACADEMIC). CI passes --only cgprobe academic "
                         "with BFF_ACADEMIC=1: the 229-file rotamer_library "
                         "is the bulk of the download and the lane reads it "
                         "through the packaged copy, not the checkout.")
    a = ap.parse_args(argv)
    if a.fetch:
        n = fetch(a.fetch, a.cache, a.quiet, a.only)
        print("%d files in %s" % (n, a.fetch))
        return 0
    with open(REGISTRY) as fh:
        stored = json.load(fh)
    # academic/ is one opt-in set: either the whole directory is fetched and
    # checked, or it is absent and its entries are simply not on this machine.
    if os.path.isdir(os.path.join(DATA, "academic")):
        current = tree()
    else:
        current = {k: v for k, v in tree().items() if not is_academic(k)}
        stored = {k: v for k, v in stored.items() if not is_academic(k)}
    if a.verify:
        missing = sorted(set(stored) - set(current))
        extra = sorted(set(current) - set(stored))
        changed = sorted(k for k in set(stored) & set(current) if stored[k] != current[k])
        for label, names in (("missing", missing), ("unlisted", extra), ("changed", changed)):
            for n in names:
                print(label, n)
        return 1 if (missing or extra or changed) else 0
    if not os.path.isdir(os.path.join(DATA, "academic")):
        # keep the opt-in set's sums even when regenerating from a checkout
        # that has not fetched them
        current.update({k: v for k, v in stored.items() if is_academic(k)})
    with open(REGISTRY, "w") as fh:
        json.dump(current, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print("%d files, %s" % (len(current), REGISTRY))
    return 0


if __name__ == "__main__":
    sys.exit(main())
