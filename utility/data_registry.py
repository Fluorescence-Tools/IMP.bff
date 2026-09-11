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
"""

import argparse
import hashlib
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FETCHED = ("rotamer_library", "cgprobe")
REGISTRY = os.path.join(DATA, "registry.json")
DATA_URL = os.environ.get("IMP_BFF_DATA_URL", "https://www.peulen.xyz/downloads/imp.bff/")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return "sha256:" + h.hexdigest()


def tree():
    out = {}
    for top in FETCHED:
        for dirpath, _, files in os.walk(os.path.join(DATA, top)):
            for f in sorted(files):
                if f.startswith("."):
                    continue
                p = os.path.join(dirpath, f)
                out[os.path.relpath(p, DATA).replace(os.sep, "/")] = sha256(p)
    return out


def fetch(target, cache=None, quiet=False, only=None):
    """Bring registry files into `target` (a data/ directory), through
    `cache` when given (pooch keeps its copy there; the file is then copied
    into place), checked against the registry's sums.

    `only` restricts the fetch to files under the given top-level data/
    directories (e.g. `["cgprobe"]`) -- CI fetches just what its lane
    reads, not the 229-file rotamer_library too."""
    import shutil
    import pooch
    with open(REGISTRY) as fh:
        registry = json.load(fh)
    if only:
        prefixes = tuple(f"{d}/" for d in only)
        registry = {k: v for k, v in registry.items() if k.startswith(prefixes)}
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
                         "(default: all). CI passes --only cgprobe: the "
                         "229-file rotamer_library is the bulk of the "
                         "download and the lane reads it through the "
                         "packaged copy, not the checkout.")
    a = ap.parse_args(argv)
    if a.fetch:
        n = fetch(a.fetch, a.cache, a.quiet, a.only)
        print("%d files in %s" % (n, a.fetch))
        return 0
    current = tree()
    if a.verify:
        with open(REGISTRY) as fh:
            stored = json.load(fh)
        missing = sorted(set(stored) - set(current))
        extra = sorted(set(current) - set(stored))
        changed = sorted(k for k in set(stored) & set(current) if stored[k] != current[k])
        for label, names in (("missing", missing), ("unlisted", extra), ("changed", changed)):
            for n in names:
                print(label, n)
        return 1 if (missing or extra or changed) else 0
    with open(REGISTRY, "w") as fh:
        json.dump(current, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print("%d files, %s" % (len(current), REGISTRY))
    return 0


if __name__ == "__main__":
    sys.exit(main())
