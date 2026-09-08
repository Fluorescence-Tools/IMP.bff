#!/usr/bin/env python3
"""data/registry.json: the data files a wheel does not carry, with sha256 sums.

The wheel (package `bff`, PRD-137) ships the small part of data/ and fetches
`rotamer_library/` and `cgprobe/` on first use from IMP.bff.DATA_URL, checked
against this registry. Regenerate it after changing either directory:

    python utility/data_registry.py            # writes data/registry.json
    python utility/data_registry.py --verify   # checks the tree against it

The files themselves are uploaded to the download host by hand (rsync the
two directories, paths kept). test/test_data_registry.py keeps the registry
and the tree in step.
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


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--verify", action="store_true")
    a = ap.parse_args(argv)
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
