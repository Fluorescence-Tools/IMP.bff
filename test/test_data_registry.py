"""data/registry.json and the two data directories a wheel does not carry.

The pip package `bff` ships the small part of data/ and fetches
`rotamer_library/` and `cgprobe/` from IMP.bff.DATA_URL, checked against the
sha256 sums in data/registry.json (written by utility/data_registry.py). The
registry has to name exactly the files in the tree, with their current
sums, or a fetched file would be refused -- or a new file never fetched.

Also the search path: the standalone build reads IMP_BFF_DATA as a PATH-like
list, so a caller (the wheel's __init__) can put the shipped set and the
fetch cache side by side. The IMP module build reads IMP's own
get_data_path and is not asked.
"""

import hashlib
import json
import os
import random
import subprocess
import sys
import tempfile
import unittest

import IMP.bff

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_DATA = os.path.join(_ROOT, "data")
_REGISTRY = os.path.join(_DATA, "registry.json")


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return "sha256:" + h.hexdigest()


class TestDataRegistry(unittest.TestCase):

    def setUp(self):
        if not os.path.isdir(os.path.join(_DATA, "rotamer_library")):
            self.skipTest("not a checkout with data/rotamer_library")
        with open(_REGISTRY) as fh:
            self.registry = json.load(fh)

    def test_the_registry_names_the_tree(self):
        tree = set()
        for top in ("rotamer_library", "cgprobe"):
            for dirpath, _, files in os.walk(os.path.join(_DATA, top)):
                for f in files:
                    if not f.startswith("."):
                        tree.add(os.path.relpath(os.path.join(dirpath, f), _DATA).replace(os.sep, "/"))
        self.assertEqual(set(self.registry), tree,
                         "run `python utility/data_registry.py` after changing data/")

    def test_a_sample_of_sums_is_current(self):
        # all 284 is 62 MB of hashing; a fixed sample catches an edited file
        # in a handful of runs and the generator's --verify does the rest
        rng = random.Random(20260908)
        for name in rng.sample(sorted(self.registry), 8):
            self.assertEqual(self.registry[name], _sha256(os.path.join(_DATA, name)), name)

    def test_the_two_largest_files_are_current(self):
        for name in ("rotamer_library/dyes.drot.pto", "rotamer_library/sidechains.drot.pto"):
            self.assertEqual(self.registry[name], _sha256(os.path.join(_DATA, name)), name)


class TestSearchPath(unittest.TestCase):

    def test_imp_bff_data_is_a_path_list(self):
        if IMP.bff.get_build() != "core":
            self.skipTest("the IMP module build resolves data through IMP")
        # a fresh interpreter, so the environment is read cold
        with tempfile.TemporaryDirectory() as tmp:
            a, b = os.path.join(tmp, "a"), os.path.join(tmp, "b")
            os.makedirs(a)
            os.makedirs(b)
            with open(os.path.join(b, "only_in_b.txt"), "w") as fh:
                fh.write("x")
            env = dict(os.environ, IMP_BFF_DATA=os.pathsep.join([a, b]))
            out = subprocess.run(
                [sys.executable, "-c",
                 "import IMP.bff; print(IMP.bff.get_data_path('only_in_b.txt'))"],
                env=env, capture_output=True, text=True)
            self.assertEqual(out.returncode, 0, out.stderr[-2000:])
            self.assertEqual(out.stdout.strip(), os.path.join(b, "only_in_b.txt"))


class TestBuildIdentity(unittest.TestCase):

    def test_get_build_names_this_build(self):
        # "core+imp" is the standalone build that links IMP as a private
        # library: no IMP in Python, but the connection layer is there
        # (PRD-139).
        build = IMP.bff.get_build()
        self.assertIn(build, ("core", "core+imp", "imp"))
        self.assertEqual(build.startswith("core"),
                         bool(getattr(IMP.bff, "IMPBFF_STANDALONE", False)))
        self.assertEqual(build in ("core+imp", "imp"),
                         hasattr(IMP.bff, "MolecularProbeSimulation"))

    def test_fetch_data_and_registry_exist_in_both_builds(self):
        self.assertIsInstance(IMP.bff.get_data_registry(), dict)
        self.assertEqual(IMP.bff.fetch_data([]), [])


if __name__ == "__main__":
    unittest.main()
