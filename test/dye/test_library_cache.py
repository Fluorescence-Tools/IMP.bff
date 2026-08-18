"""The dye library is parsed once per file, not once per lookup.

``read_dye_library`` parses a CIF. Profiling one FRETpredict comparison found it
called **21 times for the same file** -- 550 000 calls to the row cleaner and
1.65 million to the float converter, 2.2 s of a 13 s run -- because every
consumer that wanted one spectrum read the whole library.

It is now cached on ``(resolved path, mtime, size)``. Not on the path alone:
editing the file during a session must be picked up, and only an edit changing
neither timestamp nor size would be missed, which is not something that happens
to a data file.

Each call gets a **fresh mapping** over the **same** ``Dye`` values, so the
identity to test is the values'. A caller that adds or drops an entry edits only
its own dict; the species are shared, which is safe because ``Dye`` is frozen.
Reaching around that with ``object.__setattr__`` corrupts the library for every
later reader -- which is what a test in ``test_photophysics_terms`` was doing on
the day this cache landed, and it only surfaced because the cache made the
corruption persist.
"""

import shutil
import time
from pathlib import Path

import pytest

import IMP.bff
from IMP.bff.dye.cif import DYE_LIBRARY_CIF, _LIBRARY_CACHE, read_dye_library


def _bundled_path():
    return Path(IMP.bff.get_data_path("rotamer_library")) / "R0" / DYE_LIBRARY_CIF


def test_the_bundled_library_is_parsed_once():
    first = read_dye_library()
    second = read_dye_library()
    assert len(first) > 0
    name = next(iter(first))
    assert first[name] is second[name], "a second read must not re-parse"


def test_each_caller_gets_its_own_mapping():
    """So dropping an entry is a local act, not an edit to the library."""
    first = read_dye_library()
    second = read_dye_library()
    assert first is not second
    first.pop(next(iter(first)))
    assert len(read_dye_library()) == len(second)


def test_repeat_reads_are_much_faster():
    read_dye_library()                      # warm
    t0 = time.perf_counter()
    for _ in range(20):
        read_dye_library()
    cached = time.perf_counter() - t0
    assert cached < 0.5, f"20 cached reads took {cached:.3f}s"


def test_an_edited_file_is_reread(tmp_path):
    copy = tmp_path / "dye_library.cif"
    shutil.copy(_bundled_path(), copy)

    first = read_dye_library(copy)
    name = next(iter(first))
    assert read_dye_library(copy)[name] is first[name]

    with copy.open("a") as handle:      # changes both mtime and size
        handle.write("\n#\n")
    assert read_dye_library(copy)[name] is not first[name], \
        "an edited file must be parsed again"


def test_two_files_do_not_share_a_cache_entry(tmp_path):
    a, b = tmp_path / "a.cif", tmp_path / "b.cif"
    shutil.copy(_bundled_path(), a)
    shutil.copy(_bundled_path(), b)
    name = next(iter(read_dye_library(a)))
    assert read_dye_library(a)[name] is not read_dye_library(b)[name]
    assert read_dye_library(a)[name] is read_dye_library(a)[name]


def test_a_missing_file_is_not_cached(tmp_path):
    before = len(_LIBRARY_CACHE)
    with pytest.raises((FileNotFoundError, OSError)):
        read_dye_library(tmp_path / "nope.cif")
    assert len(_LIBRARY_CACHE) == before


def test_the_atom_type_lookup_is_memoised():
    """Thousands of atoms, a few dozen distinct names, five possible answers."""
    from IMP.bff.scoring.rotamer import _atom_type
    assert hasattr(_atom_type, "cache_info")
    _atom_type.cache_clear()
    for _ in range(500):
        _atom_type("CB")
    info = _atom_type.cache_info()
    assert info.hits == 499 and info.misses == 1
    assert _atom_type("SG") == "S" and _atom_type("XX") == "C"


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
