"""The approved hard migration: paths, exports and direct container ownership."""

import json
from pathlib import Path

import pytest
import IMP.bff as bff

ROOT = Path(__file__).resolve().parents[1]
MAP = json.loads((Path(__file__).parent / "api_taxonomy.json").read_text())


@pytest.mark.parametrize("old,new_paths", MAP["paths"].items())
def test_canonical_paths_replace_retired_paths(old, new_paths):
    old_path = ROOT / old
    # Directory entries also test canonical acronym casing on macOS.
    assert old_path.name not in {p.name for p in old_path.parent.iterdir()}
    for new in new_paths:
        assert (ROOT / new).is_file(), new


@pytest.mark.parametrize("old,new", MAP["symbols"].items())
def test_no_retired_export_survives(old, new):
    assert not hasattr(bff, old), old


@pytest.mark.parametrize("old", MAP["exported_before"])
def test_migrated_exports_resolve_and_have_documentation(old):
    new = MAP["symbols"][old]
    value = getattr(bff, new)
    if callable(value):
        assert value.__doc__, new


def test_container_implementation_is_owned_by_ptolib():
    for name in ("PtoReader", "PtoWriter", "PtoObject"):
        assert not hasattr(bff, name)
    assert (ROOT / "src/internal/Ptolib.cpp").is_file()
    assert not (ROOT / "include/Pto.h").exists()
    assert not (ROOT / "src/Pto.cpp").exists()
