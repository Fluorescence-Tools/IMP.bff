"""Utility helpers for cgdye."""

from __future__ import annotations

import os
from pathlib import Path


def _data_root() -> Path:
    """Return the cgdye data directory.

    Templates, input structures and restraint files are IMP module *data*, not
    package sources: they live in ``imp.bff/data/cgdye`` and are reached through
    ``IMP.bff.get_data_path``. Deriving them from ``__file__`` instead would tie
    them to where the package happens to sit, which is exactly what broke when
    cgdye moved out of imp-tricks.
    """
    import IMP.bff
    return Path(IMP.bff.get_data_path("cgdye"))


def _join_parts(parts: tuple[str, ...]) -> Path:
    """Join path parts, handling nested path strings."""
    if not parts:
        return Path()
    result = Path()
    for part in parts:
        result = result / Path(part)
    return result


def get_template_dir(*parts: str) -> Path:
    """Return the template directory, optionally with subpath parts."""
    return _data_root() / "templates" / _join_parts(parts) if parts else _data_root() / "templates"


def get_structure_dir(*parts: str) -> Path:
    """Return the structures directory, optionally with subpath parts."""
    return _data_root() / "inputs" / "structures" / _join_parts(parts) if parts else _data_root() / "inputs" / "structures"


def get_output_dir(*parts: str) -> Path:
    """Return the output directory, optionally with subpath parts.

    Relative to the working directory, not to the package: module data is
    read-only and installed, so nothing may be written beside it.
    """
    return Path.cwd() / "output" / _join_parts(parts) if parts else Path.cwd() / "output"


def ensure_dir(path: Path) -> Path:
    """Ensure a directory exists and return it."""
    path.mkdir(parents=True, exist_ok=True)
    return path
