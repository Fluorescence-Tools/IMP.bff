"""Look a dye species up by name from the bundled tables.

One place that answers "what is AlexaFluor 488": the extinction coefficient and
quantum yield and the excitation/emission curves from the bundled
``dye_library.cif``, and -- when a template CIF is supplied -- the transition
dipole atoms, the chromophore centre and the formal charges.

Before PRD-113 those two halves were read by unrelated code
(``fret/forster.py`` and ``cgdye/io/template_cif.py``) and never met, and the
spectral half was CSV. **CIF is the data format now** -- see
:mod:`IMP.bff.dye.cif` for the categories and why they are bff-native.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Optional

from .cif import read_dye_library
from .species import Dye, Spectrum
from .spectra import normalize_dye_name

__all__ = ["find_dye", "available_dyes"]


def find_dye(
    name: str,
    library_cif: str | Path | None = None,
    template_cif: str | Path | None = None,
) -> Dye:
    """Assemble a :class:`IMP.bff.dye.Dye` from the bundled data.

    :param name: e.g. ``"AlexaFluor 488"`` or ``"AlexaFluor488"``.
    :param library_cif: optional path to a dye-library CIF.
    :param template_cif: optional dye template, for the molecular half --
        transition-dipole atoms, centre atom, formal charges.
    :raises KeyError: if the dye is not in the table.
    """
    _dye_type, _number, compact = normalize_dye_name(name)
    library = read_dye_library(library_cif)
    dye = library.get(compact)
    if dye is None:
        raise KeyError(
            f"{name!r} is not in the bundled dye library; "
            f"known names look like {sorted(library)[:3]}")

    if template_cif is None:
        return dye

    from IMP.bff.cgdye.io.template_cif import read_dye_template_cif
    template = read_dye_template_cif(str(template_cif))
    d1, d2 = template.get("dipole_atom_1"), template.get("dipole_atom_2")
    return replace(
        dye,
        dipole_atoms=(d1, d2) if d1 and d2 else None,
        chromophore_center_atom=template.get("center_atom"),
        positive_atoms=tuple(template.get("positive_atoms") or ()),
        negative_atoms=tuple(template.get("negative_atoms") or ()),
    )


def available_dyes(library_cif: str | Path | None = None) -> list[str]:
    """Chromophore names of every dye in the bundled library."""
    return sorted(read_dye_library(library_cif))
