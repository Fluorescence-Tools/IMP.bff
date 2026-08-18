"""Look a dye species up by name from the bundled tables.

One place that answers "what is AlexaFluor 488": the extinction coefficient and
quantum yield from ``Dyes_extinction_QD.csv``, the excitation/emission curves
from ``<name>.csv``, and -- when a template CIF is supplied -- the transition
dipole atoms, the chromophore centre and the formal charges.

Before PRD-113 those two halves were read by unrelated code
(``fret/forster.py`` and ``cgdye/io/template_cif.py``) and never met.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from .species import Dye, Spectrum
from .spectra import find_r0_file, normalize_dye_name, read_dye_table, read_spectrum

__all__ = ["find_dye", "available_dyes"]


def find_dye(
    name: str,
    r0_dir: str | Path | None = None,
    template_cif: str | Path | None = None,
) -> Dye:
    """Assemble a :class:`IMP.bff.dye.Dye` from the bundled data.

    :param name: e.g. ``"AlexaFluor 488"`` or ``"AlexaFluor488"``.
    :param r0_dir: optional directory holding the spectra and the dye table.
    :param template_cif: optional dye template, for the molecular half --
        transition-dipole atoms, centre atom, formal charges.
    :raises KeyError: if the dye is not in the table.
    """
    dye_type, number, compact = normalize_dye_name(name)
    table = read_dye_table(None if r0_dir is None else Path(r0_dir) / "Dyes_extinction_QD.csv")
    entry = table.get((dye_type, number))
    if entry is None:
        raise KeyError(
            f"{name!r} is not in the bundled dye table; "
            f"known names look like {sorted(table)[:3]}")

    spectrum: Optional[Spectrum] = None
    try:
        spectrum = read_spectrum(find_r0_file(f"{compact}.csv", r0_dir))
    except FileNotFoundError:
        # The table carries more dyes than the spectra directory does. A dye
        # without curves is still a usable species -- it just cannot derive R0.
        pass

    molecular = {}
    if template_cif is not None:
        from IMP.bff.cgdye.io.template_cif import read_dye_template_cif
        template = read_dye_template_cif(str(template_cif))
        d1, d2 = template.get("dipole_atom_1"), template.get("dipole_atom_2")
        molecular = {
            "dipole_atoms": (d1, d2) if d1 and d2 else None,
            "chromophore_center_atom": template.get("center_atom"),
            "positive_atoms": tuple(template.get("positive_atoms") or ()),
            "negative_atoms": tuple(template.get("negative_atoms") or ()),
        }

    return Dye(
        name=compact,
        spectrum=spectrum,
        extinction_coefficient=entry["Ext_coeff"],
        quantum_yield=entry["QD"],
        **molecular,
    )


def available_dyes(r0_dir: str | Path | None = None) -> list[str]:
    """Compact names of every dye in the bundled table."""
    table = read_dye_table(None if r0_dir is None else Path(r0_dir) / "Dyes_extinction_QD.csv")
    return sorted(f"{t}{n}" for t, n in table)
