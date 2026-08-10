"""Förster radius calculations for rotamer FRET."""

from __future__ import annotations

import csv
import re
from collections.abc import Iterable
from importlib import resources
from pathlib import Path

import numpy as np

_DYE_NAME_RE = re.compile(r"^(?P<type>.+?)\s+(?P<number>[A-Za-z0-9]+)$")


# numpy 2 removed np.trapz in favour of np.trapezoid, and this module is the
# Foerster-radius calculation -- so under the numpy the stack actually runs
# (2.4) calculate_r0 raised AttributeError and R0 could not be computed at all.
try:  # numpy >= 2
    _trapezoid = np.trapezoid
except AttributeError:  # pragma: no cover - numpy < 2
    _trapezoid = np.trapz

def _package_resource_path(*parts: str) -> Path:
    """Return a path to a bundled resource file.

    Parameters
    ----------
    *parts : str
        Relative path parts below the FRETpredict package data directory.

    Returns
    -------
    pathlib.Path
        Existing resource path.
    """
    # The rotamer library ships as IMP.bff module data (data/rotamer_library),
    # not as a vendored package: it is data the rotamer route loads at run
    # time, so it moved into the module rather than staying in junk/.
    import IMP.bff
    root = Path(IMP.bff.get_data_path("rotamer_library"))
    parts = tuple(p for p in parts if p != "lib")
    return root.joinpath(*parts)


def find_r0_file(filename: str, r0_dir: str | Path | None = None) -> Path:
    """Find a bundled R0 data file.

    Parameters
    ----------
    filename : str
        File name, for example ``AlexaFluor488.csv``.
    r0_dir : pathlib.Path or str, optional
        Optional directory containing R0 CSV files.

    Returns
    -------
    pathlib.Path
        Path to the requested file.
    """
    if r0_dir is not None:
        path = Path(r0_dir) / filename
        if path.exists():
            return path
    return _package_resource_path("lib", "R0", filename)


def normalize_dye_name(dye_name: str) -> tuple[str, str, str]:
    """Normalize a FRETpredict-style dye name.

    Parameters
    ----------
    dye_name : str
        Dye name such as ``AlexaFluor 488`` or ``ATTO Thio12``.

    Returns
    -------
    tuple[str, str, str]
        ``(type_name, dye_number, compact_name)``.
    """
    match = _DYE_NAME_RE.match(str(dye_name).strip())
    if match is None:
        compact = str(dye_name).replace(" ", "")
        return compact, compact, compact
    dye_type = match.group("type").strip()
    number = match.group("number").strip()
    return dye_type, number, f"{dye_type}{number}"


def read_dye_table(path: str | Path | None = None) -> dict[tuple[str, str], dict[str, float]]:
    """Read the bundled dye extinction and quantum-yield table.

    Parameters
    ----------
    path : pathlib.Path or str, optional
        Optional path to ``Dyes_extinction_QD.csv``.

    Returns
    -------
    dict
        Mapping ``(type, chromophore)`` to ``Ext_coeff`` and ``QD`` values.
    """
    table_path = Path(path) if path is not None else find_r0_file("Dyes_extinction_QD.csv")
    table: dict[tuple[str, str], dict[str, float]] = {}
    with table_path.open(newline="") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if len(row) < 4:
                continue
            dye_type, chromophore, ext_coeff, qd = row[:4]
            table[(dye_type.strip(), chromophore.strip())] = {
                "Ext_coeff": float(ext_coeff),
                "QD": float(qd),
            }
    return table


def _read_spectrum(path: Path) -> np.ndarray:
    """Read a normalized donor or acceptor spectrum.

    Parameters
    ----------
    path : pathlib.Path
        Spectrum CSV path.

    Returns
    -------
    numpy.ndarray
        Structured array with ``Wavelength``, ``Excitation``, and ``Emission``.
    """
    data: list[tuple[float, float, float]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            data.append(
                (
                    float(row["Wavelength"]),
                    float(row["Excitation"]) / 100.0,
                    float(row["Emission"]) / 100.0,
                )
            )
    return np.array(data, dtype=[("Wavelength", float), ("Excitation", float), ("Emission", float)])


def calculate_r0(
    donor: str,
    acceptor: str,
    k2: float,
    r0_dir: str | Path | None = None,
) -> float:
    """Calculate the Förster radius for a dye pair.

    Parameters
    ----------
    donor : str
        Donor dye name, for example ``AlexaFluor 488``.
    acceptor : str
        Acceptor dye name, for example ``AlexaFluor 594``.
    k2 : float
        Orientation factor.
    r0_dir : pathlib.Path or str, optional
        Optional directory containing R0 CSV files.

    Returns
    -------
    float
        Förster radius in nm.
    """
    donor_type, donor_number, donor_name = normalize_dye_name(donor)
    _acceptor_type, acceptor_number, acceptor_name = normalize_dye_name(acceptor)

    dye_table = read_dye_table()
    donor_data = dye_table.get((donor_type, donor_number))
    acceptor_data = dye_table.get((_acceptor_type, acceptor_number))
    if donor_data is None or acceptor_data is None:
        raise ValueError(f"No R0 dye data found for {donor!r} / {acceptor!r}")

    donor_spectrum = _read_spectrum(find_r0_file(f"{donor_name}.csv", r0_dir))
    acceptor_spectrum = _read_spectrum(find_r0_file(f"{acceptor_name}.csv", r0_dir))

    wavelengths = donor_spectrum["Wavelength"].astype(float)
    donor_emission = donor_spectrum["Emission"].astype(float)
    acceptor_excitation = acceptor_spectrum["Excitation"].astype(float)

    if donor_emission.size != wavelengths.size or acceptor_excitation.size != wavelengths.size:
        raise ValueError("Donor and acceptor spectra must share the same wavelength grid")

    emission_integral = _trapezoid(donor_emission, x=wavelengths)
    if emission_integral == 0:
        return 0.0

    ext_coeff_max = acceptor_data["Ext_coeff"]
    ext_coeff_acceptor = ext_coeff_max * acceptor_excitation
    overlap = _trapezoid(donor_emission * ext_coeff_acceptor * np.power(wavelengths, 4), x=wavelengths)
    overlap /= emission_integral

    factor = 0.02108
    refractive_index_4 = 1.4**4
    quantum_yield = donor_data["QD"]
    return float(factor * np.power(k2 * quantum_yield / refractive_index_4 * overlap, 1.0 / 6.0))


def iter_r0_pairs(r0_dir: str | Path | None = None) -> Iterable[dict[str, object]]:
    """Iterate precomputed R0 pairs if available.

    Parameters
    ----------
    r0_dir : pathlib.Path or str, optional
        Optional directory containing ``R0_pairs.csv``.

    Yields
    ------
    dict
        Precomputed R0 row.
    """
    path = Path(r0_dir) / "R0_pairs.csv" if r0_dir is not None else find_r0_file("R0_pairs.csv")
    if not path.exists():
        return
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        yield from reader
