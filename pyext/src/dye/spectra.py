"""Spectra, the overlap integral, and the Förster radius derived from them.

**R0 is derived, not supplied.** It follows from two dyes' spectra, the donor's
quantum yield, the solvent's refractive index and kappa^2 -- so it is a property
of a pair in a medium, not a number to be passed around. Before PRD-113 this
module lived in ``fret/`` and every consumer took ``forster_radius=52.0`` as a
default argument instead; ``refractive_index`` was not even a parameter, it was
the literal ``1.4**4`` inside the calculation.

Lives in ``dye`` because the inputs are dye species properties. The bundled
spectra are module data (``data/rotamer_library/R0``, the FRETpredict tables).
"""

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
# (2.4) forster_radius_from_spectra raised AttributeError and R0 could not be computed at all.
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


def read_spectrum(path: Path) -> "Spectrum":
    """Read a normalized donor or acceptor spectrum.

    Parameters
    ----------
    path : pathlib.Path
        Spectrum CSV path.

    Returns
    -------
    IMP.bff.dye.Spectrum
        Excitation and emission on the file's shared wavelength grid.
    """
    from .species import Spectrum
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
    array = np.array(
        data, dtype=[("Wavelength", float), ("Excitation", float), ("Emission", float)])
    return Spectrum(
        wavelength=array["Wavelength"].astype(float),
        excitation=array["Excitation"].astype(float),
        emission=array["Emission"].astype(float),
    )


#: Numerical factor of the Foerster expression with R0 in nm, the overlap
#: integral in M^-1 cm^-1 nm^4 and wavelengths in nm.
_R0_FACTOR = 0.02108

#: Refractive index of the medium when the caller does not say. 1.4 is the
#: conventional value for a dye on a protein surface -- between water (1.33) and
#: protein interior (~1.6). It was hard-coded as `1.4**4` before PRD-113, which
#: is why nothing could ask what R0 would be in a different solvent.
DEFAULT_REFRACTIVE_INDEX = 1.4


def _r0_from_overlap(
    overlap: float, quantum_yield: float, k2: float, refractive_index: float
) -> float:
    """R0 in nm from the overlap integral, quantum yield, kappa^2 and n."""
    return float(_R0_FACTOR * np.power(
        k2 * quantum_yield / refractive_index ** 4 * overlap, 1.0 / 6.0))


def forster_radius(
    donor: "Dye",
    acceptor: "Dye",
    k2: float = 2.0 / 3.0,
    refractive_index: float = DEFAULT_REFRACTIVE_INDEX,
) -> float:
    """R0 in nm for a pair of :class:`IMP.bff.dye.Dye`, in a medium.

    The species-level entry point: everything it needs is a property of the two
    dyes and the solvent, so nothing has to be threaded in from a call site.

    :param donor: needs a spectrum and a quantum yield.
    :param acceptor: needs a spectrum and an extinction coefficient.
    :param k2: orientation factor; the isotropic 2/3 by default.
    :param refractive_index: of the medium between the dyes.
    """
    if not donor.has_spectrum or not acceptor.has_spectrum:
        raise ValueError(
            f"both dyes need a spectrum to derive R0 "
            f"({donor.name}: {donor.has_spectrum}, "
            f"{acceptor.name}: {acceptor.has_spectrum})")
    if donor.quantum_yield is None:
        raise ValueError(f"{donor.name} has no quantum yield")
    if acceptor.extinction_coefficient is None:
        raise ValueError(f"{acceptor.name} has no extinction coefficient")
    overlap = spectral_overlap(donor, acceptor)
    return _r0_from_overlap(
        overlap, donor.quantum_yield, k2, refractive_index)


def spectral_overlap(donor: "Dye", acceptor: "Dye") -> float:
    """The overlap integral J, in M^-1 cm^-1 nm^4.

    Donor emission against acceptor extinction, weighted by lambda^4 and
    normalised by the donor's emission integral.
    """
    wavelengths = np.asarray(donor.spectrum.wavelength, dtype=float)
    donor_emission = np.asarray(donor.spectrum.emission, dtype=float)
    acceptor_excitation = np.asarray(acceptor.spectrum.excitation, dtype=float)
    if acceptor_excitation.size != wavelengths.size:
        raise ValueError("donor and acceptor spectra must share one wavelength grid")
    emission_integral = _trapezoid(donor_emission, x=wavelengths)
    if emission_integral == 0:
        return 0.0
    extinction = acceptor.extinction_coefficient * acceptor_excitation
    overlap = _trapezoid(
        donor_emission * extinction * np.power(wavelengths, 4), x=wavelengths)
    return float(overlap / emission_integral)


def forster_radius_from_spectra(
    donor: str,
    acceptor: str,
    k2: float,
    r0_dir: str | Path | None = None,
    refractive_index: float = DEFAULT_REFRACTIVE_INDEX,
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
    refractive_index : float, optional
        Of the medium between the dyes. Was hard-coded before PRD-113.

    Returns
    -------
    float
        Förster radius in nm.

    Notes
    -----
    The name-based route, kept because the rotamer code addresses dyes by
    string. :func:`forster_radius` is the same calculation over two
    :class:`IMP.bff.dye.Dye` objects and is what new code should use.
    """
    donor_type, donor_number, donor_name = normalize_dye_name(donor)
    _acceptor_type, acceptor_number, acceptor_name = normalize_dye_name(acceptor)

    dye_table = read_dye_table()
    donor_data = dye_table.get((donor_type, donor_number))
    acceptor_data = dye_table.get((_acceptor_type, acceptor_number))
    if donor_data is None or acceptor_data is None:
        raise ValueError(f"No R0 dye data found for {donor!r} / {acceptor!r}")

    donor_spectrum = read_spectrum(find_r0_file(f"{donor_name}.csv", r0_dir))
    acceptor_spectrum = read_spectrum(find_r0_file(f"{acceptor_name}.csv", r0_dir))

    wavelengths = donor_spectrum.wavelength
    donor_emission = donor_spectrum.emission
    acceptor_excitation = acceptor_spectrum.excitation

    if donor_emission.size != wavelengths.size or acceptor_excitation.size != wavelengths.size:
        raise ValueError("Donor and acceptor spectra must share the same wavelength grid")

    emission_integral = _trapezoid(donor_emission, x=wavelengths)
    if emission_integral == 0:
        return 0.0

    ext_coeff_max = acceptor_data["Ext_coeff"]
    ext_coeff_acceptor = ext_coeff_max * acceptor_excitation
    overlap = _trapezoid(donor_emission * ext_coeff_acceptor * np.power(wavelengths, 4), x=wavelengths)
    overlap /= emission_integral

    return _r0_from_overlap(overlap, donor_data["QD"], k2, refractive_index)


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
