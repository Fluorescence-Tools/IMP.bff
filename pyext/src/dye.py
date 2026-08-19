"""``IMP.bff.dye`` -- the dye as a species, independent of how it is modelled.

A dye's spectra, quantum yield, extinction coefficient, transition dipole and
formal charges are what it *is*; an accessible volume, a rotamer library, a
coarse-grained conformer set and an MD trajectory are ways of *representing* it.
Only the former lives here -- see ``IMP.bff.representation`` for the latter.

The practical consequence is that the Förster radius is derived rather than
supplied: :func:`forster_radius` takes two dyes, kappa^2 and the medium's
refractive index, so no call site needs to carry ``forster_radius=52.0``.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple
import ihm.format
import re

import numpy as np

__all__ = [
    'DEFAULT_REFRACTIVE_INDEX',
    'DYE_FLRCIF_ITEMS',
    'DYE_LIBRARY_CIF',
    'Dye',
    'Spectrum',
    'available_dyes',
    'find_dye',
    'forster_radius',
    'forster_radius_from_spectra',
    'normalize_dye_name',
    'read_dye_library',
    'spectral_overlap',
]

# --------------------------------------------------------------------------
# species
# --------------------------------------------------------------------------
"""The dye as a species: what a dye *is*, independent of how it is modelled.

Until PRD-113 a dye had two unrelated descriptions in this package. The
*spectral* one lived in ``fret/forster.py`` -- names, extinction coefficients,
quantum yields and emission/excitation curves read from the bundled tables. The
*molecular* one lived in ``cgdye/io/template_cif.py`` -- atoms, bonds, the two
atoms defining the transition dipole, formal charges. Nothing connected them, so
"what is AlexaFluor 488" had two answers, and a Förster radius could not be
derived from a dye pair without going through a name string twice.

A :class:`Dye` is that one answer. It is deliberately **model-independent**: an
accessible volume, a rotamer library, a coarse-grained conformer set and an MD
trajectory are all ways of *representing* the same dye, and none of their
parameters belong here. ``linker_length`` and ``allowed_sphere_radius`` are
properties of the AV representation; a rotamer library has neither.
"""

#: flrCIF item for each :class:`Dye` field, or ``None`` where the dictionary has
#: none. Same convention as ``fret/fps_schema.py``, and checked against
#: ``python-ihm``'s FLR model (``ihm.flr``) rather than guessed.
#:
#: **flrCIF's word for a dye is "probe", and for the fluorescent moiety
#: "chromophore".** A *probe* is the labelling reagent (``reactive_probe_name``,
#: e.g. a maleimide) and its *chromophore* is what fluoresces
#: (``chromophore_name``). ``Dye`` here is the chromophore plus the photophysics.
#:
#: **Note what is missing**: **no dictionary in the stack has an item for
#: quantum yield, extinction coefficient or a spectrum.** Checked across all ten
#: ``.dic`` files in ``../mmfdb/src/mmfdb/data`` (mmCIF std, PDBx v50 and
#: v5_next, DDL, MA, IHM, IHM-FLR, mmfdb FLR and workflow extensions): the only
#: matches are ``_em_detector.detective_quantum_efficiency`` and the NMR
#: spectral categories. The closest FLR item,
#: ``_flr_fret_calibration_parameters.phi_acceptor``, is an analysis calibration
#: value rather than a property of the species. Those three fields are therefore
#: bff-native and marked ``None`` -- deliberately, not by omission.
DYE_FLRCIF_ITEMS = {
    "name": "_flr_probe_list.chromophore_name",
    "reactive_probe_name": "_flr_probe_list.reactive_probe_name",
    "probe_origin": "_flr_probe_list.probe_origin",
    "probe_link_type": "_flr_probe_list.probe_link_type",
    "chromophore_center_atom": "_flr_probe_descriptor.chromophore_center_atom",
    "lifetime": "_flr_reference_measurement_lifetime.lifetime",
    # no flrCIF item exists for these
    "spectrum": None,
    "extinction_coefficient": None,
    "quantum_yield": None,
    "dipole_atoms": None,
    "radius": None,
    "hydrodynamic_radius": None,
    "positive_atoms": None,
    "negative_atoms": None,
}


@dataclass(frozen=True)
class Spectrum:
    """Excitation and emission on a shared wavelength grid.

    :param wavelength: nm.
    :param excitation: normalised to a maximum of 1.
    :param emission: normalised to a maximum of 1.

    Both curves share one grid because the overlap integral needs them
    point-for-point; :func:`IMP.bff.dye.spectra.read_spectrum` enforces it.
    """

    wavelength: np.ndarray
    excitation: np.ndarray
    emission: np.ndarray

    def __post_init__(self):
        n = len(self.wavelength)
        if len(self.excitation) != n or len(self.emission) != n:
            raise ValueError(
                f"excitation ({len(self.excitation)}) and emission "
                f"({len(self.emission)}) must share the wavelength grid ({n})")


@dataclass(frozen=True)
class Dye:
    """A dye species.

    :param name: the compact name, e.g. ``AlexaFluor488``.
    :param spectrum: excitation/emission curves, or ``None`` if unknown.
    :param extinction_coefficient: molar extinction at the absorption maximum,
        M^-1 cm^-1. With the acceptor's excitation curve this gives the
        wavelength-resolved extinction the overlap integral needs.
    :param quantum_yield: fluorescence quantum yield of the free dye.
    :param lifetime: unquenched fluorescence lifetime in ns, if known.
    :param dipole_atoms: the two atom names whose separation defines the
        transition dipole. Needed for kappa^2; ``None`` for an isotropic model.
    :param chromophore_center_atom: the atom taken as the chromophore centre.
        flrCIF ``_flr_probe_descriptor.chromophore_center_atom``.
    :param reactive_probe_name: the labelling reagent, e.g. the maleimide form.
        flrCIF distinguishes the reactive probe from its chromophore; ``name``
        is the chromophore.
    :param radius: steric radius in Angstrom, for representations that need one.
    :param hydrodynamic_radius: Angstrom, for the rotational correlation time
        and the translational diffusion coefficient.

    Nothing here says how the dye's positions are enumerated. That is the
    representation's business.

    Field names follow flrCIF where the dictionary has an item; see
    :data:`DYE_FLRCIF_ITEMS`, which also records the three that it does not cover.
    """

    name: str
    spectrum: Optional[Spectrum] = None
    extinction_coefficient: Optional[float] = None
    quantum_yield: Optional[float] = None
    lifetime: Optional[float] = None
    dipole_atoms: Optional[Tuple[str, str]] = None
    chromophore_center_atom: Optional[str] = None
    radius: Optional[float] = None
    hydrodynamic_radius: Optional[float] = None
    reactive_probe_name: Optional[str] = None
    probe_origin: Optional[str] = None
    probe_link_type: Optional[str] = None
    positive_atoms: Sequence[str] = field(default_factory=tuple)
    negative_atoms: Sequence[str] = field(default_factory=tuple)

    @property
    def has_spectrum(self) -> bool:
        return self.spectrum is not None

    def __repr__(self) -> str:
        bits = [f"name={self.name!r}"]
        if self.quantum_yield is not None:
            bits.append(f"QY={self.quantum_yield:.3g}")
        if self.extinction_coefficient is not None:
            bits.append(f"eps={self.extinction_coefficient:.4g}")
        bits.append(f"spectrum={'yes' if self.has_spectrum else 'no'}")
        return f"Dye({', '.join(bits)})"


# --------------------------------------------------------------------------
# cif
# --------------------------------------------------------------------------
"""Read the bundled dye library from mmCIF.

**CIF is the data format for this package** (PRD-113). The library was CSV --
one table of extinction coefficients and quantum yields plus one file of
excitation/emission curves per dye -- which carried no category names, no units
and no way to say where a number came from.

Two categories, both bff-native:

``_bff_dye``
    ``chromophore_name``, ``probe_type``, ``chromophore_number``,
    ``extinction_coefficient``, ``quantum_yield``
``_bff_dye_spectrum``
    ``chromophore_name``, ``wavelength``, ``excitation``, ``emission``

They are bff-native because **no dictionary in the stack defines an item for a
quantum yield, an extinction coefficient or a spectrum** -- checked across all
ten ``.dic`` files in ``../mmfdb/src/mmfdb/data`` (mmCIF std, PDBx v50 and
v5_next, DDL, MA, IHM, IHM-FLR, mmfdb FLR and workflow extensions); the only
matches are ``_em_detector.detective_quantum_efficiency`` and the NMR spectral
categories. They should be proposed for ``mmfdb_flr_ext.dic``.

Identifiers follow flrCIF: ``chromophore_name`` is
``_flr_probe_list.chromophore_name``.
"""

#: Name of the bundled library inside the rotamer-library data directory.
DYE_LIBRARY_CIF = "dye_library.cif"


class _BaseHandler:
    """``ihm.format.CifReader`` reads the keywords from ``__call__``'s signature.

    So the parameter names below *are* the CIF item names -- a ``**kwargs``
    handler is handed nothing at all, silently.
    """

    not_in_file = object()
    omitted = object()
    unknown = object()

    def __init__(self):
        self.rows = []

    def end_save_frame(self):
        pass

    def _clean(self, **kwargs):
        self.rows.append({
            k: (None if v in (self.not_in_file, self.omitted, self.unknown) else v)
            for k, v in kwargs.items()})


class _DyeHandler(_BaseHandler):
    def __call__(self, chromophore_name, probe_type, chromophore_number,
                 extinction_coefficient, quantum_yield):
        self._clean(chromophore_name=chromophore_name, probe_type=probe_type,
                    chromophore_number=chromophore_number,
                    extinction_coefficient=extinction_coefficient,
                    quantum_yield=quantum_yield)


class _SpectrumHandler(_BaseHandler):
    def __call__(self, chromophore_name, wavelength, excitation, emission):
        self._clean(chromophore_name=chromophore_name, wavelength=wavelength,
                    excitation=excitation, emission=emission)


def _float(value) -> Optional[float]:
    if value is None or value in ("", ".", "?"):
        return None
    return float(value)


#: Parsed libraries, keyed by resolved path and modification time.
_LIBRARY_CACHE: Dict[tuple, Dict[str, "Dye"]] = {}


def read_dye_library(path: str | Path | None = None) -> Dict[str, Dye]:
    """Every dye in the bundled library, keyed by chromophore name.

    :param path: the CIF file; defaults to the bundled one.
    :returns: ``{chromophore_name: Dye}``. A dye with a table entry but no
        curves is returned without a spectrum -- it is still a usable species,
        it just cannot derive R0.

    **Cached on (path, mtime, size).** The bundled library is a shipped
    read-only file, and a spectrum lookup does not change it -- but the parse is
    not cheap: profiling one FRETpredict comparison found this called 21 times
    for the same file, 550 000 calls to the row cleaner and 1.65 million to the
    float converter, for 2.2 s of a 13 s run.

    The key includes mtime and size so editing the file during a session is
    picked up; only an edit that changes neither would be missed, which is not a
    thing that happens to a data file. Callers that mutate the returned ``Dye``
    objects would now be mutating the cached ones -- nothing does, and nothing
    should: a dye is a species, not a scratch pad.
    """
    if path is None:
        import IMP.bff
        path = Path(IMP.bff.get_data_path("rotamer_library")) / "R0" / DYE_LIBRARY_CIF
    resolved = Path(path)
    key = None
    if resolved.exists():
        stat = resolved.stat()
        key = (str(resolved.resolve()), stat.st_mtime_ns, stat.st_size)
        hit = _LIBRARY_CACHE.get(key)
        if hit is not None:
            # A fresh mapping each time, so a caller that adds or drops an entry
            # does not edit the library for everyone. The `Dye` values are
            # shared, and that is safe because `Dye` is frozen -- reaching round
            # that with `object.__setattr__` corrupts the species for every
            # later reader, which is exactly what a test did on the day this
            # cache landed.
            return dict(hit)
    path = resolved

    dyes = _DyeHandler()
    spectra = _SpectrumHandler()
    with path.open() as handle:
        reader = ihm.format.CifReader(handle, {"_bff_dye": dyes,
                                               "_bff_dye_spectrum": spectra})
        reader.read_file()

    curves: Dict[str, list] = {}
    for row in spectra.rows:
        name = row.get("chromophore_name")
        if name is None:
            continue
        curves.setdefault(name, []).append(
            (_float(row.get("wavelength")),
             _float(row.get("excitation")),
             _float(row.get("emission"))))

    out: Dict[str, Dye] = {}
    for row in dyes.rows:
        name = row.get("chromophore_name")
        if name is None:
            continue
        points = sorted(curves.get(name, []), key=lambda r: r[0])
        spectrum = None
        if points:
            spectrum = Spectrum(
                wavelength=np.array([p[0] for p in points], dtype=float),
                excitation=np.array([p[1] for p in points], dtype=float),
                emission=np.array([p[2] for p in points], dtype=float),
            )
        out[name] = Dye(
            name=name,
            spectrum=spectrum,
            extinction_coefficient=_float(row.get("extinction_coefficient")),
            quantum_yield=_float(row.get("quantum_yield")),
        )
    if key is not None:
        _LIBRARY_CACHE[key] = out
        return dict(out)
    return out


# --------------------------------------------------------------------------
# spectra
# --------------------------------------------------------------------------
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

_DYE_NAME_RE = re.compile(r"^(?P<type>.+?)\s+(?P<number>[A-Za-z0-9]+)$")


# numpy 2 removed np.trapz in favour of np.trapezoid, and this module is the
# Foerster-radius calculation -- so under the numpy the stack actually runs
# (2.4) forster_radius_from_spectra raised AttributeError and R0 could not be computed at all.
try:  # numpy >= 2
    _trapezoid = np.trapezoid
except AttributeError:  # pragma: no cover - numpy < 2
    _trapezoid = np.trapz

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


#: flrCIF items for the derivation's inputs and output. ``index_of_refraction``
#: and ``kappa_squared`` are **not** in the upstream IHM-FLR dictionary -- they
#: are added by ``mmfdb_flr_ext.dic`` (``../mmfdb/src/mmfdb/data``), which is
#: what makes a stored R0 reproducible rather than a bare number.
FORSTER_RADIUS_FLRCIF_ITEMS = {
    "forster_radius": "_flr_fret_forster_radius.forster_radius",
    "k2": "_flr_fret_forster_radius.kappa_squared",
    "refractive_index": "_flr_fret_forster_radius.index_of_refraction",
    "donor": "_flr_fret_forster_radius.donor_probe_id",
    "acceptor": "_flr_fret_forster_radius.acceptor_probe_id",
    # no dictionary in the stack has an item for these
    "spectral_overlap": None,
}

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


    # (was: from .cif import ...) -- merged into this module
    library = read_dye_library()
    donor_dye, acceptor_dye = library.get(donor_name), library.get(acceptor_name)
    if donor_dye is None or acceptor_dye is None or not (
            donor_dye.has_spectrum and acceptor_dye.has_spectrum):
        raise ValueError(f"No spectra for {donor!r} / {acceptor!r}")

    wavelengths = donor_dye.spectrum.wavelength
    donor_emission = donor_dye.spectrum.emission
    acceptor_excitation = acceptor_dye.spectrum.excitation

    if donor_emission.size != wavelengths.size or acceptor_excitation.size != wavelengths.size:
        raise ValueError("Donor and acceptor spectra must share the same wavelength grid")

    emission_integral = _trapezoid(donor_emission, x=wavelengths)
    if emission_integral == 0:
        return 0.0

    ext_coeff_max = acceptor_dye.extinction_coefficient
    ext_coeff_acceptor = ext_coeff_max * acceptor_excitation
    overlap = _trapezoid(donor_emission * ext_coeff_acceptor * np.power(wavelengths, 4), x=wavelengths)
    overlap /= emission_integral

    return _r0_from_overlap(
        overlap, donor_dye.quantum_yield, k2, refractive_index)


# --------------------------------------------------------------------------
# library
# --------------------------------------------------------------------------
"""Look a dye species up by name from the bundled tables.

One place that answers "what is AlexaFluor 488": the extinction coefficient and
quantum yield and the excitation/emission curves from the bundled
``dye_library.cif``, and -- when a template CIF is supplied -- the transition
dipole atoms, the chromophore centre and the formal charges.

Before PRD-113 those two halves were read by unrelated code
(``fret/forster.py`` and ``cgdye/io/template_cif.py``) and never met, and the
spectral half was CSV. **CIF is the data format now** -- see
:mod:`IMP.bff.dye` for the categories and why they are bff-native.
"""

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

    from IMP.bff.io.cif import read_dye_template_cif
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
