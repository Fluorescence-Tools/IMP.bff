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

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence, Tuple

import numpy as np

__all__ = ["Spectrum", "Dye", "FLRCIF_ITEMS"]

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
FLRCIF_ITEMS = {
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
    :data:`FLRCIF_ITEMS`, which also records the three that it does not cover.
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
