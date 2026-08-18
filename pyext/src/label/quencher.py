"""What quenches a dye, and with what -- two different things.

A :class:`Quencher` is an *identity*: which residue type quenches and which of
its atoms is the redox-active moiety. That is a property of the protein.

:class:`PETParameters` is a **pair property**: how fast a *particular dye* is
quenched by a *particular residue type*. Photoinduced electron transfer depends
on the redox potentials of both partners, so ``kQ`` is not a property of the
tryptophan -- a rhodamine, an oxazine and a cyanine see the same tryptophan
differently, and an oxazine can be quenched where a xanthene is not.

The bundled table says so itself: ``PET_QUENCHING_REFERENCE`` is documented as
*"reference PET parameters for a xanthene dye (Alexa488-like)"* and is then
applied to every dye in the package. :func:`reference_pet_parameters` keeps that
transfer but makes it visible -- each returned parameter records the dye it was
measured for, so a set applied to a different dye can be recognised as
transferred rather than measured.

**No dictionary in the stack has a quenching item.** Checked across all ten
``.dic`` files in ``../mmfdb/src/mmfdb/data``: nothing matches "quench". These
names are bff-native; the residue and atom identifiers still follow flrCIF's
``comp_id`` / ``atom_id`` so a quencher is located the same way a label is.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

__all__ = ["Quencher", "PETParameters", "reference_quenchers",
           "reference_pet_parameters", "REFERENCE_DYE", "FLRCIF_ITEMS"]

#: The dye the bundled PET table was measured for -- a xanthene, Alexa488-like.
#: Peulen et al., J. Phys. Chem. B 2017, 121, 8211.
REFERENCE_DYE = "AlexaFluor488"

#: flrCIF items where they exist. The identifiers do; the photophysics does not.
FLRCIF_ITEMS = {
    "comp_id": "_flr_poly_probe_position.comp_id",
    "asym_id": "_flr_poly_probe_position.asym_id",
    "seq_id": "_flr_poly_probe_position.seq_id",
    "atom_ids": "_flr_poly_probe_position.atom_id",
    # bff-native: no quenching category anywhere in the stack's dictionaries
    "rate_constant": None,
    "attenuation_length": None,
    "contact_distance": None,
    "dye": None,
}


@dataclass(frozen=True)
class Quencher:
    """A quenching moiety: which residue, which atoms. **No rate.**

    :param comp_id: residue name, e.g. ``"TRP"``.
    :param atom_ids: the redox-active atoms. Deliberately not CB -- electron
        transfer happens at the indole ring, the phenol, the thioether or the
        thiol, and stamping a rate on CB puts it up to 4 A from the chemistry.
    :param asym_id, seq_id: set when this is a *specific* residue in a structure
        rather than a type.

    The rate lives in :class:`PETParameters` because it takes two partners to
    define one.
    """

    comp_id: str
    atom_ids: Tuple[str, ...] = ()
    asym_id: Optional[str] = None
    seq_id: Optional[int] = None

    def __post_init__(self):
        object.__setattr__(self, "atom_ids", tuple(self.atom_ids))

    @property
    def is_typed(self) -> bool:
        """True when this describes a residue *type* rather than one residue."""
        return self.asym_id is None and self.seq_id is None

    def at(self, asym_id: str, seq_id: int) -> "Quencher":
        """The same moiety, located at one residue of a structure."""
        return Quencher(comp_id=self.comp_id, atom_ids=self.atom_ids,
                        asym_id=asym_id, seq_id=int(seq_id))

    def __repr__(self) -> str:
        where = "" if self.is_typed else f" @{self.asym_id}{self.seq_id}"
        return f"Quencher({self.comp_id}{where}, atoms={list(self.atom_ids)})"


@dataclass(frozen=True)
class PETParameters:
    """PET between one dye and one quencher type -- a **pair** property.

    :param dye: chromophore name the parameters apply to.
    :param comp_id: quencher residue type.
    :param rate_constant: ``kQ`` in 1/ns at contact.
    :param contact_distance: from the *dye surface* to the quenching centre, in
        Angstrom -- roughly van der Waals contact plus the offset from the
        moiety centroid to its outer atoms. Surface-relative on purpose, so the
        geometry transfers between dyes of different size even when ``kQ`` does
        not.
    :param attenuation_length: ``rC``, the exponential decay length of the
        through-space rate. ``None`` for a hard contact-sphere model.
    :param measured_for: the dye the values were actually measured with. When
        this differs from ``dye`` the parameters were **transferred**, which is
        an assumption about redox chemistry rather than a measurement.
    """

    dye: str
    comp_id: str
    rate_constant: float
    contact_distance: float
    attenuation_length: Optional[float] = None
    measured_for: Optional[str] = None

    def __post_init__(self):
        if self.rate_constant < 0.0:
            raise ValueError(f"rate_constant must be >= 0, not {self.rate_constant}")
        if self.measured_for is None:
            object.__setattr__(self, "measured_for", self.dye)

    @property
    def is_transferred(self) -> bool:
        """True when these values were measured with a *different* dye."""
        return self.measured_for != self.dye

    def scaled(self, rate_scale: float) -> "PETParameters":
        """The same pair with ``kQ`` multiplied -- what a calibration turns."""
        return PETParameters(
            dye=self.dye, comp_id=self.comp_id,
            rate_constant=self.rate_constant * float(rate_scale),
            contact_distance=self.contact_distance,
            attenuation_length=self.attenuation_length,
            measured_for=self.measured_for,
        )

    def __repr__(self) -> str:
        tag = f", transferred from {self.measured_for}" if self.is_transferred else ""
        return (f"PETParameters({self.dye} x {self.comp_id}, "
                f"kQ={self.rate_constant:g}{tag})")


def reference_quenchers() -> Dict[str, Quencher]:
    """The redox-active moieties, read from the one table that defines them.

    Identity only -- see :func:`reference_pet_parameters` for the rates, which
    depend on the dye.
    """
    from IMP.bff.quenching.pet import PET_QUENCHING_REFERENCE, QUENCHER_ATOMS

    return {
        comp_id: Quencher(comp_id=comp_id,
                          atom_ids=tuple(QUENCHER_ATOMS.get(comp_id, ())))
        for comp_id in PET_QUENCHING_REFERENCE
    }


def reference_pet_parameters(
    dye: str = REFERENCE_DYE,
    rate_scale: float = 1.0,
    attenuation_length: Optional[float] = None,
) -> Dict[str, PETParameters]:
    """The published PET chemistry for one dye against each quencher type.

    Reads ``IMP.bff.quenching.pet.PET_QUENCHING_REFERENCE`` rather than
    restating it. That table was measured for a xanthene dye
    (:data:`REFERENCE_DYE`); asking for any other dye returns the same numbers
    with ``measured_for`` set to the reference, so
    :attr:`PETParameters.is_transferred` reports that an assumption was made.

    :param dye: chromophore the parameters are wanted for.
    :param rate_scale: multiplies every ``kQ``. The published values are
        *starting points to be calibrated against measured lifetimes*, and this
        is the knob a calibration turns -- PRD-111 found it recoverable to a few
        percent from six sites jointly and not at all from one.
    :param attenuation_length: ``rC`` for the exponential through-space form.
        ``None`` keeps the hard contact-sphere model.
    """
    from IMP.bff.quenching.pet import PET_QUENCHING_REFERENCE

    return {
        comp_id: PETParameters(
            dye=str(dye),
            comp_id=comp_id,
            rate_constant=float(entry["kQ"]) * float(rate_scale),
            contact_distance=float(entry["contact_distance"]),
            attenuation_length=attenuation_length,
            measured_for=REFERENCE_DYE,
        )
        for comp_id, entry in PET_QUENCHING_REFERENCE.items()
    }
