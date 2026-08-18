"""Interaction terms: the processes that deactivate a dye, and their rate constants.

Shaped the way a force field's terms are: a functional form, an **arity**, and
parameters looked up by type. Each term answers one question -- *given the
states of the participants, what is the rate constant per state* -- and the
answers **add**, because the channels are parallel. That additivity is the whole
reason these are objects rather than numbers computed inside an observable: the
field solver already relies on summing a quenching map and a FRET map, and
burying either inside its own function is what left ``fret_rate_trace`` and
``fret_rate_map`` as two incompatible calls.

Terms are **representation-agnostic**. They consume
:class:`~IMP.bff.representation.States` -- positions, weights, orientations -- so
one implementation serves an accessible volume, a rotamer library, a
coarse-grained model and an MD trajectory alike. A term that needs orientations
says so, and a representation that cannot supply them (an AV) is told, rather
than silently averaged.

Arity is what distinguishes state-dependent from emergent:

* **1-body** -- :class:`RadiativeTerm`. Depends on the dye and the solvent only.
* **2-body** -- :class:`PETTerm` (a dye against quencher atoms),
  :class:`FRETTerm` (a dye against a dye).
* **N-body** -- homo-FRET and multi-chromophore transfer, where the rate is not
  a sum over pairs. Not implemented; the interface is shaped to take it.

These currently delegate to the kernels in :mod:`IMP.bff.quenching`, and the
tests assert they reproduce them exactly. The kernels move underneath the terms
in a later stage; wrapping first and moving second is what keeps the change
reviewable.
"""

from __future__ import annotations

import abc
from dataclasses import dataclass, field
from typing import Mapping, Optional, Sequence

import numpy as np

__all__ = [
    "InteractionTerm", "RadiativeTerm", "PETTerm", "FRETTerm", "total_rate",
]


class InteractionTerm(abc.ABC):
    """One deactivation channel.

    :cvar arity: how many participants the term is defined over. 1 for a
        process of the dye alone, 2 for a pair, more for a genuinely many-body
        process.
    """

    arity: int = 1
    #: True when the term needs transition dipoles, so a caller can tell whether
    #: a representation without them (an accessible volume) forces an isotropic
    #: assumption instead of resolving the orientation.
    needs_orientations: bool = False

    @abc.abstractmethod
    def rate_constants(self, *participants, **kwargs) -> np.ndarray:
        """Rate constants in 1/ns, one per state of the first participant."""

    def __repr__(self) -> str:
        return f"{type(self).__name__}(arity={self.arity})"


@dataclass(repr=False)
class RadiativeTerm(InteractionTerm):
    """1-body: the dye's own decay, ``1/tau0``.

    The floor every other channel adds to. Depends on the dye and its medium,
    not on where it is or what is near it -- which is exactly what "1-body"
    means here.
    """

    lifetime: float
    arity: int = field(default=1, init=False)

    def rate_constants(self, states, **_) -> np.ndarray:
        if self.lifetime <= 0.0:
            raise ValueError(f"lifetime must be > 0, not {self.lifetime}")
        return np.full(int(states.n_points), 1.0 / float(self.lifetime),
                       dtype=np.float64)


@dataclass(repr=False)
class PETTerm(InteractionTerm):
    """2-body: photoinduced electron transfer between a dye and quencher atoms.

    ``k(r) = sum_a kQ_a * exp(-(|r - r_a| - r_dye) / rC_a)``

    The rate is a **pair** property -- ``kQ`` depends on the redox potentials of
    both partners -- so the parameters come in keyed by ``(dye, comp_id)`` as
    :class:`~IMP.bff.label.PETParameters` rather than as a table of the
    quencher alone.

    :param parameters: ``{comp_id: PETParameters}`` for **one** dye.
    :param dye_radius: subtracted from the centre-to-centre distance, because
        the tabulated contact distances are measured from the dye *surface*.
    """

    parameters: Mapping[str, "PETParameters"]
    dye_radius: float = 3.5
    arity: int = field(default=2, init=False)

    def rate_constants(self, states, atoms, **_) -> np.ndarray:
        """:param atoms: the structured atom array ``IMP.bff.quenching`` uses."""
        from IMP.bff.quenching import maps

        table = {
            comp_id: {
                atom: (p.rate_constant,
                       p.attenuation_length if p.attenuation_length else 1.0)
                for atom in _quencher_atoms(comp_id)
            }
            for comp_id, p in self.parameters.items()
        }
        kQ, rC = maps.atomic_quenching_parameters(atoms, table)
        xyz = np.ascontiguousarray(np.asarray(atoms["coord"], dtype=np.float64))
        positions = np.ascontiguousarray(
            np.asarray(states.positions, dtype=np.float64))
        return _pet_rates(positions, xyz, kQ, rC, float(self.dye_radius))


@dataclass(repr=False)
class FRETTerm(InteractionTerm):
    """2-body: Förster transfer between two dyes.

    ``k(r) = (1/tau0) * (R0/r)^6 * kappa^2 / (2/3)``

    **R0 is derived, not supplied** -- from the two dyes' spectra, the donor's
    quantum yield, the medium's refractive index and kappa^2. Passing it in was
    how ``forster_radius=52.0`` came to be a default in a dozen signatures.

    :param donor, acceptor: :class:`~IMP.bff.dye.Dye` species.
    :param refractive_index: of the medium between them.
    :param kappa2: orientation factor. ``None`` resolves it from the
        participants' orientations when both have them, and falls back to the
        isotropic 2/3 otherwise -- reporting which, through
        :attr:`used_isotropic_kappa2`.
    """

    donor: "Dye"
    acceptor: "Dye"
    refractive_index: float = 1.4
    kappa2: Optional[float] = None
    arity: int = field(default=2, init=False)
    needs_orientations: bool = field(default=True, init=False)

    def __post_init__(self):
        self.used_isotropic_kappa2 = self.kappa2 is None

    @property
    def forster_radius(self) -> float:
        """R0 in Angstrom, derived from the pair and the medium."""
        from IMP.bff.dye.spectra import forster_radius
        k2 = 2.0 / 3.0 if self.kappa2 is None else float(self.kappa2)
        return 10.0 * forster_radius(
            self.donor, self.acceptor, k2, self.refractive_index)

    def rate_constants(self, donor_states, acceptor_states, r_min: float = 7.0,
                       **_) -> np.ndarray:
        from IMP.bff.quenching.fret_trace import fret_rate_trace

        tau0 = self.donor.lifetime
        if tau0 is None or tau0 <= 0.0:
            raise ValueError(
                f"{self.donor.name} has no lifetime, so a FRET rate cannot be "
                "expressed as 1/tau0 * (R0/r)^6")
        return fret_rate_trace(
            np.ascontiguousarray(np.asarray(donor_states.positions, dtype=np.float64)),
            np.ascontiguousarray(np.asarray(acceptor_states.positions, dtype=np.float64)),
            R0=self.forster_radius, tau0=float(tau0), r_min=r_min,
            kappa2=self.kappa2,
        )


def total_rate(terms: Sequence[InteractionTerm], *participants, **kwargs) -> np.ndarray:
    """Sum the channels.

    Parallel deactivation channels add, which is the property that makes a
    quenching map and a FRET map summable on the same grid. Doing it here rather
    than inside each observable is what keeps that true when a channel is added.
    """
    if not terms:
        raise ValueError("no interaction terms to sum")
    out = None
    for term in terms:
        r = np.asarray(term.rate_constants(*participants[:term.arity], **kwargs),
                       dtype=np.float64)
        out = r.copy() if out is None else out + r
    return out


def _quencher_atoms(comp_id: str):
    from IMP.bff.quenching.pet import QUENCHER_ATOMS
    return QUENCHER_ATOMS.get(comp_id, ())


def _pet_rates(positions, atoms_xyz, kQ, rC, dye_radius) -> np.ndarray:
    """Per-state PET rate, the same law ``quenching_rate_map`` stamps on a grid.

    Written over *states* rather than voxels so it serves a rotamer library and
    an MD trajectory as well as an accessible volume.
    """
    active = (kQ > 0.0) & (rC > 0.0)
    if not np.any(active):
        return np.zeros(positions.shape[0], dtype=np.float64)
    xyz, k, r = atoms_xyz[active], kQ[active], rC[active]
    out = np.zeros(positions.shape[0], dtype=np.float64)
    for i in range(positions.shape[0]):
        d = np.sqrt(np.sum((xyz - positions[i]) ** 2, axis=1)) - dye_radius
        out[i] = float(np.sum(k * np.exp(-d / r)))
    return out
