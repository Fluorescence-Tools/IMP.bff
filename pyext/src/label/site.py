"""Where a dye sits on a structure, in flrCIF's vocabulary.

Before PRD-113 a labelling site was an untyped ``source_info`` dict threaded
through the AV builder, the quenching model and every benchmark -- and it mixed
two different things: *where the dye is attached* (chain, residue, atom) and
*how its accessible volume is computed* (``linker_length``, ``radius1..3``,
``allowed_sphere_radius``, ``simulation_grid_resolution``). Those are a system
property and a representation parameter respectively, and conflating them is why
`AV` ended up doubling as "the dye".

:class:`Label` is the first half only. The second half belongs to whichever
representation is used, and a rotamer library has none of those fields.

Field names follow the FLR dictionaries: ``_flr_poly_probe_position`` for the
position and ``_flr_sample_probe_details`` for what the probe is doing there.
Checked against ``../mmfdb/src/mmfdb/data/*.dic``, not recalled.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

__all__ = ["Label", "FLUOROPHORE_TYPES", "FLRCIF_ITEMS"]

#: ``_flr_sample_probe_details.fluorophore_type`` enumeration, verbatim.
FLUOROPHORE_TYPES = ("donor", "acceptor", "unspecified")

#: flrCIF item per field. ``asym_id``/``seq_id``/``comp_id``/``atom_id`` are the
#: dictionary's names for what the fps dialect calls ``chain_identifier``,
#: ``residue_seq_number``, the residue name and ``atom_name``.
FLRCIF_ITEMS = {
    "asym_id": "_flr_poly_probe_position.asym_id",
    "seq_id": "_flr_poly_probe_position.seq_id",
    "comp_id": "_flr_poly_probe_position.comp_id",
    "atom_id": "_flr_poly_probe_position.atom_id",
    "auth_name": "_flr_poly_probe_position.auth_name",
    "mutation_flag": "_flr_poly_probe_position.mutation_flag",
    "modification_flag": "_flr_poly_probe_position.modification_flag",
    "fluorophore_type": "_flr_sample_probe_details.fluorophore_type",
    "description": "_flr_sample_probe_details.description",
    "dye": "_flr_sample_probe_details.probe_id",
}


@dataclass(frozen=True)
class Label:
    """A dye attached at one position on a structure.

    :param asym_id: chain identifier.
    :param seq_id: residue sequence number.
    :param atom_id: attachment atom, e.g. ``"CB"``.
    :param comp_id: residue name, when known.
    :param dye: the :class:`IMP.bff.dye.Dye` species, when known. A label whose
        dye is unknown is still a usable site -- it just cannot derive R0 or a
        rotational correlation time.
    :param fluorophore_type: ``donor``, ``acceptor`` or ``unspecified``.
    :param auth_name: author-assigned name for the position.
    :param mutation_flag: the residue was mutated to carry the label.
    :param modification_flag: the residue was chemically modified.
    :param description: free text.
    """

    asym_id: str
    seq_id: int
    atom_id: str = "CB"
    comp_id: Optional[str] = None
    dye: Optional["Dye"] = None
    fluorophore_type: str = "unspecified"
    auth_name: Optional[str] = None
    mutation_flag: bool = False
    modification_flag: bool = False
    description: Optional[str] = None

    def __post_init__(self):
        if self.fluorophore_type not in FLUOROPHORE_TYPES:
            raise ValueError(
                f"fluorophore_type must be one of {FLUOROPHORE_TYPES}, "
                f"not {self.fluorophore_type!r}")

    @property
    def key(self) -> tuple:
        """``(asym_id, seq_id, atom_id)`` -- what identifies the position."""
        return (self.asym_id, int(self.seq_id), self.atom_id)

    @classmethod
    def from_source_info(cls, source_info: dict, dye=None) -> "Label":
        """Build a Label from the fps-dialect dict the AV builder still takes.

        The bridge while the representation layer is migrated: it reads only the
        *position* fields and ignores the AV parameters in the same dict, which
        is the split this class exists to make.
        """
        return cls(
            asym_id=str(source_info.get("chain_identifier", "")),
            seq_id=int(source_info.get("residue_seq_number", 0)),
            atom_id=str(source_info.get("atom_name", "CB")),
            dye=dye,
        )

    def to_source_info(self) -> dict:
        """The position half of an fps-dialect dict, for builders that want one."""
        return {
            "chain_identifier": self.asym_id,
            "residue_seq_number": int(self.seq_id),
            "atom_name": self.atom_id,
        }

    def __repr__(self) -> str:
        dye = self.dye.name if self.dye is not None else "?"
        return (f"Label({self.asym_id}{self.seq_id}.{self.atom_id}, "
                f"dye={dye}, {self.fluorophore_type})")
