"""One strip mechanism: fps ``strip_mask`` selections over structures (PRD-106).

The fps.json ``strip_mask`` field names atoms that must not be obstacles
(``chain A and resid 36 and not name N+CA+C+O+CB``). This module is the one
place that grammar is parsed and applied. It offers the same selection over
three representations so every consumer strips the same way:

* PDB text lines (:func:`strip_pdb_lines`) — the AV build (``fret.av``);
* an ``IMP.atom.Hierarchy`` (:func:`select_atoms`, :func:`strip_hierarchy`) —
  cgdye labelling, which removes the labelled residue's side chain before it
  attaches an explicit dye;
* an obstacle array with parsed records (:func:`strip_obstacles`).

Only the mask grammar lives here; the *defaults* belong to the consumers,
because they legitimately differ: the AV convention keeps ``N CA C O`` plus
the attachment atom (usually CB), cgdye keeps ``N CA C O OXT`` and strips CB
(the explicit linker is built off CA and replaces the whole side chain).
``fret`` never imports ``cgdye``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterable, Optional, Sequence, Tuple

import numpy as np

#: The atoms of a peptide backbone (the AV strip keeps these plus the
#: attachment atom; see :func:`default_strip_mask`).
BACKBONE_ATOM_NAMES = ("N", "CA", "C", "O")


@dataclass(frozen=True)
class StripSelection:
    """A parsed ``strip_mask``: which (chain, resseq, atom name) it selects."""

    chain: Optional[str]          #: chain id, ``None`` = any chain
    resids: Tuple[int, ...]       #: residue numbers, empty = any residue
    names: Optional[frozenset]    #: atom-name set, ``None`` = no name term
    negate: bool                  #: the name term was ``not name ...``

    def matches(self, chain: str, resseq: Optional[int], name: str) -> bool:
        """True when an atom (chain, resseq, name) is selected by the mask."""
        if resseq is None or not name:
            return False
        if self.chain is not None and (chain or "") != self.chain:
            return False
        if self.resids and resseq not in self.resids:
            return False
        if self.names is not None and (name.upper() not in self.names) != self.negate:
            return False
        return True


def parse_strip_mask(mask: str) -> StripSelection:
    """Parse an fps ``strip_mask`` into a :class:`StripSelection`.

    The accepted grammar is the dialect fps documents actually carry: an
    ``and``-chain of ``chain <id>``, ``resid <n>`` (``resi`` accepted), and
    one name term, ``name A+B+...`` or ``not name A+B+...``. ``+`` is the
    list separator (PyMOL's own; space-separated lists are a parse error in
    PyMOL and here). Anything else -- ``or``, parentheses, other keywords --
    is refused loudly: a mask this function cannot read must not be
    silently ignored or approximated.
    """
    chain = None
    resids: tuple[int, ...] = ()
    names = None
    negate = False
    seen_name_term = False
    for term in re.split(r"\s+and\s+", mask.strip(), flags=re.IGNORECASE):
        term = term.strip()
        m = re.fullmatch(r"chain\s+(\S+)", term, flags=re.IGNORECASE)
        if m:
            if chain is not None:
                raise ValueError(f"strip_mask: repeated 'chain' term in {mask!r}")
            chain = m.group(1)
            continue
        m = re.fullmatch(r"resi(?:d)?\s+(-?\d+)", term, flags=re.IGNORECASE)
        if m:
            if resids:
                raise ValueError(f"strip_mask: repeated 'resid' term in {mask!r}")
            resids = (int(m.group(1)),)
            continue
        m = re.fullmatch(r"(not\s+)?name\s+(\S+(?:\+\S+)*)", term, flags=re.IGNORECASE)
        if m:
            if seen_name_term:
                raise ValueError(f"strip_mask: repeated name term in {mask!r}")
            seen_name_term = True
            negate = bool(m.group(1))
            parts = [p for p in m.group(2).split("+") if p]
            if not parts:
                raise ValueError(f"strip_mask: empty name list in {mask!r}")
            names = frozenset(p.upper() for p in parts)
            continue
        raise ValueError(
            f"strip_mask {mask!r} is outside the fps dialect "
            "(chain <id> and resid <n> and [not] name A+B+...); "
            "evaluate the selection externally and hand the stripped "
            "structure to the consumer"
        )
    return StripSelection(chain, resids, names, negate)


def site_strip_mask(chain: str, resseq: int, keep_atom_names: Iterable[str]) -> str:
    """Mask that strips residue ``(chain, resseq)`` except ``keep_atom_names``.

    ``chain`` empty matches any chain. Spelled in the fps dialect so a
    declared mask and a consumer default round-trip through the same parser.
    """
    keep = "+".join(dict.fromkeys(str(n).upper() for n in keep_atom_names))
    parts = []
    if chain:
        parts.append(f"chain {chain}")
    parts.append(f"resid {int(resseq)}")
    if keep:
        parts.append(f"not name {keep}")
    return " and ".join(parts)


def default_strip_mask(chain: str, resseq: int, atom_name: str) -> str:
    """The AV default strip for an attachment site, as an fps ``strip_mask``.

    The attachment residue's side chain minus the attachment atom; the
    backbone stays -- the FPS convention every fps.json
    ``allowed_sphere_radius`` is calibrated against. E.g.
    ``chain A and resid 36 and not name N+CA+C+O+CB``.
    """
    return site_strip_mask(chain, resseq, BACKBONE_ATOM_NAMES + (atom_name,))


# ---------------------------------------------------------------------------
# PDB text
# ---------------------------------------------------------------------------

def _pdb_atom_fields(line: str):
    """(chain, resseq, name) of an ATOM/HETATM line, or None."""
    if not line.startswith(("ATOM  ", "HETATM")):
        return None
    try:
        resseq = int(line[22:26].strip())
    except ValueError:
        return None
    return line[21].strip(), resseq, line[12:16].strip().upper()


def strip_pdb_lines(
    lines: Sequence[str],
    mask: str | StripSelection,
    keep_attachment: Optional[Tuple[str, int, str]] = None,
) -> list[str]:
    """Return the PDB lines with the atoms selected by ``mask`` removed.

    ``keep_attachment=(chain, resseq, atom_name)`` protects the attachment
    atom whatever the mask selects (chain empty = any chain).
    """
    sel = parse_strip_mask(mask) if isinstance(mask, str) else mask
    kept = []
    for line in lines:
        fields = _pdb_atom_fields(line)
        if fields is not None:
            chain, resseq, name = fields
            if sel.matches(chain, resseq, name):
                is_attachment = (
                    keep_attachment is not None
                    and name == keep_attachment[2].strip().upper()
                    and resseq == int(keep_attachment[1])
                    and (not keep_attachment[0] or chain == keep_attachment[0])
                )
                if not is_attachment:
                    continue
        kept.append(line)
    return kept


# ---------------------------------------------------------------------------
# Obstacle arrays
# ---------------------------------------------------------------------------

def strip_obstacles(
    atoms: np.ndarray,
    records: Sequence,
    mask: str | StripSelection,
) -> np.ndarray:
    """Return ``atoms`` (N, k) with the rows selected by ``mask`` removed.

    ``records`` are the parsed PDB records aligned with ``atoms`` (as
    ``fret.av._cached_pdb_records`` returns them: ``(chain, resseq, name, x, y,
    z, r)``); a length mismatch raises rather than guessing.
    """
    sel = parse_strip_mask(mask) if isinstance(mask, str) else mask
    if len(records) != atoms.shape[0]:
        raise ValueError(
            f"strip_obstacles: {len(records)} records for {atoms.shape[0]} atoms")
    keep = np.fromiter(
        (not sel.matches(rec[0], rec[1], rec[2]) for rec in records),
        dtype=bool, count=len(records))
    return atoms[keep]


# ---------------------------------------------------------------------------
# IMP hierarchies
# ---------------------------------------------------------------------------

def _hierarchy_atom_fields(atom_particle):
    """(chain, resseq, name) of an IMP atom particle (chain '' if none)."""
    import IMP.atom
    at = IMP.atom.Atom(atom_particle)
    name = at.get_atom_type().get_string()
    if name.startswith("HET:"):
        name = name[4:].strip()
    res_p = at.get_parent()
    resseq = None
    chain = ""
    if IMP.atom.Residue.get_is_setup(res_p):
        resseq = int(IMP.atom.Residue(res_p).get_index())
        chain_p = res_p.get_parent()
        if IMP.atom.Chain.get_is_setup(chain_p):
            chain = IMP.atom.Chain(chain_p).get_id()
    return chain, resseq, name.upper()


def select_atoms(hierarchy, mask: str | StripSelection):
    """The ``IMP.atom.Atom`` leaves of ``hierarchy`` selected by ``mask``."""
    import IMP.atom
    sel = parse_strip_mask(mask) if isinstance(mask, str) else mask
    out = []
    for a in IMP.atom.get_by_type(hierarchy, IMP.atom.ATOM_TYPE):
        chain, resseq, name = _hierarchy_atom_fields(a)
        if sel.matches(chain, resseq, name):
            out.append(IMP.atom.Atom(a))
    return out


def strip_hierarchy(hierarchy, mask: str | StripSelection, *, inplace: bool = False):
    """Remove the atoms selected by ``mask`` from an IMP hierarchy.

    Returns ``(hierarchy, n_removed)``. With ``inplace=False`` (default) the
    input is left untouched and a clone (``IMP.atom.create_clone``) is
    stripped and returned; with ``inplace=True`` the given hierarchy is
    modified and returned. Removed atoms are detached from their residue and
    destroyed.
    """
    import IMP.atom
    target = hierarchy if inplace else IMP.atom.create_clone(hierarchy)
    victims = select_atoms(target, mask)
    for atom in victims:
        parent = IMP.atom.Hierarchy(atom.get_parent())
        parent.remove_child(IMP.atom.Hierarchy(atom))
        IMP.atom.destroy(IMP.atom.Hierarchy(atom))
    return target, len(victims)


__all__ = [
    "BACKBONE_ATOM_NAMES",
    "StripSelection",
    "parse_strip_mask",
    "site_strip_mask",
    "default_strip_mask",
    "strip_pdb_lines",
    "strip_obstacles",
    "select_atoms",
    "strip_hierarchy",
]
