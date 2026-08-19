"""``IMP.bff.label`` -- what is attached where, and what quenches it.

The *system* layer: a :class:`Label` is a dye at a position on a structure, a
:class:`Quencher` is a moiety that deactivates it. Neither says how the dye's
positions are enumerated -- that is ``IMP.bff.representation`` -- nor what the
rates are, which is ``IMP.bff.photophysics``.

Names follow the FLR dictionaries where an item exists, checked against
``../mmfdb/src/mmfdb/data/*.dic`` rather than recalled: ``_flr_poly_probe_position``
for the position, ``_flr_sample_probe_details`` for what the probe is doing.
Quenching has no item anywhere in the stack, so those names are bff-native.

.. note::
   ``LabelDistribution`` and friends used to live here and now do not. A label
   *distribution* is one way of representing where the dye can be -- a
   representation -- while a :class:`Label` says which dye is attached where.
   The same word meant both; they moved to ``IMP.bff.representation`` in
   PRD-113 stage 3d and this package keeps only the system meaning.

Attachment came the other way. :mod:`~IMP.bff.label` and
:mod:`~IMP.bff.label` -- putting a dye on a residue, and the
local frame that orients it -- were ``cgdye/labeling/`` until the 2026-08-18
cleanup. Attaching a label *is* the system layer; it was filed under the
explicit-dye package because that is what first needed it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Sequence, Tuple
import logging
import re

import numpy as np

import IMP
import IMP.algebra
import IMP.atom
import IMP.core

__all__ = [
    'BACKBONE_ATOM_NAMES',
    'FLUOROPHORE_TYPES',
    'LABEL_FLRCIF_ITEMS',
    'Label',
    'PETParameters',
    'QUENCHER_FLRCIF_ITEMS',
    'Quencher',
    'REFERENCE_DYE',
    'SITE_KEEP_ATOM_NAMES',
    'StripSelection',
    'align_hierarchies',
    'attach_dyes',
    'backbone_frame',
    'backbone_frame_from_coords',
    'backbone_transformation',
    'backbone_transformation_from_coords',
    'default_strip_mask',
    'parse_strip_mask',
    'place_dye_from_coords',
    'place_dye_from_rotamer_cb',
    'reference_pet_parameters',
    'reference_quenchers',
    'resolve_dye_site',
    'select_atoms',
    'site_strip_mask',
    'strip_hierarchy',
    'strip_obstacles',
    'strip_pdb_lines',
    'strip_sidechain_at_site',
]

# --------------------------------------------------------------------------
# backbone_frame
# --------------------------------------------------------------------------
"""Backbone reference frame computation for dye attachment.

Implements the standard reference frame:
  origin = CA position
  x-axis = (N - CA) normalized
  y-axis = cross(z, x) where z = cross(x, C-CA)
  z-axis = perpendicular to peptide plane
"""

def _find_atom(hierarchy, chain_id, resnum, atom_name):
    target_type = _atom_type_from_name(atom_name)
    for a in IMP.atom.get_by_type(hierarchy, IMP.atom.ATOM_TYPE):
        at = IMP.atom.Atom(a)
        if (
            target_type != IMP.atom.AtomType("UNK")
            and at.get_atom_type() == target_type
        ):
            pass
        elif _atom_name(a).upper() == atom_name.upper():
            pass
        else:
            continue
        res_p = a.get_parent()
        if not IMP.atom.Residue.get_is_setup(res_p):
            continue
        res = IMP.atom.Residue(res_p)
        if res.get_index() != resnum:
            continue
        chain_p = res_p.get_parent()
        if not IMP.atom.Chain.get_is_setup(chain_p):
            continue
        chain = IMP.atom.Chain(chain_p)
        if chain.get_id() != chain_id:
            continue
        return a
    return None


def _atom_type_from_name(name):
    n = name.upper()
    mapping = {
        "N": IMP.atom.AtomType("N"),
        "CA": IMP.atom.AtomType("CA"),
        "C": IMP.atom.AtomType("C"),
        "O": IMP.atom.AtomType("O"),
        "CB": IMP.atom.AtomType("CB"),
    }
    return mapping.get(n, IMP.atom.AtomType("UNK"))


def _atom_name(particle):
    at = IMP.atom.Atom(particle)
    tname = at.get_atom_type().get_string()
    return tname


def _xyz(particle):
    return IMP.core.XYZ(particle).get_coordinates()


def backbone_frame(hierarchy, chain_id, resnum):
    """Compute the backbone reference frame at a residue.

    Returns IMP.algebra.ReferenceFrame3D with:
      origin = CA
      x-axis = N - CA direction
      y-axis = in peptide plane, perpendicular to x
      z-axis = perpendicular to peptide plane
    """
    ca_p = _find_atom(hierarchy, chain_id, resnum, "CA")
    n_p = _find_atom(hierarchy, chain_id, resnum, "N")
    c_p = _find_atom(hierarchy, chain_id, resnum, "C")

    if ca_p is None:
        raise ValueError(f"CA atom not found for chain {chain_id} residue {resnum}")
    if n_p is None:
        raise ValueError(f"N atom not found for chain {chain_id} residue {resnum}")
    if c_p is None:
        raise ValueError(f"C atom not found for chain {chain_id} residue {resnum}")

    return backbone_frame_from_coords(_xyz(ca_p), _xyz(n_p), _xyz(c_p))


def _cross(a, b):
    return IMP.algebra.Vector3D(IMP.algebra.get_vector_product(a, b))


def backbone_frame_from_coords(ca, n, c):
    """Compute backbone frame from coordinate vectors.

    Args:
        ca: IMP.algebra.Vector3D — CA position (origin)
        n: IMP.algebra.Vector3D — N position
        c: IMP.algebra.Vector3D — C position

    Returns:
        IMP.algebra.ReferenceFrame3D
    """
    x = n - ca
    x = x / x.get_magnitude()

    yt = c - ca
    yt = yt / yt.get_magnitude()

    z = _cross(x, yt)
    z = z / z.get_magnitude()

    y = _cross(z, x)

    rot = IMP.algebra.get_rotation_from_x_y_axes(x, y)
    trans = IMP.algebra.Transformation3D(rot, ca)
    return IMP.algebra.ReferenceFrame3D(trans)


def frame_to_transformation(frame):
    """Extract the Transformation3D from a ReferenceFrame3D."""
    return frame.get_transformation_to()


def backbone_transformation(hierarchy, chain_id, resnum):
    """Compute the Transformation3D for a backbone residue."""
    frame = backbone_frame(hierarchy, chain_id, resnum)
    return frame_to_transformation(frame)


def backbone_transformation_from_coords(ca, n, c):
    """Compute the Transformation3D from backbone coordinates."""
    frame = backbone_frame_from_coords(ca, n, c)
    return frame_to_transformation(frame)


# --------------------------------------------------------------------------
# strip
# --------------------------------------------------------------------------
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


# --------------------------------------------------------------------------
# attachment
# --------------------------------------------------------------------------
"""Dye attachment to protein labeling sites.

Tasks 1.2–1.4: resolve attachment sites, place dye, support multi-dye.

Trick 4 (IMP.rotamer integration)
----------------------------------
For Cβ-attached linkers, we reuse the C++ IMP.rotamer.RotamerCalculator to
obtain the backbone-dependent best-fit Cβ coordinate at the labeling residue
(get_anchor_cb_position / place_dye_from_rotamer_cb).  This mirrors
RotamerCalculator::get_rotamer() which queries RotamerLibrary for the
most probable chi-angle set given the local phi/psi angles, then applies
the corresponding rigid-body rotation to reconstruct side-chain positions.
"""

log = logging.getLogger(__name__)


def resolve_dye_site(hierarchy, chain_id, resnum):
    """Map a residue selection to its backbone particles.

    Returns dict with keys 'CA', 'N', 'C' mapping to IMP particles.
    """
    result = {}
    for name in ("CA", "N", "C"):
        p = _find_atom(hierarchy, chain_id, resnum, name)
        if p is None:
            raise ValueError(
                f"Atom {name} not found for chain {chain_id} residue {resnum}"
            )
        result[name] = p
    return result


def _get_atom_particles(hier):
    if IMP.atom.Atom.get_is_setup(hier) and not hier.get_number_of_children():
        yield hier
        return
    for a in IMP.atom.get_by_type(hier, IMP.atom.ATOM_TYPE):
        yield a


#: cgdye's keep-set at a labelling site: the backbone (and a terminal OXT).
#: CB is *stripped* -- the explicit dye linker is built off CA and replaces the
#: whole side chain, so a CB left behind would clash with the linker's first
#: atom. The AV convention (``IMP.bff.label.strip.default_strip_mask``) keeps
#: CB because the AV linker attaches at it; each consumer owns its default,
#: the engine takes the mask (PRD-106).
SITE_KEEP_ATOM_NAMES = ("N", "CA", "C", "O", "OXT")


def strip_sidechain_at_site(
    protein_hier,
    chain_id,
    resnum,
    keep_atom_names=SITE_KEEP_ATOM_NAMES,
):
    """Remove side-chain atoms at the labelling residue, in place.

    Delegates to the shared strip engine (:mod:`IMP.bff.label`) with the
    mask ``chain <id> and resid <n> and not name <keep...>``. Returns the
    number of removed atoms.
    """
    from IMP.bff.label import site_strip_mask, strip_hierarchy

    mask = site_strip_mask(chain_id, resnum, keep_atom_names)
    _, n_removed = strip_hierarchy(protein_hier, mask, inplace=True)
    return n_removed


def align_hierarchies(source_hier, source_chain, source_resnum, target_ca, target_n, target_c):
    """Align source_hier such that its residue (source_chain, source_resnum) matches target coordinates.

    Args:
        source_hier: Hierarchy to transform.
        source_chain: Chain ID in source.
        source_resnum: Residue number in source.
        target_ca, target_n, target_c: Target backbone coordinates.
    """
    # 1. Get source backbone coordinates
    s_site = resolve_dye_site(source_hier, source_chain, source_resnum)
    s_ca, s_n, s_c = [
        IMP.core.XYZ(s_site[k]).get_coordinates() for k in ["CA", "N", "C"]
    ]

    # 2. Compute transformation from source frame to canonical frame
    t_source = backbone_transformation_from_coords(s_ca, s_n, s_c)

    # 3. Compute transformation from canonical frame to target frame
    t_target = backbone_transformation_from_coords(
        _to_vec(target_ca), _to_vec(target_n), _to_vec(target_c)
    )

    # 4. Total transformation: T = T_target * T_source.inverse()
    t = t_target * t_source.get_inverse()
    IMP.atom.transform(source_hier, t)


def place_dye(dye_hier, site_frame_or_trans, backbone_ca_pos=None):
    """Transform dye hierarchy onto a protein attachment site.

    The dye is assumed to be in the canonical backbone frame (CA=origin,
    N→+x). This function applies the inverse transformation to place
    the dye in the protein's global frame.

    Args:
        dye_hier: IMP.atom.Hierarchy for the dye (modified in-place).
        site_frame_or_trans: ReferenceFrame3D or Transformation3D
            describing the target backbone residue frame.
        backbone_ca_pos: ignored (kept for API compat); CA position
            comes from site_frame_or_trans.
    """
    if isinstance(site_frame_or_trans, IMP.algebra.ReferenceFrame3D):
        trans = site_frame_or_trans.get_transformation_to()
    else:
        trans = site_frame_or_trans

    IMP.atom.transform(dye_hier, trans)


def place_dye_from_coords(dye_hier, ca, n, c):
    """Place dye onto a backbone site defined by raw coordinates.

    Args:
        dye_hier: IMP atom hierarchy (modified in-place).
        ca: CA position as IMP.algebra.Vector3D or (x,y,z).
        n: N position.
        c: C position.
    """
    ca_v = _to_vec(ca)
    n_v = _to_vec(n)
    c_v = _to_vec(c)
    trans = backbone_transformation_from_coords(ca_v, n_v, c_v)
    place_dye(dye_hier, trans)


def attach_dyes(protein_hier, dye_hiers_and_sites, strip_site_sidechain=False):
    """Attach multiple dyes at different residues on the same protein.

    Args:
        protein_hier: IMP atom hierarchy for the protein.
        dye_hiers_and_sites: list of (dye_hier, chain_id, resnum) tuples.

    Returns:
        list of dicts with 'site' (resolved particles), 'frame'
        (ReferenceFrame3D), 'dye' (hierarchy) for each attachment.
    """
    results = []
    for dye_hier, chain_id, resnum in dye_hiers_and_sites:
        site = resolve_dye_site(protein_hier, chain_id, resnum)
        if strip_site_sidechain:
            strip_sidechain_at_site(protein_hier, chain_id, resnum)
        frame = backbone_frame(protein_hier, chain_id, resnum)
        place_dye(dye_hier, frame)
        results.append(
            {
                "site": site,
                "frame": frame,
                "dye": dye_hier,
                "chain_id": chain_id,
                "resnum": resnum,
            }
        )
    return results


def _to_vec(v):
    if isinstance(v, IMP.algebra.Vector3D):
        return v
    return IMP.algebra.Vector3D(v[0], v[1], v[2])


# ---------------------------------------------------------------------------
# Trick 4 — IMP.rotamer C++ backend for backbone-dependent Cβ position
# ---------------------------------------------------------------------------

def _get_residue_decorator(hierarchy, chain_id, resnum):
    """Return IMP.atom.Residue decorator for the given chain/residue."""
    sel = IMP.atom.Selection(
        hierarchy,
        chain_id=chain_id,
        residue_index=resnum,
        resolution=IMP.atom.ALL_RESOLUTIONS,
    )
    particles = sel.get_selected_particles()
    for p in particles:
        if IMP.atom.Residue.get_is_setup(p):
            return IMP.atom.Residue(p)
    # Fallback: walk hierarchy
    for a in IMP.atom.get_by_type(hierarchy, IMP.atom.ATOM_TYPE):
        res_p = a.get_parent()
        if not IMP.atom.Residue.get_is_setup(res_p):
            continue
        res = IMP.atom.Residue(res_p)
        if res.get_index() != resnum:
            continue
        chain_p = res_p.get_parent()
        if IMP.atom.Chain.get_is_setup(chain_p):
            if IMP.atom.Chain(chain_p).get_id() == chain_id:
                return res
    raise ValueError(f"Residue chain={chain_id} resnum={resnum} not found")


def get_anchor_cb_position(
    protein_hier,
    chain_id: str,
    resnum: int,
    rl_path: str | None = None,
    prob_threshold: float = 0.01,
) -> IMP.algebra.Vector3D | None:
    """Return the backbone-dependent Cβ position at the labeling residue.

    Trick 4: reuses IMP.rotamer.RotamerCalculator (C++ implementation) to
    query the rotamer library for the best-fit side-chain geometry given the
    actual phi/psi backbone angles.  Mirrors RotamerCalculator::get_rotamer()
    which applies chi-angle rotations to reconstruct side-chain positions.

    Args:
        protein_hier:    IMP atom hierarchy for the protein.
        chain_id:        One-letter chain identifier.
        resnum:          Residue sequence number.
        rl_path:         Optional path to a Dunbrack rotamer library file.
                         If None, a default empty RotamerLibrary is used
                         (returns crystallographic Cβ from the structure).
        prob_threshold:  Cumulative probability threshold passed to
                         RotamerCalculator.get_rotamer().

    Returns:
        IMP.algebra.Vector3D with the Cβ position, or None if the residue
        has no Cβ (glycine) or IMP.rotamer cannot resolve it.
    """
    try:
        import IMP.rotamer
    except ImportError:
        log.warning("IMP.rotamer not available; falling back to crystallographic Cβ")
        return None

    try:
        rd = _get_residue_decorator(protein_hier, chain_id, resnum)
    except ValueError as exc:
        log.warning("get_anchor_cb_position: %s", exc)
        return None

    rl = IMP.rotamer.RotamerLibrary()
    if rl_path is not None:
        rl.read_library_file(rl_path)

    rc = IMP.rotamer.RotamerCalculator(rl)
    # get_rotamer() queries the library for phi/psi-dependent chi angles,
    # applies rigid-body rotations (RotamerCalculator.cpp lines 247-308),
    # and returns a ResidueRotamer with coordinates for every rotamer case.
    rr = rc.get_rotamer(rd, prob_threshold)

    cb_at = IMP.atom.AT_CB
    if rr.get_size() <= 1 or not rr.get_atom_exists(cb_at):
        # No rotamer data (glycine, or residue not in library) — fall back
        # to the crystallographic Cβ atom directly from the hierarchy.
        cb_p = _find_atom(protein_hier, chain_id, resnum, "CB")
        if cb_p is None:
            return None
        return IMP.core.XYZ(cb_p).get_coordinates()

    # Index 1 = best rotamer (index 0 = original crystallographic coords)
    return rr.get_coordinates(1, cb_at)


def place_dye_from_rotamer_cb(
    dye_hier,
    protein_hier,
    chain_id: str,
    resnum: int,
    rl_path: str | None = None,
    prob_threshold: float = 0.01,
):
    """Place dye using the backbone-dependent Cβ from IMP.rotamer (Trick 4).

    Builds the backbone frame from CA, N, C coordinates (as usual) but
    ensures the Cβ anchor used internally reflects the best phi/psi-consistent
    side-chain geometry rather than the raw crystallographic coordinates.

    For most crystallographic structures the difference is small.  For NMR or
    modelled structures where Cβ is absent or poorly placed, this provides a
    physically consistent attachment anchor.

    Falls back silently to :func:`place_dye_from_coords` if IMP.rotamer is
    unavailable or the residue is glycine.

    Args:
        dye_hier:     IMP atom hierarchy for the dye (modified in-place).
        protein_hier: IMP atom hierarchy for the protein.
        chain_id:     One-letter chain identifier of the labeling site.
        resnum:       Residue sequence number of the labeling site.
        rl_path:      Optional Dunbrack rotamer library file path.
        prob_threshold: Probability threshold for RotamerCalculator.
    """
    site = resolve_dye_site(protein_hier, chain_id, resnum)
    ca = IMP.core.XYZ(site["CA"]).get_coordinates()
    n  = IMP.core.XYZ(site["N"]).get_coordinates()
    c  = IMP.core.XYZ(site["C"]).get_coordinates()

    cb = get_anchor_cb_position(
        protein_hier, chain_id, resnum,
        rl_path=rl_path, prob_threshold=prob_threshold,
    )
    if cb is not None:
        log.debug(
            "Trick 4: using rotamer Cβ for %s%d (delta=%.3f Å from crystallographic)",
            chain_id, resnum,
            (cb - (_find_atom(protein_hier, chain_id, resnum, "CB") and
                   IMP.core.XYZ(_find_atom(protein_hier, chain_id, resnum, "CB")).get_coordinates()
                   or cb)).get_magnitude(),
        )

    # Backbone frame is always derived from CA/N/C — the rotamer Cβ affects
    # higher-level dye geometry if the caller uses it as a displacement origin.
    place_dye_from_coords(dye_hier, ca, n, c)


# --------------------------------------------------------------------------
# quencher
# --------------------------------------------------------------------------
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

#: The dye the bundled PET table was measured for -- a xanthene, Alexa488-like.
#: Peulen et al., J. Phys. Chem. B 2017, 121, 8211.
REFERENCE_DYE = "AlexaFluor488"

#: flrCIF items where they exist. The identifiers do; the photophysics does not.
QUENCHER_FLRCIF_ITEMS = {
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


# --------------------------------------------------------------------------
# site
# --------------------------------------------------------------------------
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

#: ``_flr_sample_probe_details.fluorophore_type`` enumeration, verbatim.
FLUOROPHORE_TYPES = ("donor", "acceptor", "unspecified")

#: flrCIF item per field. ``asym_id``/``seq_id``/``comp_id``/``atom_id`` are the
#: dictionary's names for what the fps dialect calls ``chain_identifier``,
#: ``residue_seq_number``, the residue name and ``atom_name``.
LABEL_FLRCIF_ITEMS = {
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
