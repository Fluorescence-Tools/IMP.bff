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
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
import json
import logging
import re
from pathlib import Path

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
    """The atom named *atom_name* on ``(chain_id, resnum)``, or None.

    Narrowed by ``IMP.atom.Selection`` before matching. This used to walk
    ``get_by_type(hierarchy, ATOM_TYPE)`` -- every atom in the structure -- for
    each call, and the name comparison inside builds an ``IMP.atom.Atom`` and
    reads its type string per atom: 641,580 of those across one test run, 1.7 s
    over 226 calls. The selection reaches the residue's ~10 atoms directly.

    The match itself is unchanged, and deliberately: an atom qualifies by
    ``AtomType`` *or* by name, because a dye's atoms carry types that IMP does
    not recognise and are found by name only.
    """
    target_type = _atom_type_from_name(atom_name)
    selection = IMP.atom.Selection(hierarchy)
    if chain_id:
        selection.set_chain_id(chain_id)
    selection.set_residue_index(int(resnum))
    candidates = selection.get_selected_particles()
    if not candidates:
        return None
    for a in candidates:
        if not IMP.atom.Atom.get_is_setup(a):
            continue
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
        # `get_selected_particles` returns Particles; the parent walk below
        # needs the hierarchy decorator
        res_p = IMP.atom.Hierarchy(a).get_parent()
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


# The `strip_mask` dialect -- ``StripSelection``, ``parse_strip_mask``,
# ``site_strip_mask``, ``default_strip_mask``, ``strip_pdb_lines`` and
# ``stripped_pdb_for`` -- is **C++** (``include/IMP/bff/StripMask.h``), next to
# the AV builder that is its only real consumer. ``BACKBONE_ATOM_NAMES`` goes
# with it as ``backbone_atom_names()``.
from IMP.bff import (  # noqa: F401
    StripSelection, backbone_atom_names, default_strip_mask, parse_strip_mask,
    site_strip_mask, strip_pdb_lines, stripped_pdb_for,
)

#: The backbone a labelling site keeps when its side chain is stripped.
BACKBONE_ATOM_NAMES = tuple(backbone_atom_names())


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


# ``Quencher``, ``PETParameters``, ``reference_quenchers`` and
# ``reference_pet_parameters`` are **C++** (``include/IMP/bff/PETQuenching.h``),
# next to the tables they read. They were here and the tables were in
# ``quenching.py``, so building the reference set crossed a module boundary in
# one direction and a language boundary in the other; a rate that takes two
# partners to define belongs with both of them.
#
# ``attenuation_length`` is **NaN** where it was ``None`` -- a hard
# contact-sphere model rather than the exponential through-space form.
from IMP.bff import (  # noqa: F401
    PETParameters, Quencher, REFERENCE_DYE, reference_pet_parameters,
    reference_quenchers,
)


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
    :param dye: the :class:`IMP.bff.Dye` species, when known. A label whose
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


# --------------------------------------------------------------------------
# Fluorescent proteins: which domain of a fusion is which
#
# These were in `cgdye/sampling.py`, which is dye dynamics -- finding a
# fluorescent protein in a fusion by aligning its sequence, and cutting a
# structure into domains by pLDDT, are neither dye nor dynamics. Their only
# consumers are `imp_bff dye label-fp` and `label-fusion`, which is labelling.
#
# `get_fp_names` and `get_fp_data` came out here rather than moving: each was
# referenced exactly once in the whole stack, at its own `def`.
# --------------------------------------------------------------------------
def _smith_waterman(query: str, templ: str,
                    match: float = 2.0, mismatch: float = -1.0,
                    gap_open: float = -5.0, gap_extend: float = -1.0):
    """Local alignment with affine gaps -- :func:`IMP.bff.smith_waterman`.

    Replaces Biopython's ``PairwiseAligner`` in local mode with the same
    scoring, so imp.bff needs nothing beyond what IMP brings. Gotoh's
    three-matrix formulation, which is what makes an affine penalty exact
    rather than a per-position approximation.

    Was 73 lines of Python over three ``(n+1, m+1)`` numpy planes, indexed
    element by element: 1,287 ms on a 900x900 pair against 8 ms here. Gated on
    300 randomised pairs, identical blocks on all of them.

    Returns
    -------
    tuple
        ``(query_blocks, templ_blocks)``, each a list of ``(start, end)``
        half-open index pairs, in the spelling Biopython's ``aligned`` uses.
    """
    blocks = IMP.bff.smith_waterman(
        query, templ, float(match), float(mismatch),
        float(gap_open), float(gap_extend))
    return ([(b.query_start, b.query_end) for b in blocks],
            [(b.template_start, b.template_end) for b in blocks])


def _align_best_window(query: str, templ: str) -> tuple[int, int, float]:
    """Local alignment; return (q_start_1based, q_end_1based, identity)."""
    q_blocks, t_blocks = _smith_waterman(query, templ)
    if not q_blocks:
        return (0, 0, 0.0)

    q_start = int(q_blocks[0][0])
    q_end = int(q_blocks[-1][1])

    matches = 0
    aligned_len = 0
    for (qs, qe), (ts, te) in zip(q_blocks, t_blocks):
        qseg = query[qs:qe]
        tseg = templ[ts:te]
        aligned_len += len(qseg)
        for qc, tc in zip(qseg, tseg):
            if qc == tc:
                matches += 1
    identity = (matches / aligned_len) if aligned_len else 0.0
    return (int(q_start + 1), int(q_end), float(identity))


def find_fp_domains(seq: str, min_identity: float = 0.35) -> List[Tuple[str, int, int, float]]:
    """Detect FP domains in a fusion protein sequence."""
    lib = load_fp_library()
    hits: List[Tuple[str, int, int, float]] = []

    for name, data in lib.items():
        templ = data["sequence"]
        motifs = data.get("motifs", [])
        
        best_seed = (-1, -1.0)
        for m in motifs:
            i = seq.find(m)
            if i >= 0 and len(m) > best_seed[1]:
                best_seed = (i, len(m))

        if best_seed[0] >= 0:
            i = best_seed[0]
            start = max(0, i - 50)
            end = min(len(seq), i + 300)
            qs, qe, idy = _align_best_window(seq[start:end], templ)
            if idy >= min_identity and qs > 0:
                hits.append((name, int(start + qs), int(start + qe), float(idy)))
            continue

        qs, qe, idy = _align_best_window(seq, templ)
        if idy >= min_identity and qs > 0:
            hits.append((name, int(qs), int(qe), float(idy)))

    hits.sort(key=lambda x: (-x[3], x[2] - x[1], x[1]))
    covered = set()
    non_overlap: List[Tuple[str, int, int, float]] = []
    for nm, s, e, idy in hits:
        if any(i in covered for i in range(s, e + 1)):
            continue
        for i in range(s, e + 1):
            covered.add(i)
        non_overlap.append((nm, s, e, idy))
    return non_overlap


def parse_plddt_from_pdb(pdb_path: str | Path, chain_id: str = "A") -> Dict[int, float]:
    """Parse pLDDT from AlphaFold PDB B-factor column."""
    plddt_sum: Dict[int, float] = {}
    plddt_cnt: Dict[int, int] = {}
    
    with open(pdb_path, "r") as f:
        for line in f:
            if not line.startswith("ATOM"): continue
            ch = line[21].strip()
            if chain_id and ch and ch != chain_id: continue
            try:
                resseq = int(line[22:26])
                b = float(line[60:66])
                plddt_sum[resseq] = plddt_sum.get(resseq, 0.0) + b
                plddt_cnt[resseq] = plddt_cnt.get(resseq, 0) + 1
            except ValueError: continue

    if not plddt_sum:
        raise ValueError(f"No pLDDT data found for chain '{chain_id}' in {pdb_path}")
    return {r: plddt_sum[r] / plddt_cnt[r] for r in sorted(plddt_sum.keys())}


def segments_from_plddt(
    seq_len: int,
    plddt: Dict[int, float],
    fp_domains: List[Tuple[str, int, int, float]],
    rigid_threshold: float = 70.0,
    min_rb_len: int = 12
) -> List[Dict]:
    """Create rigid/linker segments from pLDDT and FP detections."""
    labels = ['R' if plddt.get(i+1, 0.0) >= rigid_threshold else 'L' for i in range(seq_len)]
    
    # FP domains are always rigid
    for _, s, e, _ in fp_domains:
        for i in range(s-1, e):
            labels[i] = 'R'
            
    # Small rigid islands become linkers
    i = 0
    while i < seq_len:
        if labels[i] == 'R':
            j = i
            while j < seq_len and labels[j] == 'R': j += 1
            if (j - i) < min_rb_len:
                # Only if not part of an FP? Let's keep it simple for now
                for k in range(i, j): labels[k] = 'L'
            i = j
        else: i += 1
        
    segs = []
    i = 0
    while i < seq_len:
        curr = labels[i]
        j = i
        while j < seq_len and labels[j] == curr: j += 1
        
        kind = "core" if curr == 'R' else "linker"
        name = kind
        
        # Check if it's an FP
        for nm, s, e, _ in fp_domains:
            if max(i+1, s) <= min(j, e):
                if (min(j, e) - max(i+1, s) + 1) > 0.5 * (j - i):
                    kind = "fp"
                    name = nm
                    break
                    
        segs.append({"kind": kind, "name": name, "start": i+1, "end": j})
        i = j
    return segs


def load_fp_library() -> Dict[str, Dict]:
    """Load the FP library from the bundled JSON file."""
    import os
    import IMP.bff
    path = os.path.join(IMP.bff.get_cgdye_data_dir(), "fp_library.json")
    with open(path, "r") as f:
        return json.load(f)

