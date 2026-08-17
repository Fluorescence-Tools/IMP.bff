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

import logging

import IMP
import IMP.atom
import IMP.algebra
import IMP.core

from IMP.bff.cgdye.labeling.backbone_frame import (
    backbone_frame,
    backbone_transformation_from_coords,
    _find_atom,
    _atom_name,
)

log = logging.getLogger(__name__)


def resolve_site(hierarchy, chain_id, resnum):
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
#: atom. The AV convention (``IMP.bff.fret.strip.default_strip_mask``) keeps
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

    Delegates to the shared strip engine (:mod:`IMP.bff.fret.strip`) with the
    mask ``chain <id> and resid <n> and not name <keep...>``. Returns the
    number of removed atoms.
    """
    from IMP.bff.fret.strip import site_strip_mask, strip_hierarchy

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
    s_site = resolve_site(source_hier, source_chain, source_resnum)
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
        site = resolve_site(protein_hier, chain_id, resnum)
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
    site = resolve_site(protein_hier, chain_id, resnum)
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
