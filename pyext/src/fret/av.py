"""Accessible-volume computation on IMP.bff's own AV machinery.

Self-contained on ``IMP.bff.AV``/``PathMap`` — there is no alternative backend
and no optional dependency. (The pre-move ChiSurf implementation carried a
LabelLib backend and selection logic; that was dropped deliberately when the
module moved here: nothing in the IMP stack depends on LabelLib.)
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Dict, Optional, Tuple

import numpy as np

from . import io

logger = logging.getLogger(__name__)

import IMP
import IMP.algebra
import IMP.atom
import IMP.bff
import IMP.core
import IMP.em


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class AccessibleVolume:
    """Represents the 3D positional distribution of a dye label."""

    points: np.ndarray  # (N, 4) float64 — xyz + weight
    density: np.ndarray  # (nx, ny, nz) float32 — 3D voxel density
    grid_origin: np.ndarray  # (3,) float64
    grid_step: float
    grid_shape: Tuple[int, int, int]
    attachment_point: np.ndarray  # (3,) float64
    position_name: str = ""
    params: Dict = field(default_factory=dict)

    @property
    def n_points(self) -> int:
        return self.points.shape[0] if self.points.ndim == 2 else 0

    @property
    def mean_position(self) -> np.ndarray:
        if self.n_points == 0:
            return self.attachment_point.copy()
        w = self.points[:, 3]
        if w.sum() == 0:
            return self.attachment_point.copy()
        return np.average(self.points[:, :3], axis=0, weights=w)

    @property
    def has_volume(self) -> bool:
        return self.n_points > 0


# ---------------------------------------------------------------------------
# Backend interface (compatibility)
# ---------------------------------------------------------------------------

def select_backend(name: str) -> None:
    """Accept the AV backend selection of the pre-move API.

    ``'auto'`` and ``'imp-bff'`` are the only valid values — IMP.bff is the
    one backend. Kept so callers written against the pre-move API keep
    working; any other name raises.

    Raises
    ------
    ValueError
        If *name* names a backend that does not exist here (including
        ``'labellib'``, which was removed in the move).
    """
    if name not in ("auto", "imp-bff"):
        raise ValueError(
            f"Unknown backend name: '{name}' — IMP.bff is the only AV "
            "backend (the LabelLib backend was removed).")


def _active_backend_name() -> str:
    """Return the active backend name for diagnostics."""
    return "imp-bff"


# ---------------------------------------------------------------------------
# AV computation
# ---------------------------------------------------------------------------

def _av_imp_bff(
    pdb_path: str,
    source_info: Dict,
    linker_length: float,
    linker_width: float,
    radii: Tuple[float, float, float],
    disc_step: float,
) -> AccessibleVolume:
    """Compute an AV with ``IMP.bff.AV``."""
    chain = source_info.get("chain_identifier", "")
    resseq = source_info.get("residue_seq_number", 0)
    aname = source_info.get("atom_name", "CA")

    # The FPS strip: positions are calibrated for a structure whose
    # attachment residue does not wall in its own dye, so the side chain
    # goes (minus the attachment atom) before anything is measured. A
    # declared ``strip_mask`` outside the fps dialect raises here -- loud,
    # because the alternative is computing against obstacles the document
    # said to remove.
    strip_mask = str(source_info.get("strip_mask") or "").strip()
    pdb_path = _stripped_pdb_for(pdb_path, chain, resseq, aname, strip_mask)

    model = IMP.Model()
    hierarchy = IMP.atom.read_pdb(pdb_path, model, IMP.atom.NonWaterPDBSelector())

    sel = IMP.atom.Selection(hierarchy)
    if chain:
        sel.set_chain_id(chain)
    sel.set_residue_index(int(resseq))
    sel.set_atom_type(IMP.atom.AtomType(str(aname)))
    particles = sel.get_selected_particles()
    if not particles:
        raise ValueError(f"Attachment site {chain}:{resseq}:{aname} not found")
    attachment_particle = particles[0]

    # Source clearance. The path search inflates obstacles by half the linker
    # width, so the free sphere around the attachment atom has to clear that
    # inflation (plus a grid step of slack) or the source tile is walled in
    # and the AV comes back empty. The strip above already removes the
    # attachment residue's side chain, which is what lets FPS-calibrated
    # small clearances (``allowed_sphere_radius: 1``) compute a real cloud.
    # A position that declares allowed_sphere_radius (schema field) keeps
    # its own value -- honoured exactly, never escalated: an empty AV at
    # the declared parameters is the answer, not a signal to retry at
    # invented ones.
    default_clearance = max(1.5, 0.5 * linker_width + 0.5 * disc_step)
    allowed_sphere_radius = float(
        source_info.get("allowed_sphere_radius", default_clearance))

    av_particle = IMP.Particle(model)
    r1, r2, r3 = radii
    IMP.bff.AV.do_setup_particle(
        model,
        av_particle,
        attachment_particle,
        linker_length=linker_length,
        linker_width=linker_width,
        radii=[r1, r2, r3],
        allowed_sphere_radius=allowed_sphere_radius,
        contact_volume_thickness=float(
            source_info.get("contact_volume_thickness", 0.0)),
        contact_volume_trapped_fraction=float(
            source_info.get("contact_volume_trapped_fraction", -1)),
        simulation_grid_resolution=disc_step,
    )
    av = IMP.bff.AV(model, av_particle)
    av.resample()

    path_map = av.get_map()
    header = path_map.get_header()
    nx, ny, nz = header.get_nx(), header.get_ny(), header.get_nz()
    tile_values = path_map.get_tile_values(
        IMP.bff.PM_TILE_ACCESSIBLE_DENSITY, (0.0, path_map.get_path_map_header().get_max_path_length())
    )
    density = np.asarray(tile_values, dtype=np.float32).reshape((nx, ny, nz), order="C")

    # get_xyz_density() returns an (N, 4) array; a bare truth test on it
    # raises under numpy, which the pre-move code never saw because the
    # LabelLib fallback swallowed the exception. Weights are uniform: the
    # density column is not normalised the way the samplers expect.
    xyz_density = np.asarray(path_map.get_xyz_density(), dtype=np.float64)
    if xyz_density.size:
        points = np.column_stack(
            [xyz_density[:, :3], np.ones(xyz_density.shape[0])])
    else:
        points = np.zeros((0, 4), dtype=np.float64)

    # IMP.em.DensityHeader.get_origin takes the axis index — the argument-less
    # call was another line the LabelLib fallback had been hiding.
    origin = np.array(
        [header.get_origin(0), header.get_origin(1), header.get_origin(2)],
        dtype=np.float64)
    step = header.get_spacing()
    att_xyz = np.array(av.get_source_coordinates(), dtype=np.float64)

    return AccessibleVolume(
        points=points,
        density=density,
        grid_origin=origin,
        grid_step=float(step),
        grid_shape=(nx, ny, nz),
        attachment_point=att_xyz,
    )


def compute_av(
    atoms: np.ndarray,
    source_xyz: np.ndarray,
    linker_length: float,
    linker_width: float,
    radii: Tuple[float, float, float],
    disc_step: float = 1.5,
    pdb_path: Optional[str] = None,
    source_info: Optional[Dict] = None,
) -> AccessibleVolume:
    """Compute an accessible volume with IMP.bff.

    Parameters
    ----------
    atoms : (N, 4) float64
        Columns: x, y, z, vdw_radius. Kept for API compatibility — the AV is
        computed from ``pdb_path``, which carries atom identity.
    source_xyz : (3,) float64
        Attachment point coordinates (informational; the attachment atom is
        resolved from ``source_info``).
    pdb_path : str
        Structure the AV is computed on. Required.
    source_info : dict
        Position definition (fps.json position fields). Required.
    """
    if pdb_path is None or source_info is None:
        raise ValueError("compute_av requires pdb_path and source_info")
    return _av_imp_bff(
        pdb_path,
        source_info,
        linker_length,
        linker_width,
        radii,
        disc_step,
    )


def compute_avs_for_structure(
    atoms: np.ndarray,
    positions: Dict,
    pdb_path: str | list[str] | None = None,
    disc_step: Optional[float] = None,
) -> Dict[str, AccessibleVolume]:
    """Compute AVs for all positions in an fps.json ``Positions`` dict.

    Parameters
    ----------
    atoms : (N, 4) float64
        xyzr from :func:`load_structure_with_vdw`. Used as fallback if pdb_path is not given.
    positions : dict
        fps.json Positions section.
    """
    if isinstance(pdb_path, (list, tuple)):
        pdb_paths = list(pdb_path)
    elif isinstance(pdb_path, str) and "," in pdb_path:
        pdb_paths = [p.strip() for p in pdb_path.split(",")]
    elif isinstance(pdb_path, str):
        pdb_paths = [pdb_path]
    else:
        pdb_paths = []

    avs: Dict[str, AccessibleVolume] = {}
    for pname, pdef in positions.items():
        bi = int(pdef.get("body_id", 0))
        curr_pdb = pdb_paths[bi] if bi < len(pdb_paths) else (pdb_paths[0] if pdb_paths else None)

        if curr_pdb is not None:
            curr_atoms = load_structure_with_vdw(curr_pdb)
        else:
            curr_atoms = atoms

        ll = float(pdef.get("linker_length", 20.0))
        lw = float(pdef.get("linker_width", 1.0))
        r1 = float(pdef.get("radius1", 3.5))
        r2 = float(pdef.get("radius2", 0.0))
        r3 = float(pdef.get("radius3", 0.0))
        ds = float(disc_step or pdef.get("simulation_grid_resolution", 1.5))

        chain = pdef.get("chain_identifier", "")
        resseq = pdef.get("residue_seq_number", 0)
        aname = pdef.get("atom_name", "CA")

        source_xyz = _find_attachment_point(curr_atoms, chain, resseq, aname, pdb_path=curr_pdb)
        if source_xyz is None:
            avs[pname] = AccessibleVolume(
                points=np.zeros((0, 4), dtype=np.float64),
                density=np.zeros((1, 1, 1), dtype=np.float32),
                grid_origin=np.zeros(3),
                grid_step=ds,
                grid_shape=(1, 1, 1),
                attachment_point=np.zeros(3),
                position_name=pname,
            )
            continue

        av = compute_av(
            atoms=curr_atoms,
            source_xyz=source_xyz,
            linker_length=ll,
            linker_width=lw,
            radii=(r1, r2, r3),
            disc_step=ds,
            pdb_path=curr_pdb,
            source_info=pdef,
        )
        av.position_name = pname
        av.params = pdef
        avs[pname] = av
    return avs


# ---------------------------------------------------------------------------
# Helper: vdW radii
# ---------------------------------------------------------------------------

# From FPS data/vdW.txt (selected common elements)
VDW_RADII = {
    1: 1.20, 2: 1.40, 3: 1.82, 4: 1.53, 5: 1.92, 6: 1.70, 7: 1.55,
    8: 1.52, 9: 1.47, 12: 1.73, 14: 2.10, 15: 1.80, 16: 1.80,
    17: 1.75, 19: 2.27, 20: 1.97, 26: 1.56, 30: 1.39,
}
_DEFAULT_VDW = 1.70
_ELEMENT_NUMBERS = {
    "H": 1,
    "HE": 2,
    "LI": 3,
    "BE": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "MG": 12,
    "SI": 14,
    "P": 15,
    "S": 16,
    "CL": 17,
    "K": 19,
    "CA": 20,
    "FE": 26,
    "ZN": 30,
}


def _element_symbol_from_pdb_line(line: str) -> str:
    """Return an element symbol parsed from a PDB ATOM/HETATM line.

    Columns 77-78 carry the element in a modern PDB and are used verbatim when
    present. Older files (the shipped FPS screening structures are 66-character
    records) leave them empty, so the element is read from the atom-name field
    instead. There the element is *right-justified in columns 13-14*: a blank or
    numeric column 13 means a one-letter element, so ``" CA "`` is an
    α-carbon while ``"CA  "`` is calcium. A two-letter reading is additionally
    rejected when columns 15-16 contain a digit, which is how a four-character
    hydrogen name such as ``"HE21"`` is written — helium would otherwise win.

    Parameters
    ----------
    line : str
        PDB ATOM or HETATM record.

    Returns
    -------
    str
        Uppercase element symbol, or an empty string when it cannot be parsed.
    """
    symbol = line[76:78].strip().upper() if len(line) >= 78 else ""
    if symbol:
        return symbol

    name_field = line[12:16].ljust(4).upper()
    candidate = name_field[:2].strip()
    if (
        len(candidate) == 2
        and candidate in _ELEMENT_NUMBERS
        and not any(ch.isdigit() for ch in name_field[2:])
    ):
        return candidate
    letters = "".join(ch for ch in name_field if ch.isalpha())
    if letters[:1] not in _ELEMENT_NUMBERS and letters[:2] in _ELEMENT_NUMBERS:
        # Left-padding a two-letter element (" ZN ") breaks the column rule, but
        # here the strict reading is not an element at all, so take the pair.
        return letters[:2]
    return letters[:1]


def _pdb_cache_token(pdb_path: str) -> tuple[str, int, int]:
    """Return a cache token that changes when a PDB file changes.

    Parameters
    ----------
    pdb_path : str
        Path to a PDB file.

    Returns
    -------
    tuple
        Absolute path, modification time in ns, and file size.
    """
    path = os.path.abspath(pdb_path)
    stat = os.stat(path)
    return path, stat.st_mtime_ns, stat.st_size


@lru_cache(maxsize=32)
def _load_pdb_records_cached(
    pdb_path: str,
    mtime_ns: int,
    size: int,
) -> tuple[tuple[str, int, str, float, float, float, float], ...]:
    """Load ATOM/HETATM records from a PDB file.

    Parameters
    ----------
    pdb_path : str
        Absolute path to a PDB file.
    mtime_ns : int
        File modification timestamp used as part of the cache key.
    size : int
        File size used as part of the cache key.

    Returns
    -------
    tuple
        Records containing chain, residue number, atom name, xyz, and vdW radius.
    """
    del mtime_ns, size
    rows = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                resseq = int(line[22:26].strip())
            except ValueError:
                continue
            chain = line[21].strip()
            atom_name = line[12:16].strip()
            element = _element_symbol_from_pdb_line(line)
            atomic_number = _ELEMENT_NUMBERS.get(element, 0)
            rows.append((
                chain,
                resseq,
                atom_name,
                xyz[0],
                xyz[1],
                xyz[2],
                VDW_RADII.get(atomic_number, _DEFAULT_VDW),
            ))
    if not rows:
        raise ValueError(f"No ATOM/HETATM coordinates found in '{pdb_path}'")
    return tuple(rows)


def _cached_pdb_records(pdb_path: str) -> tuple[tuple[str, int, str, float, float, float, float], ...]:
    """Return cached PDB records for a path.

    Parameters
    ----------
    pdb_path : str
        Path to a PDB file.

    Returns
    -------
    tuple
        Cached ATOM/HETATM records.
    """
    return _load_pdb_records_cached(*_pdb_cache_token(pdb_path))


def _load_pdb_xyzr_direct(pdb_path: str) -> np.ndarray:
    """Load PDB ATOM/HETATM coordinates and vdW radii without IMP.

    Parameters
    ----------
    pdb_path : str
        Path to a PDB file.

    Returns
    -------
    numpy.ndarray
        ``(N, 4)`` array with ``x, y, z, vdw_radius`` columns.
    """
    records = _cached_pdb_records(pdb_path)
    return np.asarray(
        [(x, y, z, radius) for _, _, _, x, y, z, radius in records],
        dtype=np.float64,
    )


def load_structure_with_vdw(pdb_path: str) -> np.ndarray:
    """Load a PDB and return (N, 4) array: x, y, z, vdw_radius.

    Parses PDB records directly to avoid IMP/CHARMM warnings for unsupported
    HETATM residues, then falls back to IMP.atom if direct parsing fails.
    """
    try:
        return _load_pdb_xyzr_direct(pdb_path)
    except Exception:
        pass

    coords, particles, _model, _hier = io.load_structure_with_particles(pdb_path)
    vdw = np.full(coords.shape[0], _DEFAULT_VDW, dtype=np.float64)
    for i, p in enumerate(particles):
        try:
            at = IMP.atom.Atom(p)
            elem = at.get_element()
            vdw[i] = VDW_RADII.get(elem, _DEFAULT_VDW)
        except Exception:
            pass
    return np.column_stack([coords, vdw])


def _find_attachment_point(
    atoms: np.ndarray,
    chain: str,
    resseq: int,
    atom_name: str,
    pdb_path: Optional[str] = None,
) -> Optional[np.ndarray]:
    """Find the coordinates of an attachment atom.

    If pdb_path is provided, the atom is resolved by identity — chain, residue
    number and atom name — and a miss stays a miss: an unresolvable site
    returns ``None`` rather than a positional guess, because a dye attached to
    an unrelated atom yields a plausible and entirely wrong accessible volume.
    Without a PDB file the atoms array carries no identity at all, so the
    residue sequence number is used as a proxy index into it.

    Parameters
    ----------
    atoms : (N, 4) ndarray
        The atoms array (x, y, z, vdw_radius).
    chain : str
        The chain identifier.
    resseq : int
        The residue sequence number.
    atom_name : str
        The attachment atom name (e.g., 'CA', 'CB').
    pdb_path : str, optional
        Path to the PDB file for exact matching.

    Returns
    -------
    ndarray or None
        The (3,) coordinates of the attachment atom, or None if not found.
    """
    if pdb_path and os.path.exists(pdb_path):
        try:
            records = _cached_pdb_records(pdb_path)
        except Exception:
            logger.warning("Could not read atom records from %s", pdb_path, exc_info=True)
            return None
        for line_chain, line_resseq, line_atom_name, x, y, z, _ in records:
            if line_resseq == resseq and line_atom_name == atom_name:
                if not chain or line_chain == chain:
                    return np.array([x, y, z], dtype=np.float64)
        logger.warning(
            "Attachment atom '%s:%s:%s' does not exist in %s",
            chain, resseq, atom_name, pdb_path,
        )
        return None

    return atoms[resseq - 1, :3] if resseq > 0 and resseq <= atoms.shape[0] else None


def _strip_residue_atoms(
    atoms: np.ndarray,
    chain: str,
    resseq: int,
    pdb_path: Optional[str] = None,
) -> np.ndarray:
    """Remove the attachment residue's atoms from the coordinate array.

    Parameters
    ----------
    atoms : (N, 4) ndarray
        The atoms array (x, y, z, vdw_radius).
    chain : str
        The chain identifier of the residue to exclude.
    resseq : int
        The residue sequence number of the residue to exclude.
    pdb_path : str, optional
        Path to the PDB file for exact matching of residue atoms.

    Returns
    -------
    clean_atoms : (M, 4) ndarray
        The coordinate array with residue atoms removed.
    """
    if not pdb_path or not os.path.exists(pdb_path):
        return atoms

    try:
        records = _cached_pdb_records(pdb_path)
    except Exception:
        return atoms

    if len(records) == atoms.shape[0]:
        keep_mask = np.asarray(
            [
                not (line_resseq == resseq and (not chain or line_chain == chain))
                for line_chain, line_resseq, _, _, _, _, _ in records
            ],
            dtype=bool,
        )
        return atoms[keep_mask]

    coords_to_exclude = [
        [x, y, z]
        for line_chain, line_resseq, _, x, y, z, _ in records
        if line_resseq == resseq and (not chain or line_chain == chain)
    ]
    if not coords_to_exclude:
        return atoms

    coords_to_exclude_arr = np.array(coords_to_exclude, dtype=np.float64)
    distances = np.linalg.norm(
        atoms[:, None, :3] - coords_to_exclude_arr[None, :, :],
        axis=2,
    )
    return atoms[~np.any(distances < 0.01, axis=1)]


# ---------------------------------------------------------------------------
# The FPS strip: obstacles removed around the attachment site
# ---------------------------------------------------------------------------

#: The atoms of a peptide backbone. The default strip keeps these plus the
#: attachment atom and removes the rest of the attachment residue's side
#: chain -- the FPS convention, which every fps.json ``allowed_sphere_radius``
#: is calibrated against (FPS computes on a structure with that residue
#: stripped, so small declared clearances are meaningful, not walled in).
_BACKBONE_ATOM_NAMES = ("N", "CA", "C", "O")

#: Stripped-PDB cache, keyed like the record cache plus the site and mask.
#: One file per distinct (structure, site, mask) per process.
_STRIPPED_PDB_CACHE: Dict[tuple, str] = {}


def default_strip_mask(chain: str, resseq: int, atom_name: str) -> str:
    """The default strip for an attachment site, as an fps ``strip_mask``.

    The attachment residue's side chain minus the attachment atom; the
    backbone stays. The spelling is the PyMOL selection dialect fps
    documents use, so a declared mask and the default are the same
    vocabulary and round-trip through the same evaluation.

    Parameters
    ----------
    chain : str
        Chain identifier; empty matches any chain.
    resseq : int
        Residue sequence number of the attachment site.
    atom_name : str
        Name of the attachment atom (kept, never stripped).

    Returns
    -------
    str
        The mask expression, e.g.
        ``chain A and resid 36 and not name N+CA+C+O+CB``.
    """
    keep = "+".join(dict.fromkeys(_BACKBONE_ATOM_NAMES + (atom_name,)))
    parts = []
    if chain:
        parts.append(f"chain {chain}")
    parts.append(f"resid {int(resseq)}")
    parts.append(f"not name {keep}")
    return " and ".join(parts)


def _parse_strip_mask(mask: str) -> tuple:
    """Parse an fps ``strip_mask`` into a predicate over PDB records.

    The accepted grammar is the dialect fps documents actually carry: an
    ``and``-chain of ``chain <id>``, ``resid <n>`` (``resi`` accepted), and
    one name term, ``name A+B+...`` or ``not name A+B+...``. ``+`` is the
    list separator (PyMOL's own; space-separated lists are a parse error in
    PyMOL and here). Anything else -- ``or``, parentheses, other keywords --
    is refused loudly: a mask this function cannot read must not be
    silently ignored or approximated.

    Parameters
    ----------
    mask : str
        The strip mask expression.

    Returns
    -------
    tuple
        ``(chain, resids, names, negate)`` -- the chain id (or ``None`` for
        any chain), the residue numbers (empty for any), the name set (or
        ``None`` for no name term), and whether the name term was negated.
    """
    import re

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
        m = re.fullmatch(
            r"(not\s+)?name\s+(\S+(?:\+\S+)*)", term, flags=re.IGNORECASE
        )
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
            "structure to compute_av"
        )
    return chain, resids, names, negate


def _stripped_pdb_for(
    pdb_path: str,
    chain: str,
    resseq: int,
    atom_name: str,
    strip_mask: Optional[str] = None,
) -> str:
    """A copy of *pdb_path* with the atoms selected by *strip_mask* removed.

    The strip is the FPS convention: fps.json positions are calibrated for a
    structure whose attachment residue does not wall in its own dye, so the
    default removes the attachment residue's side chain minus the attachment
    atom (the backbone stays). A declared mask in the fps dialect is
    honoured as given; one outside the dialect raises rather than being
    ignored.

    The attachment atom itself is always kept, whatever the mask selects --
    the attachment is resolved by ``(chain, residue, atom name)`` from this
    file.

    Parameters
    ----------
    pdb_path : str
        Structure to strip.
    chain : str
        Attachment chain identifier (empty matches any chain).
    resseq : int
        Attachment residue sequence number.
    atom_name : str
        Attachment atom name.
    strip_mask : str, optional
        Declared fps ``strip_mask``; empty means the default strip.

    Returns
    -------
    str
        Path of the stripped copy (the original path on I/O failure -- a
        cloud computed against the unstripped structure beats no cloud).
    """
    mask = (strip_mask or "").strip()
    key = (*_pdb_cache_token(pdb_path), chain, int(resseq), atom_name, mask)
    cached = _STRIPPED_PDB_CACHE.get(key)
    if cached is not None and os.path.exists(cached):
        return cached
    try:
        with open(pdb_path) as source:
            lines = source.readlines()
    except OSError:
        return str(pdb_path)

    if mask:
        sel_chain, sel_resids, sel_names, negate = _parse_strip_mask(mask)
    else:
        sel_chain, sel_resids = (chain or None), (int(resseq),)
        sel_names = frozenset(
            n.upper() for n in _BACKBONE_ATOM_NAMES + (atom_name,)
        )
        negate = True

    kept = []
    for line in lines:
        stripped = False
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                line_resseq = int(line[22:26].strip())
            except ValueError:
                line_resseq = None
            line_chain = line[21].strip()
            line_name = line[12:16].strip().upper()
            if line_resseq is not None and line_name:
                selected = (
                    (sel_chain is None or line_chain == sel_chain)
                    and (not sel_resids or line_resseq in sel_resids)
                    and (
                        sel_names is None
                        or (line_name not in sel_names)
                        == negate
                    )
                )
                is_attachment = (
                    line_name == atom_name.strip().upper()
                    and line_resseq == int(resseq)
                    and (not chain or line_chain == chain)
                )
                stripped = selected and not is_attachment
        if not stripped:
            kept.append(line)

    try:
        import tempfile

        with tempfile.NamedTemporaryFile(
            "w", suffix=".pdb", prefix="bff_av_", delete=False
        ) as handle:
            handle.writelines(kept)
            result = handle.name
    except OSError:
        return str(pdb_path)
    _STRIPPED_PDB_CACHE[key] = result
    return result
