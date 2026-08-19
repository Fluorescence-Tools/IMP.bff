"""A minimal reader for CHARMM/NAMD DCD trajectories.

The rotamer libraries ship as PDB + DCD pairs, and reading them used to pull in
MDAnalysis. IMP.bff carries no dependency beyond what IMP itself brings, so the
format is read here instead: DCD is a small, fixed, well-documented layout, and
the subset the libraries use is a header, a title block, an atom count, and one
Fortran record per coordinate axis per frame.

Only what the libraries need is implemented — fixed atom counts, no velocity
blocks, no 4-dimensional trajectories. Anything outside that raises rather than
guessing, because a trajectory silently read as the wrong shape is worse than
one that refuses.

Verified frame-for-frame against MDAnalysis on the bundled library; the parity
test skips when MDAnalysis is absent (``test/test_dcd_reader.py``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple
import RMF
import ast
import json
import math
import os
import re
import struct

import numpy as np

import IMP
import IMP.atom
import IMP.bff.io.fps as fps_schema
import IMP.rmf

__all__ = [
    'compute_rmsd',
    'count_frames',
    'infer_bonds',
    'load_structure',
    'load_structure_with_particles',
    'parse_conect_bonds',
    'parse_pdb_atoms',
    'read_dcd',
    'read_dcd_header',
    'read_score_series',
    'write_mol2',
    'write_pdb',
    'write_rmf',
]

# --------------------------------------------------------------------------
# dcd
# --------------------------------------------------------------------------
"""A minimal reader for CHARMM/NAMD DCD trajectories.

The rotamer libraries ship as PDB + DCD pairs, and reading them used to pull in
MDAnalysis. IMP.bff carries no dependency beyond what IMP itself brings, so the
format is read here instead: DCD is a small, fixed, well-documented layout, and
the subset the libraries use is a header, a title block, an atom count, and one
Fortran record per coordinate axis per frame.

Only what the libraries need is implemented — fixed atom counts, no velocity
blocks, no 4-dimensional trajectories. Anything outside that raises rather than
guessing, because a trajectory silently read as the wrong shape is worse than
one that refuses.

Verified frame-for-frame against MDAnalysis on the bundled library; the parity
test skips when MDAnalysis is absent (``test/test_dcd_reader.py``).
"""

_MAGIC = b"CORD"
#: 4-byte record marker + 84-byte header payload.
_HEADER_RECORD = 84


class DCDFormatError(ValueError):
    """Raised when a file is not a DCD this reader is prepared to handle."""


def _endianness(raw: bytes) -> str:
    """Return the struct prefix matching the file's byte order.

    DCD carries no magic number for endianness; the convention is to read the
    leading Fortran record length and see which order makes it the expected 84.
    """
    for prefix in ("<", ">"):
        if struct.unpack(prefix + "i", raw[0:4])[0] == _HEADER_RECORD:
            return prefix
    raise DCDFormatError("leading record is not 84 bytes in either byte order")


def read_dcd_header(path: str | Path) -> dict:
    """Read a DCD header without loading coordinates.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the ``.dcd`` file.

    Returns
    -------
    dict
        ``n_frames``, ``n_atoms``, ``first_step``, ``step_stride``,
        ``time_step``, ``has_unit_cell``, ``charmm_version`` and the byte
        ``offset`` at which frame data begins.
    """
    path = Path(path)
    # The header, title block and atom count sit in the first few hundred
    # bytes; reading the whole trajectory to parse them made a header-only
    # call as expensive as a full load, which is what pushed the test suite
    # past its time budget. 64 KiB is far more than any title block needs.
    with path.open("rb") as fh:
        raw = fh.read(65536)
    return _parse_header(raw, path)


def _parse_header(raw: bytes, path: Path) -> dict:
    if len(raw) < 4 + _HEADER_RECORD + 8:
        raise DCDFormatError(f"{path} is too short to be a DCD file")
    end = _endianness(raw)
    if raw[4:8] != _MAGIC:
        raise DCDFormatError(f"{path} does not start with the CORD magic")

    ints = struct.unpack(end + "9i", raw[8:44])
    n_frames, first_step, step_stride, _last_step = ints[0], ints[1], ints[2], ints[3]
    # The float at offset 44 is the timestep; CHARMM writes it single precision.
    (time_step,) = struct.unpack(end + "f", raw[44:48])
    flags = struct.unpack(end + "9i", raw[48:84])
    has_unit_cell = bool(flags[0])
    four_dimensions = bool(flags[3])
    (charmm_version,) = struct.unpack(end + "i", raw[84:88])
    if four_dimensions:
        raise DCDFormatError(f"{path} is a 4-dimensional trajectory, which is not supported")

    pos = 4 + _HEADER_RECORD + 4  # closing marker of the header record

    # Title block: one Fortran record holding a count and that many 80-char lines.
    (title_len,) = struct.unpack(end + "i", raw[pos:pos + 4])
    pos += 4 + title_len + 4

    # Atom-count block.
    (natom_len,) = struct.unpack(end + "i", raw[pos:pos + 4])
    if natom_len != 4:
        raise DCDFormatError(f"{path} has an unexpected atom-count record ({natom_len} bytes)")
    (n_atoms,) = struct.unpack(end + "i", raw[pos + 4:pos + 8])
    pos += 4 + natom_len + 4

    return {
        "n_frames": n_frames,
        "n_atoms": n_atoms,
        "first_step": first_step,
        "step_stride": step_stride,
        "time_step": time_step,
        "has_unit_cell": has_unit_cell,
        "charmm_version": charmm_version,
        "endianness": end,
        "offset": pos,
    }


def read_dcd(path: str | Path, max_frames: int | None = None) -> np.ndarray:
    """Read coordinates from a DCD trajectory.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the ``.dcd`` file.
    max_frames : int, optional
        Stop after this many frames. ``None`` reads all of them.

    Returns
    -------
    numpy.ndarray
        ``(n_frames, n_atoms, 3)`` float64 array of coordinates, in the file's
        own units (Angstrom for the bundled libraries).
    """
    path = Path(path)
    with path.open("rb") as fh:
        raw = fh.read()
    head = _parse_header(raw, path)
    end, n_atoms, pos = head["endianness"], head["n_atoms"], head["offset"]

    n_frames = head["n_frames"]
    if max_frames is not None:
        n_frames = min(n_frames, max_frames)

    axis_bytes = 4 * n_atoms
    dtype = np.dtype(np.float32).newbyteorder(end)
    frames = np.empty((n_frames, n_atoms, 3), dtype=np.float64)

    for frame in range(n_frames):
        if head["has_unit_cell"]:
            # A 48-byte double record ahead of the coordinates.
            (cell_len,) = struct.unpack(end + "i", raw[pos:pos + 4])
            pos += 4 + cell_len + 4
        for axis in range(3):
            (rec_len,) = struct.unpack(end + "i", raw[pos:pos + 4])
            if rec_len != axis_bytes:
                raise DCDFormatError(
                    f"{path}: frame {frame} axis {axis} record is {rec_len} bytes, "
                    f"expected {axis_bytes} for {n_atoms} atoms"
                )
            start = pos + 4
            frames[frame, :, axis] = np.frombuffer(
                raw, dtype=dtype, count=n_atoms, offset=start
            )
            pos = start + rec_len + 4

    return frames


# --------------------------------------------------------------------------
# mol2
# --------------------------------------------------------------------------
"""Reading a PDB into mol2, and writing mol2 out.

The library could **parse** mol2 (:func:`IMP.bff.parse_dye_mol2`,
``io/cif.py``) but not write one, and mol2 is the format a dye topology has to
be in before a force field can be built for it. The writer lived in a CLI
driver, ``cgdye/scripts/pdb_to_mol2.py``, which is why it was invisible: a
capability reachable only by running a script is a capability nobody finds.

Lifted into the library in the PRD-113 cleanup (2026-08-18) when that script
moved to ``junk/``. The driver is gone; this is what it was for.

It landed in ``cgdye/io/`` first and moved here when ``cgdye`` left the domain
layout -- a capability rescued from a script should not then be buried in a
legacy package.

Bonds come from ``CONECT`` records when the PDB has them and from geometry when
it does not -- :func:`infer_bonds` uses a covalent-radius cutoff, which is a
guess and is why it is not the default.
"""

_ELEM_TO_TRIPOS = {
    "C": "C.3",
    "N": "N.3",
    "O": "O.3",
    "S": "S.3",
    "H": "H",
    "P": "P.3",
    "F": "F",
    "Cl": "Cl",
    "Br": "Br",
    "I": "I",
}


def parse_pdb_atoms(path: Path) -> dict:
    atoms = {}
    with open(path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            serial = int(line[6:11])
            name = line[12:16].strip()
            resname = line[17:20].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            elem_raw = line[76:78].strip() if len(line) > 76 else ""
            if not elem_raw:
                m = re.match(r"[A-Za-z]", name)
                elem_raw = m.group(0) if m else "C"
            atoms[serial] = dict(
                serial=serial,
                name=name,
                resname=resname,
                x=x,
                y=y,
                z=z,
                elem=elem_raw[0].upper(),
            )
    return atoms


def parse_conect_bonds(path: Path) -> set:
    bonds = set()
    with open(path) as fh:
        for line in fh:
            if line[:6].strip() != "CONECT":
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            a = int(fields[1])
            for tok in fields[2:]:
                b = int(tok)
                if a != b:
                    bonds.add(tuple(sorted((a, b))))
    return bonds


def infer_bonds(atoms: dict, scale: float = 1.22) -> set:
    """Geometry-based bond inference using covalent radii."""
    cov_r = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07}
    serials = sorted(atoms)
    bonds = set()
    for i, a in enumerate(serials):
        ea = atoms[a]["elem"]
        ra = cov_r.get(ea, 0.77)
        for b in serials[i + 1 :]:
            eb = atoms[b]["elem"]
            rb = cov_r.get(eb, 0.77)
            dx = atoms[a]["x"] - atoms[b]["x"]
            dy = atoms[a]["y"] - atoms[b]["y"]
            dz = atoms[a]["z"] - atoms[b]["z"]
            d = math.sqrt(dx * dx + dy * dy + dz * dz)
            if d <= scale * (ra + rb):
                bonds.add(tuple(sorted((a, b))))
    return bonds


def write_mol2(path: Path, atoms: dict, bonds: set, mol_name: str = "MOL") -> None:
    serials = sorted(atoms)
    idx = {s: i + 1 for i, s in enumerate(serials)}  # 1-based atom indices

    lines = []
    lines.append("@<TRIPOS>MOLECULE")
    lines.append(mol_name)
    lines.append(f"{len(atoms)} {len(bonds)} 0 0 0")
    lines.append("SMALL")
    lines.append("NO_CHARGES")
    lines.append("")

    lines.append("@<TRIPOS>ATOM")
    for s in serials:
        a = atoms[s]
        el = a["elem"]
        atype = _ELEM_TO_TRIPOS.get(el, el)
        lines.append(
            f"{idx[s]:6d} {a['name']:<8s}"
            f" {a['x']:10.4f} {a['y']:10.4f} {a['z']:10.4f}"
            f" {atype:<8s} 1 {a['resname']} 0.0000"
        )

    lines.append("@<TRIPOS>BOND")
    for bi, (a, b) in enumerate(sorted(bonds), 1):
        lines.append(f"{bi:6d} {idx[a]:6d} {idx[b]:6d} 1")

    lines.append("")
    path.write_text("\n".join(lines))


# --------------------------------------------------------------------------
# pmi_stat
# --------------------------------------------------------------------------
"""Read PMI ``stat.*.out`` files for live docking progress and score curves.

The PMI ``ReplicaExchange`` macro appends one line per frame to ``stat.0.out``:
a first header line (string keys describing each column) followed by one
``repr(dict)`` per frame whose integer keys map to the header columns. We only
need two columns — the frame index and the total score — to drive a progress
bar and a score-vs-frame convergence plot, so the parsing here is deliberately
tolerant: unreadable lines (including the giant header) are skipped.
"""

#: Integer column keys PMI uses in the per-frame stat dicts.
_TOTAL_SCORE_KEY = 1
_NFRAME_KEY = 4

# Field names in the header line, in case PMI ever renumbers the columns.
_TOTAL_SCORE_NAME = "Total_Score"
_NFRAME_NAME = "MonteCarlo_Nframe"


def count_frames(stat_path) -> int:
    """Return the number of completed frames recorded in ``stat_path``.

    One header line plus one line per frame, so ``data_lines = total - 1``.
    Returns 0 when the file is missing or unreadable.
    """
    try:
        with open(stat_path, "r", encoding="utf-8", errors="ignore") as fh:
            n = sum(1 for _ in fh)
    except OSError:
        return 0
    return max(0, n - 1)


def _resolve_keys(header_line: str) -> Tuple[int, int]:
    """Map the score/frame columns from the header line; fall back to defaults."""
    score_key, frame_key = _TOTAL_SCORE_KEY, _NFRAME_KEY
    for col, name in _safe_dict(header_line).items():
        if name == _TOTAL_SCORE_NAME and isinstance(col, int):
            score_key = col
        elif name == _NFRAME_NAME and isinstance(col, int):
            frame_key = col
    return score_key, frame_key


def _safe_dict(line: str) -> dict:
    """``ast.literal_eval`` a stat line into a dict, or ``{}`` if it cannot."""
    line = line.strip()
    if not line.startswith("{"):
        return {}
    try:
        obj = ast.literal_eval(line)
    except (ValueError, SyntaxError):  # e.g. header's environ(...) call
        return {}
    return obj if isinstance(obj, dict) else {}


def read_score_series(stat_path) -> Tuple[List[float], List[float]]:
    """Return ``(frames, scores)`` parsed from a PMI stat file.

    Parameters
    ----------
    stat_path : str
        Path to a ``stat.0.out`` file (need not be complete; partial files from
        a running job are fine).

    Returns
    -------
    (list of float, list of float)
        Frame indices and the matching total scores, in file order. Empty when
        the file is missing or has no parseable data lines.
    """
    if not os.path.exists(stat_path):
        return [], []
    frames: List[float] = []
    scores: List[float] = []
    score_key, frame_key = _TOTAL_SCORE_KEY, _NFRAME_KEY
    with open(stat_path, "r", encoding="utf-8", errors="ignore") as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if line.startswith("{"):  # PMI Monte-Carlo stat dict
                if i == 0:
                    score_key, frame_key = _resolve_keys(line)
                    continue
                d = _safe_dict(line)
                if score_key not in d:
                    continue
                try:
                    scores.append(float(d[score_key]))
                    frames.append(float(d.get(frame_key, len(frames))))
                except (TypeError, ValueError):
                    continue
            else:  # plain "frame,score" convergence CSV (minimisation)
                parts = line.split(",")
                if len(parts) < 2:
                    continue
                try:
                    frames.append(float(parts[0]))
                    scores.append(float(parts[1]))
                except ValueError:
                    continue  # header row "frame,score"
    return frames, scores


# --------------------------------------------------------------------------
# rotamer_rmf
# --------------------------------------------------------------------------
"""Read/write rotamer libraries in RMF format.

Each rotamer is stored as a frame in the RMF file.
Weights are stored as frame-level attributes (Score category).
"""

def write_rotamer_library_rmf(path, library):
    """Write a rotamer library to an RMF file.
    
    Args:
        path: Path to the .rmf3 file
        library: Dict with 'weight' (list), 'atom_names' (list), 
                 'coords' (dict {id: array})
    """
    if not path.endswith(".rmf3"):
        path += ".rmf3"
        
    model = IMP.Model()
    root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model, "rotamers"))
    
    # Create a dummy chain and residue to hold the atoms
    chain = IMP.atom.Chain.setup_particle(IMP.Particle(model, "A"), "A")
    root.add_child(chain)
    res = IMP.atom.Residue.setup_particle(IMP.Particle(model, "DYE"), IMP.atom.ResidueType("DYE"), 1)
    chain.add_child(res)

    # Create atoms once
    atom_names = library["atom_names"]
    particles = []
    for name in atom_names:
        p = IMP.Particle(model, name)
        # Use Mass instead of Atom to preserve name in RMF
        IMP.atom.Mass.setup_particle(p, 1.0)
        IMP.core.XYZR.setup_particle(p)
        IMP.core.XYZR(p).set_radius(1.0)
        root_hier = IMP.atom.Hierarchy.setup_particle(p)
        res.add_child(root_hier)
        particles.append(p)
        
    fh = RMF.create_rmf_file(path)
    IMP.rmf.add_hierarchies(fh, [root])
    
    # Setup score category for weights
    score_cat = fh.get_category("score")
    weight_key = fh.get_key(score_cat, "weight", RMF.FloatTag())
    
    # Setup kinetic category for transitions
    kinetic_cat = fh.get_category("kinetic")
    trans_key = fh.get_key(kinetic_cat, "transition_matrix", RMF.IntsTag())
    
    root_node = fh.get_root_node()
    
    # Save transition matrix as a static value on the root node
    if "transitions" in library:
        # Flatten 2D matrix to 1D for RMF IntsTag
        flat_trans = np.array(library["transitions"]).flatten().tolist()
        root_node.set_static_value(trans_key, flat_trans)
    
    ids = sorted(library["coords"].keys())
    for i, rid in enumerate(ids):
        coords = library["coords"][rid]
        for p, c in zip(particles, coords):
            IMP.core.XYZ(p).set_coordinates(IMP.algebra.Vector3D(c[0], c[1], c[2]))
        
        # Save weight to frame
        IMP.rmf.save_frame(fh, str(rid))
        root_node.set_frame_value(weight_key, float(library["weight"][i]))
        
    del fh # Close file


def read_rotamer_library_rmf(path):
    """Read a rotamer library from an RMF file.
    
    Returns dict with 'weight', 'atom_names', 'coords', 'transitions'.
    """
    if not os.path.exists(path):
        if os.path.exists(path + ".rmf3"):
            path += ".rmf3"
        else:
            raise FileNotFoundError(f"RMF library not found: {path}")
            
    model = IMP.Model()
    fh = RMF.open_rmf_file_read_only(path)
    hs = IMP.rmf.create_hierarchies(fh, model)
    root = hs[0]
    
    # Use get_leaves to get all atoms in the hierarchy
    atoms = IMP.atom.get_leaves(root)
    # Filter to only keep particles that are likely atoms (have a name and XYZ)
    atom_names = [a.get_name() for a in atoms]
    
    score_cat = fh.get_category("score")
    weight_key = fh.get_key(score_cat, "weight", RMF.FloatTag())
    
    kinetic_cat = fh.get_category("kinetic")
    trans_key = fh.get_key(kinetic_cat, "transition_matrix", RMF.IntsTag())
    
    root_node = fh.get_root_node()
    
    n_frames = fh.get_number_of_frames()
    weights = []
    coords_dict = {}
    
    for i in range(n_frames):
        IMP.rmf.load_frame(fh, RMF.FrameID(i))
        weights.append(root_node.get_frame_value(weight_key))
        
        frame_coords = np.zeros((len(atoms), 3))
        for j, a in enumerate(atoms):
            c = IMP.core.XYZ(a).get_coordinates()
            frame_coords[j] = [c[0], c[1], c[2]]
        coords_dict[i + 1] = frame_coords
        
    # Read transition matrix if present
    transitions = None
    if root_node.get_has_value(trans_key):
        flat_trans = root_node.get_static_value(trans_key)
        n_clusters = int(np.sqrt(len(flat_trans)))
        transitions = np.array(flat_trans).reshape((n_clusters, n_clusters)).tolist()
        
    return {
        "id": list(range(1, n_frames + 1)),
        "weight": weights,
        "atom_names": atom_names,
        "coords": coords_dict,
        "transitions": transitions
    }


# --------------------------------------------------------------------------
# structure
# --------------------------------------------------------------------------
"""Structure readers and writers: PDB, RMF, and coordinate comparison.

Nothing here is about fluorescence. It sits in :mod:`IMP.bff.io` because it is
serialisation, and it is a separate module from :mod:`IMP.bff.io.fps` because a
labelling file and a coordinate file are different formats that change for
different reasons -- they shared one module until PRD-113 stage 7, which is how
``fret/io.py`` came to be the thing that writes PDBs.
"""

def write_pdb(
    atoms: np.ndarray,
    path: str | os.PathLike,
    chain: str = "A",
    res_name: str = "ALA",
    transform: Optional[np.ndarray] = None,
    model_index: int = 0,
) -> None:
    """Write (N,3) or (N,4) coordinates as a multi-model PDB.

    Parameters
    ----------
    atoms : (N, 3) or (N, 4) ndarray
        If (N,4), column 4 is ignored (element hint).
    transform : (4,4), (3,3), or (3,) optional
        Homogeneous transform, rotation matrix, or translation.
    """
    if atoms.ndim != 2 or atoms.shape[1] < 3:
        raise ValueError(f"Expected (N,3) or (N,4) array, got {atoms.shape}")
    coords = atoms[:, :3].copy()
    if transform is not None:
        coords = _apply_transform(coords, transform)

    lines = [f"MODEL     {model_index:>4d}"]
    for i in range(coords.shape[0]):
        x, y, z = coords[i]
        name = "CA"
        name_fmt = f"{name:>4s}" if len(name) >= 4 else f" {name:<3s}"
        lines.append(
            f"ATOM  {i + 1:>5d} {name_fmt}{res_name:<3s} {chain}"
            f"{i + 1:>4d}    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
            f"{1.0:>6.2f}{0.0:>6.2f}"
        )
    lines.append("ENDMDL")
    lines.append("END")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _apply_transform(coords: np.ndarray, transform: np.ndarray) -> np.ndarray:
    """Apply a (4,4) homogeneous, (3,3) rotation, or (3,) translation."""
    t = np.asarray(transform, dtype=np.float64)
    if t.shape == (4, 4):
        return coords @ t[:3, :3].T + t[:3, 3]
    if t.shape == (3, 3):
        return coords @ t.T
    if t.shape == (3,):
        return coords + t
    raise ValueError(f"Unexpected transform shape {t.shape}")


def write_rmf(
    atoms: np.ndarray,
    path: str | os.PathLike,
    model_name: str = "structure",
    transform: Optional[np.ndarray] = None,
    metadata: Optional[Dict] = None,
    radius: float = 1.5,
) -> None:
    """Write ``(N, 3)`` coordinates to an RMF file via ``IMP.rmf``.

    Each coordinate becomes an ``XYZR`` ball under a single hierarchy named
    ``model_name``; ``metadata`` (when given) is stored as the RMF file
    description. This is a self-contained replacement for the application
    RMF writer the pre-move implementation borrowed.

    Parameters
    ----------
    atoms : (N, 3) or (N, 4) ndarray
        Coordinates; a 4th column is ignored.
    path : str or PathLike
        Output ``.rmf``/``.rmf3`` path.
    model_name : str
        Name of the root hierarchy node.
    transform : (4,4), (3,3), or (3,) optional
        Applied to the coordinates before writing.
    metadata : dict, optional
        JSON-serialised into the RMF file description.
    radius : float
        Ball radius (Angstrom) given to every particle.
    """
    import IMP
    import IMP.atom
    import IMP.core
    import IMP.rmf
    import RMF

    coords = np.asarray(atoms[:, :3], dtype=np.float64).copy()
    if transform is not None:
        coords = _apply_transform(coords, transform)

    model = IMP.Model()
    root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model, str(model_name)))
    for i in range(coords.shape[0]):
        p = IMP.Particle(model, f"p{i}")
        IMP.core.XYZR.setup_particle(
            p, IMP.algebra.Sphere3D(IMP.algebra.Vector3D(*coords[i]), radius))
        IMP.atom.Mass.setup_particle(p, 1.0)
        root.add_child(IMP.atom.Hierarchy.setup_particle(p))

    fh = RMF.create_rmf_file(str(path))
    try:
        if metadata:
            fh.set_description(json.dumps(metadata, default=str))
        IMP.rmf.add_hierarchy(fh, root)
        IMP.rmf.save_frame(fh, "frame_0")
    finally:
        del fh


def load_structure(pdb_path: str | os.PathLike) -> np.ndarray:
    """Load a PDB and return (N, 3) coordinates.

    Requires ``IMP.atom``.
    """
    coords, _leaves, _model, _hier = load_structure_with_particles(pdb_path)
    return coords


def load_structure_with_particles(pdb_path: str | os.PathLike):
    """Load a PDB and return (coords, leaves, model, hierarchy).

    Requires ``IMP.atom``.
    """
    import IMP
    import IMP.atom
    import IMP.core

    model = IMP.Model()
    hier = IMP.atom.read_pdb(str(pdb_path), model, IMP.atom.NonWaterPDBSelector())
    leaves = IMP.atom.get_leaves(hier)
    n = len(leaves)
    coords = np.empty((n, 3), dtype=np.float64)
    for i, p in enumerate(leaves):
        xyz = IMP.core.XYZ(p)
        coords[i, 0] = xyz.get_x()
        coords[i, 1] = xyz.get_y()
        coords[i, 2] = xyz.get_z()
    return coords, leaves, model, hier


def compute_rmsd(
    coords_a: np.ndarray,
    coords_b: np.ndarray,
    selection_mask: Optional[np.ndarray] = None,
    superpose: bool = False,
) -> float:
    """Compute RMSD between two coordinate arrays, optionally with superposition.

    Parameters
    ----------
    coords_a : (N, 3) ndarray
        First coordinate array.
    coords_b : (N, 3) ndarray
        Second coordinate array.
    selection_mask : (N,) bool ndarray, optional
        If provided, only the selected atoms are used to calculate the RMSD
        (and the superposition transform, if enabled).
    superpose : bool
        If True, apply Kabsch rotation to align coords_a onto coords_b
        before computing RMSD.

    Returns
    -------
    rmsd : float
        The root-mean-square deviation.
    """
    if coords_a.shape != coords_b.shape:
        raise ValueError(f"Shape mismatch: {coords_a.shape} vs {coords_b.shape}")

    if selection_mask is not None:
        a_sel = coords_a[selection_mask]
        b_sel = coords_b[selection_mask]
    else:
        a_sel = coords_a
        b_sel = coords_b

    if len(a_sel) == 0:
        return 0.0

    if superpose:
        # Compute centroids
        centroid_a = a_sel.mean(axis=0)
        centroid_b = b_sel.mean(axis=0)

        # Center coordinates
        a_centered = a_sel - centroid_a
        b_centered = b_sel - centroid_b

        # SVD of covariance matrix
        c = a_centered.T @ b_centered
        v, s, wt = np.linalg.svd(c)

        # Det for reflection detection
        d = np.linalg.det(v @ wt)
        s_diag = np.array([1.0, 1.0, np.sign(d) if d != 0 else 1.0])
        r = v @ np.diag(s_diag) @ wt

        # Rotate entire coords_a
        coords_a_aligned = (coords_a - centroid_a) @ r + centroid_b
        if selection_mask is not None:
            diff = coords_a_aligned[selection_mask] - coords_b[selection_mask]
        else:
            diff = coords_a_aligned - coords_b
    else:
        diff = a_sel - b_sel

    return float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))
