"""Structure readers and writers: PDB, RMF, and coordinate comparison.

Nothing here is about fluorescence. It sits in :mod:`IMP.bff.io` because it is
serialisation, and it is a separate module from :mod:`IMP.bff.io.fps` because a
labelling file and a coordinate file are different formats that change for
different reasons -- they shared one module until PRD-113 stage 7, which is how
``fret/io.py`` came to be the thing that writes PDBs.
"""

from __future__ import annotations

import json
import os
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

from IMP.bff.io import fps_schema

__all__ = [
    "load_structure",
    "load_structure_with_particles",
    "write_pdb",
    "write_rmf",
    "compute_rmsd",
]


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
