"""Loading a rotamer library, and drawing a state from it.

The mechanics of the discrete-states representation: read the library off disk
(DCD frames plus weights), put one of its conformers onto a hierarchy, and draw
an index according to the weights.

Was ``cgdye/sampling/rotamer.py`` until the 2026-08-18 cleanup. It sat with the
force-field samplers because it uses the same file formats, but a rotamer
library is not something you *sample from a potential* -- it is a list of states
with weights already attached, which is a representation.
"""

from __future__ import annotations

import random
from pathlib import Path

import numpy as np
import IMP
import IMP.atom
import IMP.core


def load_rotamer_library_dcd(pdb_path, dcd_path, weights_path=None, max_frames=None):
    """Load a reference rotamer library from PDB+DCD (+ optional weights).

    Atom names come from the PDB through IMP.atom and coordinates from the DCD
    through the in-tree reader, so this needs nothing beyond IMP and numpy.
    """
    from IMP.bff.io.dcd import read_dcd

    model = IMP.Model()
    hierarchy = IMP.atom.read_pdb(str(pdb_path), model, IMP.atom.AllPDBSelector())
    atom_names = [
        IMP.atom.Atom(leaf).get_atom_type().get_string().strip()
        for leaf in IMP.atom.get_leaves(hierarchy)
    ]

    coords = read_dcd(dcd_path, max_frames=max_frames)

    if coords.shape[0] == 0:
        raise ValueError(f"No frames found in DCD: {dcd_path}")

    if weights_path is not None and Path(weights_path).exists():
        w = []
        with open(weights_path) as fh:
            for line in fh:
                txt = line.strip()
                if not txt:
                    continue
                w.append(float(txt))
        weights = np.asarray(w, dtype=float)
        if max_frames is not None:
            weights = weights[:max_frames]
        if len(weights) != coords.shape[0]:
            n = min(len(weights), coords.shape[0])
            coords = coords[:n]
            weights = weights[:n]
    else:
        weights = np.ones(coords.shape[0], dtype=float)

    total = float(weights.sum())
    if total <= 0:
        weights = np.ones_like(weights)
        total = float(weights.sum())
    weights = weights / total

    return {
        "coords": coords,
        "weights": weights,
        "atom_names": atom_names,
    }


def find_reference_rotamer_files(lib_dir, dye_name, cutoff=30):
    """Return (pdb, dcd, weights) paths for a reference dye+linker name."""
    lib_dir = Path(lib_dir)
    pdb = lib_dir / f"{dye_name}.pdb"
    dcd = lib_dir / f"{dye_name}_cutoff{cutoff}.dcd"
    weights = lib_dir / f"{dye_name}_cutoff{cutoff}_weights.txt"
    if not pdb.exists() or not dcd.exists():
        raise FileNotFoundError(f"Missing required reference files for {dye_name}")
    return pdb, dcd, weights if weights.exists() else None


def apply_rotamer_coordinates(dye_hier, coords):
    """Apply one rotamer coordinate set to an IMP dye hierarchy in-place."""
    atoms = list(IMP.atom.get_by_type(dye_hier, IMP.atom.ATOM_TYPE))
    arr = np.asarray(coords, dtype=float)
    if arr.shape[0] != len(atoms):
        raise ValueError(
            f"Atom count mismatch: coords={arr.shape[0]} hierarchy={len(atoms)}"
        )
    for a, c in zip(atoms, arr):
        IMP.core.XYZ(a).set_coordinates(
            IMP.algebra.Vector3D(float(c[0]), float(c[1]), float(c[2]))
        )


def sample_rotamer_index(weights, rng=None):
    """Sample a rotamer index from normalized weights."""
    if rng is None:
        rng = random
    idx = list(range(len(weights)))
    return rng.choices(idx, weights=weights, k=1)[0]
