#!/usr/bin/env python
"""Read/write rotamer libraries.

Coordinates stored as numpy .npy files (shape: n_rotamers x n_atoms x 3).
Weights stored as plain text (one per line). Atom names stored as text.
"""

import os

import numpy as np


def _base_path(rmf_or_npy_path):
    if rmf_or_npy_path.endswith(".rmf"):
        return rmf_or_npy_path[:-4]
    if rmf_or_npy_path.endswith(".npy"):
        return rmf_or_npy_path[:-4]
    return rmf_or_npy_path


def read_rotamer_library(path):
    """Read a rotamer library from numpy/text files.

    ``path`` can be the .npy, .rmf, or base path (without extension).
    Files expected: ``{base}_coords.npy``, ``{base}_weights.txt``,
    ``{base}_atoms.txt``.

    Returns dict with:
        - id: list of rotamer IDs (1-indexed ints)
        - weight: list of weights (float)
        - atom_names: list of atom name strings
        - coords: dict {rotamer_id: ndarray shape (n_atoms, 3)}
    """
    base = _base_path(path)
    coords = np.load(base + "_coords.npy")
    n_rotamers = coords.shape[0]

    weights_file = base + "_weights.txt"
    if os.path.exists(weights_file):
        with open(weights_file) as fh:
            weights = [float(line.strip()) for line in fh if line.strip()]
    else:
        weights = [1.0 / n_rotamers] * n_rotamers

    atoms_file = base + "_atoms.txt"
    if os.path.exists(atoms_file):
        with open(atoms_file) as fh:
            atom_names = [line.strip() for line in fh if line.strip()]
    else:
        atom_names = [f"AT{i}" for i in range(coords.shape[1])]

    data = {
        "id": list(range(1, n_rotamers + 1)),
        "weight": weights,
        "atom_names": atom_names,
        "coords": {i + 1: coords[i] for i in range(n_rotamers)},
    }
    return data


def write_rotamer_library(path, library):
    """Write a rotamer library to numpy/text files.

    ``library`` must have:
        - weight: list of floats
        - atom_names: list of atom name strings
        - coords: dict {rotamer_id: list/array of [x, y, z] per atom}
    """
    base = _base_path(path)
    atom_names = library["atom_names"]
    n_atoms = len(atom_names)
    ids = sorted(library["coords"].keys())

    coord_array = np.zeros((len(ids), n_atoms, 3))
    for i, rid in enumerate(ids):
        c = library["coords"][rid]
        arr = np.asarray(c)
        if arr.ndim == 1 and arr.shape[0] == 3:
            arr = arr.reshape(1, 3)
        coord_array[i, : arr.shape[0]] = arr[:n_atoms]

    np.save(base + "_coords.npy", coord_array)

    with open(base + "_weights.txt", "w") as fh:
        for w in library["weight"]:
            fh.write(f"{w}\n")

    with open(base + "_atoms.txt", "w") as fh:
        for name in atom_names:
            fh.write(f"{name}\n")


def normalize_weights(library):
    """Normalize rotamer weights to sum to 1.0 in-place."""
    total = sum(library["weight"])
    if total > 0:
        library["weight"] = [w / total for w in library["weight"]]
