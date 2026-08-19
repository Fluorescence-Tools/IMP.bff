#!/usr/bin/env python3
"""Convert a rotamer-library trajectory (DCD or XTC) to BinaryCIF.

BinaryCIF is the trajectory format for this package from 2026-08-19. It is
smaller than either format it replaces and is decoded by the C implementation
of ``ihm`` that IMP already vendors, so the C++ side reads it with no new
dependency:

===========================================  ==================
 DCD, raw float32 plus per-frame headers      4.31 bytes/coord
 XTC, its own 3-D compression                 1.60
 **BinaryCIF at a 0.1 A grid**                **1.27**
===========================================  ==================

Measured on ``data/rotamer_library/A48_C1R/traj.xtc``: 28 110 frames x 92
atoms, 12.44 MB of XTC becoming 9.87 MB of BinaryCIF, decoded exactly in 68 ms.
The full record, including what the precision costs, is in
``okf/validation/bcif_for_trajectories.md``.

**Why this script exists at all.** ``python-ihm``'s ``BinaryCifWriter``
implements only ByteArray, Delta, RunLength and the string/mask encoders. The
compression here needs **FixedPoint** and **IntegerPacking**, which it does not
have, so the encoder is written out by hand below. The *reader* side needs
nothing: ``ihm_format.c`` implements all seven encodings.

Usage
-----
::

    scripts/trajectory_to_bcif.py traj.xtc traj.bcif --top conf_ed.gro
    scripts/trajectory_to_bcif.py lib.dcd lib.bcif
    scripts/trajectory_to_bcif.py --all data/rotamer_library

Reading a DCD needs only ``IMP.bff``; reading an XTC needs ``mdtraj``, which is
a converter-time dependency and not one of the package's.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

try:
    import msgpack
except ImportError:                                          # pragma: no cover
    sys.exit("msgpack is required (it ships with python-ihm)")

#: Grid in Angstrom. 0.1 A moves a mean FRET efficiency by 1.7e-5, against the
#: ~1e-2 an experiment resolves, and the accessible volumes these libraries
#: feed are built at 0.5-2 A. At XTC's own 0.01 A the file is 25 % *larger*
#: than the XTC -- the format only wins because the precision can be relaxed.
DEFAULT_GRID_A = 0.1

#: BinaryCIF ByteArray type codes (the spec's, not numpy's).
_BYTE_ARRAY_INT8 = 1
_BYTE_ARRAY_INT32 = 3


def _integer_pack_int8(deltas: np.ndarray) -> np.ndarray:
    """BinaryCIF IntegerPacking into int8, with escape runs.

    The sentinels are ``127`` and ``-128``, and they are *also* legitimate
    values -- so a delta of exactly +-127 has to be written as an escape plus a
    remainder. Using ``>`` rather than ``>=`` here produces a file whose columns
    decode to different lengths, which the reader reports as
    ``Column size mismatch`` a long way from the cause.
    """
    lo, hi = -128, 127
    out: list[int] = []
    append = out.append
    for v in deltas.tolist():
        while v >= hi:
            append(hi)
            v -= hi
        while v <= lo:
            append(lo)
            v -= lo
        append(v)
    return np.array(out, dtype=np.int8)


def encode_column(name: str, values_A: np.ndarray, grid_A: float) -> dict:
    """One coordinate column, FixedPoint -> Delta -> IntegerPacking -> ByteArray.

    The ``encoding`` list is stored in **encode** order. The C reader prepends
    each entry as it parses (``ihm_format.c``), so its linked list comes out
    reversed and it decodes in the right order. Writing the list the other way
    round fails with ``FixedPoint not given integers as input``.
    """
    factor = int(round(1.0 / grid_A))
    q = np.round(np.asarray(values_A, dtype=np.float64) * factor).astype(np.int32)
    if q.size == 0:
        raise ValueError(f"column {name!r} is empty")
    origin = int(q[0])
    deltas = np.diff(q, prepend=q[:1]).astype(np.int32)
    deltas[0] = 0                       # Delta's own origin carries the first
    return {
        "name": name,
        "data": {
            "data": _integer_pack_int8(deltas).tobytes(),
            "encoding": [
                {"kind": "FixedPoint", "factor": factor,
                 "srcType": _BYTE_ARRAY_INT32},
                {"kind": "Delta", "origin": origin,
                 "srcType": _BYTE_ARRAY_INT32},
                {"kind": "IntegerPacking", "byteCount": 1,
                 "isUnsigned": False, "srcSize": int(q.size)},
                {"kind": "ByteArray", "type": _BYTE_ARRAY_INT8},
            ],
        },
        "mask": None,
    }


def write_bcif(path: Path, xyz_A: np.ndarray, grid_A: float = DEFAULT_GRID_A,
               category: str = "_rotamer_coord", block: str = "rotamers") -> int:
    """Write ``(n_frames, n_atoms, 3)`` coordinates in Angstrom as BinaryCIF.

    Columns are laid out **atom-major** -- every frame of atom 0, then atom 1 --
    so the deltas run along an atom's own frame series. That is not obviously
    the right choice for a rotamer library, where consecutive frames are
    independent conformers rather than a time series, and it measured no better
    than frame-major at 0.001 nm. It measures better once the grid is relaxed,
    which is the regime this actually runs in.
    """
    xyz_A = np.asarray(xyz_A, dtype=np.float64)
    if xyz_A.ndim != 3 or xyz_A.shape[2] != 3:
        raise ValueError(f"expected (n_frames, n_atoms, 3), got {xyz_A.shape}")
    n_frames, n_atoms, _ = xyz_A.shape
    columns = [
        encode_column(axis, xyz_A[:, :, i].T.reshape(-1), grid_A)
        for i, axis in enumerate(("x", "y", "z"))
    ]
    doc = {
        "version": "0.3.0",
        "encoder": "IMP.bff scripts/trajectory_to_bcif.py",
        "dataBlocks": [{
            "header": block,
            "categories": [{
                "name": category,
                "rowCount": n_frames * n_atoms,
                "columns": columns,
            }],
        }],
    }
    blob = msgpack.packb(doc, use_bin_type=True)
    path.write_bytes(blob)
    return len(blob)


def read_bcif(path: Path, n_atoms: int, category: str = "_rotamer_coord") -> np.ndarray:
    """Read it back, for verification. Uses python-ihm, which *can* read this."""
    import ihm.format_bcif

    got: dict[str, list[float]] = {"x": [], "y": [], "z": []}

    class Handler:
        not_in_file = omitted = unknown = None
        _int_keys = _bool_keys = ()
        _float_keys = ("x", "y", "z")
        _keys = ("x", "y", "z")

        def __init__(self):
            self.category = category

        def __call__(self, x, y, z):
            got["x"].append(x)
            got["y"].append(y)
            got["z"].append(z)

        def end_save_frame(self):
            pass

    with open(path, "rb") as fh:
        ihm.format_bcif.BinaryCifReader(fh, {category: Handler()}).read_file()
    flat = np.array([got["x"], got["y"], got["z"]], dtype=np.float64)   # (3, A*F)
    return flat.reshape(3, n_atoms, -1).transpose(2, 1, 0)             # (F, A, 3)


def load_trajectory(path: Path, top: Path | None) -> np.ndarray:
    """Coordinates in **Angstrom** from a DCD or an XTC."""
    if path.suffix.lower() == ".dcd":
        import IMP.bff.io.structure as ios
        return np.asarray(ios.read_dcd(str(path)), dtype=np.float64)   # already A
    if path.suffix.lower() in (".xtc", ".trr"):
        try:
            import mdtraj as md
        except ImportError:                                   # pragma: no cover
            sys.exit(f"reading {path.suffix} needs mdtraj (converter-time only)")
        if top is None:
            sys.exit(f"{path.suffix} needs --top (a .gro/.pdb topology)")
        return np.asarray(md.load(str(path), top=str(top)).xyz,
                          dtype=np.float64) * 10.0            # nm -> A
    sys.exit(f"unsupported trajectory: {path}")


def convert(src: Path, dst: Path, top: Path | None, grid_A: float,
            verify: bool) -> None:
    xyz = load_trajectory(src, top)
    n_frames, n_atoms, _ = xyz.shape
    n_bytes = write_bcif(dst, xyz, grid_A)
    before = src.stat().st_size
    coords = xyz.size
    line = (f"{src.name}: {n_frames} x {n_atoms} x 3  "
            f"{before / 1e6:6.2f} MB -> {n_bytes / 1e6:6.2f} MB  "
            f"({n_bytes / coords:.2f} B/coord, {n_bytes / before * 100:.0f} %)")
    if verify:
        back = read_bcif(dst, n_atoms)
        want = np.round(xyz / grid_A) * grid_A
        if back.shape != want.shape:
            sys.exit(f"{src.name}: shape {back.shape} != {want.shape}")
        err = float(np.abs(back - want).max())
        if err > 1e-9:
            sys.exit(f"{src.name}: round trip differs by {err:.3e} A")
        line += f"  verified exact on the {grid_A} A grid"
    print(line)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("src", type=Path, nargs="?", help="a .dcd or .xtc")
    ap.add_argument("dst", type=Path, nargs="?", help="the .bcif to write")
    ap.add_argument("--top", type=Path, default=None,
                    help="topology for an XTC (.gro/.pdb)")
    ap.add_argument("--grid", type=float, default=DEFAULT_GRID_A,
                    metavar="A", help=f"grid in A (default {DEFAULT_GRID_A})")
    ap.add_argument("--all", type=Path, default=None, metavar="DIR",
                    help="convert every .dcd in DIR next to its source")
    ap.add_argument("--no-verify", action="store_true",
                    help="skip the round-trip check")
    a = ap.parse_args(argv)

    if a.all:
        for dcd in sorted(a.all.glob("*.dcd")):
            convert(dcd, dcd.with_suffix(".bcif"), None, a.grid, not a.no_verify)
        return 0
    if not a.src or not a.dst:
        ap.error("give src and dst, or --all DIR")
    convert(a.src, a.dst, a.top, a.grid, not a.no_verify)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
