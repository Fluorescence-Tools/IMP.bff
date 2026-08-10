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

import struct
from pathlib import Path

import numpy as np

__all__ = ["read_dcd", "read_dcd_header"]

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
