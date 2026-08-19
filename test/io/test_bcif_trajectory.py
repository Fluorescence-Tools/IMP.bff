"""BinaryCIF trajectories: the encoder script and the C++ reader, round-tripped.

BinaryCIF is this package's trajectory format as of 2026-08-19. It replaces DCD
and XTC and is smaller than either -- 1.27 bytes per coordinate against DCD's
4.31 and XTC's 1.60 -- because it quantises to a 0.1 A grid, delta-encodes
along each atom's own frame series, and packs the deltas into int8 with escape
runs. The reasoning and what the precision costs are in
``okf/validation/bcif_for_trajectories.md``.

Both halves are exercised here against real shipped data, because they are only
correct *together*: the encoder is in ``scripts/`` and the decoder is
``ihm_format.c``, which IMP vendors, and neither is covered by the other's
tests.
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

import IMP.bff
import IMP.bff.io.structure as ios

REPO = Path(__file__).resolve().parent.parent.parent
SCRIPT = REPO / "scripts" / "trajectory_to_bcif.py"
#: Small enough to convert in a test, real enough to be worth converting.
DCD = REPO / "data" / "rotamer_library" / "A56_C1R_cutoff30.dcd"
GRID = 0.1


@pytest.fixture(scope="module")
def encoder():
    spec = importlib.util.spec_from_file_location("traj_to_bcif", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["traj_to_bcif"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def converted(encoder, tmp_path_factory):
    if not DCD.exists():
        pytest.skip(f"{DCD} not present")
    xyz = np.asarray(ios.read_dcd(str(DCD)), dtype=np.float64)
    out = tmp_path_factory.mktemp("bcif") / "lib.bcif"
    n_bytes = encoder.write_bcif(out, xyz, GRID)
    return out, xyz, n_bytes


def test_the_reader_recovers_the_coordinates_on_the_grid(converted):
    """Exact, not approximate: the only difference allowed is the quantisation
    the encoder announces, and that is applied to the reference here too."""
    path, xyz, _ = converted
    n_frames, n_atoms, _ = xyz.shape
    flat = np.asarray(IMP.bff.read_bcif_trajectory(
        str(path), n_atoms, "_rotamer_coord"))
    got = flat.reshape(n_frames, n_atoms, 3)
    want = np.round(xyz / GRID) * GRID
    np.testing.assert_allclose(got, want, rtol=0, atol=1e-9)


def test_the_error_against_the_original_is_half_the_grid(converted):
    """A tighter statement than "close": rounding to a grid of *g* cannot be
    wrong by more than *g*/2, and if it were systematically less the grid would
    not be doing what the size table claims."""
    path, xyz, _ = converted
    n_frames, n_atoms, _ = xyz.shape
    got = np.asarray(IMP.bff.read_bcif_trajectory(
        str(path), n_atoms, "_rotamer_coord")).reshape(xyz.shape)
    err = np.abs(got - xyz).max()
    assert err <= GRID / 2 + 1e-9, err


def test_it_is_smaller_than_the_dcd_it_replaces(converted):
    path, xyz, n_bytes = converted
    assert n_bytes < DCD.stat().st_size
    assert n_bytes / xyz.size < 2.0, "should be under 2 bytes per coordinate"


def test_the_row_count_is_readable_without_decoding_everything(converted):
    path, xyz, _ = converted
    n_frames, n_atoms, _ = xyz.shape
    assert IMP.bff.bcif_trajectory_rows(str(path), "_rotamer_coord") == \
        n_frames * n_atoms


def test_a_row_count_that_is_not_a_multiple_of_n_atoms_raises(converted):
    """Silence here would reshape the trajectory wrongly and return something
    that looks like coordinates."""
    path, xyz, _ = converted
    with pytest.raises(ValueError, match="not a multiple"):
        IMP.bff.read_bcif_trajectory(str(path), xyz.shape[1] + 1,
                                     "_rotamer_coord")


def test_a_missing_file_raises_rather_than_returning_nothing(tmp_path):
    with pytest.raises(Exception):
        IMP.bff.read_bcif_trajectory(str(tmp_path / "nope.bcif"), 1,
                                     "_rotamer_coord")


def test_the_escape_sentinel_is_not_mistaken_for_a_value(encoder, tmp_path):
    """IntegerPacking's int8 sentinels, 127 and -128, are legitimate deltas.

    Emitting one as a plain value makes the decoder read it as an escape, and
    the three columns then decode to *different lengths*. Constructed here
    directly: coordinates whose successive differences land exactly on the
    sentinels at the chosen grid.
    """
    steps = np.array([127, -128, 127, 126, -127, 128, 0, 255], dtype=np.float64)
    coords = np.cumsum(steps) * GRID
    xyz = np.zeros((len(coords), 1, 3))
    xyz[:, 0, 0] = coords
    xyz[:, 0, 1] = coords[::-1]
    xyz[:, 0, 2] = -coords
    out = tmp_path / "sentinel.bcif"
    encoder.write_bcif(out, xyz, GRID)
    got = np.asarray(IMP.bff.read_bcif_trajectory(
        str(out), 1, "_rotamer_coord")).reshape(xyz.shape)
    np.testing.assert_allclose(got, np.round(xyz / GRID) * GRID,
                               rtol=0, atol=1e-9)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
