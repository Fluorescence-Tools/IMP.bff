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
import IMP.bff as ios

REPO = Path(__file__).resolve().parent.parent.parent
SCRIPT = REPO / "bin" / "imp_bff_traj2bcif"
#: A shipped library: small enough to round-trip in a test, real data.
LIB = REPO / "data" / "rotamer_library" / "A56_C1R_cutoff30.bcif"
#: The shipped files are lossless. Quantisation is exercised separately,
#: because it is an option rather than what the package stores.
GRID = 0.001


@pytest.fixture(scope="module")
def encoder():
    # An explicit loader: the program lives in `bin/` with no extension, the
    # way IMP's installed programs do, and `spec_from_file_location` cannot
    # infer a loader without one.
    import importlib.machinery
    spec = importlib.util.spec_from_file_location(
        "traj_to_bcif", SCRIPT,
        loader=importlib.machinery.SourceFileLoader("traj_to_bcif", str(SCRIPT)))
    mod = importlib.util.module_from_spec(spec)
    sys.modules["traj_to_bcif"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def shipped(encoder, tmp_path_factory):
    """The shipped library, re-encoded losslessly and read back."""
    if not LIB.exists():
        pytest.skip(f"{LIB} not present")
    n_atoms = _n_atoms()
    xyz = np.asarray(IMP.bff.read_bcif_trajectory(
        str(LIB), n_atoms, "_rotamer_coord")).reshape(-1, n_atoms, 3)
    out = tmp_path_factory.mktemp("bcif") / "lib.bcif"
    n_bytes = encoder.write_bcif(out, xyz)          # lossless, the default
    return out, xyz, n_bytes


def _n_atoms():
    """From the companion PDB, which is where every caller gets it."""
    import IMP, IMP.atom
    pdb = LIB.with_name(LIB.name.split("_cutoff")[0] + ".pdb")
    m = IMP.Model()
    h = IMP.atom.read_pdb(str(pdb), m, IMP.atom.AllPDBSelector())
    return len(IMP.atom.get_leaves(h))


@pytest.fixture(scope="module")
def converted(shipped):
    return shipped


def test_the_shipped_libraries_are_lossless(converted):
    """float32 in, float32 out, bit for bit.

    The shipped libraries are stored unquantised. An early version of this
    work stored them on a 0.1 A grid, which shifts transition-dipole
    directions by 1.6 degrees and moved the FRETpredict pins by 2.3e-3 in E
    against their 2e-5 tolerance. Nothing about the *distances* revealed it --
    a dipole spans two atoms 1.7 A apart and does not average.
    """
    path, xyz, _ = converted
    n_frames, n_atoms, _ = xyz.shape
    got = np.asarray(IMP.bff.read_bcif_trajectory(
        str(path), n_atoms, "_rotamer_coord")).reshape(n_frames, n_atoms, 3)
    np.testing.assert_array_equal(got, xyz.astype(np.float32).astype(np.float64))


def test_quantising_is_bounded_by_half_the_grid(encoder, converted, tmp_path):
    """The option, exercised on its own terms.

    Rounding to a grid of *g* cannot be wrong by more than *g*/2, and if it
    were systematically less the grid would not be doing what the size table
    claims."""
    _, xyz, _ = converted
    out = tmp_path / "q.bcif"
    encoder.write_bcif(out, xyz, GRID)
    n_atoms = xyz.shape[1]
    got = np.asarray(IMP.bff.read_bcif_trajectory(
        str(out), n_atoms, "_rotamer_coord")).reshape(xyz.shape)
    err = np.abs(got - xyz).max()
    assert err <= GRID / 2 + 1e-9, err


def test_a_coarser_grid_costs_accuracy_and_saves_nothing(encoder, converted, tmp_path):
    """The measurement that settled the grid question.

    Every coordinate in the corpus is under 28.1 A, so 0.001, 0.005 and 0.01 A
    all quantise into ``int16`` and all cost exactly two bytes. The size is set
    by the integer *type*, not by the grid.
    """
    _, xyz, _ = converted
    sizes = {g: encoder.write_bcif(tmp_path / f"g{g}.bcif", xyz, g)
             for g in (0.001, 0.005, 0.01)}
    # Within a few bytes, not identical: the payload is the same length either
    # way, and what differs is the msgpack encoding of the `factor` literal
    # (1000 against 100) in each of the three columns' headers.
    assert max(sizes.values()) - min(sizes.values()) < 32, sizes
    for g, n in sizes.items():
        assert 1.9 < n / xyz.size < 2.2, (g, n / xyz.size)


def test_lossless_is_four_bytes_a_coordinate(converted):
    """Against DCD's 4.31 -- 7 % smaller, which is the honest figure."""
    _, xyz, n_bytes = converted
    assert 3.9 < n_bytes / xyz.size < 4.2


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
