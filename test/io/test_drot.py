"""`.drot` rotamer libraries: write, read back, and what the shipped set holds.

`.drot` is the internal-coordinate library store of PRD-118 -- a brotli+tar
container carrying the template, the Z-matrix and, per conformer, its own base
coordinates, bond lengths, bond angles and dihedrals. Both halves live in this
package (`write_drot`/`read_drot`, C++, vendored brotli both ways), so they are
tested together against real shipped data: a writer is only correct with
respect to a reader, and a container format that round-trips through its own
two halves has still proved nothing until the numbers come back.

The gate that matters is **losslessness**. The FRETpredict parity pins
(`test/cgprobe/rotamer/test_fretpredict_pins.py`) are recorded at 1e-5 in E, and
a rotamer library that reconstructs to 0.02 A moves them by 1.3e-3 -- which is
how the torsion-only v8 layout was caught. The default rung stores what a
conformer actually has and comes back at ~1e-6 A.
"""

import importlib.machinery
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

import IMP.bff

REPO = Path(__file__).resolve().parent.parent.parent
DATA = REPO / "data" / "rotamer_library"
SCRIPT = REPO / "bin" / "imp_bff_traj2drot"

#: A shipped library small enough to write inside a test, real data.
STEM, CUTOFF = "A56_C1R", 30
#: float32 grids put the floor at ~1e-5 A for a 30 A molecule.
LOSSLESS_TOL_A = 1e-4


def _pdb(stem):
    """Atom names, elements, residue names of a shipped template."""
    names, elements, resnames = [], [], []
    path = DATA / f"{stem}.pdb"
    if not path.exists():
        pytest.skip(f"{path} not present")
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            names.append(name)
            resnames.append(line[17:20].strip())
            elements.append(line[76:78].strip() or name[:1])
    return names, elements, resnames


@pytest.fixture(scope="module")
def ensemble():
    """(coords, weights, names, elements, resnames) of a shipped library."""
    bcif = DATA / f"{STEM}_cutoff{CUTOFF}.bcif"
    if not bcif.exists():
        pytest.skip(f"{bcif} not present")
    names, elements, resnames = _pdb(STEM)
    xyz = np.asarray(IMP.bff.read_bcif_trajectory(
        str(bcif), len(names), "_rotamer_coord"),
        dtype=np.float64).reshape(-1, len(names), 3)
    weights = np.loadtxt(DATA / f"{STEM}_cutoff{CUTOFF}_weights.txt",
                         dtype=np.float64).reshape(-1)
    return xyz, weights, names, elements, resnames


@pytest.fixture(scope="module")
def program():
    # The program lives in `bin/` with no extension, the way IMP's installed
    # programs do, so the loader has to be spelled out.
    spec = importlib.util.spec_from_file_location(
        "traj_to_drot", SCRIPT,
        loader=importlib.machinery.SourceFileLoader("traj_to_drot",
                                                    str(SCRIPT)))
    mod = importlib.util.module_from_spec(spec)
    sys.modules["traj_to_drot"] = mod
    spec.loader.exec_module(mod)
    return mod


def _write(path, ensemble, encoding=None):
    xyz, weights, names, elements, resnames = ensemble
    IMP.bff.write_probe_rotamer_drot(str(path), np.ascontiguousarray(xyz).ravel(),
                       names, elements, resnames,
                       np.ascontiguousarray(weights),
                       encoding or IMP.bff.ProbeRotamerDrotEncoding())
    return path


def _read(path):
    lib = IMP.bff.read_probe_rotamer_drot(str(path))
    xyz = np.asarray(lib.get_coords(), dtype=np.float64).reshape(
        lib.n_rotamers, lib.n_atoms, 3)
    return lib, xyz


def test_roundtrip_is_lossless(tmp_path, ensemble):
    """The default rung reproduces the conformers, not an approximation."""
    xyz, weights, names, _, resnames = ensemble
    lib, back = _read(_write(tmp_path / "lib.drot", ensemble))

    assert back.shape == xyz.shape
    assert np.abs(back - xyz).max() < LOSSLESS_TOL_A
    assert list(lib.atom_names) == names
    assert list(lib.resnames) == resnames
    # Weights travel as written -- normalising is the loader's business.
    assert np.allclose(np.asarray(lib.get_weights()), weights)


def test_the_store_is_self_contained(tmp_path, ensemble):
    """No PDB beside it: names, residues and elements ride in the file."""
    out = tmp_path / "alone" / "lib.drot"
    out.parent.mkdir()
    _write(out, ensemble)
    assert not list(out.parent.glob("*.pdb"))

    library = IMP.bff.load_probe_rotamer_library(str(out))
    _, _, names, _, resnames = ensemble
    assert list(library.atom_names) == names
    assert list(library.resnames) == resnames
    assert np.isclose(np.asarray(library.weights).sum(), 1.0)


def test_compact_rung_is_smaller_and_says_so(tmp_path, ensemble):
    """`--grid` buys size at a cost the format states rather than hides."""
    xyz = ensemble[0]
    lossless = _write(tmp_path / "lossless.drot", ensemble)
    encoding = IMP.bff.ProbeRotamerDrotEncoding()
    encoding.lossless = False
    encoding.grid_a = 0.001
    encoding.grid_deg = 0.01
    compact = _write(tmp_path / "compact.drot", ensemble, encoding)

    assert compact.stat().st_size < lossless.stat().st_size
    err = np.abs(_read(compact)[1] - xyz).max()
    # Well outside the pins' reach, well inside a picture's.
    assert LOSSLESS_TOL_A < err < 0.1


def test_writing_twice_gives_the_same_content(tmp_path, ensemble):
    """Container UIDs vary; every decoded library value must be identical."""
    a = _write(tmp_path / "a.drot", ensemble)
    b = _write(tmp_path / "b.drot", ensemble)
    first, first_coords = _read(a)
    second, second_coords = _read(b)
    np.testing.assert_array_equal(first_coords, second_coords)
    np.testing.assert_array_equal(first.weights, second.weights)
    for column in ("atom_names", "elements", "resnames"):
        assert list(getattr(first, column)) == list(getattr(second, column))
    assert first.metadata == second.metadata


def test_shipped_libraries_reproduce_their_bcif():
    """Every shipped `.drot` still is the ensemble the `.bcif` holds.

    This is the standing guard on the corpus: the loader prefers `.drot`
    (`resolve_rotamer_library_path`), so a re-encoding that drifted would move
    every rotamer result in the package without any test naming it.
    """
    family = DATA / "dyes.drot.pto"
    if family.exists():                       # the shipped shape: one family
        libraries = list(IMP.bff.probe_rotamer_drot_catalog(str(family)))[:6]
        pairs = [f"{family}::{name}" for name in libraries]
    else:                                     # a directory of one-library files
        pairs = [str(p) for p in sorted(DATA.glob("*_cutoff*.drot.pto"))[:6]]
    if not pairs:
        pytest.skip("no shipped .drot libraries")
    for drot in pairs:
        container, _, inside = str(drot).partition("::")
        stem = inside or Path(container).name.split(".drot")[0]
        bcif = DATA / f"{stem}.bcif"
        if not bcif.exists():
            continue
        lib, back = _read(drot)
        xyz = np.asarray(IMP.bff.read_bcif_trajectory(
            str(bcif), lib.n_atoms, "_rotamer_coord"),
            dtype=np.float64).reshape(-1, lib.n_atoms, 3)
        assert back.shape == xyz.shape, stem
        assert np.abs(back - xyz).max() < LOSSLESS_TOL_A, stem


def test_a_family_container_holds_many_libraries(tmp_path, ensemble):
    """Bundling is a copy: the parts come back out unchanged.

    The shipped layout is one container per *family* -- every dye in
    `dyes.drot.pto` -- because PTO addresses objects individually, so reading
    one library out of ninety-five costs that library and not the file. What
    has to hold is that bundling changes nothing: same coordinates, same
    weights, same names, and a catalog that says what is inside.
    """
    xyz, weights, names, elements, resnames = ensemble
    singles = []
    for i, n in enumerate((3, 5)):                 # two small, different libraries
        one = tmp_path / f"lib{i}.drot.pto"
        IMP.bff.write_probe_rotamer_drot(str(one), np.ascontiguousarray(xyz[:n]).ravel(),
                           names, elements, resnames, weights[:n],
                           IMP.bff.ProbeRotamerDrotEncoding())
        singles.append(one)

    family = tmp_path / "family.drot.pto"
    IMP.bff.write_probe_rotamer_drot_bundle([str(p) for p in singles], ["alpha", "beta"],
                              str(family))
    assert list(IMP.bff.probe_rotamer_drot_catalog(str(family))) == ["alpha", "beta"]

    for one, name in zip(singles, ("alpha", "beta")):
        want, want_xyz = _read(one)
        got, got_xyz = _read(f"{family}::{name}")      # the locator form
        assert got_xyz.shape == want_xyz.shape
        assert np.array_equal(got_xyz, want_xyz), name
        assert list(got.atom_names) == list(want.atom_names)
        assert np.array_equal(np.asarray(got.get_weights()),
                              np.asarray(want.get_weights()))

    with pytest.raises(Exception):
        IMP.bff.read_probe_rotamer_drot(f"{family}::gamma")          # names what it has


def test_rejects_shapes_that_do_not_agree(tmp_path, ensemble):
    xyz, weights, names, elements, resnames = ensemble
    out = str(tmp_path / "bad.drot")
    with pytest.raises(Exception):
        IMP.bff.write_probe_rotamer_drot(out, np.ascontiguousarray(xyz).ravel()[:-3],
                           names, elements, resnames, weights,
                           IMP.bff.ProbeRotamerDrotEncoding())
    with pytest.raises(Exception):
        IMP.bff.write_probe_rotamer_drot(out, np.ascontiguousarray(xyz).ravel(),
                           names, elements, resnames, weights[:-1],
                           IMP.bff.ProbeRotamerDrotEncoding())


def test_rejects_a_file_that_is_not_a_drot(tmp_path):
    junk = tmp_path / "not.drot"
    junk.write_bytes(b"\x00\x01\x02not a brotli stream")
    with pytest.raises(Exception):
        IMP.bff.read_probe_rotamer_drot(str(junk))


def test_program_converts_a_shipped_library(tmp_path, program, capsys):
    """`imp_bff_traj2drot lib.bcif lib.drot` -- the builder end to end."""
    bcif = DATA / f"{STEM}_cutoff{CUTOFF}.bcif"
    if not bcif.exists():
        pytest.skip(f"{bcif} not present")
    out = tmp_path / "built.drot"
    program.convert(bcif, out, DATA / f"{STEM}.pdb",
                    DATA / f"{STEM}_cutoff{CUTOFF}_weights.txt",
                    None, None, 0.01, True)
    assert "verified" in capsys.readouterr().out
    lib, back = _read(out)
    xyz = np.asarray(IMP.bff.read_bcif_trajectory(
        str(bcif), lib.n_atoms, "_rotamer_coord"),
        dtype=np.float64).reshape(-1, lib.n_atoms, 3)
    assert np.abs(back - xyz).max() < LOSSLESS_TOL_A


def test_program_clusters_raw_frames(tmp_path, program, ensemble):
    """The raw-input path: leaders become rotamers, populations the weights."""
    xyz = ensemble[0]
    leaders, counts = program.cluster(xyz, 1.5)
    assert 0 < len(leaders) <= len(xyz)
    assert counts.sum() == len(xyz)
    assert len(counts) == len(leaders)


def test_a_library_without_residue_names_round_trips(tmp_path, ensemble):
    """Residue names are optional -- the template loop drops the column."""
    xyz, weights, names, elements, _ = ensemble
    out = tmp_path / "nores.drot"
    IMP.bff.write_probe_rotamer_drot(str(out), np.ascontiguousarray(xyz).ravel(),
                       names, elements, [], np.ascontiguousarray(weights),
                       IMP.bff.ProbeRotamerDrotEncoding())
    lib, back = _read(out)
    assert list(lib.atom_names) == names
    assert list(lib.resnames) == []
    assert np.abs(back - xyz).max() < LOSSLESS_TOL_A


def test_compact_rung_refuses_a_grid_it_cannot_hold(tmp_path, ensemble):
    """int16 has 32767 counts; a grid finer than the molecule is an error.

    Silently wrapping would write a library that reads back as a different
    molecule, which is the one failure mode a store must not have.
    """
    xyz, weights, names, elements, resnames = ensemble
    encoding = IMP.bff.ProbeRotamerDrotEncoding()
    encoding.lossless = False
    encoding.grid_a = 1e-6
    with pytest.raises(Exception, match="int16|lossless"):
        IMP.bff.write_probe_rotamer_drot(str(tmp_path / "toofine.drot"),
                           np.ascontiguousarray(xyz).ravel(),
                           names, elements, resnames,
                           np.ascontiguousarray(weights), encoding)


def test_a_payload_that_compresses_hugely_still_reads(tmp_path):
    """A library of near-identical conformers compresses far better than eight
    times, and used to be reported as a corrupt brotli stream.

    The reader decoded with brotli's **one-shot** API and sized its output
    buffer as `compressed * 8`. That API reports a buffer that is too small as
    `BROTLI_DECODER_RESULT_ERROR`, not as `NEEDS_MORE_OUTPUT`, so the growth
    branch never fired and the too-small buffer came back as "corrupt". Every
    shipped `.drot` is under the ratio, which is why nothing had noticed; this
    file is not. The reader streams now, so there is nothing to guess.
    """
    n_rotamers, n_atoms = 4000, 40
    # a connected chain -- the writer builds a Z-matrix, so the molecule has to
    # be one -- repeated 4000 times so fixed container indexes do not dominate
    # the compression-ratio check, with a whisper of noise. This makes
    # the payload compress the way a real near-rigid library would
    rng = np.random.default_rng(4)
    one = np.array([[1.5 * i, 0.35 * (i % 2), 0.0] for i in range(n_atoms)])
    xyz = np.repeat(one[None, :, :], n_rotamers, axis=0)
    # only the first few move, so the grid is 4000 near-identical blocks -- the
    # shape a library of a rigid dye on a short linker has
    xyz[:, :4] += rng.normal(scale=0.05, size=(n_rotamers, 4, 3))
    weights = np.full(n_rotamers, 1.0 / n_rotamers)
    names = [f"C{i}" for i in range(n_atoms)]
    elements = ["C"] * n_atoms
    resnames = ["DYE"] * n_atoms

    path = tmp_path / "flat.drot.pto"
    IMP.bff.write_probe_rotamer_drot(str(path), np.ascontiguousarray(xyz).ravel(),
                       names, elements, resnames,
                       np.ascontiguousarray(weights))
    raw = n_rotamers * n_atoms * 3 * 4          # float32 grids
    assert path.stat().st_size * 8 < raw, "this file has to beat the old guess"

    lib = IMP.bff.read_probe_rotamer_drot(str(path))
    assert lib.n_rotamers == n_rotamers
    assert lib.n_atoms == n_atoms
    back = np.asarray(lib.get_coords(), dtype=np.float64).reshape(
        n_rotamers, n_atoms, 3)
    np.testing.assert_allclose(back, xyz, atol=1e-3)
