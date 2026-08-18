"""FRETpredict parity tests (expensive: need FRETpredict + MDAnalysis installed)."""

from __future__ import annotations

from pathlib import Path
import pickle

import pytest

from IMP.bff.representation.rotamer.fret import RotamerFRET


@pytest.mark.parametrize("fixed_r0", [False, True])
@pytest.mark.parametrize("electrostatic", [False, True])
@pytest.mark.parametrize("temperature", [293, 300])
def test_rotamer_fretpredict_parity_parameter_matrix(
    tmp_path: Path,
    fixed_r0: bool,
    electrostatic: bool,
    temperature: int,
) -> None:
    """Compare IMP rotamer FRET against FRETpredict across core options."""
    FRETpredict = pytest.importorskip("FRETpredict")
    MDAnalysis = pytest.importorskip("MDAnalysis")
    pdb = _hsp90_pdb()
    if not pdb.exists():
        pytest.skip("Hsp90 test fixture is not available")

    fp_prefix = tmp_path / "fretpredict"
    imp_prefix = tmp_path / "imp"
    kwargs = _base_fretpredict_kwargs(fixed_r0=fixed_r0, electrostatic=electrostatic, temperature=temperature)

    fretpredict = FRETpredict.FRETpredict(
        protein=MDAnalysis.Universe(str(pdb)),
        output_prefix=str(fp_prefix),
        verbose=False,
        **kwargs,
    )
    fretpredict.run()
    imp = RotamerFRET(pdb, output_prefix=str(imp_prefix), **kwargs)
    imp.run()

    _assert_output_files_match(fp_prefix, imp_prefix, residues=kwargs["residues"])
    _assert_summary_files_match(fp_prefix, imp_prefix, residues=kwargs["residues"])


@pytest.mark.parametrize("fixed_r0", [False, True])
def test_rotamer_fretpredict_parity_user_reweight(tmp_path: Path, fixed_r0: bool) -> None:
    """Compare user-weighted IMP rotamer FRET against FRETpredict."""
    FRETpredict = pytest.importorskip("FRETpredict")
    MDAnalysis = pytest.importorskip("MDAnalysis")
    from MDAnalysis.coordinates.memory import MemoryReader

    pdb = _hsp90_pdb()
    if not pdb.exists():
        pytest.skip("Hsp90 test fixture is not available")

    fp_prefix = tmp_path / "fretpredict"
    imp_prefix = tmp_path / "imp"
    rmf = tmp_path / "hsp90_two_frames.rmf3"
    _write_two_frame_rmf(pdb, rmf)
    universe = MDAnalysis.Universe(str(pdb))
    positions = [universe.atoms.positions.copy(), universe.atoms.positions.copy()]
    universe.load_new(_as_numpy_positions(positions), format=MemoryReader, order="fac")
    universe.trajectory.filename = "memory"
    kwargs = _base_fretpredict_kwargs(fixed_r0=fixed_r0, electrostatic=True, temperature=293)

    fretpredict = FRETpredict.FRETpredict(
        protein=universe,
        output_prefix=str(fp_prefix),
        verbose=False,
        **kwargs,
    )
    fretpredict.run()
    imp = RotamerFRET(rmf, output_prefix=str(imp_prefix), **kwargs)
    imp.run()

    weights = [0.25, 0.75]
    fretpredict.reweight(user_weights=weights)
    imp.reweight(user_weights=weights)
    _assert_reweighted_pkl_match(fp_prefix, imp_prefix)


def test_rotamer_fretpredict_parity_boltzmann_reweight(tmp_path: Path) -> None:
    """Compare partition-weighted IMP rotamer FRET against FRETpredict."""
    FRETpredict = pytest.importorskip("FRETpredict")
    MDAnalysis = pytest.importorskip("MDAnalysis")
    pdb = _hsp90_pdb()
    if not pdb.exists():
        pytest.skip("Hsp90 test fixture is not available")

    fp_prefix = tmp_path / "fretpredict"
    imp_prefix = tmp_path / "imp"
    kwargs = _base_fretpredict_kwargs(fixed_r0=True, electrostatic=True, temperature=293)

    fretpredict = FRETpredict.FRETpredict(
        protein=MDAnalysis.Universe(str(pdb)),
        output_prefix=str(fp_prefix),
        verbose=False,
        **kwargs,
    )
    fretpredict.run()
    imp = RotamerFRET(pdb, output_prefix=str(imp_prefix), **kwargs)
    imp.run()

    fretpredict.reweight(boltzmann_weights=True)
    imp.reweight(boltzmann_weights=True)
    _assert_summary_files_match(fp_prefix, imp_prefix)


def test_rotamer_fretpredict_parity_swapped_sites(tmp_path: Path) -> None:
    """Compare swapped donor/acceptor sites against FRETpredict."""
    FRETpredict = pytest.importorskip("FRETpredict")
    MDAnalysis = pytest.importorskip("MDAnalysis")
    pdb = _hsp90_pdb()
    if not pdb.exists():
        pytest.skip("Hsp90 test fixture is not available")

    fp_prefix = tmp_path / "fretpredict"
    imp_prefix = tmp_path / "imp"
    kwargs = {
        "residues": [637, 452],
        "temperature": 293,
        "chains": ["B", "A"],
        "donor": "AlexaFluor 568",
        "acceptor": "AlexaFluor 594",
        "libname_1": "AlexaFluor 568 C1R cutoff30",
        "libname_2": "AlexaFluor 594 C1R cutoff30",
        "electrostatic": True,
        "fixed_R0": True,
        "r0": 5.5,
    }

    fretpredict = FRETpredict.FRETpredict(
        protein=MDAnalysis.Universe(str(pdb)),
        output_prefix=str(fp_prefix),
        verbose=False,
        **kwargs,
    )
    fretpredict.run()
    imp = RotamerFRET(pdb, output_prefix=str(imp_prefix), **kwargs)
    imp.run()

    _assert_output_files_match(fp_prefix, imp_prefix, residues=kwargs["residues"])
    _assert_summary_files_match(fp_prefix, imp_prefix, residues=kwargs["residues"])


def test_rotamer_fretpredict_parity_two_frame_rmf(tmp_path: Path) -> None:
    """Compare a two-frame IMP RMF trajectory against FRETpredict memory frames."""
    FRETpredict = pytest.importorskip("FRETpredict")
    MDAnalysis = pytest.importorskip("MDAnalysis")
    from MDAnalysis.coordinates.memory import MemoryReader

    pdb = _hsp90_pdb()
    if not pdb.exists():
        pytest.skip("Hsp90 test fixture is not available")

    fp_prefix = tmp_path / "fretpredict"
    imp_prefix = tmp_path / "imp"
    rmf = tmp_path / "hsp90_two_frames.rmf3"
    _write_two_frame_rmf(pdb, rmf)
    universe = MDAnalysis.Universe(str(pdb))
    positions = [universe.atoms.positions.copy(), universe.atoms.positions.copy()]
    universe.load_new(_as_numpy_positions(positions), format=MemoryReader, order="fac")
    universe.trajectory.filename = "memory"
    kwargs = _base_fretpredict_kwargs(fixed_r0=True, electrostatic=True, temperature=293)

    fretpredict = FRETpredict.FRETpredict(
        protein=universe,
        output_prefix=str(fp_prefix),
        verbose=False,
        **kwargs,
    )
    fretpredict.run()
    imp = RotamerFRET(rmf, output_prefix=str(imp_prefix), **kwargs)
    imp.run()

    _assert_output_files_match(fp_prefix, imp_prefix, residues=kwargs["residues"])
    _assert_summary_files_match(fp_prefix, imp_prefix, residues=kwargs["residues"])


def _base_fretpredict_kwargs(*, fixed_r0: bool, electrostatic: bool, temperature: int) -> dict[str, object]:
    """Return shared Hsp90 FRETpredict/IMP kwargs."""
    kwargs: dict[str, object] = {
        "residues": [452, 637],
        "temperature": temperature,
        "chains": ["A", "B"],
        "donor": "AlexaFluor 594",
        "acceptor": "AlexaFluor 568",
        "libname_1": "AlexaFluor 594 C1R cutoff30",
        "libname_2": "AlexaFluor 568 C1R cutoff30",
        "electrostatic": electrostatic,
        "fixed_R0": fixed_r0,
    }
    if fixed_r0:
        kwargs["r0"] = 5.5
    return kwargs


def _hsp90_pdb() -> Path:
    """Unpack the bundled Hsp90 fixture (test/cgdye/rotamer/data) once per process."""
    import gzip
    import shutil
    import tempfile
    src = Path(__file__).resolve().parent / "data" / "openHsp90.pdb.gz"
    dst = Path(tempfile.gettempdir()) / "imp_bff_cgdye_openHsp90.pdb"
    if not dst.exists():
        with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    return dst


def _assert_output_files_match(fp_prefix: Path, imp_prefix: Path, residues: tuple[int, int] | list[int] = (452, 637)) -> None:
    """Assert IMP and FRETpredict output arrays match."""
    r1, r2 = residues
    for suffix in ["Z", "w_s", "k2", "Es", "Ed1", "Ed2"]:
        fp = _read_dat_values(fp_prefix.with_name(fp_prefix.name + f"-{suffix}-{r1}-{r2}.dat"))
        imp = _read_dat_values(imp_prefix.with_name(imp_prefix.name + f"-{suffix}-{r1}-{r2}.dat"))
        assert imp == pytest.approx(fp, rel=1e-6, abs=1e-6)


def _assert_summary_files_match(fp_prefix: Path, imp_prefix: Path, residues: tuple[int, int] | list[int] = (452, 637)) -> None:
    """Assert IMP and FRETpredict per-frame summary files match."""
    r1, r2 = residues
    for suffix in ["k2", "Es", "Ed1", "Ed2"]:
        fp = _read_dat_values(fp_prefix.with_name(fp_prefix.name + f"-{suffix}-{r1}-{r2}.dat"))
        imp = _read_dat_values(imp_prefix.with_name(imp_prefix.name + f"-{suffix}-{r1}-{r2}.dat"))
        assert imp == pytest.approx(fp, rel=1e-6, abs=1e-6)


def _assert_reweighted_pkl_match(fp_prefix: Path, imp_prefix: Path) -> None:
    """Assert IMP and FRETpredict reweighted summary pickles match."""
    fp = _read_pkl(fp_prefix.with_name(fp_prefix.name + "-data-452-637.pkl"))
    imp = _read_pkl(imp_prefix.with_name(imp_prefix.name + "-data-452-637.pkl"))
    for key in ["k2", "Estatic", "Edynamic1", "Edynamic2"]:
        assert _pkl_value(imp, key) == pytest.approx(_pkl_value(fp, key), rel=1e-6, abs=1e-6)


def _read_pkl(path: Path) -> object:
    """Read a pickle file without tying the test to pandas."""
    return pickle.loads(path.read_bytes())


def _pkl_value(data: object, key: str) -> float:
    """Return a numeric value from a Series-like or DataFrame-like object."""
    value = data.loc[key]  # type: ignore[attr-defined]
    if hasattr(value, "iloc"):
        value = value.iloc[0]
    return float(value)


def _read_dat_values(path: Path) -> list[float]:
    """Read a FRETpredict-style numeric text file."""
    values: list[float] = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        values.extend(float(value) for value in line.split())
    return values


def _read_single_dat_value(prefix: Path, suffix: str, r1: int = 452, r2: int = 637) -> float:
    """Read a single numeric value from a FRETpredict-style output file."""
    values = _read_dat_values(prefix.with_name(prefix.name + f"-{suffix}-{r1}-{r2}.dat"))
    assert len(values) == 1
    return values[0]


def _as_numpy_positions(positions: list[object]) -> object:
    """Return positions as a frame-atom-coordinate numpy array."""
    import numpy as np

    return np.asarray(positions, dtype=np.float32)


def _write_two_frame_rmf(pdb: Path, rmf: Path) -> None:
    """Write a two-frame RMF file from a PDB hierarchy."""
    import IMP
    import IMP.atom
    import RMF

    model = IMP.Model()
    hierarchy = IMP.atom.read_pdb(str(pdb), model, IMP.atom.NonWaterPDBSelector())
    fh = RMF.create_rmf_file(str(rmf))
    IMP.rmf.add_hierarchies(fh, [hierarchy])
    IMP.rmf.save_frame(fh, "frame_0")
    IMP.rmf.save_frame(fh, "frame_1")
    del fh


# IMP runs every .py under test/ as a standalone script, and a file of bare
# pytest functions would import cleanly and exit 0 -- reporting success without
# running a single assertion. Hand the file to pytest explicitly so a failure
# here is a failure in ctest.
if __name__ == "__main__":
    import sys
    try:
        import pytest
    except ImportError:
        print("pytest not installed; skipping", __file__)
        sys.exit(0)
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
