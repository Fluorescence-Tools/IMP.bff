"""Tests for IMP.bff.cgdye.rotamer."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner
from IMP.bff.fret import io as fps
from IMP.bff.cgdye.rotamer.fps import read_rotamer_fps, rotamer_fret_from_fps
from IMP.bff.cgdye.rotamer.fret import RotamerFRET
from IMP.bff.cgdye.rotamer.io import load_protein_frames, load_rotamer_library
from IMP.bff.cgdye.rotamer.r0 import calculate_r0
# The rotamer CLI ships with cgdye now; invoking it through imp-tricks'
# aggregator would make an imp.bff test depend on the layer above it.
from IMP.bff.cgdye.rotamer.cli import rotamer as _rotamer_cli

def _fretpredict_test_system(*parts):
    """Locate a FRETpredict reference system, if a checkout is around.

    Returns None when it is not, which is the normal case: these systems are
    reference material for parity checks, not data IMP.bff ships.
    """
    from pathlib import Path as _P
    import os
    env = os.environ.get("FRETPREDICT_DIR")
    if env:
        return _P(env) / "tests" / "test_systems" / _P(*parts)
    return None



def test_calculate_r0_matches_fretpredict_pp11_reference() -> None:
    """R0 calculation matches the pp11 tutorial reference."""
    r0 = calculate_r0("AlexaFluor 488", "AlexaFluor 594", 0.684587)
    assert r0 == pytest.approx(5.712982, abs=1e-5)


def test_load_rotamer_library_from_fretpredict_name() -> None:
    """FRETpredict-style library names resolve to bundled RMF libraries."""
    lib = load_rotamer_library("AlexaFluor 488 C1R cutoff30")
    assert lib["coords"].ndim == 3
    assert lib["coords"].shape[0] > 1
    assert lib["weights"].sum() == pytest.approx(1.0)
    assert lib["metadata"]["mu"]


def test_fps_json_roundtrip(tmp_path: Path) -> None:
    """Generic fps.json read/write helpers round-trip payload data."""
    path = tmp_path / "labels.fps.json"
    fps.write_fps_json(path, {"p": {"residue_seq_number": 1}}, {"d": {"position1_name": "p", "position2_name": "p"}}, {"s": {}})
    positions, distances, score_sets, extra = fps.read_fps_json(path)
    assert positions["p"]["residue_seq_number"] == 1
    assert distances["d"]["position1_name"] == "p"
    assert score_sets == {"s": {}}
    assert extra == {}


def test_read_rotamer_fps_label_distributions(tmp_path: Path) -> None:
    """fps.json files expose rotamer label distributions."""
    fps = tmp_path / "labels.fps.json"
    payload = {
        "Positions": {
            "d1": {
                "chain_identifier": "A",
                "residue_seq_number": 10,
                "atom_name": "CA",
                "dye": "AlexaFluor 488",
                "rotamer_library": "AlexaFluor 488 C1R cutoff30",
            },
            "a1": {
                "chain": "B",
                "residue": 20,
                "dye_name": "AlexaFluor 594",
                "library": "AlexaFluor 594 C1R cutoff30",
            },
        },
        "Distances": {
            "d1_a1": {
                "position1_name": "d1",
                "position2_name": "a1",
            }
        },
    }
    fps.write_text(json.dumps(payload))
    donor, acceptor, distance, _positions, _distances, _extra = read_rotamer_fps(fps)
    assert donor.name == "d1"
    assert donor.chain == "A"
    assert donor.residue == 10
    assert donor.dye == "AlexaFluor 488"
    assert donor.library == "AlexaFluor 488 C1R cutoff30"
    assert acceptor.chain == "B"
    assert acceptor.residue == 20
    assert acceptor.dye == "AlexaFluor 594"
    assert acceptor.library == "AlexaFluor 594 C1R cutoff30"
    assert distance.name == "d1_a1"
    assert distance.donor_position == "d1"
    assert distance.acceptor_position == "a1"


def test_rotamer_fret_from_fps(tmp_path: Path) -> None:
    """fps.json files can construct a RotamerFRET object."""
    fps = tmp_path / "labels.fps.json"
    fps.write_text(
        json.dumps(
            {
                "Positions": {
                    "d1": {
                        "chain_identifier": "A",
                        "residue_seq_number": 452,
                        "dye": "AlexaFluor 594",
                        "rotamer_library": "AlexaFluor 594 C1R cutoff30",
                    },
                    "a1": {
                        "chain_identifier": "B",
                        "residue_seq_number": 637,
                        "dye": "AlexaFluor 568",
                        "rotamer_library": "AlexaFluor 568 C1R cutoff30",
                    },
                },
                "Distances": {
                    "d1_a1": {
                        "position1_name": "d1",
                        "position2_name": "a1",
                    }
                },
            }
        )
    )
    # The Hsp90 reference system belonged to the vendored FRETpredict checkout,
    # which is reference material rather than shipped data and did not follow
    # cgdye into imp.bff. Only the rotamer *library* did. Skip loudly rather
    # than silently comparing against nothing.
    pdb = _fretpredict_test_system("Hsp90", "openHsp90.pdb")
    if pdb is None or not pdb.exists():
        pytest.skip(
            "FRETpredict's Hsp90 reference system is not present; it is "
            "vendored reference material, not imp.bff module data")
    fret = rotamer_fret_from_fps(fps, pdb, fixed_R0=True, r0=5.5, output_prefix=str(tmp_path / "from_fps"))
    assert isinstance(fret, RotamerFRET)
    assert fret.residues == [452, 637]
    assert fret.chains == ["A", "B"]
    assert fret.donor == "AlexaFluor 594"
    assert fret.acceptor == "AlexaFluor 568"


def test_rotamer_cli_help_and_r0() -> None:
    """Rotamer CLI exposes help and R0 subcommand."""
    runner = CliRunner()
    result = runner.invoke(_rotamer_cli, ["--help"])
    assert result.exit_code == 0
    assert "predict" in result.output
    assert "r0" in result.output

    result = runner.invoke(_rotamer_cli, ["r0", "--donor", "AlexaFluor 488", "--acceptor", "AlexaFluor 594", "--k2", "0.684587"])
    assert result.exit_code == 0
    assert float(result.output.split()[0]) == pytest.approx(5.712982, abs=1e-5)


def test_load_protein_frames_from_rmf_trajectory() -> None:
    """RMF protein trajectories load multiple frames."""
    root = _repo_root()
    rmf = root / "examples" / "cgdye" / "hgbp1_alexa488_langevin.rmf3"
    if not rmf.exists():
        pytest.skip("RMF trajectory fixture is not available")
    frames = load_protein_frames(rmf, max_frames=2)
    assert len(frames) == 2
    assert frames[0]["coords"].ndim == 2
    assert frames[0]["residue_indices"]
    assert frames[1]["coords"].shape == frames[0]["coords"].shape


def _repo_root() -> Path:
    """Return the chisurf repository root."""
    for parent in Path(__file__).resolve().parents:
        if (parent / "modules" / "imp-tricks" / "src").exists():
            return parent
        if (parent / "src").exists() and (parent / "tests").exists():
            return parent
    return Path(__file__).resolve().parents[1]


def test_rotamer_fret_hsp90_fixed_r0_writes_outputs(tmp_path: Path) -> None:
    """Hsp90 single-frame fixed-R0 run writes FRETpredict-style files."""
    # The Hsp90 reference system belonged to the vendored FRETpredict checkout,
    # which is reference material rather than shipped data and did not follow
    # cgdye into imp.bff. Only the rotamer *library* did. Skip loudly rather
    # than silently comparing against nothing.
    pdb = _fretpredict_test_system("Hsp90", "openHsp90.pdb")
    if pdb is None or not pdb.exists():
        pytest.skip(
            "FRETpredict's Hsp90 reference system is not present; it is "
            "vendored reference material, not imp.bff module data")

    prefix = tmp_path / "E30_fixedR0"
    fret = RotamerFRET(
        pdb,
        [452, 637],
        chains=["A", "B"],
        donor="AlexaFluor 594",
        acceptor="AlexaFluor 568",
        libname_1="AlexaFluor 594 C1R cutoff30",
        libname_2="AlexaFluor 568 C1R cutoff30",
        electrostatic=True,
        fixed_R0=True,
        r0=5.5,
        output_prefix=str(prefix),
    )
    fret.run()

    for suffix in ["Z", "w_s", "k2", "Es", "Ed1", "Ed2"]:
        assert (tmp_path / f"E30_fixedR0-{suffix}-452-637.dat").exists()
    k2 = _read_single_dat_value(prefix, "k2")
    estatic = _read_single_dat_value(prefix, "Es")
    assert k2 == pytest.approx(0.970927, abs=1e-5)
    assert estatic == pytest.approx(0.539034, abs=1e-5)


def test_rotamer_fret_hsp90_calculated_r0_matches_fretpredict(tmp_path: Path) -> None:
    """Hsp90 single-frame calculated-R0 run matches FRETpredict."""
    # The Hsp90 reference system belonged to the vendored FRETpredict checkout,
    # which is reference material rather than shipped data and did not follow
    # cgdye into imp.bff. Only the rotamer *library* did. Skip loudly rather
    # than silently comparing against nothing.
    pdb = _fretpredict_test_system("Hsp90", "openHsp90.pdb")
    if pdb is None or not pdb.exists():
        pytest.skip(
            "FRETpredict's Hsp90 reference system is not present; it is "
            "vendored reference material, not imp.bff module data")

    prefix = tmp_path / "E30_calcR0"
    fret = RotamerFRET(
        pdb,
        [452, 637],
        chains=["A", "B"],
        donor="AlexaFluor 594",
        acceptor="AlexaFluor 568",
        libname_1="AlexaFluor 594 C1R cutoff30",
        libname_2="AlexaFluor 568 C1R cutoff30",
        electrostatic=True,
        fixed_R0=False,
        output_prefix=str(prefix),
    )
    fret.run()

    k2 = _read_single_dat_value(prefix, "k2")
    estatic = _read_single_dat_value(prefix, "Es")
    assert k2 == pytest.approx(0.970927, abs=1e-5)
    assert estatic == pytest.approx(0.468768, abs=1e-5)


@pytest.mark.parametrize("fixed_r0", [False, True])
def test_rotamer_fret_hsp90_matches_fretpredict(tmp_path: Path, fixed_r0: bool) -> None:
    """IMP rotamer FRET matches FRETpredict on the Hsp90 reference case."""
    FRETpredict = pytest.importorskip("FRETpredict")
    MDAnalysis = pytest.importorskip("MDAnalysis")
    pdb = _hsp90_pdb()
    if not pdb.exists():
        pytest.skip("Hsp90 test fixture is not available")

    fp_prefix = tmp_path / "fretpredict"
    imp_prefix = tmp_path / "imp"
    kwargs = {
        "residues": [452, 637],
        "temperature": 293,
        "chains": ["A", "B"],
        "donor": "AlexaFluor 594",
        "acceptor": "AlexaFluor 568",
        "libname_1": "AlexaFluor 594 C1R cutoff30",
        "libname_2": "AlexaFluor 568 C1R cutoff30",
        "electrostatic": True,
        "fixed_R0": fixed_r0,
    }
    if fixed_r0:
        kwargs["r0"] = 5.5

    fretpredict = FRETpredict.FRETpredict(protein=MDAnalysis.Universe(str(pdb)), output_prefix=str(fp_prefix), verbose=False, **kwargs)
    fretpredict.run()
    imp = RotamerFRET(pdb, output_prefix=str(imp_prefix), **kwargs)
    imp.run()

    _assert_rotamer_output_files_match(fp_prefix, imp_prefix)

    weights = [0.5]
    fretpredict.reweight(user_weights=weights)
    imp.reweight(user_weights=weights)
    _assert_rotamer_summary_files_match(fp_prefix, imp_prefix)


def _hsp90_pdb() -> Path:
    """Return the bundled Hsp90 test PDB path."""
    root = _repo_root()
    if (root / "src").exists():
        return root / "src" / "IMP" / "bff" / "cgdye" / "thirdparty" / "FRETpredict" / "tests" / "test_systems" / "Hsp90" / "openHsp90.pdb"
    return root / "modules" / "imp-tricks" / "src" / "IMP" / "bff" / "cgdye" / "thirdparty" / "FRETpredict" / "tests" / "test_systems" / "Hsp90" / "openHsp90.pdb"


def _assert_rotamer_output_files_match(fp_prefix: Path, imp_prefix: Path) -> None:
    """Assert IMP and FRETpredict output arrays match."""
    for suffix in ["Z", "w_s", "k2", "Es", "Ed1", "Ed2"]:
        fp = _read_dat_values(fp_prefix.with_name(fp_prefix.name + f"-{suffix}-452-637.dat"))
        imp = _read_dat_values(imp_prefix.with_name(imp_prefix.name + f"-{suffix}-452-637.dat"))
        assert imp == pytest.approx(fp, rel=1e-6, abs=1e-6)


def _assert_rotamer_summary_files_match(fp_prefix: Path, imp_prefix: Path) -> None:
    """Assert IMP and FRETpredict summary files match."""
    for suffix, key in [("k2", "k2"), ("Es", "Estatic"), ("Ed1", "Edynamic1"), ("Ed2", "Edynamic2")]:
        fp = _read_single_dat_value(fp_prefix, suffix)
        imp = _read_single_dat_value(imp_prefix, suffix)
        assert imp == pytest.approx(fp, rel=1e-6, abs=1e-6)


def _read_dat_values(path: Path) -> list[float]:
    """Read a FRETpredict-style numeric text file."""
    values: list[float] = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        values.extend(float(value) for value in line.split())
    return values


def _read_single_dat_value(prefix: Path, suffix: str) -> float:
    """Read a single numeric value from a FRETpredict-style output file."""
    values = _read_dat_values(prefix.with_name(prefix.name + f"-{suffix}-452-637.dat"))
    assert len(values) == 1
    return values[0]


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
