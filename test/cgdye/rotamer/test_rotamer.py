"""Tests for IMP.bff.representation.rotamer."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner
from IMP.bff.io import fps
from IMP.bff.representation.rotamer.fps import read_rotamer_fps, rotamer_fret_from_fps
from IMP.bff.representation.rotamer.fret import RotamerFRET
from IMP.bff.representation.rotamer.io import load_protein_frames, load_rotamer_library
from IMP.bff.dye.spectra import forster_radius_from_spectra
# The rotamer CLI ships with cgdye now; invoking it through imp-tricks'
# aggregator would make an imp.bff test depend on the layer above it.
from IMP.bff.representation.rotamer.cli import rotamer as _rotamer_cli

def _fixture_pdb(name: str, tmp_path: Path) -> Path:
    """Unpack a bundled (gzip'd) fixture structure into tmp_path and return it."""
    import gzip
    import shutil
    src = Path(__file__).resolve().parent / "data" / f"{name}.gz"
    dst = tmp_path / name
    with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    return dst


def test_calculate_r0_matches_fretpredict_pp11_reference() -> None:
    """R0 calculation matches the pp11 tutorial reference."""
    r0 = forster_radius_from_spectra("AlexaFluor 488", "AlexaFluor 594", 0.684587)
    assert r0 == pytest.approx(5.712982, abs=1e-5)


def test_load_rotamer_library_from_fretpredict_name() -> None:
    """FRETpredict-style library names resolve to the bundled libraries."""
    lib = load_rotamer_library("AlexaFluor 488 C1R cutoff30")
    assert lib["coords"].ndim == 3
    assert lib["coords"].shape[0] > 1
    assert lib["weights"].sum() == pytest.approx(1.0)
    assert lib["metadata"]["mu"]


def test_library_name_cutoff_selects_that_cutoff() -> None:
    """``cutoff10/20/30`` in a name loads that library, not always cutoff 30.

    Every name used to resolve to the cutoff-30 RMF template (33 rotamers for
    Alexa488 C1R), so a cutoff-10 request silently scored 33 rotamers instead
    of 711 -- the pp11 tutorial numbers were off by 0.2 in E. The FRETpredict
    DCD sets in data/rotamer_library are the canonical libraries now.
    """
    from IMP.bff.representation.rotamer.io import resolve_rotamer_library_path
    sizes = {}
    for cutoff in (10, 20, 30):
        name = f"AlexaFluor 488 C1R cutoff{cutoff}"
        path = resolve_rotamer_library_path(name)
        assert path.name == f"A48_C1R_cutoff{cutoff}.dcd"
        lib = load_rotamer_library(name)
        sizes[cutoff] = lib["coords"].shape[0]
        assert lib["weights"].shape == (sizes[cutoff],)
        assert lib["weights"].sum() == pytest.approx(1.0)
        assert lib["resnames"] and len(lib["resnames"]) == len(lib["atom_names"]) == lib["coords"].shape[1]
    assert sizes == {10: 711, 20: 123, 30: 33}
    # the default (no cutoff in the name) is FRETpredict's default, cutoff 30
    assert load_rotamer_library("AlexaFluor 488 C1R")["coords"].shape[0] == 33


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
    pdb = _fixture_pdb("openHsp90.pdb", tmp_path)
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


def test_load_protein_frames_from_rmf_trajectory(tmp_path: Path) -> None:
    """RMF protein trajectories load multiple frames."""
    import IMP
    import IMP.atom
    import IMP.core
    import IMP.algebra
    import IMP.rmf
    import RMF
    from IMP.bff.tools.paths import get_structure_dir

    # Two-frame RMF written on the fly from 148L: frame 1 is frame 0
    # translated by 1 Angstrom along x.
    model = IMP.Model()
    hier = IMP.atom.read_pdb(str(get_structure_dir("148L.pdb")), model,
                             IMP.atom.NonWaterPDBSelector())
    rmf = tmp_path / "two_frames.rmf3"
    fh = RMF.create_rmf_file(str(rmf))
    IMP.rmf.add_hierarchies(fh, [hier])
    IMP.rmf.save_frame(fh, "0")
    for a in IMP.atom.get_by_type(hier, IMP.atom.ATOM_TYPE):
        xyz = IMP.core.XYZ(a)
        xyz.set_coordinates(xyz.get_coordinates() + IMP.algebra.Vector3D(1.0, 0.0, 0.0))
    IMP.rmf.save_frame(fh, "1")
    del fh

    frames = load_protein_frames(rmf, max_frames=2)
    assert len(frames) == 2
    assert frames[0]["coords"].ndim == 2
    assert frames[0]["residue_indices"]
    assert frames[1]["coords"].shape == frames[0]["coords"].shape
    assert np.allclose(frames[1]["coords"] - frames[0]["coords"], [1.0, 0.0, 0.0])


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
