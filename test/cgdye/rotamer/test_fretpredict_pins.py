"""FRETpredict parity from bundled pins -- runs without FRETpredict installed.

The reference numbers in test/references/cgdye_fretpredict_pins.json were
recorded once by running FRETpredict itself (provenance in the file's
``_note``); the fixtures are FRETpredict's Hsp90 structure and the first 20
frames of its pp11 tutorial trajectory. This is the hard parity gate of
PRD-107: RotamerFRET must reproduce FRETpredict's Es/Ed1/Ed2/<k2>/Z (and R0
when calculated) to FRETpredict's own tolerance (1e-5; 2e-5 for pp11 whose
frames went through PDB coordinate precision).
"""

from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path

import numpy as np
import pytest

from IMP.bff.representation.rotamer import RotamerFRET

_HERE = Path(__file__).resolve().parent
_PINS = _HERE.parents[1] / "references" / "cgdye_fretpredict_pins.json"


def _pins():
    with open(_PINS) as fh:
        return json.load(fh)


def _fixture(name: str, tmp_path: Path) -> Path:
    src = _HERE / "data" / f"{name}.gz"
    dst = tmp_path / name
    with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    return dst


def _run(system: dict, case: dict, tmp_path: Path) -> RotamerFRET:
    pdb = _fixture(system["structure"], tmp_path)
    residues = case.get("residues", system["residues"])
    chains = case.get("chains", system["chains"])
    libs = case.get("libraries", system["libraries"])
    fret = RotamerFRET(
        pdb, residues, chains=chains, libname_1=libs[0], libname_2=libs[1],
        output_prefix=str(tmp_path / "res"), **case["kwargs"])
    fret.trajectory_analysis()
    return fret


def _cases(system_key: str):
    pins = _pins()
    return [pytest.param(system_key, name, id=f"{system_key}-{name}")
            for name in pins[system_key]["cases"]]


@pytest.mark.parametrize("system_key,case_name", _cases("hsp90") + _cases("pp11"))
def test_rotamer_fret_matches_fretpredict(system_key: str, case_name: str, tmp_path: Path) -> None:
    pins = _pins()
    system = pins[system_key]
    case = system["cases"][case_name]
    ref = case["fretpredict"]
    tol = system["tolerance"]
    fret = _run(system, case, tmp_path)

    got = {
        "Es": np.asarray(fret.estatic_values, dtype=float),
        "Ed1": np.asarray(fret.edynamic1_values, dtype=float),
        "Ed2": np.asarray(fret.edynamic2_values, dtype=float),
        "k2": np.asarray(fret.k2_values, dtype=float),
    }
    for key, values in got.items():
        expected = np.atleast_1d(np.asarray(ref[key], dtype=float))
        assert values.shape == expected.shape, key
        np.testing.assert_allclose(values, expected, rtol=0, atol=system.get(f"tolerance_{key}", tol), err_msg=key)

    if "Z" in ref:
        z = np.asarray(fret.z_values, dtype=float).reshape(-1)
        expected = np.asarray(ref["Z"], dtype=float).reshape(-1)
        np.testing.assert_allclose(z, expected, rtol=0, atol=system.get("tolerance_Z", tol), err_msg="Z")
    if "R0_nm" in ref:
        assert float(np.ravel(fret.r0)[0]) == pytest.approx(ref["R0_nm"], abs=1e-4)


def test_hsp90_reweight_with_user_weights_matches_fretpredict(tmp_path: Path) -> None:
    """A single frame reweighted with one user weight leaves Es unchanged (as in FRETpredict's test)."""
    pins = _pins()
    system = pins["hsp90"]
    case = system["cases"]["calcR0_T293_DH"]
    fret = _run(system, case, tmp_path)
    fret.save()
    fret.reweight(user_weights=[0.5])
    ref = case["fretpredict"]
    rows = np.loadtxt(tmp_path / "res-data-452-637.dat")  # k2, Es, Ed1, Ed2 x (avg, sd, se)
    assert rows[1, 0] == pytest.approx(ref["Estatic_reweighted_user_0.5"], abs=system["tolerance"])
    assert rows[1, 0] == pytest.approx(ref["Estatic_pkl"], abs=system["tolerance"])


def test_hsp90_writes_fretpredict_style_outputs(tmp_path: Path) -> None:
    pins = _pins()
    system = pins["hsp90"]
    case = system["cases"]["fixedR0_5.5_T293_DH"]
    fret = _run(system, case, tmp_path)
    fret.save()
    for suffix in ("Z", "w_s", "k2", "Es", "Ed1", "Ed2"):
        assert (tmp_path / f"res-{suffix}-452-637.dat").exists(), suffix
    es = float(np.loadtxt(tmp_path / "res-Es-452-637.dat").ravel()[0])
    assert es == pytest.approx(case["fretpredict"]["Es"], abs=system["tolerance"])


# IMP runs every .py under test/ as a standalone script, and a file of bare
# pytest functions would import cleanly and exit 0 -- reporting success without
# running a single assertion. Hand the file to pytest explicitly so a failure
# here is a failure in ctest.
if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
