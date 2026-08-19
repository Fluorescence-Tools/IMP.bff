"""AV ↔ rotamer-ensemble table (PRD-108 stage 2): drift pins + loose sanity bounds.

The table okf/validation/av_vs_rotamer.md is authoritative; the differences
between an AV1 cloud and a screened rotamer library are two physical models
disagreeing, recorded here, not gated. This test pins the recorded numbers
(exact -- both computations are deterministic) and asserts only that the two
models describe the same label: mean positions within 20 Å, ⟨R_DA⟩ within
20 Å, no rotamer with weight > 1e-3 clashing with the protein.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

from IMP.bff.representation.compare import compare_av_and_rotamer_pairs, compare_av_and_rotamer_positions
from IMP.bff.representation.compare import R0_A488_A594, hgbp1_case, t4l_case

_PINS = Path(__file__).resolve().parents[2] / "references" / "cgdye_av_vs_rotamer_pins.json"


def _pins():
    with open(_PINS) as fh:
        return json.load(fh)


@pytest.mark.parametrize("case_name,case", [("hgbp1_cutoff30", hgbp1_case), ("t4l_cutoff30", t4l_case)])
def test_av_vs_rotamer_pins_and_sanity(case_name, case):
    pins = _pins()
    settings = pins["settings"]
    title, pdb, positions, libraries, pairs, experimental = case(30)
    per_pos = compare_av_and_rotamer_positions(pdb, positions, libraries, n_samples=settings["n_samples"], temperature=settings["temperature"])
    rows = compare_av_and_rotamer_pairs(per_pos, pairs, R0_A488_A594, experimental=experimental, n_samples=settings["n_samples"])
    ref = pins[case_name]

    # drift pins (deterministic: lattice AV + fixed-seed sampling; full pair matrix)
    for name, rec in ref["positions"].items():
        got = per_pos[name]
        assert got["n_rotamers"] == rec["n_rotamers"] and got["av_n_points"] == rec["av_n_points"]
        assert got["d_mean_position"] == pytest.approx(rec["d_mean_position"], abs=1e-6)
        assert got["partition"] == pytest.approx(rec["partition"], abs=1e-9)
    ref_rows = {r["pair"]: r for r in ref["pairs"]}
    for row in rows:
        rec = ref_rows[row["pair"]]
        for key in ("Rmp_av", "Rmp_rot", "RDAMean_av", "RDAMean_rot", "RDAMeanE_av", "RDAMeanE_rot", "sigma_av", "sigma_rot", "kappa2_rot"):
            assert row[key] == pytest.approx(rec[key], abs=1e-6), (row["pair"], key)

    # sanity: two models of the same label (skipping sites the rotamer model
    # finds buried, Z < 0.05 -- FRETpredict's uniform fallback)
    from IMP.bff.representation.compare import Z_CUTOFF
    for name, got in per_pos.items():
        assert got["av_n_points"] > 0 and got["n_rotamers"] > 0
        if got["partition"] < Z_CUTOFF:
            continue
        assert got["d_mean_position"] < 20.0, name
        # the screening is FRETpredict's soft LJ (sigma_scaling 0.5): rotamers
        # may sit ~1.4 A from protein atoms with full weight (recorded in the
        # table as "weight within 2.5 A"); a rotamer *inside* the protein
        # (< 1.0 A) would be a scoring defect
        assert got["min_heavy_atom_distance"] > 1.0, name
        assert got["min_heavy_atom_distance"] == pytest.approx(ref["positions"][name]["min_heavy_atom_distance"], abs=1e-6)
        assert got["interpenetration_weight_2.5A"] == pytest.approx(ref["positions"][name]["interpenetration_weight_2.5A"], abs=1e-9)
    for row in rows:
        if not row["rotamer_valid"]:
            continue
        assert abs(row["RDAMean_av"] - row["RDAMean_rot"]) < 20.0, row["pair"]
        assert abs(row["RDAMeanE_av"] - row["RDAMeanE_rot"]) < 20.0, row["pair"]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
