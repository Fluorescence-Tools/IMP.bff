"""RotamerEnsemble (fps R1): a screened library that is an AccessibleVolume (PRD-108 stage 0)."""

from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path

import numpy as np
import pytest

from IMP.bff.cgdye.rotamer.ensemble import RotamerEnsemble, rotamer_ensembles_from_fps
from IMP.bff.cgdye.rotamer.fret import RotamerFRET
from IMP.bff.fret.av import AccessibleVolume
from IMP.bff.fret.distance import av_pair_statistics, histogram_rda, mean_fret_distance

_HERE = Path(__file__).resolve().parent
_PINS = _HERE.parents[1] / "references" / "cgdye_fretpredict_pins.json"


@pytest.fixture(scope="module")
def hsp90(tmp_path_factory):
    dst = tmp_path_factory.mktemp("hsp90") / "openHsp90.pdb"
    with gzip.open(_HERE / "data" / "openHsp90.pdb.gz", "rb") as fin, open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    return dst


@pytest.fixture(scope="module")
def pair(hsp90):
    d = RotamerEnsemble.from_site(hsp90, "A", 452, "AlexaFluor 594 C1R cutoff30", temperature=293, electrostatic=True)
    a = RotamerEnsemble.from_site(hsp90, "B", 637, "AlexaFluor 568 C1R cutoff30", temperature=293, electrostatic=True)
    return d, a


def test_is_an_accessible_volume_with_full_library(pair):
    d, a = pair
    assert isinstance(d, AccessibleVolume)
    assert d.n_rotamers == 37 and a.n_rotamers == 7      # cutoff-30 libraries
    assert d.points.shape == (37, 4) and d.mu.shape == (37, 3) and d.atoms.shape[0] == 37
    assert d.atoms.shape[1] == len(d.atom_names)
    assert d.weights.sum() == pytest.approx(1.0)
    np.testing.assert_allclose(np.linalg.norm(d.mu, axis=1), 1.0)
    assert d.params["simulation_type"] == "R1" and d.library == "AlexaFluor 594 C1R cutoff30"
    assert d.chain == "A" and d.residue == 452
    # centres sit within a linker length of the attachment CA
    assert np.linalg.norm(d.mean_position - d.attachment_point) < 25.0
    assert np.all(np.linalg.norm(d.centres - d.attachment_point, axis=1) < 30.0)
    assert d.has_volume


def test_av_helpers_accept_it(pair):
    d, a = pair
    rmp, rda, rda_e, sigma = av_pair_statistics(d, a, forster_radius=50.4)
    assert 30 < rmp < 80 and 30 < rda < 80 and sigma > 0
    assert mean_fret_distance(d, a, forster_radius=50.4) > 0
    hist = histogram_rda(d, a)
    assert hist is not None


def test_matches_fretpredict_and_rotamer_fret(pair, hsp90, tmp_path):
    d, a = pair
    pins = json.load(open(_PINS))["hsp90"]["cases"]["calcR0_T293_DH"]["fretpredict"]
    eff = d.fret_efficiencies(a, donor="AlexaFluor 594", acceptor="AlexaFluor 568")
    assert eff["static"] == pytest.approx(pins["Es"], abs=1e-5)
    assert eff["dynamic1"] == pytest.approx(pins["Ed1"], abs=1e-5)
    assert eff["dynamic2"] == pytest.approx(pins["Ed2"], abs=1e-5)
    assert eff["kappa2_avg"] == pytest.approx(pins["k2"], abs=1e-5)
    assert eff["forster_radius_nm"] == pytest.approx(pins["R0_nm"], abs=1e-4)
    assert d.partition == pytest.approx(pins["Z"][0], abs=1e-5) and a.partition == pytest.approx(pins["Z"][1], abs=1e-5)
    # the same numbers RotamerFRET (built on the ensembles) reports
    fret = RotamerFRET(hsp90, [452, 637], chains=["A", "B"], temperature=293, electrostatic=True,
                       donor="AlexaFluor 594", acceptor="AlexaFluor 568",
                       libname_1="AlexaFluor 594 C1R cutoff30", libname_2="AlexaFluor 568 C1R cutoff30",
                       output_prefix=str(tmp_path / "x"))
    fret.trajectory_analysis()
    assert float(fret.estatic_values[0]) == pytest.approx(eff["static"], abs=1e-9)


def test_pair_distribution_is_the_full_matrix(pair):
    d, a = pair
    out = d.pair_distribution(a, forster_radius=50.4, tau0=4.0)
    assert out["R"].shape == (37, 7) and out["kappa2"].shape == (37, 7) and out["k_fret"].shape == (37, 7)
    assert out["weight"].sum() == pytest.approx(1.0)
    np.testing.assert_allclose(np.outer(d.weights, a.weights), out["weight"], atol=1e-12)
    # kappa2 against a plain AV is isotropic
    plain = AccessibleVolume(points=a.points.copy(), density=np.zeros((0, 0, 0), np.float32),
                             grid_origin=a.attachment_point, grid_step=0.0, grid_shape=(0, 0, 0),
                             attachment_point=a.attachment_point)
    assert np.all(d.pair_geometry(plain)["kappa2"] == 2.0 / 3.0)


def test_from_fps_positions(hsp90, tmp_path):
    fps = tmp_path / "hsp90.fps.json"
    fps.write_text(json.dumps({
        "Positions": {
            "d1": {"chain_identifier": "A", "residue_seq_number": 452, "atom_name": "CA",
                   "rotamer_library": "AlexaFluor 594 C1R cutoff30"},
            "a1": {"chain_identifier": "B", "residue_seq_number": 637, "atom_name": "CA"},
        },
        "Distances": {"d1_a1": {"position1_name": "d1", "position2_name": "a1"}},
    }))
    ens = rotamer_ensembles_from_fps(fps, hsp90, temperature=293, electrostatic=True)
    assert set(ens) == {"d1"}                             # a1 names no library
    ens = rotamer_ensembles_from_fps(fps, hsp90, library_map={"a1": "AlexaFluor 568 C1R cutoff30"},
                                     temperature=293, electrostatic=True)
    assert set(ens) == {"d1", "a1"} and ens["a1"].n_rotamers == 7 and ens["a1"].position_name == "a1"


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
