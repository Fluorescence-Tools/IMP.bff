"""RotamerEnsemble (fps R1): a screened library, and a sibling of the AV.

It *was* a subclass of :class:`AccessibleVolume` (PRD-108 stage 0), which is why
distance code worked for rotamers: by inheritance rather than by design. It then
had to carry a grid it does not have, filling ``density``, ``grid_step`` and
``grid_shape`` with empty placeholders that no consumer ever read.

PRD-113 stage 3c cuts that edge. Both are now
:class:`~IMP.bff.representation.States` -- positions, weights and orientations --
which is the surface every consumer was actually using, and which an MD or
coarse-grained representation can supply just as well.
"""

from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path

import numpy as np
import pytest

from IMP.bff import RotamerEnsemble, rotamer_ensembles_from_fps
from IMP.bff import RotamerFRET
from IMP.bff import AccessibleVolume, States
from IMP.bff import av_pair_statistics, histogram_rda
from IMP.bff import states_mean_fret_distance as mean_fret_distance

_HERE = Path(__file__).resolve().parent
_PINS = _HERE.parents[1] / "references" / "cgdye_fretpredict_pins.json"


_HSP90_PATH: dict = {}


@pytest.fixture(scope="module")
def hsp90(tmp_path_factory):
    dst = tmp_path_factory.mktemp("hsp90") / "openHsp90.pdb"
    with gzip.open(_HERE / "data" / "openHsp90.pdb.gz", "rb") as fin, open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    _HSP90_PATH["path"] = dst
    return dst


@pytest.fixture(scope="module")
def pair(hsp90):
    d = RotamerEnsemble.from_site(hsp90, "A", 452, "AlexaFluor 594 C1R cutoff30", temperature=293, electrostatic=True)
    a = RotamerEnsemble.from_site(hsp90, "B", 637, "AlexaFluor 568 C1R cutoff30", temperature=293, electrostatic=True)
    return d, a


def test_is_a_sibling_of_the_av_not_a_subclass(pair):
    """The shared surface is States; the grid is not part of it."""
    d, a = pair
    assert isinstance(d, States)
    assert not isinstance(d, AccessibleVolume), (
        "a rotamer library is not an accessible volume: it has no grid, and "
        "inheriting one forced it to fake density, grid_step and grid_shape")
    for grid_field in ("density", "grid_step", "grid_shape", "grid_origin"):
        assert not hasattr(d, grid_field)


def test_the_full_library_is_screened(pair):
    d, a = pair
    assert d.n_rotamers == 37 and a.n_rotamers == 7      # cutoff-30 libraries
    assert d.points.shape == (37, 4) and d.mu.shape == (37, 3) and d.atoms.shape[0] == 37
    assert d.atoms.shape[1] == len(d.atom_names)
    assert d.weights.sum() == pytest.approx(1.0)
    np.testing.assert_allclose(np.linalg.norm(d.mu, axis=1), 1.0)
    # `mu` and the representation-agnostic `orientations` are one array -- `mu`
    # is a property over `States.orientations`, not a second copy -- so a
    # consumer written against States sees the dipoles a rotamer library has and
    # an accessible volume has not. Each read is a fresh view over the C++
    # buffer, so the check is that a write through one is seen through the
    # other, which is what "one array" has to mean here.
    np.testing.assert_array_equal(d.orientations, d.mu)
    assert d.has_orientations
    probe = -d.mu
    d.mu = probe
    np.testing.assert_array_equal(d.orientations, probe)
    d.mu = -probe                                    # put the fixture back
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
    plain = AccessibleVolume(points=a.points.copy(),
                             grid_origin=a.attachment_point, grid_step=0.0,
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


def test_r1_positions_round_trip_and_docking_filter(pair, tmp_path):
    """R1 entries write, validate, read back; the docking filter drops them."""
    from IMP.bff import (
        distances_from_ensembles, rotamer_ensemble_payload, write_rotamer_fps)
    import IMP.bff as fps_schema
    from IMP.bff import fps_positions_for_docking, read_fps_json

    d, a = pair
    positions = {"d1": rotamer_ensemble_payload(d), "a1": rotamer_ensemble_payload(a)}
    assert positions["d1"]["simulation_type"] == "R1"
    assert positions["d1"]["rotamer_library"] == "AlexaFluor 594 C1R cutoff30"
    for dtype in ("RDAMean", "RDAMeanE", "Rmp"):
        distances = distances_from_ensembles({"d1": d, "a1": a}, [("d1", "a1")], 50.4, distance_type=dtype)
        entry = distances["d1_a1"]
        assert entry["distance_type"] == dtype and 20 < entry["distance"] < 90
    distances = distances_from_ensembles({"d1": d, "a1": a}, [("d1", "a1")], 50.4)
    # <R_DA>_E (kappa2 = 2/3, the fps.json convention) is what av_pair_statistics gives
    _rmp, _rda, rda_e_av, _sig = av_pair_statistics(d, a, forster_radius=50.4, n_samples=200000)
    assert distances["d1_a1"]["distance"] == pytest.approx(rda_e_av, abs=0.3)
    # with the ensembles' dipoles it reproduces the orientation-resolved static efficiency
    dip = distances_from_ensembles({"d1": d, "a1": a}, [("d1", "a1")], 50.4, kappa2="dipoles")
    eff = d.fret_efficiencies(a, forster_radius=50.4)
    r_e = 50.4 * (1.0 / eff["static"] - 1.0) ** (1.0 / 6.0)
    assert dip["d1_a1"]["distance"] == pytest.approx(r_e, rel=1e-9)

    out = tmp_path / "rotamer.fps.json"
    write_rotamer_fps(out, positions, distances)                    # validated on write
    doc = read_fps_json(str(out))
    p = json.loads(doc.positions)
    dist = json.loads(doc.distances)
    assert set(p) == {"d1", "a1"} and set(dist) == {"d1_a1"}
    errors = fps_schema.fps_schema_validate(json.dumps(
        {"Positions": p, "Distances": dist})).errors
    assert not errors
    # what the C++ scorer may see: no R1 positions, no dangling distances
    kept_doc = fps_positions_for_docking(json.dumps(p), json.dumps(dist))
    kept_p = json.loads(kept_doc.positions)
    kept_d = json.loads(kept_doc.distances)
    assert kept_p == {} and kept_d == {}
    # merge into an AV-style file keeps the AV positions and adds R1
    av_file = tmp_path / "av.fps.json"
    av_file.write_text(json.dumps({
        "Positions": {"p1": {"chain_identifier": "A", "residue_seq_number": 10, "atom_name": "CB",
                             "simulation_type": "AV1", "linker_length": 20.0, "linker_width": 4.5, "radius1": 3.5}},
        "Distances": {}}))
    write_rotamer_fps(out, positions, distances, merge_into=av_file)
    doc = read_fps_json(str(out))
    p = json.loads(doc.positions)
    dist = json.loads(doc.distances)
    assert set(p) == {"p1", "d1", "a1"}
    kept_doc = fps_positions_for_docking(json.dumps(p), json.dumps(dist))
    kept_p = json.loads(kept_doc.positions)
    kept_d = json.loads(kept_doc.distances)
    assert set(kept_p) == {"p1"} and kept_d == {}
    # the ensembles come back from the file
    ens = rotamer_ensembles_from_fps(out, hsp90_path(pair), temperature=293, electrostatic=True)
    assert set(ens) == {"d1", "a1"} and ens["d1"].n_rotamers == 37


def test_r1_requires_a_library():
    import IMP.bff as fps_schema
    report = fps_schema.validate_position(json.dumps(
        {"chain_identifier": "A", "residue_seq_number": 1,
         "atom_name": "CA", "simulation_type": "R1"}), "x")
    errors = report.errors
    assert any("rotamer_library" in e for e in errors)
    report = fps_schema.validate_position(json.dumps(
        {"chain_identifier": "A", "residue_seq_number": 1, "atom_name": "CA",
         "simulation_type": "R1", "rotamer_library": "AlexaFluor 488 C1R cutoff30",
         "linker_length": 20.0}), "x")
    assert not report.errors and any("AV parameter" in w for w in report.warnings)


def hsp90_path(pair):
    """The PDB the module fixture unpacked (kept in the ensembles' params? no: re-derive)."""
    return _HSP90_PATH["path"]


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
