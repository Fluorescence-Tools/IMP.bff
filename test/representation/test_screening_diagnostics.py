"""FPS screening diagnostics: sigma counts, InvalidR, the pair table and the
reference-frame fit (PRD-121 G3).

The reference-frame tests are synthetic on purpose. A frame fitted onto a
structure has an answer only if the frame was *made* from that structure by a
transform we chose, and then the answer is exact: the fit must undo the
transform, and the RMSD must be the distortion we put in and nothing else. A
plausible-looking number off a real file would not distinguish a correct
Kabsch from a nearly correct one.
"""

import json
import math
from pathlib import Path

import pytest

import IMP
import IMP.algebra
import IMP.atom
import IMP.core
import IMP.bff

REPO = Path(__file__).resolve().parents[2]
T4L = REPO / "examples" / "structure" / "T4L"
PDB = T4L / "3GUN.pdb"
FPS_JSON = T4L / "fret.fps.json"


@pytest.fixture(scope="module")
def structure():
    IMP.set_log_level(IMP.SILENT)
    m = IMP.Model()
    h = IMP.atom.read_pdb(str(PDB), m,
                          IMP.atom.NonWaterNonHydrogenPDBSelector())
    return m, h


def _ca_atoms(hierarchy, residues, chain="A"):
    """(residue number, coordinate) of the CA atoms that exist."""
    out = []
    for r in residues:
        ps = IMP.atom.Selection(hierarchy, chain_id=chain, residue_index=r,
                                atom_type=IMP.atom.AT_CA
                                ).get_selected_particles()
        if ps:
            out.append((r, IMP.core.XYZ(ps[0]).get_coordinates()))
    return out


def _frame(atoms, transform=None, chain="A"):
    """A ReferenceAtoms list from (residue, coordinate) pairs."""
    frame = IMP.bff.ReferenceAtoms()
    for residue, coordinate in atoms:
        v = transform * coordinate if transform is not None else coordinate
        a = IMP.bff.ReferenceAtom()
        a.chain_identifier = chain
        a.residue_seq_number = residue
        a.atom_name = "CA"
        a.x, a.y, a.z = v[0], v[1], v[2]
        frame.append(a)
    return frame


def _frame_json(atoms, transform=None, chain="A"):
    rows = []
    for residue, coordinate in atoms:
        v = transform * coordinate if transform is not None else coordinate
        rows.append({"chain_identifier": chain, "residue_seq_number": residue,
                     "atom_name": "CA", "x": v[0], "y": v[1], "z": v[2]})
    return rows


# --------------------------------------------------------------------------
# the reference-frame kernel
# --------------------------------------------------------------------------
def test_kabsch_fit_undoes_a_known_rigid_transform(structure):
    """A frame is the structure moved by T; the fit must be exactly T^-1."""
    _, h = structure
    atoms = _ca_atoms(h, range(10, 21))
    assert len(atoms) == 11

    rotation = IMP.algebra.get_rotation_about_axis(
        IMP.algebra.Vector3D(0.3, -0.5, 0.81), 0.7)
    t = IMP.algebra.Transformation3D(rotation,
                                     IMP.algebra.Vector3D(12.0, -3.0, 5.5))

    point = IMP.algebra.Vector3D(1.0, 2.0, 3.0)
    fit = IMP.bff.fit_reference_atoms(_frame(atoms, t), h, point)

    assert fit.fitted
    assert fit.n_atoms == 11
    # A rigid frame fits its own structure exactly. The tolerance is
    # arithmetic noise, not a modelling allowance: the fit goes through an SVD
    # and a rotation-matrix-to-quaternion extraction, and over the rotations
    # tried here the residual came out between 5e-15 and 2e-7 A on coordinates
    # of tens of angstrom.
    assert fit.rmsd < 1e-6
    expected = t.get_inverse() * point
    assert (IMP.algebra.Vector3D(*fit.coordinates) - expected
            ).get_magnitude() < 1e-6


def test_reference_rmsd_is_the_distortion_that_was_put_in(structure):
    """Scale the frame about its centroid: the best rigid fit cannot follow.

    For a frame scaled by s about its own centroid the optimal rotation is the
    identity and the residual of atom i is (s - 1) * (r_i - centroid), so the
    fit RMSD is |s - 1| times the frame's root-mean-square radius. That is an
    exact number, so it is what is asserted.
    """
    _, h = structure
    atoms = _ca_atoms(h, range(30, 45))
    coordinates = [c for _, c in atoms]
    n = len(coordinates)
    assert n >= 8

    centroid = IMP.algebra.Vector3D(
        sum(c[0] for c in coordinates) / n,
        sum(c[1] for c in coordinates) / n,
        sum(c[2] for c in coordinates) / n)
    rms_radius = math.sqrt(
        sum((c - centroid).get_squared_magnitude() for c in coordinates) / n)

    s = 1.05
    frame = IMP.bff.ReferenceAtoms()
    for (residue, c) in atoms:
        v = centroid + (c - centroid) * s
        a = IMP.bff.ReferenceAtom()
        a.chain_identifier, a.residue_seq_number, a.atom_name = "A", residue, "CA"
        a.x, a.y, a.z = v[0], v[1], v[2]
        frame.append(a)

    fit = IMP.bff.fit_reference_atoms(frame, h)
    assert fit.fitted
    assert fit.rmsd == pytest.approx(abs(s - 1.0) * rms_radius, rel=1e-9)


def test_a_frame_of_two_atoms_is_refused(structure):
    """Two atoms leave the rotation about their axis free, so the point they
    would transport is arbitrary rather than uncertain."""
    _, h = structure
    fit = IMP.bff.fit_reference_atoms(_frame(_ca_atoms(h, [10, 11])), h,
                                      IMP.algebra.Vector3D(1, 2, 3))
    assert not fit.fitted
    assert fit.n_atoms == 2
    # the point comes back as given, not silently moved
    assert (IMP.algebra.Vector3D(*fit.coordinates)
            - IMP.algebra.Vector3D(1, 2, 3)).get_magnitude() == 0.0


def test_atoms_the_structure_does_not_have_are_skipped_not_guessed(structure):
    _, h = structure
    frame = _frame(_ca_atoms(h, range(10, 21)))
    ghost = IMP.bff.ReferenceAtom()
    ghost.chain_identifier, ghost.residue_seq_number = "A", 99999
    ghost.atom_name, ghost.x, ghost.y, ghost.z = "CA", 0.0, 0.0, 0.0
    frame.append(ghost)
    fit = IMP.bff.fit_reference_atoms(frame, h)
    assert fit.n_atoms == 11
    assert fit.rmsd < 1e-6


def test_an_empty_chain_identifier_means_any_chain(structure):
    """The shipped fps.json positions leave `chain_identifier` empty, so a
    reference atom must be allowed to as well; it resolves when the residue
    and atom name pick out exactly one atom across the chains."""
    _, h = structure
    fit = IMP.bff.fit_reference_atoms(
        _frame(_ca_atoms(h, range(10, 21)), chain=""), h)
    assert fit.fitted
    assert fit.n_atoms == 11
    assert fit.rmsd < 1e-6


def test_fit_reference_positions_reads_the_fps_json_spelling(structure):
    _, h = structure
    atoms = _ca_atoms(h, range(10, 21))
    rotation = IMP.algebra.get_rotation_about_axis(
        IMP.algebra.Vector3D(0.0, 0.0, 1.0), 1.1)
    t = IMP.algebra.Transformation3D(rotation,
                                     IMP.algebra.Vector3D(-4.0, 7.0, 1.0))
    fixed = IMP.algebra.Vector3D(5.0, -6.0, 7.0)

    positions = {
        "plain_av": {"residue_seq_number": 132, "atom_name": "CB"},
        "fixed_site": {"simulation_type": "XYZ",
                       "x": fixed[0], "y": fixed[1], "z": fixed[2],
                       "reference_atoms": _frame_json(atoms, t)},
    }
    fits = IMP.bff.fit_reference_positions(json.dumps(positions), h)
    # only positions that declare a frame are reported
    assert [f.position for f in fits] == ["fixed_site"]
    assert fits[0].n_atoms == 11
    assert fits[0].rmsd < 1e-6
    expected = t.get_inverse() * fixed
    assert (IMP.algebra.Vector3D(*fits[0].coordinates)
            - expected).get_magnitude() < 1e-6
    assert IMP.bff.reference_rmsd(fits) < 1e-6


def test_reference_rmsd_pools_over_positions_it_does_not_average(structure):
    """FPS divides one sum of squares by one total count, so a twenty-atom
    frame and a four-atom frame do not weigh the same."""
    _, h = structure
    fits = IMP.bff.ReferenceFits()
    big, small = IMP.bff.ReferenceFit(), IMP.bff.ReferenceFit()
    big.fitted, big.n_atoms, big.squared_deviation = True, 20, 20 * 4.0
    small.fitted, small.n_atoms, small.squared_deviation = True, 4, 4 * 100.0
    fits.append(big)
    fits.append(small)
    pooled = math.sqrt((20 * 4.0 + 4 * 100.0) / 24.0)
    assert IMP.bff.reference_rmsd(fits) == pytest.approx(pooled)
    # and it is not the mean of the two RMSDs (2 and 10)
    assert abs(IMP.bff.reference_rmsd(fits) - 6.0) > 1.0


def test_reference_rmsd_of_nothing_is_zero_and_of_an_unfittable_frame_is_nan():
    assert IMP.bff.reference_rmsd(IMP.bff.ReferenceFits()) == 0.0
    fits = IMP.bff.ReferenceFits()
    fits.append(IMP.bff.ReferenceFit())          # fitted == False
    assert math.isnan(IMP.bff.reference_rmsd(fits))


# --------------------------------------------------------------------------
# the schema
# --------------------------------------------------------------------------
def test_schema_declares_reference_atoms():
    # 1.3 added the `ATOM` position type (PRD-121 G4); `reference_atoms` came
    # in at 1.2 and every step is additive, so this asserts the field rather
    # than the number it arrived at.
    assert IMP.bff.fps_schema_version() >= "1.2"
    schema = json.loads(IMP.bff.fps_json_schema())
    position = schema["properties"]["Positions"]["additionalProperties"]
    assert "reference_atoms" in position["properties"]
    # An optional field needs a real default: an empty one means *required*.
    assert position["properties"]["reference_atoms"]["default"] == []
    assert "reference_atoms" not in position.get("required", [])


def test_schema_checks_what_a_reference_atom_has_to_carry():
    good = {"simulation_type": "XYZ", "x": 1.0, "y": 2.0, "z": 3.0,
            "reference_atoms": [
                {"chain_identifier": "A", "residue_seq_number": 10,
                 "atom_name": "CA", "x": 0.0, "y": 0.0, "z": 0.0},
                {"atom_serial": 42, "x": 1.0, "y": 0.0, "z": 0.0},
                {"chain_identifier": "A", "residue_seq_number": 12,
                 "atom_name": "CA", "x": 0.0, "y": 1.0, "z": 0.0}]}
    report = IMP.bff.validate_position(json.dumps(good), "p")
    assert not report.errors, report.errors
    assert not report.warnings, report.warnings

    # no coordinates, and no way to find the atom again
    bad = dict(good, reference_atoms=[{"atom_name": "CA"}])
    errors = IMP.bff.validate_position(json.dumps(bad), "p").errors
    assert any("'x'" in e for e in errors), errors
    assert any("atom_serial" in e for e in errors), errors

    # two atoms are well formed but cannot fix a rotation: a warning, not an
    # error -- the file is fine, the frame is not usable
    short = dict(good, reference_atoms=good["reference_atoms"][:2])
    report = IMP.bff.validate_position(json.dumps(short), "p")
    assert not report.errors
    assert any("three are needed" in w for w in report.warnings)


def test_shipped_fps_json_still_validates_under_1_2():
    payload = json.loads(FPS_JSON.read_text())
    for name, position in payload["Positions"].items():
        errors = IMP.bff.validate_position(json.dumps(position), name).errors
        assert not errors, (name, errors)


# --------------------------------------------------------------------------
# screening
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def screened():
    IMP.set_log_level(IMP.SILENT)
    return list(IMP.bff.screen_structures(
        [str(PDB)], str(FPS_JSON), "chi2_C2_33p", "", False))


def test_screen_carries_the_pair_table(screened):
    entry = screened[0]
    assert len(entry.pairs) == 33
    assert all(p.name for p in entry.pairs)
    assert all(math.isfinite(p.distance_model) for p in entry.pairs)


def test_sigma_counts_are_the_asymmetric_ones(screened):
    """Recount them in Python, from the table, with FPS's rule: a model that
    is too long is judged against error_pos, one that is too short against
    error_neg."""
    entry = screened[0]
    expected = [0, 0, 0]
    for p in entry.pairs:
        dr = p.distance_model - p.distance_exp
        err = p.error_pos if dr > 0 else p.error_neg
        if err <= 0:
            continue
        for k in (1, 2, 3):
            if abs(dr) > k * err:
                expected[k - 1] += 1
    assert [entry.sigma1, entry.sigma2, entry.sigma3] == expected
    # nested, not exclusive
    assert entry.sigma1 >= entry.sigma2 >= entry.sigma3
    # this fixture is not a perfect fit, so the counts are not all zero --
    # a test that only saw zeros would not have seen the counting
    assert entry.sigma1 > 0


def test_chi2_r_is_the_mean_over_scored_pairs(screened):
    entry = screened[0]
    valid = [p.get_chi2() for p in entry.pairs
             if math.isfinite(p.distance_model)]
    assert entry.chi2_r == pytest.approx(sum(valid) / len(valid), rel=1e-9)
    assert entry.invalid_r == len(entry.pairs) - len(valid)
    # score is half the chi2 sum (a Gaussian -log L) plus the clash term,
    # which is zero for a single rigid body
    assert entry.score == pytest.approx(0.5 * sum(valid), rel=1e-6)


def test_an_empty_volume_is_counted_not_scored(tmp_path):
    """InvalidR: a distance with no model value at one end.

    Built by making one position's dye radius larger than anything the linker
    can reach -- the volume comes out empty and the two distances that touch
    it have no model value, while the third still scores.
    """
    payload = json.loads(FPS_JSON.read_text())
    chosen = ["19-119_C2", "19-132_C2", "44-119_C2"]
    for name in chosen:
        assert name in payload["Distances"], name
    keep = {payload["Distances"][n][k]
            for n in chosen for k in ("position1_name", "position2_name")}
    payload["Distances"] = {n: payload["Distances"][n] for n in chosen}
    payload["Positions"] = {n: payload["Positions"][n] for n in keep}
    payload["χ²"] = {"small": {"distances": chosen}}
    # 19D is an end of two of the three distances
    payload["Positions"]["19D"]["radius1"] = 40.0
    path = tmp_path / "empty_volume.fps.json"
    path.write_text(json.dumps(payload))

    entry = list(IMP.bff.screen_structures([str(PDB)], str(path), "small",
                                           "", False))[0]
    assert entry.invalid_r == 2
    assert [math.isnan(p.distance_model) for p in entry.pairs].count(True) == 2
    # the surviving pair is still scored, and chi2_r is over that pair alone
    assert math.isfinite(entry.chi2_r)
    assert math.isinf(entry.score)  # score goes infinite on the first NaN


def test_screen_reports_ref_rmsd_of_a_frame_it_was_given(tmp_path, structure):
    """End to end: a position carries a frame, screening fits it and the
    pooled RMSD is the distortion put into that frame.

    The frame is attached to an ordinary AV position rather than to a fixed
    (`XYZ`) one, which is what says the two paths agree:
    `fit_reference_positions` does not read `simulation_type`, so a frame is
    fitted and reported the same way wherever it hangs. The fixed-position
    path -- where the fit also *moves* the coordinate that is then scored --
    is `test_fps_point_positions.py`.
    """
    _, h = structure
    atoms = _ca_atoms(h, range(30, 45))
    coordinates = [c for _, c in atoms]
    n = len(coordinates)
    centroid = IMP.algebra.Vector3D(
        sum(c[0] for c in coordinates) / n,
        sum(c[1] for c in coordinates) / n,
        sum(c[2] for c in coordinates) / n)
    rms_radius = math.sqrt(
        sum((c - centroid).get_squared_magnitude() for c in coordinates) / n)
    s = 1.05
    scaled = [(residue, centroid + (c - centroid) * s)
              for residue, c in atoms]

    payload = json.loads(FPS_JSON.read_text())
    chosen = ["19-119_C2"]
    keep = {payload["Distances"][n_][k]
            for n_ in chosen for k in ("position1_name", "position2_name")}
    payload["Distances"] = {n_: payload["Distances"][n_] for n_ in chosen}
    payload["Positions"] = {n_: payload["Positions"][n_] for n_ in keep}
    payload["χ²"] = {"small": {"distances": chosen}}
    payload["Positions"]["19D"]["reference_atoms"] = _frame_json(scaled)
    path = tmp_path / "with_frame.fps.json"
    path.write_text(json.dumps(payload))

    entry = list(IMP.bff.screen_structures([str(PDB)], str(path), "small",
                                           "", False))[0]
    assert entry.ref_rmsd == pytest.approx(abs(s - 1.0) * rms_radius, rel=1e-6)


def test_a_file_without_frames_reports_ref_rmsd_zero(screened):
    """FPS's own answer for a file that declares no reference atoms."""
    assert screened[0].ref_rmsd == 0.0


def test_ranked_csv_carries_the_diagnostic_columns(tmp_path):
    out = tmp_path / "ranked.csv"
    ranked = list(IMP.bff.screen_structures(
        [str(PDB), str(T4L / "3GUN_faspr_port.pdb")], str(FPS_JSON),
        "chi2_C2_33p", str(out), False))
    lines = out.read_text().strip().split("\n")
    assert lines[0] == ("pdb,score,chi2_r,sigma1,sigma2,sigma3,invalid_r,"
                        "ref_rmsd")
    assert len(lines) == 3
    first = lines[1].split(",")
    assert first[0] == ranked[0].path
    assert float(first[1]) == pytest.approx(ranked[0].score, rel=1e-4)
    assert float(first[2]) == pytest.approx(ranked[0].chi2_r, rel=1e-4)
    assert [int(x) for x in first[3:7]] == [
        ranked[0].sigma1, ranked[0].sigma2, ranked[0].sigma3,
        ranked[0].invalid_r]


def test_screening_pairs_csv_stacks_every_structures_table(tmp_path):
    ranked = list(IMP.bff.screen_structures(
        [str(PDB), str(T4L / "3GUN_faspr_port.pdb")], str(FPS_JSON),
        "chi2_C2_33p", "", False))
    out = tmp_path / "pairs.csv"
    IMP.bff.write_screening_pairs_csv(str(out), ranked)
    lines = out.read_text().strip().split("\n")
    assert lines[0].startswith("pdb,name,position1,position2")
    assert len(lines) == 1 + 2 * 33
    # every row carries the structure it belongs to, so a reader can group
    assert {line.split(",")[0] for line in lines[1:]} == {
        e.path for e in ranked}
