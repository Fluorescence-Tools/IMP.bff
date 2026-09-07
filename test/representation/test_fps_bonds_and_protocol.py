"""FPS's bond restraints (PRD-121 G4) and the protocol knobs (G7).

**G4.** `SpringEngine.cs:111-116`: a distance is a *bond* iff **both** ends are
plain atoms -- no dye, a real atom, no accessible volume. That is a crosslink or
a covalent tie between subunits, not a FRET measurement, and FPS treats it
differently in two ways: both anchor atoms drop to `vdWRNoClash = 0.4 A` so they
cannot repel each other (`:119-134`), and the energy is accumulated separately
as `Ebond` (`:445`) -- **a subset of the total, never an addition to it**.

**G7.** Three knobs that change a docked answer rather than decorating it:
`OptimizeSelected` (over *distances*, with clashes never gated), `MaxForce`
(harmonic inside `MaxForce*err^2/2` and **linear** outside), and
`ClashTolerance` (`kclash = 2/ClashTolerance^2`).

The clash cases are on a synthetic two-body pair rather than a real structure,
because the point being asserted is an exact zero: two atoms of radius 0.4 A
that are 3.5 A apart contribute nothing, and on a real complex that zero would
be buried under every other contact.
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
PDB = str(T4L / "3GUN.pdb")
FPS_JSON = T4L / "fret.fps.json"

#: The radius FPS gives a bonded atom (`StaticData.cs:33`).
NO_CLASH_RADIUS = 0.4


def _residue_pdb(path, chain, x0):
    """Two backbone atoms, so a body has an orientation as well as a position."""
    rows = []
    for i, (name, element, x) in enumerate((("N", "N", 0.0),
                                            ("CA", "C", 1.5))):
        rows.append("ATOM  %5d  %-3s ALA %s%4d    %8.3f%8.3f%8.3f  1.00  0.00"
                    "          %2s" % (i + 1, name, chain, 1, x + x0, 0.0, 0.0,
                                       element))
    rows.append("END")
    Path(path).write_text("\n".join(rows) + "\n")
    return str(path)


#: The radii `IMP.atom.read_pdb` gives the two atoms below (measured, not
#: assumed: they are united-atom radii carrying implicit hydrogens, so N is
#: 1.85 and an alpha carbon is 2.275 -- not the Bondi 1.55 / 1.70 FPS uses).
R_N, R_CA = 1.85, 2.275


def _two_bodies(tmp_path, gap=5.0):
    """Body A at 0/1.5 and body B at gap/gap+1.5 along x.

    `gap = 5.0` is chosen so that exactly **one** cross pair overlaps: A's CA
    against B's N, 3.5 A apart against radii 2.275 + 1.85 = 4.125, an overlap
    of 0.625 A. A's N against B's N is 5.0 A against 3.7, A's CA against B's CA
    is 5.0 against 4.55, and A's N against B's CA is 6.5 against 4.125 -- none
    of them touch. So the excluded-volume term *is* that one pair, and a zero
    means the pair stopped clashing rather than something cancelling.
    """
    a = _residue_pdb(tmp_path / "a.pdb", "A", 0.0)
    b = _residue_pdb(tmp_path / "b.pdb", "B", gap)
    return a, b


def _anchor(chain, residue, atom, body):
    return {"simulation_type": "ATOM", "chain_identifier": chain,
            "residue_seq_number": residue, "atom_name": atom, "body_id": body}


def _distance(p1, p2, value=2.5, err=0.5):
    return {"position1_name": p1, "position2_name": p2, "distance": value,
            "error_neg": err, "error_pos": err, "distance_type": "RDAMean",
            "Forster_radius": 52.0}


def _write(tmp_path, payload, name="bond.fps.json"):
    path = tmp_path / name
    path.write_text(json.dumps(payload))
    return str(path)


@pytest.fixture(scope="module")
def t4l():
    IMP.set_log_level(IMP.SILENT)
    return json.loads(FPS_JSON.read_text())


# --------------------------------------------------------------------------
# G4 -- what counts as a bond
# --------------------------------------------------------------------------
def test_a_bond_is_atom_to_atom_and_nothing_else(tmp_path, t4l):
    """Both ends `ATOM`. An AV end or an `XYZ` end does not qualify.

    FPS's test is `Dye == Unknown && AtomID > 0 && AVType == None` on *both*
    ends. An `XYZ` position fails it on `AtomID > 0` -- it is a coordinate, not
    an atom -- which is why the two point types are not interchangeable here.
    """
    av_name = t4l["Distances"]["19-119_C2"]["position1_name"]
    payload = {
        "Positions": {
            av_name: t4l["Positions"][av_name],
            "a1": _anchor("A", 10, "CA", 0),
            "a2": _anchor("A", 60, "CA", 0),
            "fixed": {"simulation_type": "XYZ", "body_id": 0,
                      "x": 0.0, "y": 0.0, "z": 0.0},
        },
        "Distances": {
            "bond": _distance("a1", "a2", 12.0, 1.0),
            "atom_to_av": _distance("a1", av_name, 30.0, 2.0),
            "atom_to_xyz": _distance("a1", "fixed", 30.0, 2.0),
        },
        "χ²": {"s": {"distances": ["bond", "atom_to_av", "atom_to_xyz"]}},
    }
    assembly = IMP.bff.create_docking_assembly([PDB], _write(tmp_path, payload),
                                              "s")
    network = assembly.get_network()
    assert list(network.get_bond_names()) == ["bond"]
    assert network.get_is_bond("bond")
    assert not network.get_is_bond("atom_to_av")
    assert not network.get_is_bond("atom_to_xyz")
    assert not network.get_is_bond("no such distance")


def test_the_anchors_of_a_bond_stop_clashing(tmp_path):
    """The exact zero: the one clashing pair becomes the bonded pair.

    Without the rule the two anchors are 3.5 A apart with radii 2.275 and 1.85
    and the excluded volume is 0.5*k*(4.125-3.5)^2 = 0.1953 at k = 1 (measured:
    0.19531). With it they are 0.4 and 0.4, their sum is 0.8, and there is
    nothing left to score -- and no other cross pair is close enough to
    contribute.
    """
    a, b = _two_bodies(tmp_path)
    payload = {"Positions": {"a1": _anchor("A", 1, "CA", 0),
                             "b1": _anchor("B", 1, "N", 1)},
               "Distances": {"bond": _distance("a1", "b1")},
               "χ²": {"s": {"distances": ["bond"]}}}
    assembly = IMP.bff.create_docking_assembly([a, b], _write(tmp_path, payload),
                                              "s")
    network = assembly.get_network()
    model = assembly.get_model()
    assert list(network.get_bond_names()) == ["bond"]

    radii = {}
    for name in ("a1", "b1"):
        pi = network.get_position_particle_index(name)
        radii[name] = IMP.core.XYZR(model, pi).get_radius()
        assert radii[name] == pytest.approx(NO_CLASH_RADIUS)

    clash = [r for r in assembly.get_restraints().get_restraints()
             if r.get_name() == "excluded_volume"][0]
    assert clash.evaluate(False) == 0.0

    # put a full van der Waals radius back and the same pair does clash, which
    # is what says the zero above came from the rule and not from the geometry
    n = IMP.bff.set_bond_anchor_radii(network, R_N)
    assert n == 2
    assert clash.evaluate(False) == pytest.approx(0.5 * (2 * R_N - 3.5) ** 2,
                                                  rel=1e-5)


def test_the_radius_change_does_not_leak_into_another_run(tmp_path):
    """FPS's D12, deliberately not reproduced.

    There the 0.4 A is written into session-global `Molecule` objects and never
    reverted, so it leaks across clones, modes and runs. Here each assembly
    reads its own PDBs into its own `IMP::Model`, so a second assembly built
    from the same files -- with the same atom used as a position but no bond --
    sees the radius the PDB reader gave it.
    """
    a, b = _two_bodies(tmp_path)
    bonded = {"Positions": {"a1": _anchor("A", 1, "CA", 0),
                            "b1": _anchor("B", 1, "N", 1)},
              "Distances": {"bond": _distance("a1", "b1")},
              "χ²": {"s": {"distances": ["bond"]}}}
    first = IMP.bff.create_docking_assembly([a, b], _write(tmp_path, bonded),
                                           "s")
    assert IMP.core.XYZR(
        first.get_model(),
        first.get_network().get_position_particle_index("a1")).get_radius() \
        == pytest.approx(NO_CLASH_RADIUS)

    # the same atom, now paired with a fixed point: not a bond
    plain = {"Positions": {"a1": _anchor("A", 1, "CA", 0),
                           "p": {"simulation_type": "XYZ", "body_id": 1,
                                 "x": 4.0, "y": 0.0, "z": 0.0}},
             "Distances": {"d": _distance("a1", "p")},
             "χ²": {"s": {"distances": ["d"]}}}
    second = IMP.bff.create_docking_assembly(
        [a, b], _write(tmp_path, plain, "plain.fps.json"), "s")
    assert list(second.get_network().get_bond_names()) == []
    radius = IMP.core.XYZR(
        second.get_model(),
        second.get_network().get_position_particle_index("a1")).get_radius()
    assert radius > 1.0, radius     # 1.85 from the PDB reader, not 0.4


def test_e_bond_is_a_subset_of_the_score_never_an_addition(tmp_path, t4l):
    """`SpringEngine.cs:445` accumulates `Ebond` inside the total's own loop."""
    av_name = t4l["Distances"]["19-119_C2"]["position1_name"]
    payload = {"Positions": {av_name: t4l["Positions"][av_name],
                             "a1": _anchor("A", 10, "CA", 0),
                             "a2": _anchor("A", 60, "CA", 0)},
               "Distances": {"bond": _distance("a1", "a2", 12.0, 1.0),
                             "fret": _distance("a1", av_name, 30.0, 2.0)},
               "χ²": {"s": {"distances": ["bond", "fret"]}}}
    assembly = IMP.bff.create_docking_assembly([PDB], _write(tmp_path, payload),
                                              "s")
    result = IMP.bff.score_assembly(assembly)

    assert result.n_bonds == 1
    bonds = [p for p in result.pairs if p.is_bond]
    others = [p for p in result.pairs if not p.is_bond]
    assert [p.name for p in bonds] == ["bond"]
    assert [p.name for p in others] == ["fret"]

    # each restraint is worth half a chi-square; one body has no clashes
    assert result.e_bond == pytest.approx(0.5 * bonds[0].chi2, rel=1e-9)
    assert result.score == pytest.approx(
        0.5 * sum(p.chi2 for p in result.pairs), rel=1e-9)
    # the subset relation, stated the way it is easy to get wrong
    assert result.e_bond < result.score
    assert result.score != pytest.approx(result.score + result.e_bond)
    assert result.score - result.e_bond == pytest.approx(0.5 * others[0].chi2,
                                                         rel=1e-9)


def test_a_file_without_bonds_reports_zero(tmp_path):
    result = IMP.bff.score_structures([PDB], str(FPS_JSON), "chi2_C2_33p")
    assert result.n_bonds == 0
    assert result.e_bond == 0.0
    assert not any(p.is_bond for p in result.pairs)


# --------------------------------------------------------------------------
# G7 -- MaxForce: the linear tail
# --------------------------------------------------------------------------
@pytest.mark.parametrize("err", [0.1, 1.0, 3.0])
def test_the_restraint_is_harmonic_inside_the_knee(err):
    """Inside `drmax` the capped form must equal the plain chi-square."""
    f = 400.0
    drmax = f * err * err / 2.0
    for frac in (0.0, 0.25, 0.9999):
        model = 50.0 + frac * drmax
        assert IMP.bff.chi2_score_capped(model, 50.0, err, err, f) == \
            pytest.approx(IMP.bff.chi2_score(model, 50.0, err, err), rel=1e-12)


@pytest.mark.parametrize("err", [0.1, 1.0, 3.0])
def test_the_tail_past_the_knee_is_linear_with_slope_max_force(err):
    """Two points past the knee, and the slope between them.

    A parabola would give a slope that grows with the deviation; a capped
    restraint gives `MaxForce` at every separation past `drmax`. That is the
    difference that changes a docked pose, so it is measured rather than
    assumed -- twice, at different deviations, so a coincidence cannot pass.
    """
    f = 400.0
    drmax = f * err * err / 2.0
    a, b, c = 2.0 * drmax, 3.0 * drmax, 7.0 * drmax
    ea = IMP.bff.chi2_score_capped(50.0 + a, 50.0, err, err, f)
    eb = IMP.bff.chi2_score_capped(50.0 + b, 50.0, err, err, f)
    ec = IMP.bff.chi2_score_capped(50.0 + c, 50.0, err, err, f)
    assert (eb - ea) / (b - a) == pytest.approx(f, rel=1e-9)
    assert (ec - eb) / (c - b) == pytest.approx(f, rel=1e-9)
    # and it is genuinely below the parabola out there
    assert eb < IMP.bff.chi2_score(50.0 + b, 50.0, err, err)


def test_the_two_pieces_join_smoothly_at_the_knee():
    """C^1: same value and same slope, so a minimiser sees no step."""
    f, err = 400.0, 1.0
    drmax = f * err * err / 2.0
    eps = 1e-6
    below = IMP.bff.chi2_score_capped(50.0 + drmax - eps, 50.0, err, err, f)
    above = IMP.bff.chi2_score_capped(50.0 + drmax + eps, 50.0, err, err, f)
    # they cannot be equal -- the function has slope MaxForce there -- so what
    # is asserted is that the gap is the slope times the step and nothing more
    assert above - below == pytest.approx(2.0 * eps * f, rel=1e-4)
    # the value at the knee is MaxForce * drmax / 2 -- FPS's own expression
    assert below == pytest.approx(0.5 * f * drmax, rel=1e-6)


def test_the_cap_uses_the_error_bar_of_the_side_it_is_on():
    """Asymmetric errors put the knee in a different place on each side."""
    f, neg, pos = 100.0, 1.0, 4.0
    knee_pos = f * pos * pos / 2.0      # 800
    knee_neg = f * neg * neg / 2.0      # 50
    # 100 A too long is still inside the positive knee -> parabola
    assert IMP.bff.chi2_score_capped(150.0, 50.0, neg, pos, f) == \
        pytest.approx(IMP.bff.chi2_score(150.0, 50.0, neg, pos), rel=1e-12)
    assert knee_pos > 100.0
    # 100 A too short is far outside the negative knee -> linear
    short = IMP.bff.chi2_score_capped(-50.0, 50.0, neg, pos, f)
    assert short == pytest.approx(f * (100.0 - 0.5 * knee_neg), rel=1e-9)
    assert short < IMP.bff.chi2_score(-50.0, 50.0, neg, pos)


def test_a_non_positive_max_force_is_the_plain_chi_square():
    for f in (0.0, -1.0):
        assert IMP.bff.chi2_score_capped(60.0, 50.0, 1.0, 10.0, f) == \
            pytest.approx(IMP.bff.chi2_score(60.0, 50.0, 1.0, 10.0))


def test_the_mean_distance_restraint_honours_max_force():
    """The same tail through the restraint a docking run actually minimises.

    The restraint scores *half* a chi-square, so the slope it shows is half
    `MaxForce`; the knee is where FPS puts it either way, and the knee is what
    moves a pose.
    """
    m = IMP.Model()
    p1 = IMP.Particle(m, "a")
    p2 = IMP.Particle(m, "b")
    IMP.core.XYZ.setup_particle(p1, IMP.algebra.Vector3D(0, 0, 0))
    IMP.core.XYZ.setup_particle(p2, IMP.algebra.Vector3D(0, 0, 50))
    e = IMP.bff.AVPairDistanceMeasurement()
    e.position_1, e.position_2 = "a", "b"
    e.distance, e.error_neg, e.error_pos = 50.0, 1.0, 1.0
    e.distance_type = IMP.bff.PROBE_PAIR_DISTANCE_MP

    f = 20.0                       # knee at 20 * 1 / 2 = 10 A
    capped = IMP.bff.AVMeanDistanceRestraint(m, p1, p2, e, 6.0, 1.0, f)
    plain = IMP.bff.AVMeanDistanceRestraint(m, p1, p2, e, 6.0, 1.0, 0.0)
    assert capped.get_max_force() == pytest.approx(f)

    # inside the knee the two agree
    assert capped.get_score_at(55.0) == pytest.approx(plain.get_score_at(55.0))
    # outside it the capped one is linear with slope f/2, and lower
    a, b = 80.0, 120.0
    slope = (capped.get_score_at(b) - capped.get_score_at(a)) / (b - a)
    assert slope == pytest.approx(0.5 * f, rel=1e-9)
    assert capped.get_score_at(b) < plain.get_score_at(b)
    # and it is the restraint that is evaluated, not just a helper
    assert capped.evaluate(False) == pytest.approx(capped.get_score_at(50.0))


# --------------------------------------------------------------------------
# G7 -- ClashTolerance
# --------------------------------------------------------------------------
def test_clash_tolerance_becomes_the_soft_sphere_constant(tmp_path):
    """`kclash = 2/ClashTolerance^2` (`SpringEngine.cs:147`).

    IMP's soft sphere scores `0.5*k*overlap^2`, so an overlap of exactly one
    tolerance costs one chi-square unit -- which is what `ClashTolerance`
    *means*. Measured on the one-pair fixture, where the overlap is 0.625 A:
    0.19531 at k = 1, 0.39063 at FPS's docking tolerance of 1.0 A.
    """
    a, b = _two_bodies(tmp_path)
    # One ATOM end and one fixed point, so this is *not* a bond and the
    # clashing pair keeps the radii the PDB reader gave it.
    payload = {"Positions": {"a1": _anchor("A", 1, "N", 0),
                             "p": {"simulation_type": "XYZ", "body_id": 1,
                                   "x": 5.0, "y": 0.0, "z": 0.0}},
               "Distances": {"d": _distance("a1", "p", 5.0, 1.0)},
               "χ²": {"s": {"distances": ["d"]}}}
    path = _write(tmp_path, payload, "ct.fps.json")
    overlap = (R_CA + R_N) - 3.5

    def clash_of(tolerance):
        assembly = IMP.bff.create_docking_assembly(
            [a, b], path, "s", clash_tolerance=tolerance)
        return [r for r in assembly.get_restraints().get_restraints()
                if r.get_name() == "excluded_volume"][0].evaluate(False)

    # non-positive keeps IMP's historical k = 1
    assert clash_of(0.0) == pytest.approx(0.5 * overlap ** 2, rel=1e-6)
    for tolerance in (1.0, 0.5, 2.0):
        k = 2.0 / tolerance ** 2
        assert clash_of(tolerance) == pytest.approx(0.5 * k * overlap ** 2,
                                                    rel=1e-6)
    # FPS's refinement tolerance penalises four times harder per angstrom than
    # its docking one
    assert clash_of(0.5) == pytest.approx(4.0 * clash_of(1.0), rel=1e-9)


def test_the_default_scoring_door_keeps_imps_own_constant():
    """The default is deliberately *not* FPS's, and this pins that.

    `create_docking_assembly` is the scoring door and its numbers are pinned
    elsewhere (HIV-RT `resolved` at 59.0404, of which 24.8 is the clash term).
    FPS's docking constant is opted into through `DockingParameters`, which is
    where the docking protocol lives.
    """
    assert IMP.bff.DockingParameters().clash_tolerance == pytest.approx(1.0)
    protein = IMP.bff.get_example_path("structure/HIV_RT/protein_1R0A.pdb")
    dna = IMP.bff.get_example_path("structure/HIV_RT/dna.pdb")
    fps = IMP.bff.get_example_path("structure/HIV_RT/hiv_rt.fps.json")
    result = IMP.bff.score_structures([protein, dna], fps, "resolved")
    # Olga's radii gave 59.2547 for the hours they were the default on
    # 2026-09-01 (this fixture's DNA carries the old `C1*`/`O1P` spellings,
    # which that name-keyed table does not have, so those 85 atoms took its
    # 1.50 A unknown-name fallback). IMP's own radii are the default again
    # (AV::set_radii_source) and 59.0404 returns exactly.
    assert result.score == pytest.approx(59.040389, abs=1e-3)


# --------------------------------------------------------------------------
# G7 -- OptimizeSelected
# --------------------------------------------------------------------------
@pytest.mark.parametrize("mode,expected", [
    ("Selected", (33, 0)),
    ("All", (33, 66)),
    ("SelectedThenAll", (33, 66)),
    ("selected-then-all", (33, 66)),     # spelling is forgiving
    ("nonsense", (33, 0)),               # and unknown text is "Selected"
])
def test_optimize_selected_gates_distances_not_molecules(tmp_path, mode,
                                                         expected):
    """"Selected" is about the *distances*; here that is the score set.

    The T4L file carries 99 distances and its `chi2_C2_33p` set names 33, so
    the three modes are trivially distinguishable -- which is the point.
    """
    params = IMP.bff.DockingParameters()
    params.score_set = "chi2_C2_33p"
    params.shuffle_max_translation = 0.0
    params.n_frames = 5
    params.optimize_selected = mode
    result = IMP.bff.dock_minimize([PDB], str(FPS_JSON),
                                   str(tmp_path / mode.replace("-", "_")),
                                   params=params)
    extra = json.loads(result.extra)
    assert (extra["n_selected_distances"],
            extra["n_deselected_distances"]) == expected
    assert result.n_distances == sum(expected)
    assert math.isfinite(result.score)


def test_the_protocol_knobs_are_recorded_in_the_result(tmp_path):
    """A pose that does not say what it was produced under cannot be compared."""
    params = IMP.bff.DockingParameters()
    params.score_set = "chi2_C2_33p"
    params.shuffle_max_translation = 0.0
    params.n_frames = 5
    params.max_force = 250.0
    params.clash_tolerance = 0.5
    result = IMP.bff.dock_minimize([PDB], str(FPS_JSON), str(tmp_path),
                                   params=params)
    extra = json.loads(result.extra)
    assert extra["optimize_selected"] == "selected"
    assert extra["max_force"] == pytest.approx(250.0)
    assert extra["clash_tolerance"] == pytest.approx(0.5)
    assert extra["k_clash"] == pytest.approx(2.0 / 0.25)
    assert extra["n_bond_anchors"] == 0
    assert extra["n_points"] == 0


def test_the_shipped_defaults_are_the_fps_docking_ones():
    """`ProjectData.cs:67-81`: the settings behind the published work."""
    p = IMP.bff.DockingParameters()
    assert p.optimize_selected == "Selected"
    assert p.max_force == pytest.approx(400.0)
    assert p.clash_tolerance == pytest.approx(1.0)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
