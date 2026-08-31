"""Flat-bottom FRET restraints, and the pieces an MD run needs under them.

`AVMeanDistanceRestraint` scores a chi-squared: it pulls at every separation and
never stops, which is a scoring function, not a potential to integrate. A
flat-bottom well is zero inside the experimental error bars and only pushes back
outside them, and past the error bars again it goes **linear** so the force is
capped -- which is what keeps a badly-placed start from blowing up the first MD
step. It is AMBER's `&rst` form, so the restraint here and the restraint an MD
engine would be handed are the same restraint.

What is tested: the well's shape and its gradient against finite differences;
the conversion of a measured distance into mean-position bounds; and the whole
`md_flat_bottom_restraints` system, including the one thing whose absence is
silent -- a probe particle nobody marked optimizable, which leaves every well
pushing on something that cannot move.
"""

import json

import numpy as np
import pytest

import IMP
import IMP.algebra
import IMP.atom
import IMP.core
import IMP.bff


BOUNDS = (35.0, 40.0, 46.0, 51.0)
K2, K3 = 0.5, 0.4


@pytest.fixture
def well():
    m = IMP.Model()
    a = IMP.core.XYZ.setup_particle(IMP.Particle(m), IMP.algebra.Vector3D(0, 0, 0))
    b = IMP.core.XYZ.setup_particle(IMP.Particle(m), IMP.algebra.Vector3D(50, 0, 0))
    for p in (a, b):
        p.set_coordinates_are_optimized(True)
    r = IMP.bff.AVFlatBottomRestraint(m, a, b, *BOUNDS, K2, K3)
    return m, a, b, r


def _at(b, r, d):
    b.set_coordinates(IMP.algebra.Vector3D(d, 0, 0))
    return r.unprotected_evaluate(None)


def test_the_bottom_is_flat_and_the_walls_are_harmonic(well):
    m, a, b, r = well
    r1, r2, r3, r4 = BOUNDS
    for d in (r2, 42.0, 44.0, r3):
        assert _at(b, r, d) == 0.0, f"the well should be flat at {d}"
    for d, k, edge in ((38.0, K2, r2), (49.0, K3, r3)):
        assert _at(b, r, d) == pytest.approx(k * (d - edge) ** 2)


def test_beyond_the_walls_the_force_is_capped(well):
    """The reason the linear continuation is there: a start 50 A out of the
    well must not deliver a 50-fold force to the first integration step."""
    m, a, b, r = well
    r1, r2, r3, r4 = BOUNDS
    h = 1e-4

    def slope(d):
        return (_at(b, r, d + h) - _at(b, r, d - h)) / (2 * h)

    at_r4 = slope(r4 + 1e-3)
    assert slope(100.0) == pytest.approx(at_r4, rel=1e-6), "not linear"
    assert slope(1000.0) == pytest.approx(at_r4, rel=1e-6)
    assert slope(0.5 * r1) == pytest.approx(slope(r1 - 1e-3), rel=1e-6)
    # and the value is continuous where the pieces meet: the tolerance is the
    # step times the slope there, not an arbitrary epsilon
    for edge in (r1, r4):
        step, slope_here = 1e-7, abs(slope(edge + 1e-3))
        assert _at(b, r, edge + step) == pytest.approx(
            _at(b, r, edge - step), abs=4 * step * slope_here + 1e-9)


@pytest.mark.parametrize("d", [20.0, 34.0, 37.0, 39.9, 43.0, 47.0, 50.9, 55.0, 90.0])
def test_the_gradient_is_the_derivative_of_the_score(well, d):
    """An MD engine integrates the gradient, so it is the gradient that has to
    be right -- a score that is correct and a force that is not integrates into
    a trajectory that is neither."""
    m, a, b, r = well
    sf = IMP.core.RestraintsScoringFunction([r])
    b.set_coordinates(IMP.algebra.Vector3D(d, 0, 0))
    sf.evaluate(True)
    analytic = b.get_derivatives()[0]

    h = 1e-5
    numeric = (_at(b, r, d + h) - _at(b, r, d - h)) / (2 * h)
    assert analytic == pytest.approx(numeric, abs=1e-4)
    # and equal and opposite on the other particle
    b.set_coordinates(IMP.algebra.Vector3D(d, 0, 0))
    sf.evaluate(True)
    np.testing.assert_allclose(list(a.get_derivatives()),
                               [-x for x in b.get_derivatives()], atol=1e-9)


def test_bounds_must_be_ordered():
    m = IMP.Model()
    a = IMP.core.XYZ.setup_particle(IMP.Particle(m), IMP.algebra.Vector3D(0, 0, 0))
    b = IMP.core.XYZ.setup_particle(IMP.Particle(m), IMP.algebra.Vector3D(1, 0, 0))
    with pytest.raises(Exception):
        IMP.bff.AVFlatBottomRestraint(m, a, b, 40.0, 35.0, 46.0, 51.0, 1.0, 1.0)


# ---------------------------------------------------------------------------
# the system an MD run is handed
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def t4l():
    IMP.set_log_level(IMP.SILENT)
    m = IMP.Model()
    h = IMP.atom.read_pdb(IMP.bff.get_example_path("structure/T4L/3GUN.pdb"), m,
                          IMP.atom.NonWaterNonHydrogenPDBSelector())
    return m, h, IMP.bff.get_example_path("structure/T4L/fret.fps.json")


@pytest.fixture(scope="module")
def system(t4l):
    m, h, fps = t4l
    return m, h, IMP.bff.md_flat_bottom_restraints(h, fps, "chi2_C1_33p", 5.0)


def test_one_probe_per_position_and_one_well_per_pair(system):
    m, h, s = system
    assert len(s.get_probes()) == 17
    assert len(s.get_wells()) == 33 == len(s.get_pair_names())
    # the set is the wells plus one tether per probe
    assert s.get_restraints().get_number_of_restraints() == 33 + 17
    assert len(s.get_bounds()) == 4 * 33


def test_a_probe_is_a_movable_particle_at_the_volumes_mean_position(system):
    """The failure this catches is silent: a probe with mass and coordinates
    that nobody marked optimizable leaves every well pushing on something that
    cannot move, and an MD run reports no error and no motion."""
    m, h, s = system
    for i, probe in enumerate(s.get_probes()):
        xyz = IMP.core.XYZ(m, probe.particle)
        assert xyz.get_coordinates_are_optimized(), \
            f"{probe.position_name} cannot be moved by an optimiser"
        assert IMP.atom.Mass.get_is_setup(m, probe.particle)
        assert IMP.core.XYZR(m, probe.particle).get_radius() > 0.0
        # tethered at the offset the volume itself has
        at = IMP.core.XYZ(m, probe.attachment).get_coordinates()
        assert IMP.algebra.get_distance(xyz.get_coordinates(), at) == \
            pytest.approx(s.get_tether_length(i), abs=1e-6)
        assert s.get_tether_length(i) > 0.0
        assert probe.tether_k > 0.0


def test_the_bounds_bracket_the_measurement_in_mean_position_space(system):
    m, h, s = system
    bounds = np.asarray(list(s.get_bounds())).reshape(-1, 4)
    r1, r2, r3, r4 = bounds.T
    assert np.all(r1 <= r2) and np.all(r2 <= r3) and np.all(r3 <= r4)
    assert np.all(r1 >= 0.0)
    # a well with no width would be a chi-squared with extra steps
    assert np.all(r3 - r2 > 0.0)


def test_tethering_to_another_atom_of_the_same_residue(t4l):
    """Only alpha carbons move in a coarse run, so a probe on the C-beta the
    fps.json names would be tethered to something nothing optimises."""
    m, h, fps = t4l
    ca = IMP.bff.md_flat_bottom_restraints(h, fps, "chi2_C1_33p", 5.0, 10.0,
                                           100.0, "CA")
    for probe in ca.get_probes():
        atom = IMP.atom.Atom(m, probe.attachment)
        assert atom.get_atom_type() == IMP.atom.AT_CA
    with pytest.raises(Exception):
        IMP.bff.md_flat_bottom_restraints(h, fps, "chi2_C1_33p", 5.0, 10.0,
                                          100.0, "NOSUCHATOM")


def test_molecular_dynamics_moves_the_structure_toward_the_data(t4l):
    """The whole point, end to end: restrained MD should satisfy more of the
    measured distances than the structure it started from."""
    m, h, fps = t4l
    cas = IMP.atom.Selection(h, atom_type=IMP.atom.AT_CA).get_selected_particles()
    for p in cas:
        if not IMP.core.XYZR.get_is_setup(p):
            IMP.core.XYZR.setup_particle(p, 2.5)
        if not IMP.atom.Mass.get_is_setup(p):
            IMP.atom.Mass.setup_particle(p, 110.0)
        IMP.core.XYZ(p).set_coordinates_are_optimized(True)

    xyz = [IMP.core.XYZ(p).get_coordinates() for p in cas]
    network = []
    for i in range(len(cas)):
        for j in range(i + 1, len(cas)):
            d = IMP.algebra.get_distance(xyz[i], xyz[j])
            if d < 10.0:
                k = 20.0 if j == i + 1 else 0.3
                network.append(IMP.core.DistanceRestraint(
                    m, IMP.core.Harmonic(d, k), cas[i], cas[j]))

    s = IMP.bff.md_flat_bottom_restraints(h, fps, "chi2_C2_33p", 15.0, 30.0,
                                          100.0, "CA")
    wells = s.get_wells()

    def satisfied():
        return sum(1 for w in wells
                   if w.get_bounds()[1] <= w.get_distance() <= w.get_bounds()[2])

    before = satisfied()
    sf = IMP.core.RestraintsScoringFunction(network + [s.get_restraints()])
    md = IMP.atom.MolecularDynamics(m)
    md.set_scoring_function(sf)
    md.set_maximum_time_step(2.0)
    md.add_optimizer_state(IMP.atom.VelocityScalingOptimizerState(m, cas, 300.0))
    md.optimize(6000)
    IMP.core.ConjugateGradients(m).set_scoring_function(sf)
    cg = IMP.core.ConjugateGradients(m)
    cg.set_scoring_function(sf)
    cg.optimize(200)

    after = satisfied()
    assert after > before, (
        f"restrained MD satisfied {after} of {len(wells)} distances, "
        f"no better than the {before} it started with")


# ---------------------------------------------------------------------------
# out to OpenMM
# ---------------------------------------------------------------------------

#: A -> nm and kcal/mol/A^2 -> kJ/mol/nm^2, the conversion the writer applies.
A_TO_NM = 0.1
KCAL_TO_KJ = 4.184
K_TO_OPENMM = KCAL_TO_KJ * 100.0


def openmm_eval(expression, **values):
    """Evaluate an OpenMM expression the way OpenMM would.

    `select(c, a, b)` and `step(x)` are OpenMM's; `^` is its exponent operator.
    This is what makes the exported restraint checkable without OpenMM
    installed -- and what would catch the expression drifting away from the C++
    it is supposed to mirror.
    """
    ns = dict(values)
    ns["select"] = lambda c, a, b: a if c else b
    ns["step"] = lambda x: 1.0 if x >= 0 else 0.0
    return eval(expression.replace("^", "**"), {"__builtins__": {}}, ns)


def test_the_exported_expression_is_the_restraint_this_module_evaluates(well):
    """The whole point of exporting: the well that runs in OpenMM has to be the
    well that ran here. Checked over the full range, through both walls and both
    linear tails, with the unit conversion applied."""
    m, a, b, r = well
    expression = IMP.bff.openmm_flat_bottom_energy()
    r1, r2, r3, r4 = BOUNDS
    for d in np.linspace(5.0, 120.0, 80):
        here = _at(b, r, d) * KCAL_TO_KJ
        there = openmm_eval(expression, r=d * A_TO_NM,
                            r1=r1 * A_TO_NM, r2=r2 * A_TO_NM,
                            r3=r3 * A_TO_NM, r4=r4 * A_TO_NM,
                            k2=K2 * K_TO_OPENMM, k3=K3 * K_TO_OPENMM)
        assert there == pytest.approx(here, abs=1e-9, rel=1e-12), \
            f"the two disagree at {d} A: {here} kJ/mol here, {there} there"


def test_the_expression_names_exactly_the_parameters_it_declares(system, tmp_path):
    m, h, s = system
    path = tmp_path / "openmm.json"
    IMP.bff.write_openmm_restraints(s, str(path), h)
    doc = json.loads(path.read_text())
    declared = set(doc["per_bond_parameters"])
    assert declared == {"r1", "r2", "r3", "r4", "k2", "k3"}
    # every declared parameter is used, and nothing else is referenced
    expression = doc["energy_expression"]
    import re
    used = set(re.findall(r"\b[a-z]\w*\b", expression))
    assert declared <= used
    assert used - declared - {"select", "step", "r"} == set()


def test_the_document_carries_what_an_openmm_script_has_to_add(system, tmp_path):
    m, h, s = system
    path = tmp_path / "openmm.json"
    IMP.bff.write_openmm_restraints(s, str(path), h)
    doc = json.loads(path.read_text())

    assert doc["units"] == {"length": "nanometer",
                            "energy": "kilojoule_per_mole",
                            "mass": "dalton"}
    assert len(doc["probes"]) == len(s.get_probes())
    assert len(doc["restraints"]) == len(s.get_wells())

    names = {p["name"] for p in doc["probes"]}
    for probe in doc["probes"]:
        assert len(probe["position"]) == 3
        assert probe["mass"] > 0.0
        assert probe["tether"]["length"] > 0.0
        assert probe["tether"]["k"] > 0.0
        # named, not indexed: an index depends on how the reader built its
        # topology and would silently restrain the wrong atom
        at = probe["attachment"]
        assert set(at) == {"chain", "residue", "atom"}
        assert at["chain"] and isinstance(at["residue"], int) and at["atom"]

    for rst in doc["restraints"]:
        assert rst["probe_1"] in names and rst["probe_2"] in names
        assert rst["probe_1"] != rst["probe_2"]
        assert rst["r1"] <= rst["r2"] <= rst["r3"] <= rst["r4"]
        assert rst["k2"] > 0.0 and rst["k3"] > 0.0


def test_the_document_is_in_openmm_units_and_not_this_modules(system, tmp_path):
    """A restraint exported in the wrong unit is a restraint that runs and is
    wrong, so the conversion is asserted against the values it came from."""
    m, h, s = system
    path = tmp_path / "openmm.json"
    IMP.bff.write_openmm_restraints(s, str(path), h)
    doc = json.loads(path.read_text())

    bounds = np.asarray(list(s.get_bounds())).reshape(-1, 4)
    for i, rst in enumerate(doc["restraints"]):
        for j, key in enumerate(("r1", "r2", "r3", "r4")):
            assert rst[key] == pytest.approx(bounds[i][j] * A_TO_NM)
        k2, k3 = s.get_wells()[i].get_force_constants()
        assert rst["k2"] == pytest.approx(k2 * K_TO_OPENMM)
        assert rst["k3"] == pytest.approx(k3 * K_TO_OPENMM)

    for i, probe in enumerate(doc["probes"]):
        assert probe["tether"]["length"] == pytest.approx(
            s.get_tether_length(i) * A_TO_NM)
    # a nanometre-scale document: nothing should look like angstrom
    assert max(r["r4"] for r in doc["restraints"]) < 20.0


# ---------------------------------------------------------------------------
# the generated OpenMM script
# ---------------------------------------------------------------------------

def test_the_generated_script_is_valid_python_and_self_contained(system,
                                                                 tmp_path):
    """OpenMM does not have to be installed to write it, so what can be checked
    here is that it parses, that it carries its own restraint table, and that
    the table is the one the system holds."""
    import ast

    m, h, s = system
    pdb = IMP.bff.get_example_path("structure/T4L/3GUN.pdb")
    path = tmp_path / "run.py"
    IMP.bff.write_openmm_script(s, str(path), h, pdb)

    text = path.read_text()
    tree = ast.parse(text)           # it is Python

    literals = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name):
            try:
                literals[node.targets[0].id] = ast.literal_eval(node.value)
            except (ValueError, SyntaxError):
                pass

    assert len(literals["PROBES"]) == len(s.get_probes())
    assert len(literals["WELLS"]) == len(s.get_wells())
    assert literals["PARAMETERS"] == ["r1", "r2", "r3", "r4", "k2", "k3"]
    assert literals["ENERGY"] == IMP.bff.openmm_flat_bottom_energy()

    # the forces it builds, and the fact that it reads no other file
    for needed in ("CustomBondForce", "HarmonicBondForce", "addParticle",
                   "LangevinMiddleIntegrator", "minimizeEnergy"):
        assert needed in text, f"the script never calls {needed}"
    assert "json" not in text.lower() or "PROBES" in text

    names = {p[0] for p in literals["PROBES"]}
    for well in literals["WELLS"]:
        assert well[1] in names and well[2] in names
    bounds = np.asarray(list(s.get_bounds())).reshape(-1, 4)
    for i, well in enumerate(literals["WELLS"]):
        for j in range(4):
            assert well[3 + j] == pytest.approx(bounds[i][j] * A_TO_NM, rel=1e-5)


def test_the_script_names_atoms_the_way_the_document_does(system, tmp_path):
    import ast
    m, h, s = system
    pdb = IMP.bff.get_example_path("structure/T4L/3GUN.pdb")
    script = tmp_path / "run.py"
    doc_path = tmp_path / "run.json"
    IMP.bff.write_openmm_script(s, str(script), h, pdb)
    IMP.bff.write_openmm_restraints(s, str(doc_path), h)

    literals = {}
    for node in ast.parse(script.read_text()).body:
        if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name):
            try:
                literals[node.targets[0].id] = ast.literal_eval(node.value)
            except (ValueError, SyntaxError):
                pass
    doc = json.loads(doc_path.read_text())

    by_name = {p["name"]: p for p in doc["probes"]}
    for name, chain, residue, atom, mass, length, k, xyz in literals["PROBES"]:
        other = by_name[name]["attachment"]
        assert (chain, residue, atom) == (other["chain"], other["residue"],
                                          other["atom"])
        assert length == pytest.approx(by_name[name]["tether"]["length"], rel=1e-5)


# ---------------------------------------------------------------------------
# rebuilding the volumes during a run
# ---------------------------------------------------------------------------

def test_the_volumes_can_be_rebuilt_during_a_run(t4l):
    """`md_flat_bottom_restraints` derives every well from volumes computed once
    on the starting structure. That approximation decays as the structure moves;
    this is the state that refreshes it, and the test is that refreshing
    actually changes the wells."""
    m, h, fps = t4l
    cas = IMP.atom.Selection(h, atom_type=IMP.atom.AT_CA).get_selected_particles()
    for p in cas:
        if not IMP.core.XYZR.get_is_setup(p):
            IMP.core.XYZR.setup_particle(p, 2.5)
        if not IMP.atom.Mass.get_is_setup(p):
            IMP.atom.Mass.setup_particle(p, 110.0)
        IMP.core.XYZ(p).set_coordinates_are_optimized(True)

    xyz = [IMP.core.XYZ(p).get_coordinates() for p in cas]
    network = []
    for i in range(len(cas)):
        for j in range(i + 1, len(cas)):
            d = IMP.algebra.get_distance(xyz[i], xyz[j])
            if d < 10.0:
                network.append(IMP.core.DistanceRestraint(
                    m, IMP.core.Harmonic(d, 20.0 if j == i + 1 else 0.3),
                    cas[i], cas[j]))

    s = IMP.bff.md_flat_bottom_restraints(h, fps, "chi2_C2_33p", 15.0, 30.0,
                                          100.0, "CA")
    assert s.get_network() is not None
    before = np.array([list(w.get_bounds()) for w in s.get_wells()])
    offsets_before = [s.get_tether_length(i)
                      for i in range(len(s.get_probes()))]

    state = IMP.bff.AVRebuildOptimizerState(m, s, 2000)
    sf = IMP.core.RestraintsScoringFunction(network + [s.get_restraints()])
    md = IMP.atom.MolecularDynamics(m)
    md.set_scoring_function(sf)
    md.set_maximum_time_step(2.0)
    md.add_optimizer_state(IMP.atom.VelocityScalingOptimizerState(m, cas, 300.0))
    md.add_optimizer_state(state)
    md.optimize(6000)

    assert state.get_number_of_updates() >= 2, "the state never fired"
    after = np.array([list(w.get_bounds()) for w in s.get_wells()])
    assert not np.allclose(before, after), \
        "the wells were not re-derived, so the volumes were not rebuilt"
    # still a well, not a crossed pair of bounds
    assert np.all(np.diff(after, axis=1) >= -1e-9)
    # the tethers moved with the volumes they hold
    assert [s.get_tether_length(i) for i in range(len(s.get_probes()))] \
        != offsets_before, "the tethers did not follow their volumes"


def test_a_rebuild_can_be_asked_for_directly(system):
    m, h, s = system
    state = IMP.bff.AVRebuildOptimizerState(m, s, 1000)
    assert state.get_number_of_updates() == 0
    state.update_now()
    assert state.get_number_of_updates() == 1
    once = np.array([list(w.get_bounds()) for w in s.get_wells()])
    # idempotent at fixed coordinates: a second rebuild of the same structure
    # must not drift, or every rebuild would move the wells a little
    state.update_now()
    twice = np.array([list(w.get_bounds()) for w in s.get_wells()])
    np.testing.assert_allclose(twice, once, atol=1e-9)


def test_a_well_refuses_bounds_that_are_not_ordered(well):
    m, a, b, r = well
    r.set_bounds(1.0, 2.0, 3.0, 4.0)
    assert list(r.get_bounds()) == [1.0, 2.0, 3.0, 4.0]
    with pytest.raises(Exception):
        r.set_bounds(4.0, 3.0, 2.0, 1.0)


# ---------------------------------------------------------------------------
# the commands
# ---------------------------------------------------------------------------

def test_the_openmm_command_writes_a_runnable_script(imp_bff_program, tmp_path):
    import ast
    from click.testing import CliRunner
    script = tmp_path / "fret.py"
    doc = tmp_path / "fret.json"
    result = CliRunner().invoke(imp_bff_program.openmm_command, [
        "--fps-json", IMP.bff.get_example_path("structure/T4L/fret.fps.json"),
        "--pdb", IMP.bff.get_example_path("structure/T4L/3GUN.pdb"),
        "--score-set", "chi2_C1_33p", "--tether-atom", "CA",
        "--output", str(script), "--json", str(doc), "--n-steps", "500"])
    assert result.exit_code == 0, result.output
    assert script.is_file() and doc.is_file()
    ast.parse(script.read_text())
    assert "N_STEPS = 500" in script.read_text()
    assert "33 FRET wells" in result.output


def test_the_av_export_command_writes_each_format(imp_bff_program, tmp_path):
    from click.testing import CliRunner
    for ext in ("xyz", "pqr", "dx"):
        out = tmp_path / f"site.{ext}"
        result = CliRunner().invoke(imp_bff_program.av_export, [
            "--pdb", IMP.bff.get_example_path("structure/T4L/3GUN.pdb"),
            "--chain", "A", "--residue", "132", "--output", str(out)])
        assert result.exit_code == 0, result.output
        assert out.is_file() and out.stat().st_size > 0
        assert "mean position" in result.output


def test_av_export_says_so_when_the_site_is_not_there(imp_bff_program, tmp_path):
    from click.testing import CliRunner
    result = CliRunner().invoke(imp_bff_program.av_export, [
        "--pdb", IMP.bff.get_example_path("structure/T4L/3GUN.pdb"),
        "--chain", "Z", "--residue", "9999",
        "--output", str(tmp_path / "nowhere.pqr")])
    assert result.exit_code != 0
    assert "no atom" in result.output
