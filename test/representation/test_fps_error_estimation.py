"""FPS's error estimation, the parametric bootstrap (PRD-121 G2).

`ErrorEstimation.cs` takes the *docked model's own* distances as the truth --
residuals zeroed -- perturbs each by its own asymmetric error bars, re-minimises
**from the parent pose**, and reports the spread. Three things about it are
decisions rather than details, and each is pinned here with a measured number.

**The perturbation.** FPS draws `z ~ N(0,1)` and scales by `ErrPlus` when
`z > 0` and by `ErrMinus` when `z < 0` (`:57`). That gives each side mass 1/2
whatever the two widths are, so the density *jumps* by `ErrMinus/ErrPlus` at the
target. The object the two error bars actually describe is the two-piece (split)
normal, which weights the sides by their own widths and is continuous. Both are
implemented; `BFF_SPLIT_NORMAL` is the default. Fixing the discontinuity does
**not** remove the outward drift -- it doubles it, from `(s+ - s-)/sqrt(2pi)` to
`2(s+ - s-)/sqrt(2pi)` -- and that is the point worth knowing before choosing,
so it is asserted rather than mentioned.

**D2, the zero-noise pins.** FPS rewrites the target of *every* distance but
perturbs only the selected ones, and error estimation ships
`OptimizeSelected = All`. A deselected distance therefore enters the re-fit as a
restraint that is already satisfied exactly, holding the replica at the parent
pose and biasing the estimate down. `BootstrapParameters.perturb_deselected`
defaults to fixing it.

**The RMSD sign.** `SimulationResult.RMSD` accumulates `|U*r - t|^2` where the
displacement it derived is `U*r + t` (`SimulationResult.cs:56-57,67`). The two
differ by `(4/N) sum_i (U r_i) . t`, which is exactly zero when the body-local
coordinates sum to zero -- so the size of FPS's error is set by how far each
body's frame origin sits from the unweighted centroid of its atoms. That
identity is checked here to machine precision, which is stronger than checking
that the two numbers merely differ.
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
HIV = REPO / "examples" / "structure" / "HIV_RT"
PDBS = [str(HIV / "protein_1R0A.pdb"), str(HIV / "dna.pdb")]
FPS_JSON = str(HIV / "hiv_rt.fps.json")

#: An ordinary asymmetric FRET error bar, and the one the reference page uses.
ERR_NEG, ERR_POS = 3.0, 10.0
#: Enough draws that a 1 % check on a mean of ~5.6 A is not a coin toss.
N_DRAWS = 400_000

INV_SQRT_2PI = 1.0 / math.sqrt(2.0 * math.pi)


def _draws(model, seed=7, n=N_DRAWS, err_neg=ERR_NEG, err_pos=ERR_POS):
    return list(IMP.bff.sample_distance_perturbations(err_neg, err_pos, n,
                                                      seed, model))


def _mean(values):
    return sum(values) / len(values)


def _sd(values):
    mean = _mean(values)
    return math.sqrt(sum((v - mean) ** 2 for v in values) / len(values))


# --------------------------------------------------------------------------
# The perturbation: the closed forms, and the draws that must match them
# --------------------------------------------------------------------------

def test_the_closed_forms_are_the_published_ones():
    """The numbers the reference page states, to four decimals."""
    fps = IMP.bff.FPS_SIGN_SPLIT_NORMAL
    split = IMP.bff.BFF_SPLIT_NORMAL
    # (ErrPlus - ErrMinus)/sqrt(2 pi) = 7/2.5066 = 2.7926 -- the +2.79 A the
    # reference page quotes for 10/3.
    assert IMP.bff.perturbation_mean(ERR_NEG, ERR_POS, fps) == pytest.approx(
        2.7926, abs=1e-4)
    assert IMP.bff.perturbation_sd(ERR_NEG, ERR_POS, fps) == pytest.approx(
        6.8338, abs=1e-4)
    assert IMP.bff.perturbation_upper_mass(ERR_NEG, ERR_POS, fps) == 0.5

    # Twice as much, because the wide side is now *chosen* more often as well
    # as being wider. The correct density is not the unbiased one; there is no
    # unbiased one.
    assert IMP.bff.perturbation_mean(ERR_NEG, ERR_POS, split) == pytest.approx(
        5.5852, abs=1e-4)
    assert IMP.bff.perturbation_sd(ERR_NEG, ERR_POS, split) == pytest.approx(
        6.9142, abs=1e-4)
    assert IMP.bff.perturbation_upper_mass(
        ERR_NEG, ERR_POS, split) == pytest.approx(10.0 / 13.0)


def test_the_split_normal_mean_shift_is_exactly_twice_fps():
    fps = IMP.bff.perturbation_mean(ERR_NEG, ERR_POS,
                                    IMP.bff.FPS_SIGN_SPLIT_NORMAL)
    split = IMP.bff.perturbation_mean(ERR_NEG, ERR_POS,
                                      IMP.bff.BFF_SPLIT_NORMAL)
    assert split == pytest.approx(2.0 * fps, rel=1e-12)
    assert split - fps == pytest.approx(2.7926, abs=1e-4)


@pytest.mark.parametrize("model", [IMP.bff.BFF_SPLIT_NORMAL,
                                   IMP.bff.FPS_SIGN_SPLIT_NORMAL])
def test_the_draws_match_their_own_closed_forms(model):
    """The engine's own RNG stream, measured against the algebra."""
    draws = _draws(model)
    assert _mean(draws) == pytest.approx(
        IMP.bff.perturbation_mean(ERR_NEG, ERR_POS, model), abs=0.05)
    assert _sd(draws) == pytest.approx(
        IMP.bff.perturbation_sd(ERR_NEG, ERR_POS, model), abs=0.05)
    upper = sum(1 for d in draws if d > 0) / len(draws)
    assert upper == pytest.approx(
        IMP.bff.perturbation_upper_mass(ERR_NEG, ERR_POS, model), abs=0.005)


@pytest.mark.parametrize("model,expected", [
    # FPS: each half gets mass 1/2, so the density at 0+ is (1/2)(2/s+)phi(0)
    # and at 0- is (1/2)(2/s-)phi(0). Their ratio is s-/s+ = 0.3.
    (IMP.bff.FPS_SIGN_SPLIT_NORMAL, 0.3),
    # The split normal is normalised once over both halves, so the two
    # one-sided densities at the mode are equal. That is what "continuous"
    # means, and it is the whole difference between the two models.
    (IMP.bff.BFF_SPLIT_NORMAL, 1.0),
])
def test_the_density_at_the_target_jumps_only_for_fps(model, expected):
    draws = _draws(model)
    width = 0.25
    above = sum(1 for d in draws if 0.0 < d < width)
    below = sum(1 for d in draws if -width < d < 0.0)
    assert above > 1000 and below > 1000  # enough to divide
    assert above / below == pytest.approx(expected, rel=0.06)


def test_symmetric_errors_make_the_two_models_agree():
    """With one error bar there is nothing to disagree about."""
    for model in (IMP.bff.BFF_SPLIT_NORMAL, IMP.bff.FPS_SIGN_SPLIT_NORMAL):
        assert IMP.bff.perturbation_mean(5.0, 5.0, model) == 0.0
        assert IMP.bff.perturbation_sd(5.0, 5.0, model) == pytest.approx(5.0)
        assert IMP.bff.perturbation_upper_mass(5.0, 5.0, model) == 0.5
    a = _draws(IMP.bff.BFF_SPLIT_NORMAL, n=20000, err_neg=5.0, err_pos=5.0)
    b = _draws(IMP.bff.FPS_SIGN_SPLIT_NORMAL, n=20000, err_neg=5.0,
               err_pos=5.0)
    assert _mean(a) == pytest.approx(0.0, abs=0.1)
    assert _mean(b) == pytest.approx(0.0, abs=0.1)


def test_a_seed_reproduces_a_stream_and_a_different_seed_does_not():
    a = _draws(IMP.bff.BFF_SPLIT_NORMAL, seed=11, n=200)
    b = _draws(IMP.bff.BFF_SPLIT_NORMAL, seed=11, n=200)
    c = _draws(IMP.bff.BFF_SPLIT_NORMAL, seed=12, n=200)
    assert a == b
    assert a != c


def test_a_zero_error_bar_gives_that_side_no_width():
    """A one-sided error bar is legal input; it must not become a NaN."""
    draws = _draws(IMP.bff.BFF_SPLIT_NORMAL, n=5000, err_neg=0.0, err_pos=4.0)
    assert all(d >= 0.0 for d in draws)
    assert max(draws) > 0.0


# --------------------------------------------------------------------------
# FPS's shipped parameters
# --------------------------------------------------------------------------

def test_the_shipped_error_estimation_parameters_are_fpss():
    """`ProjectData.cs:101-115`, and the one number that is deliberately not."""
    p = IMP.bff.fps_error_estimation_parameters()
    assert p.max_force == 10000.0        # the Huber tail is effectively off
    assert p.clash_tolerance == 0.5      # four times harder than docking
    assert p.optimize_selected == "All"
    # Error estimation re-optimises from the parent pose; there is no random
    # restart in this mode (`SimulationJobManager.cs:96-104`).
    assert p.shuffle_max_translation == 0.0
    # Not FPS's 100000: that counts damped-Verlet steps of an integrator this
    # module does not have, and `n_frames` counts conjugate-gradient steps.
    assert p.n_frames == IMP.bff.DockingParameters().n_frames
    # And not FPS's either: FPS's clash is all-atom over Bondi radii, the
    # coarse term is one 2.5 A bead per residue bought for docking's step cost,
    # and a mode whose pose moves by fractions of an angstrom cannot spend it.
    assert p.coarse_clash is False
    assert IMP.bff.DockingParameters().coarse_clash is True


def test_the_bootstrap_defaults_choose_the_correct_density_and_fix_d2():
    b = IMP.bff.BootstrapParameters()
    assert b.perturbation == IMP.bff.BFF_SPLIT_NORMAL
    assert b.perturb_deselected is True
    assert b.n_replicas == 10   # FPS's repetitions spinner


# --------------------------------------------------------------------------
# The RMSD, and FPS's sign error in it
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def assembly():
    return IMP.bff.create_docking_assembly(PDBS, FPS_JSON)


def _poses(assembly):
    return json.loads(IMP.bff.capture_poses(assembly))


def _local_sums(assembly):
    """Per body: the sum of its atoms' body-local coordinates, and how many.

    Zero only if the body's frame origin is the *unweighted* centroid. IMP puts
    it at the mass-weighted centre, exactly as FPS's `Molecule.CM` is
    mass-weighted, so it is not zero -- which is why FPS's sign error is a real
    difference and not an algebraic no-op.
    """
    model = assembly.model
    bodies = list(assembly.rigid_bodies)
    index = {b.get_particle_index(): i for i, b in enumerate(bodies)}
    sums = [IMP.algebra.Vector3D(0, 0, 0) for _ in bodies]
    counts = [0 for _ in bodies]
    for leaf in IMP.atom.get_leaves(assembly.root):
        p = leaf.get_particle_index()
        if not IMP.core.RigidBodyMember.get_is_setup(model, p):
            continue
        member = IMP.core.RigidBodyMember(model, p)
        k = index[member.get_rigid_body().get_particle_index()]
        sums[k] = sums[k] + member.get_internal_coordinates()
        counts[k] += 1
    return sums, counts


def test_a_pose_against_itself_is_zero(assembly):
    poses = IMP.bff.capture_poses(assembly)
    assert IMP.bff.pose_rmsd(assembly, poses, poses) == 0.0
    assert IMP.bff.pose_rmsd(assembly, poses, poses, True) == 0.0


def test_translating_one_body_gives_the_analytic_rmsd(assembly):
    """A pure translation of one body: sqrt(n_body/n_total) * |d|.

    No superposition is done and none should be: the two poses are already in
    one frame and the displacement between them is the whole quantity.
    """
    poses = _poses(assembly)
    moved = json.loads(json.dumps(poses))
    moved[1]["t"][0] += 10.0
    _, counts = _local_sums(assembly)
    expected = 10.0 * math.sqrt(counts[1] / float(sum(counts)))
    got = IMP.bff.pose_rmsd(assembly, json.dumps(moved), json.dumps(poses))
    assert got == pytest.approx(expected, rel=1e-9)
    # U = R_a - R_b is zero here, so `|Ur + t|` and `|Ur - t|` agree: a pure
    # translation cannot expose the sign error.
    assert IMP.bff.pose_rmsd(assembly, json.dumps(moved), json.dumps(poses),
                             True) == pytest.approx(got, rel=1e-12)


def test_fpss_rmsd_sign_error_is_exactly_four_t_dot_u_rbar(assembly):
    """rmsd_correct^2 - rmsd_fps^2 = (4/N) sum_i (U r_i) . t, measured.

    Both a rotation *and* a translation have to differ for the cross term to
    survive, which is why FPS's error hides in the common case of a pure
    displacement.
    """
    poses = _poses(assembly)
    turned = json.loads(json.dumps(poses))
    rotation = IMP.algebra.get_rotation_about_axis(
        IMP.algebra.Vector3D(0.0, 0.0, 1.0), 0.35)
    original = IMP.algebra.Rotation3D(*turned[1]["q"])
    combined = rotation * original
    turned[1]["q"] = list(combined.get_quaternion())
    turned[1]["t"][0] += 7.0
    turned[1]["t"][2] -= 3.0

    a, b = json.dumps(turned), json.dumps(poses)
    correct = IMP.bff.pose_rmsd(assembly, a, b, False)
    fps = IMP.bff.pose_rmsd(assembly, a, b, True)
    assert correct != fps

    sums, counts = _local_sums(assembly)
    total = float(sum(counts))
    cross = 0.0
    for k, pose_a in enumerate(turned):
        r_a = IMP.algebra.Rotation3D(*pose_a["q"])
        r_b = IMP.algebra.Rotation3D(*poses[k]["q"])
        u_s = r_a.get_rotated(sums[k]) - r_b.get_rotated(sums[k])
        t = (IMP.algebra.Vector3D(*pose_a["t"])
             - IMP.algebra.Vector3D(*poses[k]["t"]))
        cross += t * u_s
    assert correct ** 2 - fps ** 2 == pytest.approx(4.0 * cross / total,
                                                    rel=1e-8)


def test_a_superposition_is_computed_not_remembered(assembly):
    """`pose_superposition` puts one pose onto another to machine precision.

    FPS never computes this in any `Save*` method -- it reads whatever the GUI
    last left in `BestFitRotation`, which under the default "RMSD vs previous"
    is a chain of pairwise fits and, with the best-fit box unticked, an all-zero
    matrix that `AngleAndAxis` reads as a 180 degree rotation.
    """
    poses = _poses(assembly)
    moved = json.loads(json.dumps(poses))
    rotation = IMP.algebra.get_rotation_about_axis(
        IMP.algebra.Vector3D(1.0, 1.0, 0.0).get_unit_vector(), 0.6)
    for pose in moved:
        original = IMP.algebra.Rotation3D(*pose["q"])
        pose["q"] = list((rotation * original).get_quaternion())
        pose["t"] = list(rotation.get_rotated(IMP.algebra.Vector3D(*pose["t"]))
                         + IMP.algebra.Vector3D(2.0, -1.0, 0.5))
    fit = IMP.bff.pose_superposition(assembly, json.dumps(moved),
                                     json.dumps(poses))
    # A rigid motion of the whole assembly is undone exactly.
    recovered = fit.get_rotation() * rotation
    assert IMP.algebra.get_axis_and_angle(recovered)[1] == pytest.approx(
        0.0, abs=1e-6)
    assert IMP.bff.pose_rmsd(assembly, json.dumps(poses),
                             json.dumps(poses)) == 0.0
    # Fitting a pose onto itself is the identity, never a 180 degree rotation.
    same = IMP.bff.pose_superposition(assembly, json.dumps(poses),
                                      json.dumps(poses))
    assert IMP.algebra.get_axis_and_angle(same.get_rotation())[1] == \
        pytest.approx(0.0, abs=1e-9)
    # 1e-4 A, not zero: the fit is a Kabsch through an SVD over 9034 atoms
    # spread across ~150 A, and the residual translation of a pose onto itself
    # measures 3.4e-6 A -- floating-point noise, two parts in 10^8 of the
    # structure's own size, and nothing a downstream reader can see.
    assert same.get_translation().get_magnitude() == pytest.approx(0.0,
                                                                   abs=1e-4)


# --------------------------------------------------------------------------
# The bootstrap end to end
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def bootstrap(tmp_path_factory):
    """One short run on HIV-RT: the parent pose, three perturbed re-fits.

    `steps` is small on purpose. What is being asserted is the *protocol* --
    that the truth is the parent's own model distances, that the replicas start
    from the parent pose and end somewhere else, and that the run is
    reproducible from its seed -- none of which needs a converged minimisation.
    """
    out = tmp_path_factory.mktemp("bootstrap")
    params = IMP.bff.fps_error_estimation_parameters()
    params.score_set = "resolved"
    params.n_frames = 40
    boot = IMP.bff.BootstrapParameters()
    boot.n_replicas = 3
    boot.seed = 5
    assembly = IMP.bff.create_docking_assembly(PDBS, FPS_JSON,
                                              score_set="resolved")
    parent = IMP.bff.capture_poses(assembly)
    return IMP.bff.fps_bootstrap(PDBS, FPS_JSON, str(out), parent, params,
                                 boot)


def test_the_truth_is_the_parents_own_model_distances(bootstrap):
    """Residuals zeroed: what the replicas are drawn around is the model."""
    result = bootstrap
    truth = list(result.truth)
    # Every distance in the file, not only the score set's: FPS rewrites `R`
    # for all of them (`ErrorEstimation.cs:23-29`) and gates only the *force*
    # on the selection.
    assert len(truth) == 20
    finite = [p for p in truth if p.distance_model == p.distance_model]
    assert len(finite) >= 18
    # And it is *not* the experiment -- if it were, this would be a re-fit of
    # the data rather than a bootstrap.
    assert any(abs(p.distance_model - p.distance_exp) > 1.0 for p in finite)


def test_every_replica_moved_away_from_the_parent(bootstrap):
    result = bootstrap
    assert len(result.replicas) == 3
    rmsds = [r.rmsd_to_parent for r in result.replicas]
    assert all(r > 0.0 for r in rmsds)
    assert result.rmsd_mean == pytest.approx(sum(rmsds) / len(rmsds))
    assert result.rmsd_max == pytest.approx(max(rmsds))
    # The sample standard deviation, n-1: three replicas are draws, not a
    # population.
    mean = result.rmsd_mean
    expected = math.sqrt(sum((r - mean) ** 2 for r in rmsds) / 2.0)
    assert result.rmsd_sd == pytest.approx(expected, rel=1e-9)


def test_the_run_is_reproducible_from_its_seed(bootstrap, tmp_path):
    params = IMP.bff.fps_error_estimation_parameters()
    params.score_set = "resolved"
    params.n_frames = 40
    boot = IMP.bff.BootstrapParameters()
    boot.n_replicas = 3
    boot.seed = 5
    again = IMP.bff.fps_bootstrap(PDBS, FPS_JSON, str(tmp_path / "again"),
                                  bootstrap.parent_poses, params, boot)
    for a, b in zip(bootstrap.replicas, again.replicas):
        assert a.rmsd_to_parent == pytest.approx(b.rmsd_to_parent, rel=1e-9)


def test_the_settings_are_recorded_in_the_result(bootstrap):
    extra = json.loads(bootstrap.extra)
    assert extra["perturbation"] == "split_normal"
    assert extra["perturb_deselected"] is True
    assert extra["seed"] == 5
    assert extra["optimize_selected"] == "All"
    assert extra["max_force"] == 10000.0


def test_d2_leaves_the_deselected_distances_as_zero_noise_pins(tmp_path):
    """`resolved` leaves two of the twenty out; under FPS they become pins.

    The count is the assertion, not the resulting spread: with only three
    replicas the spread is too noisy to test, but *which distances carry noise*
    is exact.
    """
    params = IMP.bff.fps_error_estimation_parameters()
    params.score_set = "resolved"
    params.n_frames = 20
    assembly = IMP.bff.create_docking_assembly(PDBS, FPS_JSON,
                                              score_set="resolved")
    parent = IMP.bff.capture_poses(assembly)

    fixed = IMP.bff.BootstrapParameters()
    fixed.n_replicas = 1
    fixed.perturb_deselected = True
    a = IMP.bff.fps_bootstrap(PDBS, FPS_JSON, str(tmp_path / "fixed"), parent,
                              params, fixed)
    # A distance with no model value has no truth to perturb around and is
    # neither noisy nor a pin, so the totals are counted from the data rather
    # than from the file's line count.
    scorable = {p.name for p in a.truth
                if p.distance_model == p.distance_model}
    assert (a.n_perturbed, a.n_pinned) == (len(scorable), 0)

    like_fps = IMP.bff.BootstrapParameters()
    like_fps.n_replicas = 1
    like_fps.perturb_deselected = False
    b = IMP.bff.fps_bootstrap(PDBS, FPS_JSON, str(tmp_path / "fps"), parent,
                              params, like_fps)
    with open(FPS_JSON) as handle:
        resolved = set(json.load(handle)["\u03c7\u00b2"]["resolved"]["distances"])
    assert b.n_perturbed == len(scorable & resolved)
    assert b.n_pinned == len(scorable - resolved)
    assert b.n_pinned > 0, "the score set must leave something out to pin"
    assert json.loads(b.extra)["perturb_deselected"] is False


def test_the_clash_term_can_pin_the_replicas_and_the_result_says_so(tmp_path):
    """The number to read before the spread: how much of the score is clash.

    FPS's error-estimation `ClashTolerance = 0.5` gives
    `k_clash = 2/0.5^2 = 8`, and IMP's united-atom radii are larger than the
    Bondi set FPS uses, so a genuine protein-DNA interface registers as a large
    overlap. The pose is then held by geometry and the data cannot move it --
    the bootstrap reports ~0 A, correctly, and a reader who takes that for
    "the FRET network locates this subunit to a thousandth of an angstrom" has
    read it backwards.

    Measured on HIV-RT `resolved`, parent = the docked pose, 10 replicas,
    300 steps: clash 102.92 of a parent score of 134.80 gives 0.000 +- 0.000 A;
    with the excluded volume off, 14.83 of which 0 is clash, gives
    1.935 +- 1.133 A. This runs the two ends of that table cheaply.

    **Switching the accessible volume to Olga's radii (2026-09-01,
    `AV::set_radii_source`) does not fix this**, and it was worth checking
    rather than assuming: the radii set is the *volume's*, and the clash term
    is `clash_container`, which reads `IMP::core::XYZR` -- the radii the
    particles carry -- so it never saw the change. Measured on the cheap
    configuration below: clash **110.98 of 150.69** with the volume on IMP's
    radii, **110.85 of 140.58** with it on Olga's, and 0.000 +- 0.000 A either
    way. The free run does move, 2.130 +- 1.458 -> 1.216 +- 0.283 A, because
    the volumes moved.

    That asymmetry is the reason the volume's default went back to IMP's radii
    the same day: a volume on Olga's table beside a clash term on IMP's is one
    score whose two halves disagree about how big an atom is, and the half
    that dominates here is the one Olga's table never reached.

    That the clash *is* a radii artefact is nonetheless measurable: statically,
    at the input pose over the protein-DNA interface with k = 8, united-atom
    radii give 268 overlapping atom pairs / 90.5 A of overlap / energy 198.4,
    and Olga's give 30 / 7.6 A / 11.1 -- a factor of 17.9. Giving
    `clash_container` a radii source of its own is PRD-121 open item 3's
    remaining half; it is a decision about the docking protocol, not about the
    volume, so it is not made here -- and it is the *only* way the two could
    both move to Olga's numbers without disagreeing.
    """
    def spread(tag, ev_weight):
        params = IMP.bff.fps_error_estimation_parameters()
        params.score_set = "resolved"
        params.n_frames = 60
        params.ev_weight = ev_weight
        parent = IMP.bff.dock_minimize(PDBS, FPS_JSON,
                                       str(tmp_path / (tag + "_parent")),
                                       params)
        boot = IMP.bff.BootstrapParameters()
        boot.n_replicas = 3
        boot.seed = 3
        return IMP.bff.fps_bootstrap(PDBS, FPS_JSON, str(tmp_path / tag),
                                     parent.poses, params, boot)

    pinned = spread("pinned", 1.0)
    free = spread("free", 0.0)

    # The clash term carries most of the score in one and none in the other,
    # and `parent_e_clash` is what says which is which.
    assert pinned.parent_e_clash > 0.5 * pinned.parent_score
    assert free.parent_e_clash == 0.0
    # And that is what decides whether the data can move the pose at all.
    assert pinned.rmsd_mean < 0.05
    assert free.rmsd_mean > 10.0 * max(pinned.rmsd_mean, 1e-3)
    assert json.loads(free.get_json())["parent_e_clash"] == 0.0
