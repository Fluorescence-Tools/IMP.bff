"""The coarse-grained potentials, against the Python they were ported from.

The kernels were `numba.njit` functions in `IMP.cgmol.{statpot,sterics,
solvation}` (imp-tricks) and the classes around them lived in chisurf's
`structure/potential/potentials.py`. The numbers below were produced by running
**those** kernels, unmodified, on the deterministic structure `_structure()`
builds -- a loose helix of twelve residues with seeded jitter, six atoms each.

Two of them differ from the Python on purpose, and say so where they are
asserted.
"""

import math

import numpy as np
import pytest

import IMP
import IMP.atom
import IMP.container
import IMP.core
import IMP.bff

N_RES = 12
N_SITES = 6                      # N, CA, C, O, CB, H
SITE_NAMES = ["N", "CA", "C", "O", "CB", "H"]
#: chisurf's `res2id` order, which is Miyazawa and Jernigan's.
MJ_ORDER = ["CYS", "MET", "PHE", "ILE", "LEU", "VAL", "TRP", "TYR", "ALA",
            "GLY", "THR", "SER", "GLN", "ASN", "GLU", "ASP", "HIS", "ARG",
            "LYS", "PRO"]
VDW = [1.55, 1.7, 1.7, 1.52, 1.7, 1.2]
CHARGES = [-0.4, 0.1, 0.5, -0.5, 0.0, 0.3]


def _coords():
    """The reference coordinates: seeded, so the pins mean something."""
    rng = np.random.default_rng(20260826)
    xyz = np.zeros((N_RES * N_SITES, 3))
    for i in range(N_RES):
        t = i * 1.6
        base = np.array([2.3 * math.cos(t), 2.3 * math.sin(t), 1.5 * i])
        for s in range(N_SITES):
            xyz[i * N_SITES + s] = base + rng.normal(scale=1.1, size=3)
    return np.ascontiguousarray(xyz)


#: The residue types the reference used: `(i * 7) % 20` into the MJ order.
RES_TYPES = [(i * 7) % 20 for i in range(N_RES)]


def _structure(missing=True):
    """The same atoms as an IMP hierarchy.

    With `missing`, residue 3 has no amide H (a proline) and residue 7 no
    C-beta (a glycine) -- the two holes the reference punched in its lookup
    table, and the two a real structure has.
    """
    xyz = _coords()
    m = IMP.Model()
    root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(m, "root"))
    chain = IMP.atom.Chain.setup_particle(IMP.Particle(m, "A"), "A")
    root.add_child(chain)
    atoms = []
    for i in range(N_RES):
        name = MJ_ORDER[RES_TYPES[i]]
        res = IMP.atom.Residue.setup_particle(
            IMP.Particle(m, f"{name}{i}"), IMP.atom.ResidueType(name), i + 1)
        chain.add_child(res)
        for s, site in enumerate(SITE_NAMES):
            if missing and i == 3 and site == "H":
                atoms.append(None)
                continue
            if missing and i == 7 and site == "CB":
                atoms.append(None)
                continue
            p = IMP.Particle(m, site)
            a = IMP.atom.Atom.setup_particle(p, IMP.atom.AtomType(site))
            IMP.core.XYZR.setup_particle(
                p, IMP.algebra.Sphere3D(
                    IMP.algebra.Vector3D(*xyz[i * N_SITES + s]), VDW[s]))
            # IMP.atom.Atom.setup_particle already gave it a Mass.
            IMP.atom.Charged.setup_particle(p, CHARGES[s])
            res.add_child(a)
            atoms.append(p)
    return m, root, atoms


# ---------------------------------------------------------------------------
# The kernels
# ---------------------------------------------------------------------------

def test_clash_energy_matches_the_python():
    xyz = _coords().ravel()
    vdw = np.array(VDW * N_RES)
    assert IMP.bff.clash_energy(xyz, vdw) == pytest.approx(
        35.961679525439735, rel=1e-12)
    assert IMP.bff.clash_energy(xyz, vdw, 6.0, 1.0) == pytest.approx(
        6.521641123970415, rel=1e-12)


def test_go_energy_matches_the_python():
    """A Go model scored on the structure it was built from is at its floor.

    Every well's minimum *is* the distance that structure has, so no pair is
    in the `rm < r < 2.5 rm` branch and the energy is minus the sum of the
    depths. The Python in chisurf got a different number here, because its
    accumulator carried across pairs; imp-tricks' kernel is the one ported.
    """
    ca = np.ascontiguousarray(_coords()[1::N_SITES]).ravel()
    contacts = IMP.bff.go_native_contacts(ca, 1.0, 0.1, 6.5)
    assert contacts.get_n_contacts() == N_RES * (N_RES - 1) // 2
    assert IMP.bff.go_energy(ca, contacts) == pytest.approx(
        -6.676356712583448, rel=1e-12)
    assert contacts.well_depths.sum() == pytest.approx(6.676356712583449,
                                                       rel=1e-12)


def test_lennard_jones_bead_energy_matches_the_python():
    ca = np.ascontiguousarray(_coords()[1::N_SITES]).ravel()
    assert IMP.bff.lennard_jones_bead_energy(ca) == pytest.approx(
        91590.40855879632, rel=1e-12)


def test_residue_asa_matches_the_python():
    ca = np.ascontiguousarray(_coords()[1::N_SITES]).ravel()
    assert IMP.bff.residue_asa(ca, 64, 1.0, 2.5) == pytest.approx(
        288.3883881225005, rel=1e-12)


def test_generalized_born_matches_the_python():
    xyz = _coords().ravel()
    radii = np.array(VDW * N_RES)
    charges = np.array(CHARGES * N_RES)
    assert IMP.bff.generalized_born_energy(xyz, radii, charges) == \
        pytest.approx(-0.024505222258838378, rel=1e-11)


def test_ramachandran_reads_a_map_as_minus_log_p():
    n = 8
    grid = np.full((n, n), 0.5)
    # a uniform map is no information, so every angle scores zero
    assert IMP.bff.ramachandran_energy(0.0, 0.0, grid.ravel(), n) == \
        pytest.approx(0.0)
    grid[4, 4] = 1.0                       # the maximum
    assert IMP.bff.ramachandran_energy(0.0, 0.0, grid.ravel(), n) == \
        pytest.approx(0.0)                 # (0, 0) rad lands in bin (4, 4)
    grid[4, 4] = 0.0                       # an empty cell
    assert IMP.bff.ramachandran_energy(0.0, 0.0, grid.ravel(), n) == \
        pytest.approx(10.0)
    # and a map that never goes positive is pseudo-energies already
    assert IMP.bff.ramachandran_energy(0.0, 0.0, -grid.ravel(), n) == \
        pytest.approx(0.0)


def test_a_chain_end_has_no_dihedral_and_costs_nothing():
    grid = np.full(16, 0.5)
    assert IMP.bff.ramachandran_energy(float("nan"), 0.0, grid, 4) == 0.0


# ---------------------------------------------------------------------------
# The IMP layer: the same numbers, through particles
# ---------------------------------------------------------------------------

def test_go_restraint_scores_what_the_kernel_does():
    m, root, atoms = _structure()
    cas = [a.get_index() for i, a in enumerate(atoms)
           if a is not None and i % N_SITES == 1]
    r = IMP.bff.GoRestraint(m, cas, 1.0, 6.5, 0.1)
    assert r.get_n_native() + r.get_n_non_native() == \
        r.get_contacts().get_n_contacts()
    assert r.unprotected_evaluate(None) == pytest.approx(
        -6.676356712583448, rel=1e-12)


def test_generalized_born_restraint_scores_what_the_kernel_does():
    """Charges come off `IMP.atom.Charged` and radii off `IMP.core.XYZR`,
    which is where `IMP.atom.CoulombPairScore` reads them too."""
    m, root, atoms = _structure(missing=False)
    pis = [a.get_index() for a in atoms]
    r = IMP.bff.GeneralizedBornRestraint(m, pis)
    assert r.unprotected_evaluate(None) == pytest.approx(
        -0.024505222258838378, rel=1e-11)


def test_lennard_jones_bead_pair_score_sums_to_the_kernel():
    m, root, atoms = _structure(missing=False)
    cas = [a.get_index() for i, a in enumerate(atoms) if i % N_SITES == 1]
    score = IMP.bff.LennardJonesBeadPairScore()
    total = 0.0
    for i in range(len(cas)):
        for j in range(i + 1, len(cas)):
            total += score.evaluate_index(m, (cas[i], cas[j]), None)
    ca = np.ascontiguousarray(_coords()[1::N_SITES]).ravel()
    assert total == pytest.approx(IMP.bff.lennard_jones_bead_energy(ca),
                                  rel=1e-12)


def test_the_bead_score_has_a_gradient_that_matches_a_difference():
    m, root, atoms = _structure(missing=False)
    a, b = atoms[1].get_index(), atoms[7].get_index()
    score = IMP.bff.LennardJonesBeadPairScore()
    # Through a scoring function, so the model zeroes the derivatives first --
    # reading them after a bare `evaluate_index` reads whatever was there.
    lsc = IMP.container.ListSingletonContainer(m, [a, b])
    r = IMP.container.PairsRestraint(
        score, IMP.container.ClosePairContainer(lsc, 1e6, 0.0))
    IMP.core.RestraintsScoringFunction([r]).evaluate(True)
    analytic = IMP.core.XYZ(m, a).get_derivatives()[0]

    eps = 1e-6
    xyz = IMP.core.XYZ(m, a)
    x0 = xyz.get_coordinate(0)
    xyz.set_coordinate(0, x0 + eps)
    plus = score.evaluate_index(m, (a, b), None)
    xyz.set_coordinate(0, x0 - eps)
    minus = score.evaluate_index(m, (a, b), None)
    xyz.set_coordinate(0, x0)
    assert analytic == pytest.approx((plus - minus) / (2 * eps), rel=1e-4)


# ---------------------------------------------------------------------------
# Miyazawa-Jernigan, through IMP's statistical pair score
# ---------------------------------------------------------------------------

def _mj_matrix():
    """The reference's stand-in contact matrix, symmetric, 20 x 20."""
    a = np.array([[((i * 31 + j * 17) % 41 - 20) / 10.0 for j in range(20)]
                  for i in range(20)])
    return (a + a.T) / 2.0


def _write_mj_table(path, matrix, cutoff=6.5):
    """The matrix as IMP's PMF format: two bins of half the cutoff each.

    A contact potential is one value out to the cutoff and nothing beyond, so
    the table is flat. Two bins rather than one because the reader builds a
    spline over the values and a spline needs two.

    The type count is the *potential's*, not IMP's: the score reads
    `IMP.bff.ResidueContactType`, whose names are the ones this file names.
    """
    with open(path, "w") as f:
        f.write(f"{cutoff / 2.0} {len(MJ_ORDER)}\n")
        for i, a in enumerate(MJ_ORDER):
            for j, b in enumerate(MJ_ORDER):
                if j < i:
                    continue
                v = matrix[i, j]
                f.write(f"{a} {b} {v} {v}\n")


def test_miyazawa_jernigan_through_imp_matches_the_python(tmp_path):
    """21 contacts summing to 3.95 -- the same as the numba kernel.

    The Python looped over residue pairs itself, prefiltering on C-alpha
    distance and testing C-beta separation; here the loop is an
    `IMP.container.ClosePairContainer` over the C-betas and the arithmetic is
    `IMP.core.StatisticalPairScore`. The set of contacts is the same set.
    """
    table = tmp_path / "mj.txt"
    _write_mj_table(table, _mj_matrix())

    m, root, atoms = _structure()
    typed = IMP.bff.add_residue_type_score_data(root, IMP.atom.AT_CB)
    # residue 7 has no C-beta, so it takes no part
    assert len(typed) == N_RES - 1

    score = IMP.bff.MiyazawaJerniganPairScore(6.5, str(table))
    indices = [p.get_index() for p in typed]
    total = 0.0
    scored = 0
    for i in range(len(indices)):
        for j in range(i + 1, len(indices)):
            v = score.evaluate_index(m, (indices[i], indices[j]), None)
            if v != 0.0:
                scored += 1
            total += v
    assert scored == 21
    assert total == pytest.approx(3.95, abs=1e-9)

    # and the same through IMP's own pair loop, which is what a caller writes
    lsc = IMP.container.ListSingletonContainer(m, indices)
    cpc = IMP.container.ClosePairContainer(lsc, 6.5, 0.0)
    r = IMP.container.PairsRestraint(score, cpc)
    sf = IMP.core.RestraintsScoringFunction([r])
    assert sf.evaluate(False) == pytest.approx(3.95, abs=1e-9)


def test_the_residue_type_key_is_the_potentials_own_index():
    """The tables name residues and the reader maps names to indices, so
    nothing in this module carries a residue order of its own -- and the index
    space is the potential's, not IMP's, so it cannot grow past the table."""
    m, root, atoms = _structure()
    typed = IMP.bff.add_residue_type_score_data(root, IMP.atom.AT_CB)
    key = IMP.bff.get_residue_type_key()
    for p in typed:
        name = IMP.atom.Residue(
            IMP.atom.Atom(p).get_parent()).get_residue_type().get_string()
        assert p.get_value(key) == IMP.bff.ResidueContactType(name).get_index()
        assert p.get_value(key) < 20


# ---------------------------------------------------------------------------
# Hydrogen bonds
# ---------------------------------------------------------------------------

def _hbond_table(n_bins=400):
    """The reference's stand-in lookup: four channels of a sine."""
    return np.array([[math.sin(0.02 * b + 0.5 * c) for b in range(n_bins)]
                     for c in range(4)]).ravel()


def test_hydrogen_bond_restraint_matches_the_python():
    m, root, atoms = _structure()
    r = IMP.bff.HydrogenBondRestraint(m, root, _hbond_table(), 400, 8.0, 3.0,
                                      0.01)
    assert r.unprotected_evaluate(None) == pytest.approx(
        -0.1324432499910165, rel=1e-12)
    assert r.get_n_hbonds() == 3


def test_the_c_alpha_cutoff_is_applied():
    """The kernel this was ported from compared a plain C-alpha distance
    against the *square* of this cutoff, so it never excluded anything. At 8 A
    that changes nothing on this structure -- every bonded pair is closer than
    that -- and at 4 A it drops one bond, which is the proof it engages.

    See `okf/validation/hbond_ca_cutoff.md`.
    """
    m, root, atoms = _structure()
    tight = IMP.bff.HydrogenBondRestraint(m, root, _hbond_table(), 400, 4.0,
                                          3.0, 0.01)
    assert tight.unprotected_evaluate(None) == pytest.approx(
        0.32062226315756825, rel=1e-12)
    assert tight.get_n_hbonds() == 2


def test_switching_a_channel_off_changes_the_energy():
    m, root, atoms = _structure()
    r = IMP.bff.HydrogenBondRestraint(m, root, _hbond_table(), 400)
    everything = r.unprotected_evaluate(None)
    r.set_channels(True, True, False, True)      # no O-H channel
    assert r.unprotected_evaluate(None) != pytest.approx(everything)
    r.set_channels(False, False, False, False)
    assert r.unprotected_evaluate(None) == pytest.approx(0.0)
    assert r.get_n_hbonds() == 3                 # still found, just not scored


# ---------------------------------------------------------------------------
# Ramachandran
# ---------------------------------------------------------------------------

def test_ramachandran_restraint_skips_the_chain_ends():
    m, root, atoms = _structure()
    n_bins = 36
    grids = np.full((3, n_bins, n_bins), 0.5).ravel()
    r = IMP.bff.RamachandranRestraint(m, root, grids, 3, n_bins)
    phi, psi = r.phi, r.psi
    assert len(phi) == N_RES and len(psi) == N_RES
    assert math.isnan(phi[0]), "the first residue has no phi"
    assert math.isnan(psi[N_RES - 1]), "the last has no psi"
    assert np.isfinite(phi[1:]).all() and np.isfinite(psi[:-1]).all()
    # a uniform map is no information
    assert r.unprotected_evaluate(None) == pytest.approx(0.0)


def test_ramachandran_costs_the_penalty_for_an_empty_map():
    m, root, atoms = _structure()
    n_bins = 36
    grids = np.zeros((3, n_bins, n_bins))
    grids[:, 0, 0] = 1.0            # one populated cell, far from any angle
    r = IMP.bff.RamachandranRestraint(m, root, grids.ravel(), 3, n_bins, 10.0)
    # every residue with both dihedrals lands in an empty cell
    assert r.unprotected_evaluate(None) == pytest.approx(10.0 * (N_RES - 2))


# ---------------------------------------------------------------------------
# Built out of IMP's own scores
# ---------------------------------------------------------------------------

def test_the_clash_restraint_is_imps_soft_sphere_with_the_ported_constant():
    """`SoftSpherePairScore` is 0.5 k (sigma - r)^2 and the ported term is
    ((sigma - r)/t)^2, so k = 2/t^2 -- and with no bonds to filter and no pair
    below the covalent floor, the two agree exactly."""
    m, root, atoms = _structure(missing=False)
    r = IMP.bff.build_clash_restraint(root, 2.0)
    xyz = _coords().ravel()
    vdw = np.array(VDW * N_RES)
    kernel = IMP.bff.clash_energy(xyz, vdw, 2.0, 0.0)
    # through a scoring function: a container restraint needs the model to
    # bring its container up to date before it means anything
    got = IMP.core.RestraintsScoringFunction([r]).evaluate(False)
    assert got == pytest.approx(kernel, rel=1e-9)


def test_the_ca_internal_restraints_vanish_at_the_reference():
    m, root, atoms = _structure(missing=False)
    cas = [a.get_index() for i, a in enumerate(atoms) if i % N_SITES == 1]
    rs = IMP.bff.build_ca_internal_restraints(m, cas, cas, 1.0, 0.2, 0.1)
    assert len(rs) == (N_RES - 1) + (N_RES - 2) + (N_RES - 3)
    assert sum(r.unprotected_evaluate(None) for r in rs) == pytest.approx(0.0)


def test_the_ca_internal_restraints_grow_when_the_trace_moves():
    m, root, atoms = _structure(missing=False)
    cas = [a.get_index() for i, a in enumerate(atoms) if i % N_SITES == 1]
    rs = IMP.bff.build_ca_internal_restraints(m, cas, cas, 1.0, 0.0, 0.0)
    xyz = IMP.core.XYZ(m, cas[0])
    xyz.set_coordinate(0, xyz.get_coordinate(0) + 1.0)
    moved = sum(r.unprotected_evaluate(None) for r in rs)
    # one bond stretched by (up to) 1 A, scored 0.5 k d^2
    assert 0.0 < moved <= 0.5 * 1.0 * 1.0 + 1e-9


def test_residue_surface_of_a_real_structure_is_positive():
    m = IMP.Model()
    h = IMP.atom.read_pdb(IMP.bff.get_example_path("structure/T4L/3GUN.pdb"),
                          m, IMP.atom.CAlphaPDBSelector())
    area = IMP.bff.residue_solvent_accessible_surface(h, IMP.atom.AT_CA, 64)
    assert area > 0.0
    # buried residues have less surface than the whole set would if isolated
    n_res = len(IMP.atom.get_by_type(h, IMP.atom.RESIDUE_TYPE))
    assert area < n_res * 4.0 * math.pi * 2.5 ** 2
