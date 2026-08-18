"""The rotamer-protein interaction energy: the all-pairs inner loop, in C++.

Scoring a rotamer library is every dye atom against every protein atom, for
every conformer. The Python it replaced built an ``(n_dye, n_protein)`` distance
matrix per conformer -- so it allocated once per rotamer and took a square root
for every pair before discovering that almost all of them were beyond the 10 A
cutoff.

The kernel checks the squared distance first and takes the root only for pairs
that pass. 200 conformers x 30 dye atoms x 3000 protein atoms: 621 ms -> 78 ms.

These tests are the equality gate. The reference below is the pre-port numpy,
kept verbatim, because "same answer" is the only claim a port gets to make.
"""

import numpy as np
import pytest

import IMP.bff
from IMP.bff.scoring.lennard_jones import lj_energy


def _reference(rot, prot, rmin_ij, eps_ij, q_rot, q_prot, potential, electrostatic):
    """The loop as it was written in Python, unchanged."""
    n = rot.shape[0]
    pot = np.zeros(n)
    dh = np.zeros(n)
    for i, coords in enumerate(rot):
        d = np.linalg.norm(coords[:, None, :] - prot[None, :, :], axis=2)
        if electrostatic:
            m = (q_rot[:, None] * q_prot[None, :] != 0.0) & (d < 20.0)
            if np.any(m):
                ri, pj = np.nonzero(m)
                dh[i] += np.sum(q_rot[ri] * q_prot[pj] * 7.0 / d[ri, pj]
                                * np.exp(-d[ri, pj] / 10.0))
        m = d < 10.0
        if np.any(m):
            if potential == IMP.bff.ROTAMER_POTENTIAL_LJ:
                pot[i] += np.sum(lj_energy(d[m], rmin_ij[m], eps_ij[m], r_floor=0.0))
            else:
                pot[i] += np.sum(eps_ij[m] * np.exp(-0.5 * np.power(d[m] / rmin_ij[m], 2)))
    return pot, dh


def _call(rot, prot, rmin_ij, eps_ij, q_rot, q_prot, potential, electrostatic):
    n_rot, n_dye, _ = rot.shape
    out = IMP.bff.rotamer_interaction_energies(
        np.ascontiguousarray(rot).ravel(), np.ascontiguousarray(prot).ravel(),
        np.ascontiguousarray(rmin_ij).ravel(), np.ascontiguousarray(eps_ij).ravel(),
        np.ascontiguousarray(q_rot) if electrostatic else np.empty(0),
        np.ascontiguousarray(q_prot) if electrostatic else np.empty(0),
        n_rot, n_dye, prot.shape[0], potential)
    return np.asarray(out, dtype=np.float64).reshape(-1, 2)


def _case(seed, n_rot=25, n_dye=20, n_prot=400):
    rng = np.random.default_rng(seed)
    return (rng.normal(0, 6, (n_rot, n_dye, 3)),
            rng.normal(0, 9, (n_prot, 3)),
            rng.uniform(2.0, 4.5, (n_dye, n_prot)),
            rng.uniform(0.02, 0.3, (n_dye, n_prot)),
            rng.choice([-1.0, 0.0, 1.0], n_dye),
            rng.choice([-1.0, 0.0, 1.0], n_prot))


@pytest.mark.parametrize("potential", [IMP.bff.ROTAMER_POTENTIAL_LJ,
                                       IMP.bff.ROTAMER_POTENTIAL_GAUSS])
@pytest.mark.parametrize("electrostatic", [False, True])
def test_matches_the_python_it_replaced(potential, electrostatic):
    rot, prot, rmin, eps, qr, qp = _case(3)
    want_pot, want_dh = _reference(rot, prot, rmin, eps, qr, qp, potential, electrostatic)
    got = _call(rot, prot, rmin, eps, qr, qp, potential, electrostatic)
    np.testing.assert_allclose(got[:, 0], want_pot, rtol=1e-11, atol=1e-12)
    np.testing.assert_allclose(got[:, 1], want_dh, rtol=1e-11, atol=1e-12)


def test_electrostatics_are_off_when_no_charges_are_given():
    """An empty charge array is the switch, not a flag that can disagree with it."""
    rot, prot, rmin, eps, qr, qp = _case(5)
    got = _call(rot, prot, rmin, eps, qr, qp, IMP.bff.ROTAMER_POTENTIAL_LJ, False)
    assert np.all(got[:, 1] == 0.0)


def test_a_pair_beyond_the_cutoff_contributes_nothing():
    """The cutoff is on distance, and the kernel tests it squared to avoid a root."""
    rot = np.array([[[0.0, 0.0, 0.0]]])
    near = np.array([[9.9, 0.0, 0.0]])
    far = np.array([[10.1, 0.0, 0.0]])
    rmin = np.array([[3.0]])
    eps = np.array([[0.2]])
    n = _call(rot, near, rmin, eps, None, None, IMP.bff.ROTAMER_POTENTIAL_LJ, False)
    f = _call(rot, far, rmin, eps, None, None, IMP.bff.ROTAMER_POTENTIAL_LJ, False)
    assert n[0, 0] != 0.0
    assert f[0, 0] == 0.0


def test_the_lj_well_has_its_minimum_at_rmin():
    """-eps at r = rmin, which is the convention lj_cross produces."""
    rmin, eps = 3.5, 0.25
    rot = np.array([[[0.0, 0.0, 0.0]]])
    at_min = _call(rot, np.array([[rmin, 0.0, 0.0]]), np.array([[rmin]]),
                   np.array([[eps]]), None, None, IMP.bff.ROTAMER_POTENTIAL_LJ, False)
    assert at_min[0, 0] == pytest.approx(-eps)
    closer = _call(rot, np.array([[rmin * 0.8, 0.0, 0.0]]), np.array([[rmin]]),
                   np.array([[eps]]), None, None, IMP.bff.ROTAMER_POTENTIAL_LJ, False)
    assert closer[0, 0] > at_min[0, 0], "the repulsive wall must rise"


def test_screened_coulomb_falls_off_faster_than_bare():
    """The Debye-Huckel factor exp(-r/10) is what makes it screened."""
    rot = np.array([[[0.0, 0.0, 0.0]]])
    q = np.array([1.0])
    vals = []
    for r in (3.0, 6.0, 12.0):
        got = _call(rot, np.array([[r, 0.0, 0.0]]), np.array([[3.0]]), np.array([[0.0]]),
                    q, q, IMP.bff.ROTAMER_POTENTIAL_LJ, True)
        vals.append(got[0, 1])
        assert got[0, 1] == pytest.approx(7.0 / r * np.exp(-r / 10.0))
    assert vals[0] > vals[1] > vals[2] > 0.0


def test_empty_inputs_do_not_crash():
    empty = np.empty(0)
    assert list(IMP.bff.rotamer_interaction_energies(
        empty, empty, empty, empty, empty, empty, 0, 0, 0,
        IMP.bff.ROTAMER_POTENTIAL_LJ)) == []


def test_the_scorer_still_produces_normalised_weights():
    """The kernel feeds compute_rotamer_score; its contract is unchanged."""
    from IMP.bff.scoring.rotamer import compute_rotamer_score
    rng = np.random.default_rng(11)
    n_rot, n_dye, n_prot = 12, 6, 40
    result = compute_rotamer_score(
        rotamer_coords=rng.normal(0, 5, (n_rot, n_dye, 3)),
        protein_coords=rng.normal(0, 8, (n_prot, 3)),
        protein_atom_names=["CA"] * n_prot,
        protein_resnames=["ALA"] * n_prot,
        rotamer_atom_names=["C"] * n_dye,
    )
    assert result.weights.shape == (n_rot,)
    assert result.weights.sum() == pytest.approx(1.0)
    assert np.all(result.weights >= 0.0)


def test_an_unknown_potential_is_refused():
    from IMP.bff.scoring.rotamer import compute_rotamer_score
    rng = np.random.default_rng(1)
    with pytest.raises(ValueError, match="Unknown potential"):
        compute_rotamer_score(
            rotamer_coords=rng.normal(0, 5, (3, 4, 3)),
            protein_coords=rng.normal(0, 8, (10, 3)),
            protein_atom_names=["CA"] * 10, protein_resnames=["ALA"] * 10,
            rotamer_atom_names=["C"] * 4, potential="quartic")


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
