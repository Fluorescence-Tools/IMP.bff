"""One CHARMM36 table, one LJ kernel: every scorer derives from Scoring.h.

Two tables used to live side by side (``rmin_half``/``epsilon`` in
topology.dye and ``p_Rmin2``/``eps`` in rotamer.scoring); the numbers were
equal only by luck. These tests pin the single source and prove the scorers
compute the same energy on identical explicit pairs. The table itself is a
C++ function now -- `charmm36_lj(element)` returning `(rmin_half, epsilon)` --
so there is no second dict to drift.
"""

import numpy as np
import pytest

from IMP.bff import (
    charmm36_lj, lj_cross, lj_energy, lj_parameter_arrays, lj_pairs_sum,
    lj_score, scaled_parameters)


def test_table_is_the_only_one():
    import IMP.bff as rs
    assert not hasattr(rs, "LJ_PARAMETERS")
    for elem in ("C", "N", "O", "S", "H"):
        arrs = scaled_parameters([elem], 1.0, 1.0)
        assert arrs.rmin_half[0] == charmm36_lj(elem)[0]
        assert arrs.epsilon[0] == charmm36_lj(elem)[1]
    arrays = lj_parameter_arrays(["C", "X", "S"])
    assert arrays.rmin_half[1] == charmm36_lj("C")[0]   # unknown -> carbon
    assert arrays.epsilon[2] == charmm36_lj("S")[1]


def test_lj_energy_analytic():
    rmin, eps = 3.0, 0.2
    assert lj_energy([rmin], [rmin], [eps])[0] == pytest.approx(-eps)   # well
    assert lj_energy([rmin * 2 ** (-1 / 6)], [rmin], [eps])[0] == \
        pytest.approx(0.0, abs=1e-12)                                  # zero
    assert lj_energy([10.0], [rmin], [eps], repulsive_only=True)[0] == 0.0
    assert lj_energy([2.0], [rmin], [eps], cutoff=1.5)[0] == 0.0
    r = np.array([1.0, 2.0, 3.0, 4.0])
    e = np.asarray(lj_energy(r, np.full(4, rmin), np.full(4, eps)))
    assert e.shape == r.shape
    assert e[0] > 0 > e[3]


def test_scalar_and_pairwise_kernels_agree():
    rng = np.random.default_rng(3)
    elems_a = ["C", "N", "O", "S", "H", "C"]
    elems_b = ["C", "C", "N", "O", "H"]
    xa = rng.uniform(0, 6, size=(len(elems_a), 3))
    xb = rng.uniform(0, 6, size=(len(elems_b), 3))
    rmin = np.zeros(len(elems_a) * len(elems_b))
    eps = np.zeros_like(rmin)
    k = 0
    for ea in elems_a:
        for eb in elems_b:
            rmin[k], eps[k] = lj_cross(ea, eb)
            k += 1
    # scalar, repulsive-only, pair by pair
    e_scalar = 0.0
    k = 0
    for i in range(len(elems_a)):
        for j in range(len(elems_b)):
            r = float(np.linalg.norm(xa[i] - xb[j]))
            e_scalar += lj_score(r, rmin[k], eps[k])
            k += 1
    e_mf = lj_pairs_sum(xa.ravel(), xb.ravel(), rmin, eps, r_cutoff=1e9)
    assert e_mf == pytest.approx(e_scalar, rel=1e-12, abs=1e-12)


def test_rotamer_scoring_uses_the_shared_kernel():
    """The end-to-end score's 12-6 energy equals lj_energy on the same pairs."""
    from IMP.bff import compute_rotamer_score
    rng = np.random.default_rng(5)
    n_prot, n_rot = 12, 5
    protein_coords = rng.uniform(0, 8, size=(n_prot, 3))
    protein_names = ["CA", "N", "O", "SD", "CB", "C", "CG", "NZ", "OG", "SG", "CD", "C"]
    protein_resnames = ["ALA"] * n_prot
    rotamer_coords = rng.uniform(0, 8, size=(1, n_rot, 3))
    rotamer_names = ["C1", "N1", "O1", "C2", "S1"]
    res = compute_rotamer_score(
        rotamer_coords, protein_coords, protein_names, protein_resnames, rotamer_names,
        rotamer_weights=np.ones(1), temperature=300.0, electrostatic=False,
        potential="lj", sigma_scaling=1.0, epsilon_scaling=1.0, ignore_h=True,
    )
    # compute_rotamer_score masks hydrogens on both sides and backbone names
    # (CA, C, N, O) on the rotamer side only (the protein keeps its backbone;
    # only the labelled residue is masked, and no site is given here).
    keep_p = list(range(len(protein_names)))
    keep_r = [i for i, n in enumerate(rotamer_names) if n.upper() not in {"CA", "C", "N", "O"}]
    p_arrs = scaled_parameters([protein_names[i][0] for i in keep_p], 1.0, 1.0)
    r_arrs = scaled_parameters([rotamer_names[i][0] for i in keep_r], 1.0, 1.0)
    p_rmin, p_eps = np.asarray(p_arrs.rmin_half), np.asarray(p_arrs.epsilon)
    r_rmin, r_eps = np.asarray(r_arrs.rmin_half), np.asarray(r_arrs.epsilon)
    d = np.linalg.norm(rotamer_coords[0][keep_r][:, None, :] - protein_coords[keep_p][None, :, :], axis=2)
    expected = np.asarray(lj_energy(
        d.ravel(), np.add.outer(r_rmin, p_rmin).ravel(),
        np.sqrt(np.multiply.outer(r_eps, p_eps)).ravel(), cutoff=10.0)).sum()
    assert res.energies[0] == pytest.approx(expected, rel=1e-12, abs=1e-12)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p no:cacheprovider"]))
