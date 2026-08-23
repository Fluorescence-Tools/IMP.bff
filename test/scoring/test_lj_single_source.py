"""One CHARMM36 table, one LJ kernel: every cgdye scorer derives from topology.dye.

Two tables used to live side by side (``rmin_half``/``epsilon`` in
topology.dye and ``p_Rmin2``/``eps`` in rotamer.scoring); the numbers were
equal only by luck. These tests pin the single source and prove the three
scorers (sampling.scoring, sampling.mean_field, rotamer.scoring) compute the
same energy on identical explicit pairs.
"""

import numpy as np
import pytest

from IMP.bff import CHARMM36_LJ, lj_cross, lj_energy, lj_parameter_arrays
from IMP.bff import lj_score
from IMP.bff import _lj_energy_pairs
from IMP.bff import _scaled_parameters


def test_table_is_the_only_one():
    import IMP.bff as rs
    assert not hasattr(rs, "LJ_PARAMETERS")
    for elem in ("C", "N", "O", "S", "H"):
        rmin2, eps = _scaled_parameters([elem], 1.0, 1.0)
        assert rmin2[0] == CHARMM36_LJ[elem]["rmin_half"]
        assert eps[0] == CHARMM36_LJ[elem]["epsilon"]
    rmin_half, epsilon = lj_parameter_arrays(["C", "X", "S"])
    assert rmin_half[1] == CHARMM36_LJ["C"]["rmin_half"]  # unknown -> carbon
    assert epsilon[2] == CHARMM36_LJ["S"]["epsilon"]


def test_lj_energy_analytic():
    rmin, eps = 3.0, 0.2
    assert lj_energy(rmin, rmin, eps) == pytest.approx(-eps)          # well depth
    assert lj_energy(rmin * 2 ** (-1 / 6), rmin, eps) == pytest.approx(0.0, abs=1e-12)  # zero crossing
    assert lj_energy(10.0, rmin, eps, repulsive_only=True) == 0.0
    assert lj_energy(2.0, rmin, eps, cutoff=1.5) == 0.0
    r = np.array([1.0, 2.0, 3.0, 4.0])
    e = lj_energy(r, rmin, eps)
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
    e_mf = _lj_energy_pairs(xa, xb, rmin, eps, r_cutoff=1e9)
    assert e_mf == pytest.approx(e_scalar, rel=1e-12, abs=1e-12)


def test_rotamer_scoring_uses_the_shared_kernel():
    """rotamer.scoring's full 12-6 energy equals lj_energy on the same pairs."""
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
        {"weights": np.ones(1)}, temperature=300.0, electrostatic=False, potential="lj",
        sigma_scaling=1.0, epsilon_scaling=1.0, ignore_h=True,
    )
    # compute_rotamer_score masks hydrogens on both sides and backbone names
    # (CA, C, N, O) on the rotamer side only (the protein keeps its backbone;
    # only the labelled residue is masked, and no site is given here).
    keep_p = list(range(len(protein_names)))
    keep_r = [i for i, n in enumerate(rotamer_names) if n.upper() not in {"CA", "C", "N", "O"}]
    p_rmin, p_eps = _scaled_parameters([protein_names[i][0] for i in keep_p], 1.0, 1.0)
    r_rmin, r_eps = _scaled_parameters([rotamer_names[i][0] for i in keep_r], 1.0, 1.0)
    d = np.linalg.norm(rotamer_coords[0][keep_r][:, None, :] - protein_coords[keep_p][None, :, :], axis=2)
    expected = lj_energy(d, np.add.outer(r_rmin, p_rmin), np.sqrt(np.multiply.outer(r_eps, p_eps)), cutoff=10.0).sum()
    assert res.energies[0] == pytest.approx(expected, rel=1e-12, abs=1e-12)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
