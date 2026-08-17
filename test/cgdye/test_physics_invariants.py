"""Physics invariants of cgdye's FRET/rotamer machinery (PRD-107 stage 4).

Cheap analytic checks that pin the formulas rather than a reference run:
orientation factor limits, R0 scaling, the single-rotamer FRET case, the
Boltzmann / mean-field / kinetic weights, and the master-equation limits.
"""

import numpy as np
import pytest

from IMP.bff.cgdye.rotamer.scoring import kappa2_from_vectors
from IMP.bff.cgdye.rotamer.r0 import calculate_r0
from IMP.bff.cgdye.sampling.boltzmann import compute_boltzmann_weights
from IMP.bff.cgdye.sampling.kinetic import build_transition_probability_matrix
from IMP.bff.cgdye.sampling.mean_field import mean_field_weights
from IMP.bff.cgdye.analysis.fret import (
    calculate_fret_exact,
    calculate_fret_regimes,
    compute_exact_efficiency,
)

KB_KCAL = 0.0019872041


# --------------------------------------------------------------------------- κ²
def _unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def test_kappa2_limits():
    r = np.zeros((1, 1, 3))
    r[0, 0] = [0.0, 0.0, 10.0]
    z = np.array([[0.0, 0.0, 1.0]])
    x = np.array([[1.0, 0.0, 0.0]])
    y = np.array([[0.0, 1.0, 0.0]])
    assert kappa2_from_vectors(z, z, r)[0, 0] == pytest.approx(4.0)   # collinear head-to-tail
    assert kappa2_from_vectors(x, y, r)[0, 0] == pytest.approx(0.0)   # perpendicular, both ⟂ r
    assert kappa2_from_vectors(x, x, r)[0, 0] == pytest.approx(1.0)   # parallel, ⟂ r
    assert kappa2_from_vectors(x, z, r)[0, 0] == pytest.approx(0.0)   # one along r, one ⟂


def test_kappa2_isotropic_average_is_two_thirds():
    rng = np.random.default_rng(7)
    n = 400
    mu_d = _unit(rng.normal(size=(n, 3)))
    mu_a = _unit(rng.normal(size=(n, 3)))
    r = rng.normal(size=(n, n, 3))
    k2 = kappa2_from_vectors(mu_d, mu_a, r)
    assert k2.shape == (n, n)
    assert k2.min() >= 0.0 and k2.max() <= 4.0 + 1e-12
    assert k2.mean() == pytest.approx(2.0 / 3.0, abs=0.01)


# --------------------------------------------------------------------------- R0
def test_forster_radius_scales_with_kappa2_to_the_sixth():
    r0_iso = calculate_r0("AlexaFluor 488", "AlexaFluor 594", 2.0 / 3.0)
    r0_col = calculate_r0("AlexaFluor 488", "AlexaFluor 594", 4.0)
    assert r0_col / r0_iso == pytest.approx(6.0 ** (1.0 / 6.0), rel=1e-12)
    assert calculate_r0("AlexaFluor 488", "AlexaFluor 594", 0.0) == 0.0


@pytest.mark.parametrize("donor,acceptor,expected_nm", [
    ("AlexaFluor 488", "AlexaFluor 594", 5.687781310385014),
    ("AlexaFluor 594", "AlexaFluor 568", 4.737655494805884),
    ("AlexaFluor 488", "AlexaFluor 647", 5.378558552645125),
    ("AlexaFluor 546", "AlexaFluor 647", 6.858934378443619),
])
def test_forster_radius_pins(donor, acceptor, expected_nm):
    """R0 (nm) at κ² = 2/3 from the bundled spectra (recorded 2026-08-17)."""
    assert calculate_r0(donor, acceptor, 2.0 / 3.0) == pytest.approx(expected_nm, abs=1e-9)


def test_forster_radius_matches_fretpredict_hsp90():
    # FRETpredict's Hsp90 case: <κ²> = 0.9692561424158548 -> R0 5.0425682 nm
    assert calculate_r0("AlexaFluor 594", "AlexaFluor 568", 0.9692561424158548) == pytest.approx(5.0425682, abs=1e-6)


# --------------------------------------------------------------- single rotamer
def test_single_rotamer_pair_all_regimes_agree():
    d = np.array([[45.0]])
    k2 = np.array([[1.2]])
    r0 = 52.0
    out = calculate_fret_regimes(d, k2, np.array([1.0]), np.array([1.0]), R0=r0)
    ratio = (r0 / 45.0) ** 6 * 1.5 * 1.2
    expected = ratio / (1.0 + ratio)      # = 1 / (1 + (2/3/κ²) (r/R0)^6)
    assert out["static"] == pytest.approx(expected, rel=1e-12)
    assert out["dynamic"] == pytest.approx(expected, rel=1e-12)
    assert out["dynamic_plus"] == pytest.approx(expected, rel=1e-12)
    assert out["kappa2_avg"] == pytest.approx(1.2)
    # the same number through the master equation with no exchange (P = I)
    p = np.eye(1)
    e = compute_exact_efficiency(p, np.array([ratio / 4.0]), tau0=4.0, dt=1.0)
    assert e == pytest.approx(expected, rel=1e-12)


def test_regime_ordering_static_le_dynamic_plus():
    """Averaging rates (dynamic+) always gives E >= averaging efficiencies (static)
    -- Jensen on the concave x/(1+x)."""
    rng = np.random.default_rng(3)
    d = rng.uniform(30, 80, size=(6, 5))
    k2 = rng.uniform(0.0, 4.0, size=(6, 5))
    wd = rng.dirichlet(np.ones(6)); wa = rng.dirichlet(np.ones(5))
    out = calculate_fret_regimes(d, k2, wd, wa, R0=52.0)
    assert out["static"] <= out["dynamic_plus"] + 1e-12


# --------------------------------------------------------------- weights
def test_boltzmann_weights_are_softmax_of_minus_e_over_kt():
    e = np.array([0.0, 1.0, 2.5, 10.0])
    t = 298.15
    w = compute_boltzmann_weights(e, temperature=t)
    ref = np.exp(-e / (KB_KCAL * t)); ref /= ref.sum()
    np.testing.assert_allclose(w, ref, rtol=1e-12)
    assert w.sum() == pytest.approx(1.0)
    # invariant to an energy offset
    np.testing.assert_allclose(compute_boltzmann_weights(e + 123.4, t), w, rtol=1e-12)
    # huge energies do not overflow
    assert np.isfinite(compute_boltzmann_weights(np.array([0.0, 5000.0]), t)).all()


def test_mean_field_weights_prior_limits():
    rng = np.random.default_rng(1)
    rot = rng.normal(size=(4, 6, 3))            # 4 clusters, 6 atoms each
    prior = rng.dirichlet(np.ones(4))
    far_protein = rng.normal(size=(20, 3)) + 500.0
    near_protein = rot[0] + rng.normal(scale=0.3, size=rot[0].shape)  # sits on cluster 0
    elems_d = ["C"] * 6
    elems_p = ["C"] * near_protein.shape[0]
    # K = 0: the prior is returned
    np.testing.assert_allclose(mean_field_weights(rot, prior, near_protein, elems_d, elems_p, K=0.0), prior, rtol=1e-12)
    # protein far away: no interaction, prior returned for any K
    np.testing.assert_allclose(mean_field_weights(rot, prior, far_protein, elems_d, ["C"] * 20, K=5.0), prior, rtol=1e-12)
    # protein on top of cluster 0: cluster 0 loses weight, result stays a distribution
    q = mean_field_weights(rot, prior, near_protein, elems_d, elems_p, K=1.0)
    assert q.shape == (4,) and (q >= 0).all() and q.sum() == pytest.approx(1.0)
    assert q[0] < prior[0]


def test_transition_matrix_rows_and_symmetric_counts():
    counts = [[5, 2, 1], [2, 7, 3], [1, 3, 9]]  # symmetric -> detailed balance w.r.t. uniform? no: w ∝ row sums
    p = build_transition_probability_matrix(counts)
    np.testing.assert_allclose(p.sum(axis=1), 1.0)
    # stationary distribution of a reversible chain from symmetric counts is ∝ row sums
    evals, evecs = np.linalg.eig(p.T)
    stat = np.real(evecs[:, np.argmin(np.abs(evals - 1.0))]); stat /= stat.sum()
    row = np.array(counts, float).sum(axis=1); row /= row.sum()
    np.testing.assert_allclose(stat, row, atol=1e-12)
    # detailed balance: w_i P_ij == w_j P_ji
    np.testing.assert_allclose(stat[:, None] * p, (stat[:, None] * p).T, atol=1e-12)
    # a state with no outgoing counts becomes absorbing (jump to self)
    p2 = build_transition_probability_matrix([[0, 0], [1, 1]])
    assert p2[0, 0] == 1.0


# --------------------------------------------------------------- master equation
def test_master_equation_slow_and_fast_exchange_limits():
    tau0 = 4.0
    k_rad = 1.0 / tau0
    w = np.array([0.2, 0.5, 0.3])
    k = np.array([0.05, 0.4, 2.0])          # FRET rates 1/ns
    # slow exchange (P = I): E = Σ w_i k_i / (k_i + k_rad)  (static average)
    e_slow = compute_exact_efficiency(np.eye(3), k, tau0, dt=1.0, weights=w)
    assert e_slow == pytest.approx(np.sum(w * k / (k + k_rad)), rel=1e-12)
    # fast exchange (jump to the stationary distribution every tiny dt):
    # E -> <k> / (<k> + k_rad)  (dynamic average)
    p_fast = np.tile(w, (3, 1))
    e_fast = compute_exact_efficiency(p_fast, k, tau0, dt=1e-7, weights=w)
    k_avg = np.sum(w * k)
    assert e_fast == pytest.approx(k_avg / (k_avg + k_rad), abs=1e-6)
    # in between the two limits
    p_mid = 0.5 * np.eye(3) + 0.5 * p_fast
    e_mid = compute_exact_efficiency(p_mid, k, tau0, dt=1.0, weights=w)
    assert min(e_slow, e_fast) - 1e-12 <= e_mid <= max(e_slow, e_fast) + 1e-12
    # the stationary distribution is found when weights are not given
    assert compute_exact_efficiency(p_fast, k, tau0, dt=1.0) == pytest.approx(
        compute_exact_efficiency(p_fast, k, tau0, dt=1.0, weights=w), rel=1e-10)


def test_kronecker_pair_equals_explicit_product_space():
    p_d = build_transition_probability_matrix([[8, 2], [3, 7]])
    p_a = build_transition_probability_matrix([[5, 5], [1, 9]])
    wd = np.array([0.6, 0.4]); wa = np.array([0.3, 0.7])
    d = np.array([[40.0, 55.0], [62.0, 48.0]])
    k2 = np.array([[0.4, 1.1], [0.9, 2.5]])
    r0, tau0, dt = 52.0, 4.0, 0.1
    e_pair = calculate_fret_exact(d, k2, p_d, p_a, wd, wa, R0=r0, tau0=tau0, dt=dt)
    # explicit 4-state chain
    p_tot = np.kron(p_d, p_a)
    rates = ((r0 / d) ** 6 * 1.5 * k2 / tau0).ravel()
    e_ref = compute_exact_efficiency(p_tot, rates, tau0, dt, weights=np.outer(wd, wa).ravel())
    assert e_pair == pytest.approx(e_ref, rel=1e-12)
    # no exchange at all -> static average of the 2x2 grid
    e_static = calculate_fret_exact(d, k2, np.eye(2), np.eye(2), wd, wa, R0=r0, tau0=tau0, dt=dt)
    reg = calculate_fret_regimes(d, k2, wd, wa, R0=r0)
    assert e_static == pytest.approx(reg["static"], rel=1e-12)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
