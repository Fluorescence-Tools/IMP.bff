"""Pair kernels over weighted point sets with dipoles (fret.distance, PRD-108)."""

import numpy as np
import pytest

from IMP.bff import fret_pair_distribution, fret_pair_efficiencies, fret_pair_geometry


def _cloud(a):
    """A cloud's first three columns, contiguous -- an (N,4) AV cloud carries its
    weight in column 4, and the pair kernel reads coordinates only."""
    return np.ascontiguousarray(np.asarray(a, dtype=np.float64)[:, :3])


def test_isotropic_fallback_and_shapes():
    rng = np.random.default_rng(1)
    p1 = rng.normal(size=(5, 3)); p2 = rng.normal(size=(7, 3)) + 40.0
    w1 = rng.dirichlet(np.ones(5)); w2 = np.ones(7)          # unnormalised is fine
    g = fret_pair_geometry(p1, w1, p2, w2)
    assert g.get_R().reshape(5, 7).shape == (5, 7)
    assert g.get_kappa2().reshape(5, 7).shape == (5, 7)
    assert g.get_weight().reshape(5, 7).shape == (5, 7)
    assert g.get_weight().sum() == pytest.approx(1.0)
    assert np.all(g.get_kappa2() == 2.0 / 3.0) and g.kappa2_avg == pytest.approx(2.0 / 3.0)
    # (N,4) inputs (AV points) are accepted: the 4th column is ignored
    g4 = fret_pair_geometry(_cloud(np.hstack([p1, w1[:, None]])), w1,
                            _cloud(np.hstack([p2, np.ones((7, 1))])), w2)
    np.testing.assert_allclose(g4.get_R(), g.get_R())


def test_single_pair_analytic():
    p1 = np.array([[0.0, 0.0, 0.0]]); p2 = np.array([[0.0, 0.0, 45.0]])
    mu = np.array([[0.0, 0.0, 1.0]])                          # collinear: kappa2 = 4
    out = fret_pair_distribution(p1, [1.0], p2, [1.0], forster_radius=52.0, mu1=mu, mu2=mu, tau0=4.0)
    ratio = (52.0 / 45.0) ** 6 * 1.5 * 4.0
    e = ratio / (1 + ratio)
    assert out.get_kappa2().reshape(1, 1)[0, 0] == pytest.approx(4.0)
    assert out.static_efficiency == pytest.approx(e, rel=1e-12)
    assert out.dynamic1 == pytest.approx(e, rel=1e-12)
    assert out.dynamic2 == pytest.approx(e, rel=1e-12)
    assert out.get_k_fret().reshape(1, 1)[0, 0] == pytest.approx(ratio / 4.0, rel=1e-12)
    assert out.get_E().reshape(1, 1)[0, 0] == pytest.approx(e, rel=1e-12)


def test_regime_relations_and_weights():
    rng = np.random.default_rng(2)
    p1 = rng.normal(size=(6, 3)) * 5; p2 = rng.normal(size=(9, 3)) * 5 + 50.0
    mu1 = rng.normal(size=(6, 3)); mu2 = rng.normal(size=(9, 3))
    w1 = rng.dirichlet(np.ones(6)); w2 = rng.dirichlet(np.ones(9))
    g = fret_pair_geometry(p1, w1, p2, w2, mu1, mu2)
    k2 = g.get_kappa2().reshape(6, 9)
    assert k2.min() >= 0 and k2.max() <= 4 + 1e-12
    eff = fret_pair_efficiencies(g, 52.0)
    assert 0 <= eff.static_efficiency <= 1 and 0 <= eff.dynamic2 <= 1
    assert eff.static_efficiency <= eff.dynamic2 + 1e-12      # Jensen on x/(1+x)
    # E matrix averaged with the pair weights is the static value
    assert np.sum(eff.get_E().reshape(6, 9) * g.get_weight().reshape(6, 9)) \
           == pytest.approx(eff.static_efficiency, rel=1e-12)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))