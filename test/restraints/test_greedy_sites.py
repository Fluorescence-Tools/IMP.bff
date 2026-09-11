"""Greedy Olga over labelling *sites*: the unit a homo-oligomer buys.

In a homo-oligomer the labelling mix is statistical: each site is mutated once
and every cross-protomer combination of the chosen sites is measurable. A dimer
labelled at sites {1, 2} measures 1:1, 2:2, 1:2 and 2:1 -- four FRET species
for two mutations. So the greedy must add *sites* and score each step by the
pair set the enlarged site set implies, not add pairs one double mutant at a
time (``select_informative_pairs``, the monomer case).

Checked here against a numpy reference of the same algorithm, against the pair
selector (they must agree exactly when every site owns a single pair), and on
the combinatorics the oligomer adds.
"""

import numpy as np
import pytest

import IMP.bff

scipy_special = pytest.importorskip("scipy.special")


def _reference_select_sites(effs, rmsds, pair_sites, err, max_sites,
                            diag_weight=0.99):
    """The same algorithm in numpy, adding pairs in the kernel's order."""
    effs = np.asarray(effs, dtype=np.float64)
    n, m = effs.shape
    inv = 1.0 / err ** 2
    pair_sites = np.asarray(pair_sites)

    def gain(p):
        d = effs[:, p]
        return inv * np.subtract.outer(d, d) ** 2

    def expected(chi2, ndof):
        w = scipy_special.gammaincc(0.5 * ndof, 0.5 * chi2)
        w = np.where(chi2 <= 0.0, 1.0, w)
        sw = w.sum(axis=0)
        return ((w * rmsds).sum(axis=0) / (sw - 1.0 + diag_weight)).mean()

    candidates = sorted(set(pair_sites.ravel().tolist()))
    pairs_of = {s: [p for p in range(m) if s in pair_sites[p]]
                for s in candidates}

    def implied(i, chosen, active):
        """Pairs site `i` would activate: its own, and those to `chosen`."""
        out = []
        for p in pairs_of[i]:
            other = pair_sites[p][1] if pair_sites[p][0] == i else pair_sites[p][0]
            if not active[p] and (other == i or other in chosen):
                out.append(p)
        return out

    chi2 = np.zeros((n, n))
    chosen, active = set(), np.zeros(m, bool)
    sites, decay = [], []
    for _ in range(max_sites):
        best, best_score, best_new = None, np.inf, None
        for i in candidates:
            if i in chosen:
                continue
            new = implied(i, chosen, active)
            scratch = chi2.copy()
            for p in new:
                scratch += gain(p)
            score = expected(scratch, max(int(active.sum()) + len(new) - 2, 1))
            if score < best_score:
                best, best_score, best_new = i, score, new
        if best is None:
            break
        for p in best_new:
            chi2 += gain(p)
            active[p] = True
        chosen.add(best)
        sites.append(best)
        decay.append(expected(chi2, max(int(active.sum()), 1)))
    return np.asarray(sites), np.asarray(decay)


def _symmetric(rng, n, hi):
    m = rng.uniform(0, hi, (n, n))
    m = 0.5 * (m + m.T)
    np.fill_diagonal(m, 0.0)
    return m


def _dimer_pair_sites(sites):
    """Every inter-protomer combination: homotypic i:i and both i:j, j:i."""
    return np.array([[i, j] for i in sites for j in sites], dtype=np.int32)


def test_selection_matches_the_numpy_reference():
    """The kernel and the reference must pick the same sites in the same order.

    Both scan candidates in ascending site order and take the first minimum,
    and both accumulate the pair gains in ascending pair order -- so the
    agreement is exact, not statistical.
    """
    rng = np.random.default_rng(21)
    n, k = 50, 5
    effs = rng.random((n, k * k))
    rmsds = _symmetric(rng, n, 12.0)
    pair_sites = _dimer_pair_sites(range(k))
    max_sites = 4

    got_sites, got_decay = IMP.bff.select_informative_sites(
        effs, rmsds, pair_sites, err=0.06, max_sites=max_sites)
    ref_sites, ref_decay = _reference_select_sites(
        effs, rmsds, pair_sites, 0.06, max_sites)

    np.testing.assert_array_equal(np.asarray(got_sites), ref_sites)
    np.testing.assert_allclose(np.asarray(got_decay), ref_decay,
                               rtol=1e-10, atol=1e-12)


def test_one_pair_per_site_reduces_to_the_pair_selector():
    """A site that owns a single pair is a double mutant by another name.

    Every scoring decision then faces the same one pair the pair selector
    adds, with the same clamped ndof -- so the two selectors must return the
    same order and the same decay, exactly.
    """
    rng = np.random.default_rng(4)
    n, m = 40, 10
    effs = rng.random((n, m))
    rmsds = _symmetric(rng, n, 9.0)
    pair_sites = np.column_stack([np.arange(m), np.arange(m)]).astype(np.int32)

    pair_idx, pair_decay = IMP.bff.select_informative_pairs(
        effs, rmsds, err=0.05, max_pairs=6)
    site_idx, site_decay = IMP.bff.select_informative_sites(
        effs, rmsds, pair_sites, err=0.05, max_sites=6)

    np.testing.assert_array_equal(np.asarray(site_idx), np.asarray(pair_idx))
    np.testing.assert_allclose(np.asarray(site_decay), np.asarray(pair_decay),
                               rtol=1e-10, atol=1e-12)


def test_two_sites_imply_all_four_dimer_measurements():
    """The set is what gets evaluated: {0, 1} must leave the RMSD that all of
    0:0, 1:1, 0:1 and 1:0 leave together -- not the best of them, and not an
    average over double mutants.
    """
    rng = np.random.default_rng(9)
    n = 35
    effs = rng.random((n, 4))
    rmsds = _symmetric(rng, n, 11.0)
    pair_sites = _dimer_pair_sites([0, 1])

    sites, decay = IMP.bff.select_informative_sites(
        effs, rmsds, pair_sites, err=0.06, max_sites=2)
    assert sorted(np.asarray(sites).tolist()) == [0, 1]

    inv = 1.0 / 0.06 ** 2
    chi2 = np.zeros((n, n))
    for p in range(4):
        chi2 += inv * np.subtract.outer(effs[:, p], effs[:, p]) ** 2
    manual = IMP.bff.expected_rmsd(
        np.ascontiguousarray(rmsds).ravel(),
        np.ascontiguousarray(chi2).ravel(), 4, 0.99, n)
    np.testing.assert_allclose(np.asarray(decay)[-1], manual, rtol=1e-10)


def test_the_first_site_is_the_one_whose_own_pair_separates():
    """A site whose homotypic pair is informative must come first, ahead of a
    site that separates nothing on its own however well it pairs with others."""
    rng = np.random.default_rng(2)
    n = 30
    rmsds = _symmetric(rng, n, 10.0)
    # columns: 0:0 spans the ensemble, 1:1 is flat, 0:1 and 1:0 are flat-ish
    effs = np.column_stack([
        np.linspace(0.05, 0.95, n),
        np.full(n, 0.5),
        rng.uniform(0.49, 0.51, n),
        rng.uniform(0.49, 0.51, n),
    ])
    pair_sites = _dimer_pair_sites([0, 1])

    sites, decay = IMP.bff.select_informative_sites(
        effs, rmsds, pair_sites, err=0.06, max_sites=1)
    assert np.asarray(sites)[0] == 0

    start = IMP.bff.expected_rmsd(
        np.ascontiguousarray(rmsds).ravel(), np.zeros(n * n), 1, 0.99, n)
    assert np.asarray(decay)[0] < start


def test_selection_is_capped_unique_and_empty_safe():
    rng = np.random.default_rng(6)
    n, k = 20, 4
    effs = rng.random((n, k * k))
    rmsds = _symmetric(rng, n, 8.0)
    pair_sites = _dimer_pair_sites(range(k))

    sites, decay = IMP.bff.select_informative_sites(
        effs, rmsds, pair_sites, err=0.06, max_sites=k + 10)
    idx = np.asarray(sites)
    assert idx.size <= k
    assert len(set(idx.tolist())) == idx.size, "a site is never re-selected"

    empty_e = np.empty((n, 0))
    empty_p = np.empty((0, 2), dtype=np.int32)
    sites, decay = IMP.bff.select_informative_sites(
        empty_e, rmsds, empty_p, err=0.06, max_sites=3)
    assert np.asarray(sites).size == 0
    assert np.asarray(decay).size == 0


def test_bad_inputs_are_named():
    rng = np.random.default_rng(1)
    n, m = 12, 3
    effs = rng.random((n, m))
    rmsds = _symmetric(rng, n, 5.0)

    negative = np.array([[0, -1]] * m, dtype=np.int32)
    with pytest.raises(ValueError):
        IMP.bff.select_informative_sites(effs, rmsds, negative,
                                         err=0.06, max_sites=2)

    triplet = np.zeros((m, 3), dtype=np.int32)
    with pytest.raises(ValueError):
        IMP.bff.select_informative_sites(effs, rmsds, triplet,
                                         err=0.06, max_sites=2)

    short = np.zeros((n - 1, n - 1))
    with pytest.raises(ValueError):
        IMP.bff.select_informative_sites(effs, short,
                                         _dimer_pair_sites([0, 1, 2]),
                                         err=0.06, max_sites=2)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
