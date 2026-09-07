"""Boltzmann sampling of a torsion by the Langevin (md) and Brownian (bd) integrators (PRD-108).

Medium: ~25 s (interwell exchange of a 3-fold torsion is slow, a few 1e6 steps
are needed for the histogram). The cheap pins (equipartition, MSD, the attached
dye) are in test_langevin_sampler.py.
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pytest

import IMP
import IMP.algebra
import IMP.atom
import IMP.core

from IMP.bff import KB_KCAL, make_langevin_simulator
from test_langevin_sampler import _particle, _run


def _torsion_system(m, k_tor=1.0, n=3):
    """Four particles: stiff bonds/angles, one cosine torsion V = k (1 - cos(n phi)) (IMP's Cosine)."""
    xyz = [(-1.5, 1.0, 0.0), (0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (3.0, 1.0, 0.0)]
    ps = [_particle(m, c) for c in xyz]
    rs = []
    for a, b in ((0, 1), (1, 2), (2, 3)):
        d0 = (IMP.core.XYZ(ps[a]).get_coordinates() - IMP.core.XYZ(ps[b]).get_coordinates()).get_magnitude()
        rs.append(IMP.core.DistanceRestraint(m, IMP.core.Harmonic(d0, 500.0), ps[a], ps[b]))
    for a, b, c in ((0, 1, 2), (1, 2, 3)):
        va = IMP.core.XYZ(ps[a]).get_coordinates() - IMP.core.XYZ(ps[b]).get_coordinates()
        vc = IMP.core.XYZ(ps[c]).get_coordinates() - IMP.core.XYZ(ps[b]).get_coordinates()
        th = math.acos(va * vc / (va.get_magnitude() * vc.get_magnitude()))
        rs.append(IMP.core.AngleRestraint(m, IMP.core.Harmonic(th, 100.0), ps[a], ps[b], ps[c]))
    rs.append(IMP.core.DihedralRestraint(m, IMP.core.Cosine(k_tor, n, 0.0), *ps))
    # keep the middle bond in place (fixed frame) so nothing drifts away
    IMP.core.XYZ(ps[1]).set_coordinates_are_optimized(False)
    IMP.core.XYZ(ps[2]).set_coordinates_are_optimized(False)
    return ps, IMP.core.RestraintsScoringFunction(rs)


def _phi(ps):
    return IMP.core.get_dihedral(*[IMP.core.XYZ(p) for p in ps])


def _boltzmann_torsion_hist(k_tor, n, T, edges):
    # IMP.core.Cosine(k, n, 0) scores k (1 - cos(n phi)): minima at 0, +-120 deg
    phi = 0.5 * (edges[1:] + edges[:-1])
    w = np.exp(-k_tor * (1.0 - np.cos(n * phi)) / (KB_KCAL * T))
    return w / w.sum()


_TORSION = {"T": 300.0, "k": 0.3, "n": 3, "edges": np.linspace(-math.pi, math.pi, 13)}
_HISTS = {}


def _torsion_hist(integrator):
    """Sampled torsion histogram (cached per process; md ~10 s, bd ~6 s)."""
    if integrator not in _HISTS:
        m = IMP.Model()
        ps, sf = _torsion_system(m, _TORSION["k"], _TORSION["n"])
        dt, n_steps, seed = (2.0, 2000000, 11) if integrator == "md" else (1.0, 1500000, 5)
        sim = make_langevin_simulator(m, [ps[0], ps[3]], sf, integrator=integrator, temperature=_TORSION["T"],
                             timestep_fs=dt, friction_ps=5.0, seed=seed)
        sim.optimize(5000)
        phis = _run(sim, n_steps, 20, lambda: _phi(ps))
        h, _ = np.histogram(phis, bins=_TORSION["edges"])
        _HISTS[integrator] = h / h.sum()
    return _HISTS[integrator]


@pytest.mark.parametrize("integrator", ["md", "bd"])
def test_torsion_samples_the_boltzmann_distribution(integrator):
    """A 3-fold cosine torsion: histogram follows exp(-V/kT) -- three equal wells,
    valley/peak ratio exp(-2k/kT); interwell exchange is slow, so the bounds are
    statistical (a few 1e6 steps)."""
    T, k, n, edges = _TORSION["T"], _TORSION["k"], _TORSION["n"], _TORSION["edges"]
    hist = _torsion_hist(integrator)
    ref = _boltzmann_torsion_hist(k, n, T, edges)
    assert np.abs(hist - ref).max() < 0.06, (integrator, np.round(hist, 3), np.round(ref, 3))
    # the Boltzmann shape: bin-averaged valley/peak ratio (exp(-2k/kT) point-wise)
    peak_bins, valley_bins = [1, 2, 5, 6, 9, 10], [0, 3, 4, 7, 8, 11]
    ratio = hist[valley_bins].mean() / hist[peak_bins].mean()
    ratio_ref = ref[valley_bins].mean() / ref[peak_bins].mean()
    assert ratio == pytest.approx(ratio_ref, rel=0.25), (ratio, ratio_ref)
    # three equal wells -- interwell exchange is slow (residence ~100s of ps),
    # so with a few ns the populations still scatter by ~0.1
    well_pop = np.array([hist[[10, 11, 0, 1]].sum(), hist[[2, 3, 4, 5]].sum(), hist[[6, 7, 8, 9]].sum()])
    assert np.abs(well_pop - 1.0 / 3.0).max() < 0.12, well_pop
    assert (well_pop > 0.2).all()


def test_md_and_bd_agree_on_the_torsion_distribution():
    """Both integrators sample the same Boltzmann distribution (shape; well
    populations carry the slow-exchange scatter of either run)."""
    peak_bins, valley_bins = [1, 2, 5, 6, 9, 10], [0, 3, 4, 7, 8, 11]
    r_md = _torsion_hist("md")[valley_bins].mean() / _torsion_hist("md")[peak_bins].mean()
    r_bd = _torsion_hist("bd")[valley_bins].mean() / _torsion_hist("bd")[peak_bins].mean()
    assert r_md == pytest.approx(r_bd, rel=0.2)
    assert np.abs(_torsion_hist("md") - _torsion_hist("bd")).max() < 0.08



if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
