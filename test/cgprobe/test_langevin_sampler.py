"""Thermodynamic pins of the Langevin (md) and Brownian (bd) dye samplers (PRD-108 stage 3).

Real stochastic dynamics must sample the Boltzmann distribution of its
energy function: a harmonic oscillator at the right temperature, a free
Brownian particle with MSD = 6 D t, a torsion whose histogram is exp(-V/kT),
and md and bd agreeing on it. Plus the attached-dye sampler itself on hGBP1
481 (anchor fixed, dye moves, kinetic temperature right).
"""

import math

import numpy as np
import pytest

import IMP
import IMP.algebra
import IMP.atom
import IMP.core

from IMP.bff import AttachedProbeDynamics, kb_kcal, make_langevin_simulator

KB_KCAL = kb_kcal()      # kcal/mol/K, from the one place it is spelled


def _particle(m, xyz, radius=1.7, mass=12.0):
    p = IMP.Particle(m)
    IMP.core.XYZR.setup_particle(p, IMP.algebra.Sphere3D(IMP.algebra.Vector3D(*xyz), radius))
    IMP.atom.Mass.setup_particle(p, mass)
    IMP.core.XYZ(p).set_coordinates_are_optimized(True)
    return p


def _run(sim, n_steps, every, fn):
    out = []
    for _ in range(n_steps // every):
        sim.optimize(every)
        out.append(fn())
    return np.asarray(out)


def test_md_harmonic_oscillator_equipartition():
    """<V> = 3/2 kT for an isotropic 3D harmonic well; kinetic temperature = T."""
    m = IMP.Model()
    T, k = 300.0, 5.0
    p = _particle(m, (0.3, 0.0, 0.0))
    anchor = _particle(m, (0.0, 0.0, 0.0))
    IMP.core.XYZ(anchor).set_coordinates_are_optimized(False)
    r = IMP.core.DistanceRestraint(m, IMP.core.Harmonic(0.0, k), p, anchor)
    sf = IMP.core.RestraintsScoringFunction([r])
    sim = make_langevin_simulator(m, [p], sf, integrator="md", temperature=T, timestep_fs=2.0, friction_ps=20.0, seed=7)
    sim.optimize(2000)   # equilibrate
    # `make_langevin_simulator` is declared to return the base `Simulator`,
    # so Python hands back that interface; the kinetic energy is Molecular
    # Dynamics' own.
    md = IMP.atom.MolecularDynamics.get_from(sim)
    samples = _run(sim, 60000, 20, lambda: (sf.evaluate(False), md.get_kinetic_energy()))
    v_mean = samples[:, 0].mean()
    t_kin = 2.0 * samples[:, 1].mean() / (3.0 * KB_KCAL)
    assert v_mean == pytest.approx(1.5 * KB_KCAL * T, rel=0.10)
    assert t_kin == pytest.approx(T, rel=0.10)


def test_bd_free_particle_msd_is_6Dt():
    m = IMP.Model()
    T = 300.0
    p = _particle(m, (0.0, 0.0, 0.0), radius=2.0)
    # a zero-strength restraint: the scoring function must not be empty
    far = _particle(m, (1000.0, 0.0, 0.0))
    IMP.core.XYZ(far).set_coordinates_are_optimized(False)
    sf = IMP.core.RestraintsScoringFunction([IMP.core.DistanceRestraint(m, IMP.core.Harmonic(0.0, 0.0), p, far)])
    dt = 10.0
    sim = make_langevin_simulator(m, [p], sf, integrator="bd", temperature=T, timestep_fs=dt, seed=3)
    d_coef = IMP.atom.Diffusion(p).get_diffusion_coefficient()
    assert d_coef == pytest.approx(IMP.atom.get_einstein_diffusion_coefficient(2.0, T))
    n_rep, n_steps = 400, 50
    msd = []
    for _ in range(n_rep):
        IMP.core.XYZ(p).set_coordinates(IMP.algebra.Vector3D(0, 0, 0))
        sim.optimize(n_steps)
        msd.append(IMP.core.XYZ(p).get_coordinates().get_squared_magnitude())
    assert np.mean(msd) == pytest.approx(6.0 * d_coef * n_steps * dt, rel=0.15)


def test_attached_dye_sampler_hgbp1_481():
    from IMP.bff import ProbeAttachment, attach_probes
    from IMP.bff import get_structure_dir
    for integrator in ("md", "bd"):
        m = IMP.Model()
        prot = IMP.atom.read_pdb(str(get_structure_dir("1DG3.pdb")), m, IMP.atom.NonWaterPDBSelector())
        dye = IMP.atom.read_mol2(str(get_structure_dir("alexa488_r48.mol2")), m)
        attach_probes(prot, [ProbeAttachment(dye, "A", 481)], strip_site_sidechain=True)
        s = AttachedProbeDynamics(prot, dye, str(get_structure_dir("alexa488_r48.mol2")), "A", 481, integrator=integrator, seed=1)
        assert len(s.fixed) == 4 and len(s.mobile) == len(s.probe_particles) - 4
        assert len(s.obstacles) > 100
        s.minimize(200)
        c0 = s.coordinates
        traj = s.run(3000, write_every=100)
        assert traj.n_frames == 30 and np.isfinite(traj.coordinates).all()
        anchor = [i for i, n in enumerate(s.atom_names) if n.upper() in ("N", "CA", "C", "O")]
        np.testing.assert_allclose(traj.coordinates[-1][anchor], c0[anchor], atol=1e-9)   # anchor never moves
        moved = np.linalg.norm(traj.coordinates[-1] - c0, axis=1)
        assert moved.max() > 0.5                                                          # the dye does move
        if integrator == "md":
            t_kin = np.mean([s.kinetic_temperature(k) for k in traj.kinetic_energy[10:]])
            assert t_kin == pytest.approx(300.0, rel=0.15)
        # no mobile dye heavy atom ends up inside the protein (soft spheres hold);
        # the fixed anchor atoms are bonded to the neighbouring residues' backbone
        prot_xyz = np.array([IMP.core.XYZ(p).get_coordinates() for p in s.obstacles])
        heavy = [i for i, n in enumerate(s.atom_names) if not n.upper().startswith("H") and i not in anchor]
        d = np.linalg.norm(traj.coordinates[-1][heavy][:, None, :] - prot_xyz[None, :, :], axis=2)
        assert d.min() > 1.5


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
