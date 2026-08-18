"""Interaction terms reproduce the kernels they wrap, and their rates add.

These terms are the photophysics abstraction: a functional form with an arity
and parameters looked up by type, exactly as a force field has them. The point
of the tests below is that they are a *consolidation* rather than a seventh
parallel implementation -- each is checked against the kernel it wraps, so the
kernels can later move underneath the terms without changing an answer.
"""

import dataclasses

import numpy as np
import pytest

import IMP
import IMP.atom
import IMP.bff
import IMP.core
from IMP.bff.dye import find_dye
from IMP.bff.label import reference_pet_parameters
from IMP.bff.photophysics import (
    FRETTerm, PETTerm, RadiativeTerm, total_rate,
)
from IMP.bff.quenching import maps
from IMP.bff.representation import States


@pytest.fixture(scope="module")
def atoms():
    pdb = IMP.bff.get_example_path("structure/T4L/3GUN.pdb")
    model = IMP.Model()
    hier = IMP.atom.read_pdb(pdb, model, IMP.atom.NonWaterNonHydrogenPDBSelector())
    rows = []
    for leaf in IMP.atom.get_leaves(hier):
        atom = IMP.atom.Atom(leaf)
        residue = IMP.atom.get_residue(atom)
        rows.append((IMP.atom.get_chain(residue).get_id(), residue.get_index(),
                     residue.get_residue_type().get_string(),
                     atom.get_atom_type().get_string().strip(),
                     IMP.core.XYZ(leaf).get_coordinates(),
                     IMP.core.XYZR(leaf).get_radius()))
    return np.array(rows, dtype=[("chain", "U4"), ("res_id", "i8"),
                                 ("res_name", "U4"), ("atom_name", "U4"),
                                 ("coord", "f8", 3), ("radius", "f8")])


def _states(n=64, seed=0):
    rng = np.random.default_rng(seed)
    centre = np.array([30.0, 20.0, 10.0])
    pts = centre + rng.normal(0.0, 8.0, size=(n, 3))
    return States(points=np.column_stack([pts, np.ones(n)]),
                  attachment_point=centre)


class TestRadiative:

    def test_it_is_one_over_tau0_everywhere(self):
        states = _states(16)
        rates = RadiativeTerm(lifetime=4.0).rate_constants(states)
        assert rates.shape == (16,)
        np.testing.assert_allclose(rates, 0.25)

    def test_a_nonpositive_lifetime_is_refused(self):
        with pytest.raises(ValueError, match="lifetime"):
            RadiativeTerm(lifetime=0.0).rate_constants(_states(4))


class TestPET:

    def test_it_reproduces_the_grid_kernel(self, atoms):
        """The term and ``quenching_rate_map`` must be the same law.

        The map stamps ``k(r) = 1/tau0 + sum_a kQ_a exp(-(|r-r_a| - r_dye)/rC_a)``
        onto voxels; the term evaluates the sum at arbitrary states. Compared
        here at the voxel centres, minus the ``1/tau0`` floor the map adds and
        the term does not -- the radiative channel is its own term, which is the
        point of having terms at all.
        """
        params = reference_pet_parameters("AlexaFluor488", attenuation_length=1.5)
        term = PETTerm(parameters=params, dye_radius=3.5)

        ng, dg = 9, 2.0
        r0 = np.array([30.0, 20.0, 10.0])
        density = np.ones((ng, ng, ng))
        axis = maps.grid_axis(ng, dg)
        table = {c: {a: (p.rate_constant, p.attenuation_length)
                     for a in IMP.bff.QUENCHER_ATOMS[c]}
                 for c, p in params.items()}
        kQ, rC = maps.atomic_quenching_parameters(atoms, table)
        tau0 = 4.0
        grid = maps.quenching_rate_map(
            density, r0, dg, np.ascontiguousarray(atoms["coord"]), kQ, rC,
            tau0=tau0, dye_radius=3.5)

        centre = (ng - 1) // 2
        pts, expect = [], []
        for i in range(ng):
            for j in range(0, ng, 3):
                for k in range(0, ng, 3):
                    pts.append(r0 + np.array([axis[i], axis[j], axis[k]]))
                    expect.append(grid[i, j, k] - 1.0 / tau0)
        states = States(points=np.column_stack([np.array(pts), np.ones(len(pts))]),
                        attachment_point=r0)
        got = term.rate_constants(states, atoms)
        np.testing.assert_allclose(got, np.array(expect), rtol=1e-9, atol=1e-12)

    def test_it_is_a_pair_property(self, atoms):
        """Scaling the pair parameters scales the rate."""
        base = PETTerm(reference_pet_parameters("AlexaFluor488", attenuation_length=1.5))
        half = PETTerm({c: p.scaled(0.5) for c, p in
                        reference_pet_parameters("AlexaFluor488", attenuation_length=1.5).items()})
        s = _states(24)
        np.testing.assert_allclose(half.rate_constants(s, atoms),
                                   0.5 * base.rate_constants(s, atoms), rtol=1e-12)


class TestFRET:

    def test_r0_is_derived_from_the_pair_and_the_medium(self):
        d, a = find_dye("AlexaFluor 488"), find_dye("AlexaFluor 594")
        term = FRETTerm(donor=d, acceptor=a, refractive_index=1.4)
        water = FRETTerm(donor=d, acceptor=a, refractive_index=1.33)
        assert 40.0 < term.forster_radius < 70.0            # Angstrom
        assert water.forster_radius > term.forster_radius   # lower n, larger R0
        assert term.used_isotropic_kappa2

    def test_it_reproduces_the_trace_kernel(self):
        from IMP.bff.quenching.fret_trace import fret_rate_trace
        # `replace`, not `object.__setattr__`: Dye is frozen, and the library
        # is cached, so writing through the freeze edits the species for every
        # later reader. That is how this test used to leak a lifetime into
        # `test_a_dye_without_a_lifetime_is_refused`.
        d = dataclasses.replace(find_dye("AlexaFluor 488"), lifetime=4.0)
        a = find_dye("AlexaFluor 594")
        term = FRETTerm(donor=d, acceptor=a)
        donor, acceptor = _states(32, 1), _states(48, 2)
        got = term.rate_constants(donor, acceptor)
        want = fret_rate_trace(
            np.ascontiguousarray(donor.positions),
            np.ascontiguousarray(acceptor.positions),
            R0=term.forster_radius, tau0=4.0, r_min=7.0, kappa2=None)
        np.testing.assert_allclose(got, want, rtol=1e-12)

    def test_a_dye_without_a_lifetime_is_refused(self):
        d, a = find_dye("AlexaFluor 488"), find_dye("AlexaFluor 594")
        with pytest.raises(ValueError, match="lifetime"):
            FRETTerm(donor=d, acceptor=a).rate_constants(_states(4), _states(4))


class TestAdditivity:

    def test_rates_add(self, atoms):
        """Parallel channels add -- the property the field solver relies on."""
        states = _states(40)
        radiative = RadiativeTerm(lifetime=4.0)
        pet = PETTerm(reference_pet_parameters("AlexaFluor488", attenuation_length=1.5))
        total = total_rate([radiative, pet], states, atoms)
        np.testing.assert_allclose(
            total, radiative.rate_constants(states) + pet.rate_constants(states, atoms),
            rtol=1e-12)

    def test_the_sum_is_at_least_the_radiative_floor(self, atoms):
        states = _states(40)
        total = total_rate(
            [RadiativeTerm(lifetime=4.0),
             PETTerm(reference_pet_parameters("AlexaFluor488", attenuation_length=1.5))],
            states, atoms)
        assert np.all(total >= 0.25 - 1e-12)

    def test_no_terms_is_an_error_not_a_zero(self):
        with pytest.raises(ValueError, match="no interaction terms"):
            total_rate([], _states(4))


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
