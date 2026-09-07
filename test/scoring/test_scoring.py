import numpy as np
from IMP.bff import (
    BoundingBoxFilter, IntramolecularEnergy, forcefield_system_from_json,
    lj_score)

def test_bounding_box_filter():
    bbf = BoundingBoxFilter(pad=3.5)
    
    # Create rotamer coords: 2 frames, 1 atom each
    coords = np.array([
        [[0.0, 0.0, 0.0]], # Frame 0 at origin
        [[100.0, 100.0, 100.0]] # Frame 1 far away
    ])
    
    # Create reference coords: 1 atom at origin
    ref_coords = np.array([[0.0, 0.0, 0.0]])
    
    out = bbf.filter_frames(
        coords.ravel(), 2, 1, ref_coords.ravel(), 1)

    # Frame 0 should survive, Frame 1 should be filtered out
    assert len(out.mask) == 2
    assert out.mask[0] == True
    assert out.mask[1] == False
    # flat: one surviving frame of one atom
    assert np.asarray(out.coords).shape == (3,)

def test_internal_energy_evaluator_vectorized():
    import json
    system = forcefield_system_from_json(json.dumps({
        'sites': [
            {'id': 'a', 'atom_name': 'C1'},
            {'id': 'b', 'atom_name': 'C2'}
        ],
        'bonds': []
    }))
    evaluator = IntramolecularEnergy(system)
    
    # 2 frames, 2 atoms
    coords = np.array([
        [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]], # Frame 0: distance 3.0
        [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]]  # Frame 1: distance 5.0 (should be 0 energy if > rmin)
    ])
    
    energies = np.asarray(evaluator.evaluate_batch(
        coords.ravel(), 2, 2))
    assert energies.shape == (2,)
    
    # Manually calculate expected energy for Frame 0
    # Assuming C-C interaction uses default params from CHARMM36_LJ
    expected_energy_0 = evaluator.evaluate(coords[0].ravel(), 2)
    expected_energy_1 = evaluator.evaluate(coords[1].ravel(), 2)
    
    np.testing.assert_allclose(energies[0], expected_energy_0)
    np.testing.assert_allclose(energies[1], expected_energy_1)

def test_evaluate_batch_filtered():
    import json
    system = forcefield_system_from_json(json.dumps({
        'sites': [
            {'id': 'a', 'atom_name': 'C1'},
            {'id': 'b', 'atom_name': 'C2'}
        ],
        'bonds': []
    }))
    evaluator = IntramolecularEnergy(system)
    
    # Frame 0 is close, Frame 1 is far
    coords = np.array([
        [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
        [[100.0, 0.0, 0.0], [103.0, 0.0, 0.0]]
    ])
    
    ref_coords = np.array([[0.0, 0.0, 0.0]])
    
    out = evaluator.evaluate_batch_filtered(
        coords.ravel(), 2, 2, ref_coords.ravel(), 1)

    assert len(out.mask) == 2
    assert out.mask[0] == True
    assert out.mask[1] == False
    assert out.energies[1] == 0.0  # no overlap -> not scored
    assert out.energies[0] != 0.0


# IMP runs every .py under test/ as a standalone script, and a file of bare
# pytest functions would import cleanly and exit 0 -- reporting success without
# running a single assertion. Hand the file to pytest explicitly so a failure
# here is a failure in ctest.
if __name__ == "__main__":
    import sys
    try:
        import pytest
    except ImportError:
        print("pytest not installed; skipping", __file__)
        sys.exit(0)
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))


def test_improper_restraints_are_built_from_a_typed_system():
    """The improper branch of ``build_probe_restraints``, which had never run.

    No builder of a combined system fills ``impropers`` -- ``build_probe_protein_system``
    hard-codes an empty list -- so the loop that turns them into
    ``IMP.core.DihedralRestraint``s had no iterations, and the ``t["k"]`` inside
    it (``t`` is an ``FFTorsionType``, not a dict) could not raise. Injecting one
    improper is enough to hold the branch open.
    """
    import IMP
    import IMP.bff
    import IMP.atom
    from IMP.bff import get_template_dir, get_structure_dir
    from IMP.bff import build_probe_protein_system
    from IMP.bff import build_probe_restraints

    system = build_probe_protein_system(
        str(get_structure_dir("cx4.mol2")), str(get_structure_dir("atto655.mol2")),
        "CX4", "atto655",
        protein_template=str(get_template_dir("cx4.template.cif")),
        probe_template=str(get_template_dir("atto655.template.cif")),
    )
    assert system.improper_types, "the templates declare improper types"
    assert system.impropers, "and the builder produces them"

    model = IMP.Model()
    site_particles = {}
    for site in system.sites:
        p = IMP.Particle(model)
        IMP.core.XYZR.setup_particle(p)
        IMP.core.XYZ(p).set_coordinates(IMP.algebra.Vector3D(
            float(site.site_no), float(site.site_no % 7), float(site.site_no % 5)))
        site_particles[site.id] = p

    def n_dihedral(restraints):
        # By name, not isinstance: the builder is C++ now and hands back an
        # `IMP::Restraints`, whose elements SWIG presents as the base
        # `IMP.Restraint` -- IMP's own `RestraintSet.get_restraints()` does
        # exactly the same. The Python builder could be checked with
        # isinstance only because it had constructed the objects here.
        return sum(r.get_name().startswith("DihedralRestraint")
                   for r in restraints)

    before = build_probe_restraints(model, system, list(site_particles),
                                 [p.get_index() for p in site_particles.values()])
    without = len(before)

    # the builder produces 81; dropping them is what the count must show
    kept = list(system.impropers)
    system.impropers = []
    stripped = build_probe_restraints(model, system, list(site_particles),
                                 [p.get_index() for p in site_particles.values()])
    system.impropers = kept

    assert n_dihedral(before) == n_dihedral(stripped) + len(kept)
    assert len(before) == len(stripped) + len(kept)


def _dye_system_and_particles(spacing=2.0, with_radii=True):
    """The bundled dye+protein system, decorated on a line.

    A line is not a structure, which is the point: every non-excluded pair
    overlaps, so a repulsion term that is present has something to say and one
    that is missing is obvious.
    """
    import IMP
    import IMP.algebra
    import IMP.core
    from IMP.bff import build_probe_protein_system, get_structure_dir, get_template_dir

    system = build_probe_protein_system(
        str(get_structure_dir("cx4.mol2")), str(get_structure_dir("atto655.mol2")),
        "CX4", "atto655",
        protein_template=str(get_template_dir("cx4.template.cif")),
        probe_template=str(get_template_dir("atto655.template.cif")))
    model = IMP.Model()
    site_particles = {}
    for i, site in enumerate(system.sites):
        p = IMP.Particle(model)
        IMP.core.XYZ.setup_particle(
            p, IMP.algebra.Vector3D(float(i) * spacing, 0.0, 0.0))
        if with_radii:
            IMP.core.XYZR(p).set_radius(float(site.radius))
        site_particles[site.id] = p
    return system, model, site_particles


def test_the_repulsion_is_one_term_shared_by_dynamics_and_monte_carlo():
    """`build_probe_restraints` ends with the same restraint MC is scored on.

    One clash term, not two. Per-pair Lennard-Jones lower bounds *and* soft
    spheres over a container would be two implementations of one piece of
    physics, and thousands of restraints where one does. Asked
    for `nonbonded=False` the builder leaves it out entirely, which is what a
    run wants when it scores rigid moves against the repulsion *alone*: a
    rigid move cannot change a bond, an angle or a torsion.
    """
    from IMP.bff import build_probe_restraints, build_steric_restraint

    system, model, site_particles = _dye_system_and_particles()
    ids = list(site_particles)
    idx = [site_particles[s].get_index() for s in ids]

    with_repulsion = list(build_probe_restraints(model, system, ids, idx))
    bonded_only = list(build_probe_restraints(model, system, ids, idx,
                                            nonbonded=False))

    assert len(with_repulsion) == len(bonded_only) + 1     # one, not thousands
    assert with_repulsion[-1].get_name() == "steric"
    assert not any(r.get_name() == "steric" for r in bonded_only)

    steric = build_steric_restraint(model, system, ids, idx)
    assert steric is not None
    assert steric.unprotected_evaluate(None) > 0.0


def test_the_repulsion_rises_when_two_atoms_are_pushed_together():
    import IMP.core
    from IMP.bff import build_steric_restraint

    system, model, site_particles = _dye_system_and_particles(spacing=8.0)
    ids = list(site_particles)
    idx = [site_particles[s].get_index() for s in ids]
    steric = build_steric_restraint(model, system, ids, idx)

    apart = steric.unprotected_evaluate(None)
    far = IMP.core.XYZ(site_particles[ids[-1]]).get_coordinates()
    IMP.core.XYZ(site_particles[ids[0]]).set_coordinates(far)
    together = steric.unprotected_evaluate(None)
    assert together > apart


def test_a_site_without_a_radius_still_repels():
    """A soft sphere is a sphere.

    A site the caller decorated without `XYZR` would contribute nothing at
    all, and silently -- the run would look fine and the molecule would pass
    through itself. The system says what radius each site has, so the builder
    gives it that one.
    """
    import IMP.core
    from IMP.bff import build_steric_restraint

    system, model, site_particles = _dye_system_and_particles(with_radii=False)
    ids = list(site_particles)
    idx = [site_particles[s].get_index() for s in ids]
    assert not IMP.core.XYZR.get_is_setup(site_particles[ids[0]])

    steric = build_steric_restraint(model, system, ids, idx)
    assert IMP.core.XYZR.get_is_setup(site_particles[ids[0]])
    assert steric.unprotected_evaluate(None) > 0.0
