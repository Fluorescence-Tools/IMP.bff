import numpy as np
from IMP.bff.scoring import BoundingBoxFilter, DyeInternalEnergyEvaluator, lj_score

def test_bounding_box_filter():
    bbf = BoundingBoxFilter(pad=3.5)
    
    # Create rotamer coords: 2 frames, 1 atom each
    coords = np.array([
        [[0.0, 0.0, 0.0]], # Frame 0 at origin
        [[100.0, 100.0, 100.0]] # Frame 1 far away
    ])
    
    # Create reference coords: 1 atom at origin
    ref_coords = np.array([[0.0, 0.0, 0.0]])
    
    surviving, mask = bbf.filter_frames(coords, ref_coords)
    
    # Frame 0 should survive, Frame 1 should be filtered out
    assert mask.shape == (2,)
    assert mask[0] == True
    assert mask[1] == False
    assert surviving.shape == (1, 1, 3)

def test_internal_energy_evaluator_vectorized():
    system = {
        'sites': [
            {'id': 'a', 'atom_name': 'C1'},
            {'id': 'b', 'atom_name': 'C2'}
        ],
        'bonds': []
    }
    evaluator = DyeInternalEnergyEvaluator(system)
    
    # 2 frames, 2 atoms
    coords = np.array([
        [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]], # Frame 0: distance 3.0
        [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]]  # Frame 1: distance 5.0 (should be 0 energy if > rmin)
    ])
    
    energies = evaluator.evaluate_batch(coords)
    assert energies.shape == (2,)
    
    # Manually calculate expected energy for Frame 0
    # Assuming C-C interaction uses default params from CHARMM36_LJ
    expected_energy_0 = evaluator.evaluate(coords[0])
    expected_energy_1 = evaluator.evaluate(coords[1])
    
    np.testing.assert_allclose(energies[0], expected_energy_0)
    np.testing.assert_allclose(energies[1], expected_energy_1)

def test_evaluate_batch_filtered():
    system = {
        'sites': [
            {'id': 'a', 'atom_name': 'C1'},
            {'id': 'b', 'atom_name': 'C2'}
        ],
        'bonds': []
    }
    evaluator = DyeInternalEnergyEvaluator(system)
    
    # Frame 0 is close, Frame 1 is far
    coords = np.array([
        [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
        [[100.0, 0.0, 0.0], [103.0, 0.0, 0.0]]
    ])
    
    ref_coords = np.array([[0.0, 0.0, 0.0]])
    
    energies, mask = evaluator.evaluate_batch_filtered(coords, ref_coords)
    
    assert mask.shape == (2,)
    assert mask[0] == True
    assert mask[1] == False
    assert energies[1] == 0.0 # Energy should be 0 since it didn't overlap
    assert energies[0] != 0.0 # Energy should be evaluated


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
    """The improper branch of ``build_dye_restraints``, which had never run.

    No builder of a combined system fills ``impropers`` -- ``build_dye_protein_system``
    hard-codes an empty list -- so the loop that turns them into
    ``IMP.core.DihedralRestraint``s had no iterations, and the ``t["k"]`` inside
    it (``t`` is an ``FFTorsionType``, not a dict) could not raise. Injecting one
    improper is enough to hold the branch open.
    """
    import IMP
    import IMP.bff
    import IMP.atom
    from IMP.bff.tools import get_template_dir, get_structure_dir
    from IMP.bff.cgdye.topology import build_dye_protein_system
    from IMP.bff.scoring import build_dye_restraints

    system = build_dye_protein_system(
        str(get_structure_dir("cx4.mol2")), str(get_structure_dir("atto655.mol2")),
        "CX4", "atto655",
        protein_template=str(get_template_dir("cx4.template.cif")),
        dye_template=str(get_template_dir("atto655.template.cif")),
    )
    assert system.improper_types, "the templates declare improper types"

    model = IMP.Model()
    site_particles = {}
    for site in system.sites:
        p = IMP.Particle(model)
        IMP.core.XYZR.setup_particle(p)
        IMP.core.XYZ(p).set_coordinates(IMP.algebra.Vector3D(
            float(site.site_no), float(site.site_no % 7), float(site.site_no % 5)))
        site_particles[site.id] = p

    def n_dihedral(restraints):
        return sum(isinstance(r, IMP.core.DihedralRestraint) for r in restraints)

    before = build_dye_restraints(model, system, site_particles)
    without = len(before)

    first = system.dihedrals[0]
    improper = IMP.bff.FFTorsion()
    improper.site_a, improper.site_b = first.site_a, first.site_b
    improper.site_c, improper.site_d = first.site_c, first.site_d
    improper.type_id = sorted(system.improper_types)[0]
    system.impropers = [improper]

    with_one = build_dye_restraints(model, system, site_particles)
    # by count, not by position: the nonbonded pairs are appended after the
    # bonded terms, so the new restraint is not the last one in the list
    assert len(with_one) == without + 1
    assert n_dihedral(with_one) == n_dihedral(before) + 1
