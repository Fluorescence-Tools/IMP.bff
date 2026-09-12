"""Frame conversion preserves element masses on IMP's atom decorators."""

import pytest

import IMP
import IMP.bff


@pytest.mark.skipif(not hasattr(IMP.bff, "hierarchy_from_protein_frame"),
                    reason="requires IMP hierarchy bindings")
def test_frame_atoms_have_the_requested_element_mass():
    import IMP.atom
    import IMP.core

    frame = IMP.bff.ProteinFrame()
    frame.set_coords([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    frame.atom_names = ["N", "C"]
    # Explicit types must win over the atom-name defaults set up by IMP.
    frame.atom_types = ["N", "O"]
    frame.resnames = ["ALA", "ALA"]
    frame.chain_ids = ["A", "A"]
    frame.residue_indices = [1, 1]
    model = IMP.Model()
    hierarchy = IMP.bff.hierarchy_from_protein_frame(frame, model)
    atoms = IMP.atom.get_by_type(hierarchy, IMP.atom.ATOM_TYPE)
    assert len(atoms) == 2
    elements = IMP.atom.get_element_table()
    for atom, element, xyz in zip(atoms, [IMP.atom.N, IMP.atom.O],
                                  [(1., 2., 3.), (4., 5., 6.)]):
        assert IMP.atom.Mass(atom).get_mass() == elements.get_mass(element)
        assert tuple(IMP.core.XYZ(atom).get_coordinates()) == xyz
