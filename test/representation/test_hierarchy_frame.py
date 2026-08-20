"""Reading a frame out of an IMP hierarchy: one call, not a dozen per atom.

``_collect_frame`` built a decorator, read three coordinates, looked up an atom
type and walked two parents for **every leaf** — about a dozen SWIG crossings
per atom, tens of thousands per frame, to move a few kilobytes. It was the
largest single cost in loading a rotamer library: the FRETpredict pins spent
7.5 s of 13.75 s inside it.

Two C++ calls now, and the rotamer suite went from ~44 s to 6 s.

The tests compare against the Python loop, kept verbatim, because the metadata
is fiddly in a way that is easy to get subtly wrong — and was. IMP names an atom
``"Atom CB of residue 1"``: five fields, not two. A first cut took the *last*
one, so every atom in the structure came back named for its residue number, and
the failure surfaced as "missing backbone atoms CA, N, C".
"""

import numpy as np
import pytest

import IMP
import IMP.atom
import IMP.core
import IMP.bff
from IMP.bff import get_structure_dir


def _reference(hierarchy):
    """The Python loop this replaced, unchanged."""
    leaves = IMP.atom.get_leaves(hierarchy)
    particles = [p for p in leaves if IMP.core.XYZ.get_is_setup(p)]
    coords = np.zeros((len(particles), 3), dtype=np.float64)
    names, types, resnames, chains, residues = [], [], [], [], []
    for i, p in enumerate(particles):
        xyz = IMP.core.XYZ(p)
        coords[i] = [xyz.get_x(), xyz.get_y(), xyz.get_z()]
        if IMP.atom.Atom.get_is_setup(p):
            atom = IMP.atom.Atom(p)
            parts = atom.get_name().split()
            names.append(parts[1] if len(parts) > 1 else parts[0])
            types.append(atom.get_atom_type().get_string())
            res_p = p.get_parent()
            if IMP.atom.Residue.get_is_setup(res_p):
                resnames.append(str(IMP.atom.Residue(res_p).get_residue_type().get_string()))
                chain_p = res_p.get_parent()
                chains.append(str(IMP.atom.Chain(chain_p).get_id())
                              if IMP.atom.Chain.get_is_setup(chain_p) else "")
                residues.append(int(IMP.atom.Residue(res_p).get_index()))
            else:
                resnames.append("")
                chains.append("")
                residues.append(-1)
        else:
            name = p.get_name()
            names.append(name)
            types.append(name[0] if name else "C")
            resnames.append("")
            chains.append("")
            residues.append(-1)
    return coords, names, types, resnames, chains, residues


@pytest.fixture(scope="module")
def hierarchy():
    """The model is kept alive alongside the hierarchy, deliberately.

    An ``IMP.Model`` owns its particles. Returning only the hierarchy lets the
    model be collected, and every later access raises "Invalid particle
    requested" -- which is what the first version of this fixture did.
    """
    from pathlib import Path
    pdbs = sorted(Path(get_structure_dir()).rglob("*.pdb"))
    if not pdbs:
        pytest.skip("no bundled structure to read")
    model = IMP.Model()
    h = IMP.atom.read_pdb(str(pdbs[0]), model)
    h._keep_model_alive = model
    return h


def test_the_cpp_walk_matches_the_python_loop(hierarchy):
    want_xyz, want_names, want_types, want_res, want_chain, want_idx = _reference(hierarchy)

    packed = np.asarray(IMP.bff.hierarchy_atom_coordinates(hierarchy),
                        dtype=np.float64).reshape(-1, 4)
    meta = list(IMP.bff.hierarchy_atom_metadata(hierarchy))

    assert packed.shape[0] == len(want_names)
    np.testing.assert_array_equal(packed[:, :3], want_xyz)
    np.testing.assert_array_equal(packed[:, 3].astype(int), want_idx)
    assert meta[0::4] == want_names
    assert meta[1::4] == want_types
    assert meta[2::4] == want_res
    assert meta[3::4] == want_chain


def test_the_atom_name_is_the_second_field_not_the_last(hierarchy):
    """IMP's name is "Atom CB of residue 1". Taking the last field returns the
    residue number for every atom, which is how this broke the first time."""
    meta = list(IMP.bff.hierarchy_atom_metadata(hierarchy))
    names = meta[0::4]
    assert names, "the fixture produced no atoms"
    assert not any(n.isdigit() for n in names), "an atom is not named for its residue"
    assert "CA" in names and "N" in names


def test_only_xyz_leaves_are_returned(hierarchy):
    leaves = IMP.atom.get_leaves(hierarchy)
    with_xyz = [p for p in leaves if IMP.core.XYZ.get_is_setup(p)]
    packed = np.asarray(IMP.bff.hierarchy_atom_coordinates(hierarchy)).reshape(-1, 4)
    assert packed.shape[0] == len(with_xyz)


def test_the_two_calls_agree_on_length(hierarchy):
    """They are read separately, so nothing but the hierarchy keeps them in step."""
    n_coords = len(IMP.bff.hierarchy_atom_coordinates(hierarchy)) // 4
    n_meta = len(IMP.bff.hierarchy_atom_metadata(hierarchy)) // 4
    assert n_coords == n_meta


def test_an_empty_hierarchy_yields_nothing():
    model = IMP.Model()                      # held for the duration of the test
    empty = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model))
    assert list(IMP.bff.hierarchy_atom_coordinates(empty)) == []
    assert list(IMP.bff.hierarchy_atom_metadata(empty)) == []


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
