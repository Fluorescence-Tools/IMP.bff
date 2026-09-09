"""`RmfStructureWriter`: a structure as an RMF trajectory.

The writer builds a PMI-shaped tree -- root, chain, residue, atom -- through
RMF's own decorators, because the connection layer's IMP does not carry
``IMP.rmf`` and a writer that needed it could not run in a package that links
IMP privately. What these tests pin is therefore the *file*, read back with
RMF: the tree it produced, the coordinates per frame, and the per-frame
scalars in the ``stat`` category where PMI puts them.

No ``IMP.atom`` here on purpose -- the writer is core, and this file must run
in the IMP-free lane too.
"""

import json

import numpy as np
import pytest

import IMP.bff as bff

RMF = pytest.importorskip("RMF")

pytestmark = pytest.mark.skipif(
    not hasattr(bff, "RmfStructureWriter"),
    reason="this IMP.bff was built without RMF",
)


def _table(n_res=3, atoms_per_res=2, chains=("A", "B")):
    """A small structure: two chains of `n_res` residues, `atoms_per_res` each."""
    t = bff.StructureTable()
    xyz, radius, mass, res_id, chain, res_name, atom_name = [], [], [], [], [], [], []
    k = 0
    for c in chains:
        for r in range(1, n_res + 1):
            for a in range(atoms_per_res):
                xyz += [float(k), float(k) + 0.5, float(k) - 0.5]
                radius.append(1.7)
                mass.append(12.0)
                res_id.append(r)
                chain.append(c)
                res_name.append("ALA")
                atom_name.append("CA" if a == 0 else "CB")
                k += 1
    t.xyz = np.asarray(xyz)
    t.radius = np.asarray(radius)
    t.mass = np.asarray(mass)
    t.res_id = np.asarray(res_id, dtype=np.int32)
    t.chain = chain
    t.res_name = res_name
    t.atom_name = atom_name
    return t


def test_a_table_can_be_built_and_read_back():
    """The record is writable from Python, not only fillable by the reader."""
    t = _table()
    assert t.n_atoms == 12
    assert t.xyz.shape == (12, 3)
    assert list(t.chain[:2]) == ["A", "A"]
    assert t.radius[0] == pytest.approx(1.7)


def test_the_tree_is_root_chain_residue_atom(tmp_path):
    path = tmp_path / "traj.rmf3"
    t = _table()
    w = bff.RmfStructureWriter(str(path), t, root_name="test")
    w.append(t.xyz.ravel())
    w.close()

    fh = RMF.open_rmf_file_read_only(str(path))
    root = fh.get_root_node()
    tops = root.get_children()
    assert len(tops) == 1 and tops[0].get_name() == "test"
    chains = tops[0].get_children()
    assert [c.get_name() for c in chains] == ["A", "B"]
    residues = chains[0].get_children()
    assert len(residues) == 3
    assert [a.get_name() for a in residues[0].get_children()] == ["CA", "CB"]


def test_every_frame_keeps_its_own_coordinates(tmp_path):
    path = tmp_path / "traj.rmf3"
    t = _table()
    frames = [t.xyz.ravel() + offset for offset in (0.0, 10.0, 100.0)]
    w = bff.RmfStructureWriter(str(path), t)
    for f in frames:
        w.append(f)
    assert w.n_frames == 3
    w.close()

    fh = RMF.open_rmf_file_read_only(str(path))
    pf = RMF.ParticleConstFactory(fh)
    assert fh.get_number_of_frames() == 3

    def leaves(node, out):
        kids = node.get_children()
        if not kids:
            out.append(node)
        for k in kids:
            leaves(k, out)
        return out

    for i in range(3):
        fh.set_current_frame(RMF.FrameID(i))
        got = []
        for n in leaves(fh.get_root_node(), []):
            if pf.get_is(n):
                got += list(pf.get(n).get_coordinates())
        assert np.allclose(got, frames[i], atol=1e-3), i


def test_per_frame_scalars_land_in_the_stat_category(tmp_path):
    """Where PMI's own output puts them, so a PMI stat plot reads this file."""
    path = tmp_path / "traj.rmf3"
    t = _table()
    w = bff.RmfStructureWriter(str(path), t)
    w.append(t.xyz.ravel(), metadata_json=json.dumps({"score": 1.5, "step": 7}))
    w.append(t.xyz.ravel(), metadata_json=json.dumps({"score": 0.5, "step": 8}))
    w.close()

    fh = RMF.open_rmf_file_read_only(str(path))
    cat = fh.get_category("stat")
    score = fh.get_key(cat, "score", RMF.FloatTag())
    step = fh.get_key(cat, "step", RMF.IntTag())
    root = fh.get_root_node()
    fh.set_current_frame(RMF.FrameID(0))
    assert root.get_value(score) == pytest.approx(1.5)
    assert root.get_value(step) == 7
    fh.set_current_frame(RMF.FrameID(1))
    assert root.get_value(score) == pytest.approx(0.5)
    assert root.get_value(step) == 8


def test_a_key_keeps_the_type_it_was_created_with(tmp_path):
    """RMF types a key once. A name that arrives as an int and later as a
    float stays an int key rather than becoming two keys of one name."""
    path = tmp_path / "traj.rmf3"
    t = _table()
    w = bff.RmfStructureWriter(str(path), t)
    w.append(t.xyz.ravel(), metadata_json=json.dumps({"n": 3}))
    w.append(t.xyz.ravel(), metadata_json=json.dumps({"n": 4.75}))
    w.close()

    fh = RMF.open_rmf_file_read_only(str(path))
    cat = fh.get_category("stat")
    n = fh.get_key(cat, "n", RMF.IntTag())
    fh.set_current_frame(RMF.FrameID(1))
    assert fh.get_root_node().get_value(n) == 4


def test_a_frame_of_the_wrong_length_is_refused(tmp_path):
    t = _table()
    w = bff.RmfStructureWriter(str(tmp_path / "traj.rmf3"), t)
    with pytest.raises(Exception):
        w.append(np.zeros(7))
    w.close()


def test_an_empty_structure_is_refused(tmp_path):
    with pytest.raises(Exception):
        bff.RmfStructureWriter(str(tmp_path / "traj.rmf3"), bff.StructureTable())


def test_columns_that_disagree_are_refused(tmp_path):
    """A short column would silently misname atoms from that point on."""
    t = _table()
    t.res_name = ["ALA"]
    with pytest.raises(Exception):
        bff.RmfStructureWriter(str(tmp_path / "traj.rmf3"), t)


def test_the_extension_is_added(tmp_path):
    t = _table()
    w = bff.RmfStructureWriter(str(tmp_path / "traj"), t)
    w.append(t.xyz.ravel())
    w.close()
    assert (tmp_path / "traj.rmf3").is_file()


def test_the_writer_is_a_context_manager(tmp_path):
    path = tmp_path / "traj.rmf3"
    t = _table()
    with bff.RmfStructureWriter(str(path), t) as w:
        w.append(t.xyz.ravel())
        assert w.n_atoms == 12
    assert path.is_file()
