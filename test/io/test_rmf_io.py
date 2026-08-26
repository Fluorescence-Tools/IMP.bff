"""`RmfIO.h`: the four RMF readers and writers, in C++ since 2026-08-26.

They were `%pythoncode` builders in `IMP_bff.structureio.i`, registered in a
`_LAZY` table so `import IMP.bff` would not require `IMP.rmf`; `rmf` is one of
this module's modules now. What that changed for a caller, and what these
tests pin:

* the rotamer pair reads and writes an `IMP.bff.RotamerLibrary` -- the same
  value `read_drot` and `read_rotamer_library` return -- where it used to hand
  back a dict with 1-indexed `coords`;
* paths are `std::string`, so a `pathlib.Path` is `str()`-ed by the caller;
* `write_rmf` takes the description as JSON text rather than a dict.
"""

import json

import numpy as np
import pytest

import IMP.bff as fio


def _library(n_rotamers=3, n_atoms=4, with_transitions=True):
    lib = fio.RotamerLibrary()
    lib.n_rotamers = n_rotamers
    lib.n_atoms = n_atoms
    rng = np.random.default_rng(7)
    lib.set_coords(rng.normal(size=n_rotamers * n_atoms * 3).ravel().tolist())
    lib.set_weights([0.5, 0.3, 0.2][:n_rotamers])
    lib.atom_names = ["C1", "C2", "N3", "O4"][:n_atoms]
    if with_transitions:
        lib.transitions = [1, 2, 0, 0, 3, 1, 2, 0, 5][:n_rotamers ** 2]
    return lib


def test_rotamer_library_round_trips(tmp_path):
    lib = _library()
    out = str(tmp_path / "lib.rmf3")
    fio.write_rotamer_library_rmf(out, lib)

    back = fio.read_rotamer_library_rmf(out)
    assert back.n_rotamers == lib.n_rotamers
    assert back.n_atoms == lib.n_atoms
    assert list(back.atom_names) == list(lib.atom_names)
    np.testing.assert_allclose(back.weights, lib.weights, atol=1e-6)
    # RMF stores coordinates as float, so the round trip is single precision.
    np.testing.assert_allclose(back.coords, lib.coords, atol=1e-4)
    assert list(back.transitions) == list(lib.transitions)


def test_the_rmf3_suffix_is_supplied_at_both_ends(tmp_path):
    """The writer appends it and the reader looks for it, as the Python did."""
    stem = str(tmp_path / "lib")
    fio.write_rotamer_library_rmf(stem, _library())
    assert (tmp_path / "lib.rmf3").is_file()
    assert fio.read_rotamer_library_rmf(stem).n_rotamers == 3


def test_a_library_without_jump_counts_reads_back_without_them(tmp_path):
    """No transitions is not a defect: a container of static conformers has
    none, and inventing an identity matrix would hide that."""
    out = str(tmp_path / "static.rmf3")
    fio.write_rotamer_library_rmf(out, _library(with_transitions=False))
    assert len(fio.read_rotamer_library_rmf(out).transitions) == 0


def test_reading_a_library_that_is_not_there_raises(tmp_path):
    with pytest.raises(IOError):
        fio.read_rotamer_library_rmf(str(tmp_path / "absent"))


def test_write_rmf_writes_one_frame_of_spheres(tmp_path):
    import RMF
    out = str(tmp_path / "m.rmf3")
    coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 4.0, 0.0]])
    fio.write_rmf(coords, out, model_name="probe", radius=2.0,
                  metadata_json=json.dumps({"origin": "test"}))
    fh = RMF.open_rmf_file_read_only(out)
    assert fh.get_number_of_frames() == 1
    assert json.loads(fh.get_description())["origin"] == "test"


def test_write_rmf_wants_three_columns(tmp_path):
    with pytest.raises(ValueError):
        fio.write_rmf(np.zeros((2, 2)), str(tmp_path / "bad.rmf3"))


def test_protein_frames_round_trip_through_rmf(tmp_path):
    """The same `ProteinFrame` values the C++ PDB reader gives."""
    import IMP
    import IMP.atom
    import IMP.core
    import IMP.algebra
    import IMP.rmf
    import RMF

    model = IMP.Model()
    hier = IMP.atom.read_pdb(str(fio.get_structure_dir("148L.pdb")), model,
                             IMP.atom.NonWaterPDBSelector())
    path = str(tmp_path / "two.rmf3")
    fh = RMF.create_rmf_file(path)
    IMP.rmf.add_hierarchies(fh, [hier])
    IMP.rmf.save_frame(fh, "0")
    for a in IMP.atom.get_by_type(hier, IMP.atom.ATOM_TYPE):
        xyz = IMP.core.XYZ(a)
        xyz.set_coordinates(xyz.get_coordinates() +
                            IMP.algebra.Vector3D(1.0, 0.0, 0.0))
    IMP.rmf.save_frame(fh, "1")
    del fh

    frames = fio.protein_frames_from_rmf(path)
    assert len(frames) == 2
    assert frames[0].coords.shape == (frames[0].n_atoms, 3)
    assert list(frames[0].residue_indices)
    shift = frames[1].coords - frames[0].coords
    np.testing.assert_allclose(shift, np.tile([1.0, 0.0, 0.0],
                                              (shift.shape[0], 1)), atol=1e-3)

    # and `max_frames` stops early, where negative reads all
    assert len(fio.protein_frames_from_rmf(path, 1)) == 1
    assert len(fio.protein_frames_from_rmf(path, -1)) == 2


def test_reading_a_trajectory_that_is_not_there_raises(tmp_path):
    with pytest.raises(IOError):
        fio.protein_frames_from_rmf(str(tmp_path / "absent.rmf3"))
