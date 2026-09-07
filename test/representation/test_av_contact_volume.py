"""The accessible **contact** volume: what the two ACV fields do (PRD-121 G9).

`contact_volume_thickness` and `contact_volume_trapped_fraction` were accepted
by `IMP.bff.AV` and did nothing until 2026-09-01 -- read into the decorator,
stored, never used. `examples/structure/T4L/fret.fps.json` and the Zenodo
`FRET_screening.fps.json` both ask for them at every site, so files in the wild
believed they worked.

The rule is Olga's, from its source (`Olga/src/AV/fretAV.cpp`, `path2points`):
a cloud voxel is *in contact* when excluded volume lies within `thickness` of
it, and the contact voxels then carry `trapped_fraction` of the cloud's weight.
**FPS has no contact volume at all** -- this is an Olga/`IMP.bff` feature, and
`test_fps_screening_ab.py` turns it off before comparing against FPS for
exactly that reason.

The end-to-end effect is pinned in `test_fps_screening_ab.py`; what is pinned
here is the mechanism, and above all the *geometry*: that the voxels the kernel
calls "contact" really are the ones near the protein. That check is an
independent numpy re-derivation from the atom coordinates, not a second reading
of the same C++.
"""
import math

import numpy as np
import pytest

import IMP
import IMP.atom
import IMP.bff
import IMP.core


PDB = "structure/T4L/3GUN.pdb"


@pytest.fixture(scope="module")
def structure():
    IMP.set_log_level(IMP.SILENT)
    model = IMP.Model()
    hierarchy = IMP.atom.read_pdb(IMP.bff.get_example_path(PDB), model)
    atoms = np.array([list(IMP.core.XYZ(p).get_coordinates())
                      for p in IMP.atom.get_leaves(hierarchy)])
    return model, hierarchy, atoms


def _cloud(model, hierarchy, residue, thickness, trapped, step=1.0,
           space_fixed=True, atom="CB", width=4.5, clearance=-1.0):
    selection = IMP.atom.Selection(hierarchy)
    selection.set_atom_type(IMP.atom.AtomType(atom))
    selection.set_residue_index(residue)
    particle = IMP.Particle(model)
    IMP.bff.AV.do_setup_particle(
        model, particle, selection.get_selected_particles()[0],
        linker_length=20.0, radii=(3.5, 0.0, 0.0), linker_width=width,
        allowed_sphere_radius=clearance,
        simulation_grid_resolution=step,
        contact_volume_thickness=thickness,
        contact_volume_trapped_fraction=trapped)
    av = IMP.bff.AV(model, particle)
    av.set_space_fixed(space_fixed)
    av.resample()
    return np.asarray(av.get_map().get_xyz_density(), dtype=float)


def _mean(cloud):
    return (cloud[:, :3] * cloud[:, 3:4]).sum(0) / cloud[:, 3].sum()


def test_the_contact_voxels_are_the_ones_near_the_protein(structure):
    """The classification, checked against the atoms rather than against itself.

    With uniform base weights the reweighted cloud has exactly two weight
    levels, so the split the kernel made is readable from the outside. Every
    voxel of the heavier level must lie nearer the protein than the light ones:
    measured on 3GUN site 132 at a 1 A grid, 2103 contact voxels of 5029, their
    nearest-atom distance 6.81 A on average against 10.76 A for the free ones,
    the whole contact set within 8.61 A where the free set reaches 16.42 A.
    """
    model, hierarchy, atoms = structure
    cloud = _cloud(model, hierarchy, 132, thickness=3.0, trapped=0.6)
    levels = np.unique(cloud[:, 3])
    assert len(levels) == 2

    contact = cloud[cloud[:, 3] == levels[1]][:, :3]
    free = cloud[cloud[:, 3] == levels[0]][:, :3]
    assert len(contact) > 100 and len(free) > 100

    def nearest(points):
        return np.array([np.linalg.norm(atoms - p, axis=1).min()
                         for p in points])

    d_contact, d_free = nearest(contact), nearest(free)
    assert d_contact.mean() < d_free.mean() - 2.0
    # a contact voxel is *always* nearer the protein than the farthest free one
    assert d_contact.max() < d_free.max()


def test_the_trapped_share_is_the_share_that_was_asked_for(structure):
    """`trapped_fraction` is a share of the total weight, not a multiplier."""
    model, hierarchy, _ = structure
    for asked in (0.2, 0.5, 0.8):
        cloud = _cloud(model, hierarchy, 132, thickness=3.0, trapped=asked)
        weight = cloud[:, 3]
        levels = np.unique(weight)
        assert len(levels) == 2, asked
        heavy = levels[1] if asked > 0.5 else levels[0]
        # which level is the contact one flips with the fraction, so identify
        # it by the voxel *count*, which does not depend on the weighting
        counts = {level: int((weight == level).sum()) for level in levels}
        contact_level = min(counts, key=counts.get)
        got = weight[weight == contact_level].sum() / weight.sum()
        assert got == pytest.approx(asked, abs=1e-6)
        assert heavy in levels


def test_a_thicker_layer_traps_more_voxels(structure):
    """Monotone in the thickness, and quantised by the grid.

    Olga rounds the layer down to whole voxels before stamping it
    (`deltaIlist` takes an `int`), and that is reproduced, so at a 1 A grid
    thicknesses of 1.0, 2.0 and 3.0 A give three strictly growing contact sets
    while anything below 1.0 A gives none at all.
    """
    model, hierarchy, _ = structure

    def n_contact(thickness):
        cloud = _cloud(model, hierarchy, 132, thickness=thickness, trapped=0.6)
        weight = cloud[:, 3]
        levels = np.unique(weight)
        if len(levels) == 1:
            return 0
        return min(int((weight == level).sum()) for level in levels)

    sizes = [n_contact(t) for t in (1.0, 2.0, 3.0)]
    assert sizes[0] < sizes[1] < sizes[2]
    # below one grid step there is no layer: the volume comes back unweighted
    assert n_contact(0.9) == 0


def test_an_unweighted_volume_is_untouched(structure):
    """Neither field alone turns it on, and neither does an absent one."""
    model, hierarchy, _ = structure
    plain = _cloud(model, hierarchy, 132, thickness=0.0, trapped=-1.0)
    for thickness, trapped in ((3.0, -1.0), (0.0, 0.6)):
        other = _cloud(model, hierarchy, 132, thickness, trapped)
        assert np.array_equal(plain, other), (thickness, trapped)


def test_a_trapped_fraction_of_one_is_refused(structure):
    """1.0 would leave the free part weightless; Olga divides by zero here."""
    model, hierarchy, _ = structure
    plain = _cloud(model, hierarchy, 132, thickness=0.0, trapped=-1.0)
    assert np.array_equal(
        plain, _cloud(model, hierarchy, 132, thickness=3.0, trapped=1.0))


def test_the_legacy_anchoring_refuses_it(structure):
    """`space_fixed=False` ignores the ACV, and that is deliberate.

    The legacy anchoring exists for one thing: reproducing the pre-PRD-105
    numbers byte for byte (`references/prd105_legacy_pins.json`). Those were
    computed with the weighting inert, so honouring it there would change them
    and remove the only reason the path is still in the tree -- and it would
    need a second copy of the rule, written against the tiles rather than the
    SoA raster. The decorator warns once per handle instead of doing it
    silently.
    """
    model, hierarchy, _ = structure
    kwargs = dict(space_fixed=False, width=1.5, clearance=4.0)
    plain = _cloud(model, hierarchy, 132, 0.0, -1.0, **kwargs)
    contact = _cloud(model, hierarchy, 132, 3.0, 0.6, **kwargs)
    assert len(plain) > 0
    assert np.array_equal(plain, contact)


def test_the_array_door_keeps_the_contact_weights():
    """`compute_av_from_structure` forces uniform weights -- except this one.

    An ACV *is* the weighting, so flattening it there would have handed back a
    volume that silently ignored the request.
    """
    path = IMP.bff.get_example_path(PDB)
    plain = IMP.bff.compute_av_from_structure(
        path, "A", 132, "CB", 20.0, 4.5, 3.5, 0.0, 0.0, 1.0, "", -1.0,
        0.0, -1.0)
    contact = IMP.bff.compute_av_from_structure(
        path, "A", 132, "CB", 20.0, 4.5, 3.5, 0.0, 0.0, 1.0, "", -1.0,
        3.0, 0.6)
    p0 = np.asarray(plain.get_points(), dtype=float).reshape(-1, 4)
    p1 = np.asarray(contact.get_points(), dtype=float).reshape(-1, 4)
    assert len(p0) > 0 and len(p0) == len(p1)
    assert np.allclose(p0[:, 3], 1.0)
    assert len(np.unique(p1[:, 3])) == 2
    assert np.linalg.norm(_mean(p1) - _mean(p0)) > 0.2
