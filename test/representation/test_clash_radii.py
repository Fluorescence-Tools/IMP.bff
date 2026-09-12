"""The excluded-volume term's own radii: a source and a scale (PRD-121 G7).

`k_clash = 2/ClashTolerance^2` is FPS's constant and FPS applies it to **Bondi**
radii (`vdW.txt`: C 1.70). `IMP::atom::read_pdb` assigns **united-atom** radii
from the CHARMM type, which carry implicit hydrogens -- carbon is 1.85-2.275 A
on this fixture. Applying FPS's constant to IMP's radii over-penalises, and the
size of the over-penalty is measured here rather than argued: over HIV-RT's
protein-DNA interface at k = 8, IMP's radii give **268 overlapping atom pairs /
90.49 A of total overlap / energy 198.40** and Olga's Bondi-scale table gives
**30 / 7.59 / 11.11** -- a factor of 17.9.

So `clash_container` has a radii source and a radii scale, spelled and validated
exactly as `ProbeAccessibleVolumeDecorator::set_radii_source` is. Two things about them are decisions:

**Nothing is written on the structure's own radii.** A substituted radii set is
carried by *shadow spheres* -- new particles at the atoms' coordinates, added to
the same rigid bodies, holding nothing but the radius. The atoms keep theirs,
because `AV` under the default `"imp"` inflates its obstacles by exactly those
radii and `dock_minimize` resamples volumes between minimisations: rewriting
them here would silently resize every volume computed afterwards. It is also
FPS's defect D12 (`SpringEngine` mutates session-global `Molecule` objects and
never reverts), which `set_bond_anchor_radii` already declined to reproduce.

**A coarse bead is not an atom, so the combination is refused.** With
`coarse_clash` the clash term is one 2.5 A sphere per residue; a table keyed by
atom name has no entry for a residue, and a scale on the bead is `bead_radius`
spelled twice. Both are errors rather than silent approximations.
"""

import json
import math
from pathlib import Path

import pytest

import IMP
import IMP.atom
import IMP.core
import IMP.bff

REPO = Path(__file__).resolve().parents[2]
HIV = REPO / "examples" / "structure" / "HIV_RT"
PDBS = [str(HIV / "protein_1R0A.pdb"), str(HIV / "dna.pdb")]
FPS_JSON = str(HIV / "hiv_rt.fps.json")


def _read(model):
    root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model, "root"))
    for path in PDBS:
        h = IMP.atom.read_pdb(path, model,
                              IMP.atom.NonWaterNonHydrogenPDBSelector())
        root.add_child(h)
        IMP.atom.create_rigid_body(h)
    model.update()
    return root


def _radii(root):
    return [IMP.core.XYZR(p.get_particle()).get_radius()
            for p in IMP.atom.get_leaves(root)]


# --------------------------------------------------------------------------
# The default changes nothing at all
# --------------------------------------------------------------------------

def test_the_default_returns_the_atoms_themselves():
    """No shadow spheres, no new particles: today's behaviour, byte for byte.

    Asserted on the particle count rather than on a score, because the whole
    claim is that the default path does not allocate a second representation.
    """
    m = IMP.Model()
    root = _read(m)
    before = len(m.get_particle_indexes())
    container = IMP.bff.clash_container(m, root, False)
    assert len(m.get_particle_indexes()) == before
    leaves = [p.get_particle_index() for p in IMP.atom.get_leaves(root)]
    assert list(container.get_contents()) == leaves


def test_a_substituted_source_makes_shadows_and_leaves_the_atoms_alone():
    """FPS's D12, not reproduced: the structure's radii are read, never written.

    The clash term measures shadow spheres; the atoms keep the radii
    `read_pdb` gave them, which is what an accessible volume on the default
    `"imp"` source inflates its obstacles by.
    """
    m = IMP.Model()
    root = _read(m)
    original = _radii(root)
    n_leaves = len(original)
    container = IMP.bff.clash_container(m, root, False, 2.5,
                                        IMP.bff.AV_RADII_OLGA)
    assert _radii(root) == original
    spheres = list(container.get_contents())
    assert len(spheres) == n_leaves
    leaves = [p.get_particle_index() for p in IMP.atom.get_leaves(root)]
    assert not set(spheres) & set(leaves)
    # The shadows carry Olga's radii, and Olga's are smaller here: that is the
    # whole point, and the mean gap is what the 17.9x in the docstring is made
    # of.
    shadow_r = [IMP.core.XYZR(m, s).get_radius() for s in spheres]
    assert sum(shadow_r) / len(shadow_r) < sum(original) / len(original)
    for s, atom in zip(spheres, leaves):
        expect = IMP.bff.olga_vdw_particle_radius(m.get_particle(atom))
        assert IMP.core.XYZR(m, s).get_radius() == pytest.approx(expect)


def test_shadows_are_not_leaves_so_nothing_that_walks_the_structure_sees_them():
    """An RMSD or a written PDB must not grow by 9034 phantom atoms."""
    m = IMP.Model()
    root = _read(m)
    n_leaves = len(IMP.atom.get_leaves(root))
    IMP.bff.clash_container(m, root, False, 2.5, IMP.bff.AV_RADII_OLGA)
    assert len(IMP.atom.get_leaves(root)) == n_leaves


def test_a_second_call_reuses_the_shadows_rather_than_doubling_them():
    """Two containers over one model must not be two sets of spheres.

    `dock_minimize` builds the assembly (which builds one container) and then
    its own; a second set would double every overlap silently.
    """
    m = IMP.Model()
    root = _read(m)
    first = list(IMP.bff.clash_container(m, root, False, 2.5,
                                         IMP.bff.AV_RADII_OLGA).get_contents())
    n = len(m.get_particle_indexes())
    second = list(IMP.bff.clash_container(m, root, False, 2.5,
                                          IMP.bff.AV_RADII_OLGA).get_contents())
    assert first == second
    assert len(m.get_particle_indexes()) == n


def test_the_scale_multiplies_whatever_the_source_gave():
    """Source picks the table, scale multiplies it -- and they compose."""
    m = IMP.Model()
    root = _read(m)
    plain = IMP.bff.clash_container(m, root, False, 2.5,
                                    IMP.bff.AV_RADII_OLGA)
    m2 = IMP.Model()
    root2 = _read(m2)
    scaled = IMP.bff.clash_container(m2, root2, False, 2.5,
                                     IMP.bff.AV_RADII_OLGA, 0.5)
    a = [IMP.core.XYZR(m, s).get_radius() for s in plain.get_contents()]
    b = [IMP.core.XYZR(m2, s).get_radius() for s in scaled.get_contents()]
    assert len(a) == len(b)
    for x, y in zip(a, b):
        assert y == pytest.approx(0.5 * x)


def test_the_scale_alone_needs_no_source():
    """"Smaller spheres" is a request on its own, and it does not touch atoms."""
    m = IMP.Model()
    root = _read(m)
    original = _radii(root)
    container = IMP.bff.clash_container(m, root, False, 2.5,
                                        IMP.bff.AV_RADII_IMP, 0.8)
    assert _radii(root) == original
    got = [IMP.core.XYZR(m, s).get_radius() for s in container.get_contents()]
    for x, y in zip(original, got):
        assert y == pytest.approx(0.8 * x)


# --------------------------------------------------------------------------
# What is refused, and why
# --------------------------------------------------------------------------

def test_a_coarse_bead_is_not_an_atom_and_the_combination_is_refused():
    m = IMP.Model()
    root = _read(m)
    with pytest.raises(Exception) as e:
        IMP.bff.clash_container(m, root, True, 2.5, IMP.bff.AV_RADII_OLGA)
    assert "coarse" in str(e.value)
    with pytest.raises(Exception):
        IMP.bff.clash_container(m, root, True, 2.5, IMP.bff.AV_RADII_IMP, 0.8)
    # ... and the default combination is still fine.
    assert IMP.bff.clash_container(m, root, True) is not None


def test_a_non_positive_scale_is_refused():
    m = IMP.Model()
    root = _read(m)
    for bad in (0.0, -1.0):
        with pytest.raises(Exception):
            IMP.bff.clash_container(m, root, False, 2.5,
                                    IMP.bff.AV_RADII_IMP, bad)


def test_an_unknown_source_is_refused_by_name():
    params = IMP.bff.DockingParameters()
    assert params.clash_radii_source == "imp"
    assert params.clash_radii_scale == 1.0
    with pytest.raises(Exception):
        IMP.bff.create_docking_assembly(PDBS, FPS_JSON, "resolved", True, 1.0,
                                       6.0, 0.0, 0.0, "bondi")


def test_docking_parameters_refuse_a_coarse_clash_with_substituted_radii(
        tmp_path):
    """The refusal happens before the output directory exists.

    A typo must not leave half a run behind, and `coarse_clash` defaults to on
    for docking, so this is the combination a first attempt will hit.
    """
    params = IMP.bff.DockingParameters()
    params.n_frames = 1
    params.clash_radii_source = "olga"
    out = tmp_path / "never"
    with pytest.raises(Exception) as e:
        IMP.bff.dock_minimize(PDBS, FPS_JSON, str(out), params)
    assert "coarse_clash" in str(e.value)
    assert not out.exists()


# --------------------------------------------------------------------------
# What it does to a score
# --------------------------------------------------------------------------

def test_the_scoring_door_is_unmoved_by_the_defaults():
    """The pinned HIV-RT number, stated with the new arguments written out."""
    a = IMP.bff.score_structures(PDBS, FPS_JSON, "resolved", False, 6.0, "")
    b = IMP.bff.score_structures(PDBS, FPS_JSON, "resolved", False, 6.0, "",
                                 "imp", 1.0)
    assert a.score == pytest.approx(b.score, abs=1e-9)
    assert a.score == pytest.approx(59.0404, abs=5e-3)


def test_olga_radii_take_the_clash_out_of_the_score():
    """The excluded-volume part of HIV-RT's score is a radii artefact.

    The scoring door runs at IMP's historical `k = 1`, not FPS's k = 8, so the
    numbers are smaller than the static census in the module docstring -- the
    ratio is what carries over.
    """
    imp_ = IMP.bff.score_structures(PDBS, FPS_JSON, "resolved", False, 6.0, "")
    olga = IMP.bff.score_structures(PDBS, FPS_JSON, "resolved", False, 6.0, "",
                                    "olga", 1.0)
    assert olga.score < imp_.score
    # The restraint half is untouched: only the clash term moved.
    assert olga.e_clash < 0.2 * imp_.e_clash
    assert (imp_.score - imp_.e_clash) == pytest.approx(
            olga.score - olga.e_clash, abs=1e-6)


def test_the_volumes_do_not_notice_the_clash_terms_radii():
    """The two halves stay separable: this knob is the clash term's alone.

    If the clash radii were written onto the particles, every volume built
    afterwards would change size -- which is exactly the failure the shadow
    spheres exist to prevent, and the restraint half of the score is where it
    would show.
    """
    imp_ = IMP.bff.score_structures(PDBS, FPS_JSON, "resolved", False, 6.0, "")
    for source, scale in (("olga", 1.0), ("imp", 0.8), ("olga", 0.85)):
        got = IMP.bff.score_structures(PDBS, FPS_JSON, "resolved", False, 6.0,
                                       "", source, scale)
        assert (got.score - got.e_clash) == pytest.approx(
                imp_.score - imp_.e_clash, abs=1e-6), (source, scale)


def test_the_scale_is_monotone_in_the_clash_term():
    """Smaller spheres, less clash -- the property a caller reaches for it for.

    Measured statically over the same interface at k = 8, the pair count falls
    268 -> 91 -> 50 -> 20 -> 1 across these five scales; at 0.7 the last pair
    overlaps by less than 0.005 A, so the term is effectively gone.
    """
    last = float("inf")
    for scale in (1.0, 0.9, 0.85, 0.8, 0.7):
        got = IMP.bff.score_structures(PDBS, FPS_JSON, "resolved", False, 6.0,
                                       "", "imp", scale)
        assert got.e_clash <= last, scale
        last = got.e_clash
    assert last < 1e-3


# --------------------------------------------------------------------------
# The project carries all three, because one without the others cannot run
# --------------------------------------------------------------------------

def test_a_project_round_trips_the_clash_radii(tmp_path):
    p = IMP.bff.FPSProject()
    assert p.coarse_clash is True
    assert p.clash_radii_source == "imp"
    assert p.clash_radii_scale == 1.0
    p.coarse_clash = False
    p.clash_radii_source = "olga"
    p.clash_radii_scale = 0.9
    payload = json.loads(p.get_json())
    assert payload["clash_radii_source"] == "olga"
    assert payload["clash_radii_scale"] == pytest.approx(0.9)
    assert payload["coarse_clash"] is False
    back = IMP.bff.fps_project_from_json(p.get_json())
    assert back.clash_radii_source == "olga"
    assert back.clash_radii_scale == pytest.approx(0.9)
    assert back.coarse_clash is False
    params = back.get_docking_parameters("Dock")
    assert params.clash_radii_source == "olga"
    assert params.clash_radii_scale == pytest.approx(0.9)
    assert params.coarse_clash is False


def test_a_project_that_cannot_run_says_so_before_it_is_run():
    p = IMP.bff.FPSProject()
    p.clash_radii_source = "olga"  # coarse_clash still on
    errors = list(p.get_problems())
    assert any("coarse_clash" in e for e in errors), errors
    p.coarse_clash = False
    assert not any("coarse_clash" in e for e in p.get_problems())
    p.clash_radii_scale = 0.0
    assert any("clash_radii_scale" in e for e in p.get_problems())
