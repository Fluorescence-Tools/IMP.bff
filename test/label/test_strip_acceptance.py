"""PRD-106's acceptance criteria, as tests.

1. the authored masks apply -- every `strip_mask` in the shipped fps.json
   files removes exactly the atoms a hand-written per-atom list removes;
3. AV parity -- a position with its mask honoured gives the same volume as
   one built from a structure that was stripped beforehand;
5. the strip report says what a mask did, including when it did nothing.
"""

import json
from pathlib import Path

import pytest

import IMP
import IMP.atom
import IMP.bff

REPO = Path(__file__).resolve().parents[2]
BACKBONE = {"CA", "CB", "C", "N", "O"}


def _structure(pdb="structure/T4L/3GUN.pdb"):
    m = IMP.Model()
    h = IMP.atom.read_pdb(IMP.bff.get_example_path(pdb), m,
                          IMP.atom.NonWaterNonHydrogenPDBSelector())
    return m, h


#: The structure each shipped fps.json is authored against. A mask names
#: residues of *its own* structure: checking a TG2 mask against T4L selects
#: nothing, and an assertion that nothing equals nothing passes while testing
#: nothing.
STRUCTURE_OF = {
    "fret.fps.json": "structure/T4L/3GUN.pdb",
    "flex.fps.json": "structure/TG2/topology.pdb",
}


def _shipped_masks():
    """Every authored `strip_mask` in the repository, with its structure."""
    out = []
    for path in sorted(REPO.glob("examples/**/*.fps.json")):
        try:
            doc = json.loads(path.read_text())
        except ValueError:
            continue
        structure = STRUCTURE_OF.get(path.name)
        if structure is None:
            continue
        for name, position in doc.get("Positions", {}).items():
            mask = position.get("strip_mask")
            if mask:
                out.append((path.name, name, mask, structure))
    return out


def test_every_shipped_fps_json_has_a_structure_to_check_against():
    """Otherwise a file's masks are silently not covered."""
    seen = {p.name for p in REPO.glob("examples/**/*.fps.json")
            if json.loads(p.read_text()).get("Positions") and
            any(q.get("strip_mask") for q in
                json.loads(p.read_text())["Positions"].values())}
    assert seen <= set(STRUCTURE_OF), (
        "these ship masks but no structure is named for them: "
        + str(sorted(seen - set(STRUCTURE_OF))))


def test_the_repository_ships_masks_to_check():
    assert len(_shipped_masks()) >= 12


@pytest.mark.parametrize("file_name,position,mask,structure",
                         _shipped_masks())
def test_an_authored_mask_removes_exactly_what_it_names(file_name, position,
                                                        mask, structure):
    """Acceptance 1: against a hand-written per-atom list.

    Every shipped mask reads `(resid N and not name CA+CB+C+N+O) or resname
    <solvent names>`, so the atoms it must remove are: the non-backbone atoms
    of residue N, plus anything whose residue name is in the solvent list.
    """
    m, h = _structure(structure)
    report = IMP.bff.strip_report(h, mask)

    # what the document names, worked out by hand from the mask's own text
    resid = int(mask.split("resid ")[1].split(" ")[0].rstrip(")"))
    solvent = set()
    if "resname" in mask:
        solvent = {w.strip("'\"") for w in
                   mask.split("resname ")[1].rstrip(")").split()}

    expected = set()
    for i, atom in enumerate(IMP.bff.selection_atoms(h)):
        key = "{}/{}/{}".format(atom.chain, atom.resi, atom.name)
        if atom.resi == resid and atom.name not in BACKBONE:
            expected.add(key)
        elif atom.resn in solvent:
            expected.add(key)

    assert set(report.atoms) == expected, mask
    assert report.n_selected == len(expected)

    # `set() == set()` passes while testing nothing, so an empty result has to
    # be a *checked* claim: the site's residue really has no atom outside the
    # keep set. An alanine site is the honest case -- its only side-chain atom
    # is CB, which the mask keeps -- and a mask read against the wrong
    # structure is the dishonest one this catches.
    if not expected:
        residue_atoms = {a.name for a in IMP.bff.selection_atoms(h)
                         if a.resi == resid}
        assert residue_atoms, (
            "{}:{} names residue {}, which {} does not have".format(
                file_name, position, resid, structure))
        assert residue_atoms <= BACKBONE, (
            "{}:{} strips nothing, but residue {} has {} outside the keep "
            "set".format(file_name, position, resid,
                         sorted(residue_atoms - BACKBONE)))


def test_the_authored_masks_are_not_all_empty():
    """The suite above must bite on the ordinary case, not only the alanines."""
    stripping = 0
    for file_name, position, mask, structure in _shipped_masks():
        m, h = _structure(structure)
        if IMP.bff.strip_report(h, mask).n_selected > 0:
            stripping += 1
    assert stripping >= 30, (
        "only {} of {} shipped masks strip anything".format(
            stripping, len(_shipped_masks())))


def test_the_report_says_when_a_mask_selects_nothing():
    """Acceptance 5. A mask that reads and matches nothing is a fact, not an
    error -- and the report is where a caller sees it."""
    m, h = _structure()
    report = IMP.bff.strip_report(h, "resid 9999 and not name CA")
    assert report.n_selected == 0
    assert report.n_atoms > 0
    assert list(report.residues) == []


def test_av_parity_between_the_mask_and_a_pre_stripped_structure():
    """Acceptance 3: honouring the mask equals stripping beforehand.

    The decorator path used to read no `strip_mask` at all, so this is the
    test that says it does now.
    """
    mask = "resid 132 and not name CA+CB+C+N+O"

    # (a) the mask honoured, on the whole structure
    m1, h1 = _structure()
    site1 = IMP.atom.Selection(h1, residue_index=132,
                               atom_type=IMP.atom.AT_CB
                               ).get_selected_particles()[0]
    av1 = IMP.bff.AV.setup_particle(IMP.Particle(m1), site1.get_index())
    av1.set_av_parameter(json.dumps({"linker_length": 20.0, "radius1": 3.5,
                                     "linker_width": 0.5,
                                     "allowed_sphere_radius": 1.5,
                                     "simulation_grid_resolution": 1.5,
                                     "strip_mask": mask}))
    av1.resample()

    # (b) the same structure, stripped first, and no mask
    m2, h2 = _structure()
    IMP.bff.strip_hierarchy(h2, mask)
    site2 = IMP.atom.Selection(h2, residue_index=132,
                               atom_type=IMP.atom.AT_CB
                               ).get_selected_particles()[0]
    av2 = IMP.bff.AV.setup_particle(IMP.Particle(m2), site2.get_index())
    av2.set_av_parameter(json.dumps({"linker_length": 20.0, "radius1": 3.5,
                                     "linker_width": 0.5,
                                     "allowed_sphere_radius": 1.5,
                                     "simulation_grid_resolution": 1.5}))
    av2.resample()

    # Acceptance 3 says *bit-identical*, so compare the densities themselves,
    # not their lengths: two clouds of equal size can still be different
    # clouds.
    import numpy as np
    d1 = np.asarray(av1.get_map().get_xyz_density())
    d2 = np.asarray(av2.get_map().get_xyz_density())
    assert d1.size > 0
    assert d1.shape == d2.shape, "the mask and a pre-stripped structure differ"
    assert np.array_equal(d1, d2), (
        "same number of points, different points: the mask and a pre-stripped "
        "structure do not give the same volume")


def test_ignoring_the_mask_would_give_a_different_volume():
    """The parity test above proves nothing unless the mask matters here."""
    mask = "resid 132 and not name CA+CB+C+N+O"
    sizes = []
    for use_mask in (True, False):
        m, h = _structure()
        site = IMP.atom.Selection(h, residue_index=132,
                                  atom_type=IMP.atom.AT_CB
                                  ).get_selected_particles()[0]
        av = IMP.bff.AV.setup_particle(IMP.Particle(m), site.get_index())
        params = {"linker_length": 20.0, "radius1": 3.5, "linker_width": 0.5,
                  "allowed_sphere_radius": 1.5,
                  "simulation_grid_resolution": 1.5}
        if use_mask:
            params["strip_mask"] = mask
        av.set_av_parameter(json.dumps(params))
        av.resample()
        sizes.append(len(av.get_map().get_xyz_density()))
    assert sizes[0] != sizes[1], (
        "stripping the site's side chain changed nothing -- either the mask "
        "is not reaching the obstacle set, or this site has no side chain")


# --------------------------------------------------------------------------
# the document's own mask
# --------------------------------------------------------------------------

def _fps_with(tmp_path, document_mask=None, position_mask=None):
    """A one-position fps.json, with either mask, both, or neither."""
    doc = {
        "Positions": {
            "p": {"chain_identifier": "", "residue_seq_number": 132,
                  "atom_name": "CB", "simulation_type": "AV1",
                  "linker_length": 20.0, "linker_width": 0.5, "radius1": 3.5,
                  "allowed_sphere_radius": 1.5,
                  "simulation_grid_resolution": 1.5},
        },
        "Distances": {"d": {"position1_name": "p", "position2_name": "p",
                            "distance": 40.0, "error_neg": 1.0,
                            "error_pos": 1.0, "distance_type": "RDAMean",
                            "Forster_radius": 52.0}},
        "χ²": {"all": {"distances": ["d"]}},
    }
    if position_mask is not None:
        doc["Positions"]["p"]["strip_mask"] = position_mask
    if document_mask is not None:
        doc["strip_mask"] = document_mask
    path = tmp_path / "one.fps.json"
    path.write_text(json.dumps(doc))
    return str(path)


def test_the_document_mask_is_a_union_with_the_position_mask():
    """`combined_strip_mask` is `or`, not an override.

    The two answer different questions -- the document's is about the
    structure (waters, ions), the position's about its labelling site -- so a
    file that states both means both.
    """
    combined = IMP.bff.combined_strip_mask("resname HOH SOL",
                                           "resid 132 and not name CA+CB")
    assert combined == "(resname HOH SOL) or (resid 132 and not name CA+CB)"
    # either alone passes through untouched
    assert IMP.bff.combined_strip_mask("", "resid 5") == "resid 5"
    assert IMP.bff.combined_strip_mask("resname HOH", "") == "resname HOH"
    assert IMP.bff.combined_strip_mask("", "") == ""


def test_the_combination_is_exactly_the_union_of_the_two():
    """What a file that states both masks means: both, and nothing else."""
    m, h = _structure()
    document = "resid 10 or resid 11"
    position = "resid 13 and name CB"

    a = set(IMP.bff.strip_report(h, document).atoms)
    b = set(IMP.bff.strip_report(h, position).atoms)
    both = set(IMP.bff.strip_report(
        h, IMP.bff.combined_strip_mask(document, position)).atoms)

    assert both == a | b
    assert a and b and not (a & b)  # the halves are disjoint here, so this bites


def test_each_half_is_parenthesised():
    """`or`-joining two expressions is safe in this grammar -- `and` binds
    tighter than `or`, and nothing binds looser except `byres`, which is a
    prefix and reaches only rightwards. The parentheses are kept anyway: they
    cost nothing, they survive a grammar that gains a looser operator, and a
    reader of the combined text should not have to know the precedence table
    to see what was joined."""
    assert (IMP.bff.combined_strip_mask("a and b", "c or d") ==
            "(a and b) or (c or d)")


def test_a_document_mask_reaches_every_position(tmp_path):
    """Read through the document reader, a position carries the union."""
    path = _fps_with(tmp_path, document_mask="resname HOH SOL WAT",
                     position_mask="resid 132 and not name CA+CB+C+N+O")
    document = IMP.bff.read_fps_json(path, [], False)
    positions = json.loads(document.positions)
    mask = positions["p"]["strip_mask"]
    assert "resname HOH SOL WAT" in mask
    assert "resid 132" in mask


def test_a_document_mask_alone_reaches_a_position_with_none(tmp_path):
    path = _fps_with(tmp_path, document_mask="resname HOH SOL WAT")
    positions = json.loads(IMP.bff.read_fps_json(path, [], False).positions)
    assert positions["p"]["strip_mask"] == "resname HOH SOL WAT"


def test_a_file_with_no_document_mask_is_unchanged(tmp_path):
    path = _fps_with(tmp_path, position_mask="resid 132 and not name CA+CB")
    positions = json.loads(IMP.bff.read_fps_json(path, [], False).positions)
    assert positions["p"]["strip_mask"] == "resid 132 and not name CA+CB"


def test_the_document_mask_is_in_the_schema():
    schema = json.loads(IMP.bff.fps_json_schema())
    assert "strip_mask" in schema["properties"], "not a documented field"
    described = schema["properties"]["strip_mask"]["description"]
    assert "every" in described.lower()
    names = [f.name for f in IMP.bff.fps_document_fields()]
    assert "strip_mask" in names


# --------------------------------------------------------------------------
# requirement 1: the two outputs
# --------------------------------------------------------------------------

def test_the_obstacle_array_is_the_strip_as_numbers():
    """`strip_obstacles` is the second of requirement 1's two outputs, and it
    strips the same way everything else does: a **radius of zero**, every atom
    still in the array."""
    import numpy as np

    mask = "resid 132 and not name CA+CB+C+N+O"
    m, h = _structure()
    n_atoms = len(IMP.bff.selection_atoms(h))

    obstacles = np.asarray(IMP.bff.strip_obstacles(h, mask)).reshape(-1, 4)
    assert len(obstacles) == n_atoms, "rows are not dropped"

    zeroed = obstacles[:, 3] == 0.0
    assert zeroed.sum() == IMP.bff.strip_report(h, mask).n_selected

    # and they are the atoms the mask names, not merely as many
    named = {a for i, a in enumerate(IMP.bff.selection_atoms(h)) if zeroed[i]}
    assert {"{}/{}/{}".format(a.chain, a.resi, a.name) for a in named} == set(
        IMP.bff.strip_report(h, mask).atoms)


def test_the_obstacle_array_keeps_the_attachment_atom():
    """A volume must not grow through its own anchor, so the attachment atom
    keeps its size even when the mask names it."""
    mask = "resid 132"          # names the whole residue, CB included
    m, h = _structure()
    import numpy as np
    index = [i for i, a in enumerate(IMP.bff.selection_atoms(h))
             if a.resi == 132 and a.name == "CB"][0]
    without = np.asarray(IMP.bff.strip_obstacles(h, mask)).reshape(-1, 4)
    with_anchor = np.asarray(
        IMP.bff.strip_obstacles(h, mask, "A/132/CB")).reshape(-1, 4)
    assert without[index, 3] == 0.0
    assert with_anchor[index, 3] > 0.0


def test_an_empty_mask_leaves_every_radius_alone():
    import numpy as np
    m, h = _structure()
    obstacles = np.asarray(IMP.bff.strip_obstacles(h, "")).reshape(-1, 4)
    assert len(obstacles) == len(IMP.bff.selection_atoms(h))
    assert (obstacles[:, 3] > 0.0).all()


def test_a_zero_radius_row_is_transparent_to_the_array_kernel():
    """The rule, end to end through the public array API: a volume built from
    obstacles with zeroed radii equals one built from the shorter array those
    rows were dropped from."""
    import numpy as np

    mask = "resid 132 and not name CA+CB+C+N+O"
    m, h = _structure()
    atoms = IMP.bff.selection_atoms(h)
    index = [i for i, a in enumerate(atoms) if a.resi == 132 and a.name == "CB"][0]
    source = [atoms[index].x, atoms[index].y, atoms[index].z]

    zeroed = np.asarray(IMP.bff.strip_obstacles(h, mask, "A/132/CB")).reshape(-1, 4)
    dropped = zeroed[zeroed[:, 3] > 0.0]
    assert len(dropped) < len(zeroed)

    a = IMP.bff.get_av(zeroed, source, 20.0, 0.5, 3.5, 0.0, 0.0, 1.5)
    b = IMP.bff.get_av(dropped, source, 20.0, 0.5, 3.5, 0.0, 0.0, 1.5)
    assert np.array_equal(np.asarray(a.get_points()),
                          np.asarray(b.get_points())), (
        "a zero radius is not transparent to the array kernel")


def test_a_stripped_atom_keeps_its_place_and_loses_its_size():
    """How a strip reaches the obstacle set: radius zero, not a shorter list.

    A rasteriser asks `distance < radius`, so an atom of no size blocks
    nothing -- and keeping it in the list means every volume indexes the same
    particles whatever it strips, which is what the incremental window
    machinery needs, and nothing is mutated on a model another thread reads.
    """
    mask = "resid 132 and not name CA+CB+C+N+O"
    m, h = _structure()
    site = IMP.atom.Selection(h, residue_index=132,
                              atom_type=IMP.atom.AT_CB
                              ).get_selected_particles()[0]
    av = IMP.bff.AV.setup_particle(IMP.Particle(m), site.get_index())
    av.set_av_parameter(json.dumps(
        {"linker_length": 20.0, "radius1": 3.5, "linker_width": 0.5,
         "allowed_sphere_radius": 1.5, "simulation_grid_resolution": 1.5,
         "strip_mask": mask}))
    av.resample()

    n_atoms = len(IMP.bff.selection_atoms(h))
    radii = list(av.get_map().get_obstacle_radii())
    assert len(radii) == n_atoms, "the particle list is not shortened"
    zeroed = [r for r in radii if r == 0.0]
    assert len(zeroed) == IMP.bff.strip_report(h, mask).n_selected

    # the attachment atom keeps its size: a volume must not grow through its
    # own anchor
    index = [i for i, a in enumerate(IMP.bff.selection_atoms(h))
             if a.resi == 132 and a.name == "CB"][0]
    assert radii[index] > 0.0


def test_an_unmasked_volume_carries_no_override():
    """An empty mask and the model's own radii mean no override at all.

    The override vector is what a strip mask is expressed as -- a stripped
    atom keeps its place in the obstacle list and loses its *size* -- so an
    empty vector is the "nothing to say" state and must stay reachable.

    Since 2026-09-01 the radii set is the other thing that vector carries
    (AV::set_radii_source): asking for Olga's table fills it with 128 named
    radii whatever the mask says, so the empty state is `radii_source = "imp"`
    -- the default -- *and* no mask. Both halves are asserted, and the default
    is written out rather than relied on, because "the override is empty"
    passes for the wrong reason the moment the default moves.
    """
    m, h = _structure()
    site = IMP.atom.Selection(h, residue_index=132,
                              atom_type=IMP.atom.AT_CB
                              ).get_selected_particles()[0]
    parameter = json.dumps(
        {"linker_length": 20.0, "radius1": 3.5, "linker_width": 0.5,
         "allowed_sphere_radius": 1.5, "simulation_grid_resolution": 1.5})

    av = IMP.bff.AV.setup_particle(IMP.Particle(m), site.get_index())
    av.set_av_parameter(parameter)
    av.set_radii_source("imp")
    av.resample()
    assert len(list(av.get_map().get_obstacle_radii())) == 0

    # ...and Olga's table, which has to be asked for, is an override of exactly
    # one radius per leaf, none of them zero (nothing is stripped).
    olga = IMP.bff.AV.setup_particle(IMP.Particle(m), site.get_index())
    olga.set_av_parameter(parameter)
    olga.set_radii_source("olga")
    olga.resample()
    radii = list(olga.get_map().get_obstacle_radii())
    assert len(radii) == len(IMP.atom.get_leaves(h))
    assert min(radii) > 0.0
    assert set(round(r, 4) for r in radii) <= {1.0, 1.49, 1.5, 1.625, 1.7,
                                               1.782, 1.86}
