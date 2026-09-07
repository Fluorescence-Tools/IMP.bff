"""End-to-end: FPS's published screening numbers for T4L (PRD-121).

`prototypes/fast_label_score/cache/anisotropy/zenodo3376527_FRET_screening/`
holds what the real FPS produced for 421 T4L structures -- 33 donor-acceptor
pairs of <R_DA> and chi2 for three states (Sanabria *et al.*, *Nat. Commun.*
**11**, 1231, 2020, Zenodo 3376527). One of those 421 is **3GUN**, which this
module already ships, so the whole pipeline can be measured against FPS's own
answer without shipping a single new structure.

This is the complement to `test_fps_av_parity.py`. That one pins a single
volume against a single cloud, where a disagreement can be localised; this one
pins the whole path -- volumes, <R_DA>, the asymmetric chi2 -- against 33
numbers, where a systematic bias cannot hide. Both are needed, and this one is
the one that found a bias.

There are **two** references here and they are not the same yardstick.

*Against FPS itself* -- its own routine run on the same structure with the same
file (`prototypes/fps_oracle/`, pinned in
`references/fps_screening_oracle_pins.json`) -- the two agree to **-0.31 A,
rmsd 1.20 A, r = 0.9926** over the 24 pairs that resolve under the default
radii, and that comparison must be made with the contact volume **off**,
because FPS has none.

Every number here depends on which van der Waals radii the volumes inflate
their obstacles by, so every one of them is labelled with a radii source. The
default is **IMP's** -- the radii on the particles, which is what a docking
score's clash term reads (`AV::set_radii_source`). Olga's table is selectable
and reproduces the published numbers better; what that trade costs is measured
here rather than argued about.

*Against the published table*, neither reproduced it while the contact volume
the file asks for was inert: this module by +2.05 A, FPS's own routine by
+2.37 A. Honouring it (2026-09-01, PRD-121 G9) closes that: **+0.22 A bias,
0.91 A rmsd, r = 0.9960** over the 24 pairs, and **-0.02 A / 0.71 A / 0.9954
over all 33** if Olga's radii are selected too. The table was not produced by
FPS, which has no contact volume, but by an ACV-capable program reading this
file -- so the file *is* the parameter set it was computed with, which is the
opposite of what this test recorded on the first pass.

What is established, and what is not, is in
`okf/validation/fps_screening_ab.md`.
"""
import csv
import json
import math
import os
import tempfile

import numpy as np
import pytest

import IMP
import IMP.atom
import IMP.bff
import IMP.core

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
ZENODO = os.path.join(
    REPO, "prototypes", "fast_label_score", "cache", "anisotropy",
    "zenodo3376527_FRET_screening", "FRET_screening_of_PDB_structures")
TABLE = os.path.join(ZENODO, "FRET_screening_of_PDB_structures.dat")
FPS_JSON = os.path.join(ZENODO, "FRET_screening.fps.json")

pytestmark = pytest.mark.skipif(
    not os.path.exists(TABLE),
    reason="the Zenodo FPS screening reference is not in this checkout")


def _reference_row(pdb_id):
    with open(TABLE) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header = rows[1]
    for row in rows[2:]:
        if row and row[0].strip() == pdb_id:
            return header, dict(zip(header, row))
    raise AssertionError("%s is not in the reference table" % pdb_id)


def _scored(derive_clearance, contact_volume=False, radii_source=None):
    """3GUN through the whole path, with FPS's own labelling file.

    ``allowed_sphere_radius: 2`` in that file is an ``IMP.bff`` field somebody
    chose, not an FPS parameter -- FPS seeds its search from
    ``LinkerInitialSphere * linker_width`` and has no such key. Dropping it
    lets the clearance derive, which is what the file would mean if it had
    never carried the field.

    ``contact_volume_thickness``/``_trapped_fraction`` are the same kind of
    key, and are dropped for the same reason: **FPS has no accessible-contact
    volume at all**, so a run that honours them is not the run FPS made. They
    used to be safe to leave in because they did nothing (PRD-121 G9); since
    they were wired up they are not, and leaving them in silently turned this
    A/B into a comparison of two different models. ``contact_volume=True``
    keeps them, which is how the size of that difference is measured.

    ``radii_source`` is left alone by default, so the file gets the module
    default -- ``"imp"``, the radii the particles carry, chosen so that a
    volume and the clash term beside it size an atom the same way.
    ``"olga"`` selects the name-keyed table this file's fitted trapped
    fractions and the published table were actually produced with. Both are
    exercised here, always side by side: the difference between them is a real
    quantity and every claim in this file names which one it is about.
    """
    document = json.load(open(FPS_JSON))
    for position in document["Positions"].values():
        if derive_clearance:
            position.pop("allowed_sphere_radius", None)
        if not contact_volume:
            position.pop("contact_volume_thickness", None)
            position.pop("contact_volume_trapped_fraction", None)
        if radii_source is not None:
            position["radii_source"] = radii_source
    path = os.path.join(tempfile.mkdtemp(), "screening.fps.json")
    with open(path, "w") as fh:
        json.dump(document, fh)
    return IMP.bff.score_structures(
        [IMP.bff.get_example_path("structure/T4L/3GUN.pdb")],
        path, "chi2_C1_33p")


def _paired(result, reference, header):
    """FPS's <R_DA> against this module's, for the pairs that have both."""
    mine = {}
    unresolved = []
    for pair in result.pairs:
        key = pair.name.replace("_C1", "").replace("-", "_")
        if math.isfinite(pair.distance_model):
            mine[key] = pair.distance_model
        else:
            unresolved.append(pair.name)
    columns = [c for c in header if "_" in c and not c.startswith("χ")]
    theirs = np.array([float(reference[c]) for c in columns if c in mine])
    ours = np.array([mine[c] for c in columns if c in mine])
    return theirs, ours, unresolved, len(columns)


def test_the_reference_row_is_there():
    header, row = _reference_row("3GUN")
    assert float(row["χ²_C1"]) == pytest.approx(57.7)
    assert float(row["19_119"]) == pytest.approx(46.0)


ORACLE_PINS = os.path.join(HERE, "..", "references",
                           "fps_screening_oracle_pins.json")


def test_rda_agrees_with_fps_itself():
    """The real A/B: this module against FPS's own routine, pair by pair.

    Under the **default radii** -- IMP's, the ones on the particles -- 24 of
    the 33 pairs resolve and they agree to bias **-0.31 A**, rmsd **1.20 A**,
    r = **0.9926**. The residual is the two search conventions, not a
    modelling difference. Nine pairs have no model value: they are the pairs
    with site 132, which united-atom radii wall in.

    Under **Olga's** radii, which have to be asked for, two things change and
    they have to be read together:

    * **all 33 pairs resolve**. Olga's heavy atoms are smaller, so no site of
      this file comes back empty.
    * over that larger set the agreement is **+0.02 A bias, 1.62 A rmsd,
      r = 0.9777** -- a *smaller* bias and a *worse* rmsd and correlation than
      the 24-pair figure, because the nine pairs it adds are the hard ones.
    * over the *same* 24 it is strictly better: **-0.62 A / 1.13 A /
      r = 0.9940** against -0.31 / 1.20 / 0.9926.

    So the default is **not** the more accurate radii set, and that is
    asserted here rather than left implicit: it is the set the clash term of a
    docking score reads (`AV::set_radii_source`), and internal consistency was
    judged worth more than the fit. Both halves are pinned so neither can be
    quoted without the other.
    """
    with open(ORACLE_PINS) as fh:
        oracle = json.load(fh)["rda_mean"]

    # --- the default: IMP's radii ------------------------------------------
    result = _scored(derive_clearance=True)
    mine, theirs = [], []
    for pair in result.pairs:
        if pair.name in oracle and math.isfinite(pair.distance_model):
            mine.append(pair.distance_model)
            theirs.append(oracle[pair.name])
    mine, theirs = np.array(mine), np.array(theirs)

    assert len(mine) == 24, "the set of resolvable pairs moved"
    assert abs((mine - theirs).mean()) < 0.5
    assert np.sqrt(((mine - theirs) ** 2).mean()) < 1.5
    assert np.corrcoef(theirs, mine)[0, 1] > 0.99

    # --- and Olga's, which resolves all 33 and fits the shared 24 better ----
    with_olga = _scored(derive_clearance=True, radii_source="olga")
    olga = {p.name: p.distance_model for p in with_olga.pairs
            if math.isfinite(p.distance_model)}
    assert len(olga) == 33, "Olga's radii no longer open every site"

    shared = [p.name for p in result.pairs
              if p.name in olga and math.isfinite(p.distance_model)]
    a = np.array([olga[n] for n in shared])
    b = np.array([m for p, m in
                  ((p, p.distance_model) for p in result.pairs)
                  if p.name in shared])
    ref = np.array([oracle[n] for n in shared])
    assert len(shared) == 24
    assert np.sqrt(((a - ref) ** 2).mean()) < np.sqrt(((b - ref) ** 2).mean())
    assert np.corrcoef(ref, a)[0, 1] > np.corrcoef(ref, b)[0, 1]

    # ...and over its own 33 the rmsd is worse than the 24-pair figure, which
    # is the "larger, harder set" half of the finding.
    all_a = np.array([olga[n] for n in olga])
    all_ref = np.array([oracle[n] for n in olga])
    assert abs((all_a - all_ref).mean()) < 0.3          # bias shrinks
    assert np.sqrt(((all_a - all_ref) ** 2).mean()) > 1.5   # rmsd grows


def test_the_published_table_is_not_reproduced_by_fps_either():
    """FPS does not reproduce the published table, pinned so it stays known.

    Without a contact volume both implementations sit ~2 A above the published
    <R_DA>, FPS's own by *more* than this module's (+2.37 A against +2.05 A).
    Neither R_mp (-1.40 A) nor <R_DA>_E (+3.11 A) matches it either, so it is
    not a mislabelled distance type, and FPS's own AV3 radii from `linker.txt`
    make it worse (+6.36 A).

    This test was first read as "the file is not the parameter set the table
    was computed with". It is not: what the file asks for and FPS cannot do is
    the **accessible-contact volume**, and honouring it brings this module to
    +0.22 A of the table -- see
    :func:`test_the_contact_volume_explains_the_offset_against_the_published_table`.
    What survives is the narrower statement this test actually makes, and it is
    worth keeping: *FPS* does not reproduce that table, so the table cannot pin
    an FPS port. The pin for that is the oracle.
    """
    header, reference = _reference_row("3GUN")
    theirs, ours, _, _ = _paired(
        _scored(derive_clearance=True), reference, header)
    with open(ORACLE_PINS) as fh:
        oracle = json.load(fh)["rda_mean"]
    assert 1.0 < (ours - theirs).mean() < 3.0        # this module
    # and FPS's own, over the same pairs
    result = _scored(derive_clearance=True)
    fps_bias = []
    for pair in result.pairs:
        if pair.name not in oracle or not math.isfinite(pair.distance_model):
            continue
        column = "%s_%s" % (pair.position1[:-1], pair.position2[:-1])
        if column not in reference:
            column = "%s_%s" % (pair.position2[:-1], pair.position1[:-1])
        if column in reference:
            fps_bias.append(oracle[pair.name] - float(reference[column]))
    assert np.mean(fps_bias) > 1.5, (
        "FPS itself no longer misses the published table; if that changed, "
        "the file may finally match and this test should be rewritten")


def test_the_authored_clearance_buries_the_structure_and_it_is_the_radii():
    """Two findings, each asserted against the radii source it belongs to.

    **1. Under the default radii the authored clearance buries the structure.**
    With the file's ``allowed_sphere_radius: 2``, eight of its seventeen sites
    come back empty and only **nine** of 33 pairs have a model value, where
    deriving the clearance leaves one empty and gives **24**. The two
    clearances are not the same quantity: FPS seeds from
    ``LinkerInitialSphere * linker_width`` unconditionally, this module carves
    a sphere out of an obstacle map already inflated by half the linker width.
    Any fps.json carrying a small explicit clearance is asking for something
    other than it appears to.

    **2. What actually buries it is the radii, not the clearance.** Ask for
    Olga's table and the *same* authored clearance resolves **all 33** pairs,
    the same as deriving it. A carbon of 2.275 A walls a source in where one
    of 1.70 A does not.

    Those two are only compatible if the clearance convention is a real
    difference that is nonetheless *not* what empties the volumes, and that
    attribution is the point of this test. It is asserted both ways round --
    each behaviour against its own radii source -- because either one alone
    reads as the whole story and is not. The finding survived the default
    flipping to Olga's radii on 2026-09-01 and flipping back the same day; it
    is a statement about the two radii sets, not about which is default.
    """
    header, reference = _reference_row("3GUN")

    # 1. the default, IMP's radii: authored buries, derived does not
    _, as_authored, _, _ = _paired(
        _scored(derive_clearance=False), reference, header)
    _, derived, _, _ = _paired(
        _scored(derive_clearance=True), reference, header)
    assert len(as_authored) == 9
    assert len(derived) == 24

    # 2. ...and it was the radii: Olga's table resolves all 33 either way
    _, olga_authored, _, _ = _paired(
        _scored(derive_clearance=False, radii_source="olga"),
        reference, header)
    _, olga_derived, _, _ = _paired(
        _scored(derive_clearance=True, radii_source="olga"),
        reference, header)
    assert len(olga_authored) == len(olga_derived) == 33


def test_the_contact_volume_moves_the_cloud():
    """An ACV request changes the cloud it asks to reweight (PRD-121 G9).

    This was a strict ``xfail`` until 2026-09-01: ``contact_volume_thickness``
    and ``contact_volume_trapped_fraction`` were read into the decorator and
    never used, so the volume, the per-voxel weights and the mean position came
    back bit-identical with the weighting on and off.

    The cloud is deliberately *not* expected to change shape -- an ACV is a
    weighting, the same voxels with different weights -- so what is asserted is
    that the points are the same, the weights are not, and the mean moves
    towards the protein. Measured on T4L 132/CB at these parameters: two weight
    levels, 11766 contact voxels of 33174 carrying 0.55 of the total, and the
    mean 1.24 A away from the unweighted one.
    """
    IMP.set_log_level(IMP.SILENT)
    model = IMP.Model()
    hier = IMP.atom.read_pdb(
        IMP.bff.get_example_path("structure/T4L/3GUN.pdb"), model)

    def cloud(thickness, trapped):
        selection = IMP.atom.Selection(hier)
        selection.set_atom_type(IMP.atom.AtomType("CB"))
        selection.set_residue_index(132)
        particle = IMP.Particle(model)
        IMP.bff.AV.do_setup_particle(
            model, particle, selection.get_selected_particles()[0],
            linker_length=22.0, radii=(3.5, 0.0, 0.0), linker_width=1.5,
            simulation_grid_resolution=0.7,
            contact_volume_thickness=thickness,
            contact_volume_trapped_fraction=trapped)
        av = IMP.bff.AV(model, particle)
        av.resample()
        return np.asarray(av.get_map().get_xyz_density(), dtype=float)

    plain = cloud(0.0, -1.0)
    contact = cloud(3.0, 0.55)
    assert len(plain) > 0
    assert not np.array_equal(plain, contact)
    # the same voxels: the ACV re-weights the cloud, it does not carve it
    assert np.array_equal(plain[:, :3], contact[:, :3])

    weight = contact[:, 3]
    levels = np.unique(weight)
    assert len(levels) == 2, "uniform base weights give exactly two ACV levels"
    trapped = weight[weight == levels[1]].sum() / weight.sum()
    assert trapped == pytest.approx(0.55, abs=1e-6)

    def mean(cloud_):
        return (cloud_[:, :3] * cloud_[:, 3:4]).sum(0) / cloud_[:, 3].sum()

    assert np.linalg.norm(mean(contact) - mean(plain)) > 0.5


def _against(reference, contact_volume, by_name):
    """This module's <R_DA> against a reference dict, over the shared pairs."""
    result = _scored(derive_clearance=True, contact_volume=contact_volume)
    mine, theirs = [], []
    for pair in result.pairs:
        if not math.isfinite(pair.distance_model):
            continue
        key = by_name(pair)
        if key in reference:
            mine.append(pair.distance_model)
            theirs.append(float(reference[key]))
    return np.array(mine), np.array(theirs)


def test_the_contact_volume_is_not_part_of_the_fps_comparison():
    """Honouring the ACV moves this module away from FPS, as it must.

    `FRET_screening.fps.json` asks for a contact volume at every site and FPS
    has none, so a run that honours the request is not the run FPS made. While
    the fields were inert this could be ignored; it cannot be since they were
    wired up (PRD-121 G9), which is why :func:`_scored` drops them. Pinned so
    the A/B cannot quietly become a comparison of two different models.

    Measured 2026-09-01 over the same 24 pairs, under the default (IMP) radii:
    **-0.31 A** against FPS with the ACV off, **-2.15 A** with it on, every one
    of the 24 <R_DA> shorter, by 1.83 A on average.
    """
    with open(ORACLE_PINS) as fh:
        oracle = json.load(fh)["rda_mean"]
    by_name = lambda pair: pair.name
    off, fps = _against(oracle, False, by_name)
    on, _ = _against(oracle, True, by_name)
    assert len(off) == len(on) == len(fps) >= 24
    assert abs((off - fps).mean()) < 1.0            # the A/B, as pinned above
    assert (on - fps).mean() < -1.5                 # the ACV, on top of it
    assert (on < off).all()


def test_the_contact_volume_explains_the_offset_against_the_published_table():
    """What the ACV was hiding: the table *was* computed with one.

    `okf/validation/fps_screening_ab.md` recorded a +2.05 A offset against
    Zenodo 3376527's published <R_DA> that neither the distance type, nor the
    AV1/AV3 choice, nor FPS's own routine could account for -- FPS misses that
    table by +2.37 A, *more* than this module did -- and concluded the file was
    not the parameter set the table was computed with.

    Wiring up the contact volume that file asks for at every site closes it:
    bias **+2.06 -> +0.22 A**, rmsd **2.54 -> 0.91 A**, r 0.9901 -> 0.9960 over
    the 24 pairs that then resolved. The file *is* the parameter set; what was
    missing was the feature it names. The table was therefore not produced by
    FPS, which has no contact volume, but by an ACV-capable program reading
    this file.

    That program was Olga, and asking for Olga's radii (`AV::set_radii_source`)
    improves it again and widens it: **-0.02 A bias, 0.71 A rmsd, r = 0.9954
    over all 33 published pairs**, where under the default nine of them have no
    model value at all. Over the same 24 it is -0.23 A / 0.72 A / 0.9966.

    That is *not* the default, and the reason is in `AV::set_radii_source`:
    the clash term of a docking score reads the particles' radii, so the volume
    must too. What the default costs against this reference is exactly the gap
    asserted below, and `okf/validation/fps_screening_ab.md` states it as a
    price rather than hiding it. Anyone reproducing Olga-era numbers should
    select `radii_source = "olga"`.

    Kept as a test rather than a note because it is the strongest external
    check the module has: 33 published numbers, reproduced to within 0.72 A
    rmsd by a path with no fitted parameter in it.
    """
    header, reference = _reference_row("3GUN")
    by_name = lambda pair: pair.name.replace("_C1", "").replace("-", "_")
    off, published = _against(reference, False, by_name)
    on, _ = _against(reference, True, by_name)
    assert len(on) >= 24

    assert (off - published).mean() > 1.5           # the offset, as recorded
    assert abs((on - published).mean()) < 0.75      # and what closes it
    assert np.sqrt(((on - published) ** 2).mean()) < 1.5
    assert (np.sqrt(((on - published) ** 2).mean())
            < 0.5 * np.sqrt(((off - published) ** 2).mean()))
    assert np.corrcoef(published, on)[0, 1] > 0.995

    # ...and what the default costs against this reference, measured. Olga's
    # radii cover all 33 published pairs and fit them better; the default
    # covers 24 and fits them worse. Asserted so the price cannot drift out of
    # the page that quotes it.
    olga_result = _scored(derive_clearance=True, contact_volume=True,
                          radii_source="olga")
    olga_mine, olga_theirs = [], []
    for pair in olga_result.pairs:
        if not math.isfinite(pair.distance_model):
            continue
        key = by_name(pair)
        if key in reference:
            olga_mine.append(pair.distance_model)
            olga_theirs.append(float(reference[key]))
    olga_mine = np.array(olga_mine)
    olga_theirs = np.array(olga_theirs)
    assert len(olga_mine) == 33 and len(on) == 24
    assert abs((olga_mine - olga_theirs).mean()) < abs((on - published).mean())
    assert (np.sqrt(((olga_mine - olga_theirs) ** 2).mean())
            < np.sqrt(((on - published) ** 2).mean()))
