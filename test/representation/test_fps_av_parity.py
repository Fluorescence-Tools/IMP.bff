"""`IMP.bff`'s accessible volume against the FPS toolkit's own (PRD-121).

FPS (Kalinin *et al.*, *Nat. Methods* **9**, 1218, 2012) shipped two accessible
volumes it had computed itself, on a frame of an HIV-1 reverse transcriptase /
DNA complex. Both are reproduced **exactly** by compiling FPS's native routine
and running it (`prototypes/fps_oracle/`), so the numbers pinned here are FPS's
and not a retelling of them.

Two things about FPS's `.xyz` that a reader has to know, because getting either
wrong makes agreement look like disagreement:

* **The line count is not a voxel count.** `calculate3R` returns `n += dn`, the
  *sum* of densities -- one per (voxel, dye radius) that fits -- and `AVEngine`
  emits a voxel once per radius that fits. `p_1bp` is 5473 lines over 3187
  voxels, `p66` 130 531 over 55 514.
* **`Dmp` is the density-weighted mean**, i.e. the mean of all the lines.

`IMP.bff` carries the same information as unique voxels plus a weight column.
Same content, different packing.

A third: an exact voxel-set comparison reports **zero** overlap however well
the two agree, because `IMP.bff` anchors its lattice on the global absolute
grid (PRD-105) and FPS on the source atom. The volumes are compared by count
and by mean position, and the offset is bounded separately.
"""
import json
import os

import numpy as np

import IMP
import IMP.atom
import IMP.bff
import IMP.core
import IMP.test

HERE = os.path.dirname(os.path.abspath(__file__))
INPUT = os.path.join(HERE, "..", "input", "fps")
PINS = os.path.join(HERE, "..", "references", "fps_av_pins.json")
FRAME = os.path.join(INPUT, "hivrt_frame03991.pdb")

#: FPS's own van der Waals table (`Fps/data/vdW.txt`) is Bondi; IMP's
#: `read_pdb` assigns united-atom radii that carry implicit hydrogens (carbon
#: 1.85-2.275 against 1.70). Reproducing FPS means reading FPS's radii.
BONDI = {1: 1.09, 6: 1.70, 7: 1.55, 8: 1.52, 15: 1.80, 16: 1.80}


def _load_cloud(path):
    """FPS's `.xyz`: the duplicate-expanded cloud and its `Dmp` line."""
    points, mean = [], None
    with open(path) as fh:
        lines = fh.read().splitlines()
    for line in lines[2:]:
        token = line.split()
        if not token:
            continue
        if token[0].endswith("mp"):
            mean = np.array([float(v) for v in token[1:4]])
        else:
            points.append([float(v) for v in token[1:4]])
    return np.asarray(points), mean, int(lines[0].strip())


def _structure(hydrogens):
    model = IMP.Model()
    selector = (IMP.atom.NonWaterPDBSelector() if hydrogens
                else IMP.atom.NonWaterNonHydrogenPDBSelector())
    hier = IMP.atom.read_pdb(FRAME, model, selector)
    leaves = IMP.atom.get_leaves(hier)
    for particle in leaves:
        element = int(IMP.atom.Atom(particle).get_element())
        IMP.core.XYZR(particle).set_radius(BONDI.get(element, 1.80))
    serials = {IMP.atom.Atom(p).get_input_index(): p for p in leaves}
    return model, serials


def _volume(model, serials, pin):
    particle = IMP.Particle(model)
    IMP.bff.ProbeAccessibleVolumeDecorator.do_setup_particle(
        model, particle, serials[pin["atom_serial"]],
        linker_length=pin["linker_length"], linker_width=pin["linker_width"],
        radii=tuple(pin["radii"]),
        simulation_grid_resolution=pin["grid"],
        contact_volume_thickness=0.0, contact_volume_trapped_fraction=-1)
    av = IMP.bff.ProbeAccessibleVolumeDecorator(model, particle)
    # FPS's radii. `_structure()` writes Bondi onto every particle above, and
    # this says out loud that the volume must read *those* -- the particles' --
    # rather than a table of its own. It is the default again since 2026-09-01
    # (`ProbeAccessibleVolumeDecorator::set_radii_source`, spelled `"imp"`), so the call is currently
    # redundant; it is kept because it is the whole premise of this A/B. For
    # one day the default was `"olga"`, which made the Bondi substitution dead
    # code and silently turned this into a comparison against a radii set FPS
    # has never heard of -- and, the two tables being close (both Bondi-derived,
    # carbon 1.70 in each), it kept passing while it did so. A default this
    # test follows silently is a default this test cannot detect.
    av.set_radii_source("imp")
    av.resample()
    points = np.asarray(av.get_map().get_xyz_density(), dtype=float)
    return av, points


class Tests(IMP.test.TestCase):

    @classmethod
    def setUpClass(cls):
        IMP.set_log_level(IMP.SILENT)
        with open(PINS) as fh:
            cls.pins = json.load(fh)

    # -- the shipped reference cloud, read the way FPS wrote it -------------

    def test_the_shipped_cloud_is_duplicate_expanded(self):
        """The line count is the sum of densities, not the volume."""
        pin = self.pins["p_1bp"]
        points, mean, header = _load_cloud(os.path.join(INPUT, "p_1bp_D.xyz"))
        self.assertEqual(len(points), pin["lines"])
        self.assertEqual(header, pin["lines"] + 1)   # FPS writes n+1
        unique = {tuple(np.round(p, 3)) for p in points}
        self.assertEqual(len(unique), pin["unique_voxels"])
        self.assertLess(len(unique), len(points))    # the point of the test

    def test_the_reference_mean_is_density_weighted(self):
        """`Dmp` is the mean of every line, duplicates included."""
        pin = self.pins["p_1bp"]
        points, mean, _ = _load_cloud(os.path.join(INPUT, "p_1bp_D.xyz"))
        np.testing.assert_allclose(mean, points.mean(axis=0), atol=5e-4)
        np.testing.assert_allclose(mean, pin["weighted_mean"], atol=1e-9)

    # -- parity ------------------------------------------------------------

    def test_p_1bp_matches_fps(self):
        """The DNA site, on FPS's own atoms and radii.

        Reading everything FPS reads -- all 17 733 atoms, Bondi radii -- and
        letting the clearance derive, the volume is within a few percent of
        FPS's and the mean position within a quarter of an Angstrom.
        """
        pin = self.pins["p_1bp"]
        model, serials = _structure(hydrogens=True)
        av, points = _volume(model, serials, pin)
        ratio = len(points) / pin["unique_voxels"]
        self.assertGreater(ratio, 0.93, "volume lost against FPS")
        self.assertLess(ratio, 1.07, "volume gained against FPS")
        offset = np.linalg.norm(np.asarray(av.get_mean_position(False))
                                - np.asarray(pin["weighted_mean"]))
        self.assertLess(offset, 0.30)

    def test_p66_matches_fps_on_heavy_atoms(self):
        """The protein site -- and the one place the two still part.

        With explicit hydrogens this volume comes back **empty**: the free
        channels out of a buried CB are narrower than a voxel, and FPS gets
        out because its link search hops `linknodes = 3` voxels and tunnels
        through them, which is exactly the leak PRD-105 closed deliberately.
        On heavy atoms -- this module's own convention, and what
        `create_docking_assembly` reads -- the two agree to a few percent.
        """
        pin = self.pins["p66"]
        model, serials = _structure(hydrogens=False)
        av, points = _volume(model, serials, pin)
        ratio = len(points) / pin["unique_voxels"]
        self.assertGreater(ratio, 0.85)
        self.assertLess(ratio, 1.15)
        offset = np.linalg.norm(np.asarray(av.get_mean_position(False))
                                - np.asarray(pin["weighted_mean"]))
        self.assertLess(offset, 1.5)

    def test_p66_with_hydrogens_is_empty_and_says_so(self):
        """The parting recorded, so that closing it is deliberate."""
        pin = self.pins["p66"]
        model, serials = _structure(hydrogens=True)
        av, points = _volume(model, serials, pin)
        self.assertEqual(len(points), 0)
        # an empty volume reports its anchor, and nothing to average is NaN
        np.testing.assert_allclose(av.get_mean_position(),
                                   av.get_source_coordinates(), atol=1e-6)
        self.assertTrue(np.all(np.isnan(
            np.asarray(av.get_mean_position(False)))))

    def test_the_two_lattices_are_offset_not_different(self):
        """Every voxel this module finds is within one voxel of an FPS voxel.

        The volumes never share a coordinate -- different lattice anchors --
        so agreement is measured as coverage at the lattice spacing. Anything
        further away would be volume FPS does not have.
        """
        pin = self.pins["p_1bp"]
        reference, _, _ = _load_cloud(os.path.join(INPUT, "p_1bp_D.xyz"))
        reference = np.unique(np.round(reference, 3), axis=0)
        model, serials = _structure(hydrogens=True)
        _, points = _volume(model, serials, pin)
        self.assertGreater(len(points), 0)

        spacing = pin["grid"]
        cells = {}
        for point in reference:
            key = tuple(np.round(point / spacing).astype(int))
            cells.setdefault(key, []).append(point)

        def near(p, tol):
            base = np.round(p[:3] / spacing).astype(int)
            # search far enough for the tolerance; a fixed +-1 neighbourhood
            # silently under-counts as soon as tol exceeds one spacing
            reach = int(np.ceil(tol / spacing))
            span = range(-reach, reach + 1)
            for dx in span:
                for dy in span:
                    for dz in span:
                        for q in cells.get((base[0] + dx, base[1] + dy,
                                            base[2] + dz), ()):
                            if np.linalg.norm(p[:3] - q) <= tol:
                                return True
            return False

        # Measured 2026-08-31: 99.04 % within one voxel, 100 % within two,
        # furthest 1.130 A. The reverse holds too -- 99.0 % of FPS's voxels are
        # within one voxel of one of these, furthest 1.051 A -- so the two
        # volumes are the same volume on lattices that do not line up, not two
        # different answers. A voxel further than *two* spacings from anything
        # FPS found would be volume this module invented; there is none.
        within_one = sum(1 for p in points if near(p, spacing))
        within_two = sum(1 for p in points if near(p, 2.0 * spacing))
        self.assertGreater(within_one / len(points), 0.98)
        self.assertEqual(within_two, len(points),
                         "this module found volume FPS does not have")


if __name__ == "__main__":
    IMP.test.main()
