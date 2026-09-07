"""The source clearance, and what an empty accessible volume says.

Two defects found while measuring `IMP.bff` against the FPS toolkit's own AV
(PRD-121), both of which produced a plausible-looking number for a volume that
did not exist:

* the clearance derivation lived in `get_av_from_structure()` only, so the
  decorator door -- and `imp_bff av-export` behind it -- took a flat 1.5 and
  returned **nothing** at FPS's standard linker width of 4.5 A;
* `get_mean_position()` started its weight sum at 1.0, a unit of weight
  belonging to no point, so an empty volume reported exactly half the source
  coordinate and every real one was pulled toward the origin.

Neither raised, warned or set a flag. The second had been pinned as reference
data in four of the twelve `prd105_legacy_pins.json` cases.
"""
import json
import subprocess
import sys
from pathlib import Path

import IMP
import IMP.atom
import IMP.bff
import IMP.core
import IMP.test
import numpy as np

PDB = IMP.bff.get_example_path("structure/T4L/3GUN.pdb")
BIN = Path(__file__).resolve().parent.parent.parent / "bin"

# FPS's standard linker geometry: every position in its shipped test data uses
# a linker width of 4.5 A, which is exactly what used to come back empty.
FPS_WIDTH = 4.5


def _av(mdl, hier, residue_index=132, atom_name="CB", **par):
    sel = IMP.atom.Selection(hier)
    sel.set_atom_type(IMP.atom.AtomType(atom_name))
    sel.set_residue_index(residue_index)
    source = sel.get_selected_particles()[0]
    p = IMP.Particle(mdl)
    IMP.bff.AV.do_setup_particle(mdl, p, source, **par)
    return IMP.bff.AV(mdl, p)


def _n_points(av):
    return len(av.get_map().get_xyz_density())


class Tests(IMP.test.TestCase):

    def setUp(self):
        super().setUp()
        IMP.set_log_level(IMP.SILENT)
        self.mdl = IMP.Model()
        self.hier = IMP.atom.read_pdb(PDB, self.mdl)

    def test_clearance_is_derived_when_not_given(self):
        """A negative (or absent) clearance derives from width and grid."""
        av = _av(self.mdl, self.hier, linker_length=20.0, radii=(3.5, 0.0, 0.0),
                 linker_width=FPS_WIDTH, simulation_grid_resolution=0.5)
        # half the linker width (what the search inflates obstacles by) plus
        # half a grid step of slack
        self.assertAlmostEqual(av.get_effective_allowed_sphere_radius(),
                               0.5 * FPS_WIDTH + 0.5 * 0.5, places=6)

    def test_derived_clearance_floor(self):
        """A narrow linker keeps the historical 1.5 A floor."""
        av = _av(self.mdl, self.hier, linker_length=20.0, radii=(3.5, 0.0, 0.0),
                 linker_width=0.5, simulation_grid_resolution=0.5)
        self.assertAlmostEqual(av.get_effective_allowed_sphere_radius(), 1.5,
                               places=6)

    def test_explicit_clearance_is_obeyed(self):
        """An explicit value is never overridden -- including a small one."""
        av = _av(self.mdl, self.hier, linker_length=20.0, radii=(3.5, 0.0, 0.0),
                 linker_width=FPS_WIDTH, allowed_sphere_radius=1.0,
                 simulation_grid_resolution=0.5)
        self.assertAlmostEqual(av.get_effective_allowed_sphere_radius(), 1.0,
                               places=6)

    def test_fps_linker_width_is_not_empty_at_the_default_grid(self):
        """The regression `imp_bff av-export` hit: 4.5 A gave no points."""
        av = _av(self.mdl, self.hier, linker_length=20.0, radii=(3.5, 0.0, 0.0),
                 linker_width=FPS_WIDTH, simulation_grid_resolution=1.5)
        av.resample()
        self.assertGreater(_n_points(av), 0)

    def test_fine_grid_at_the_fps_width_is_still_empty(self):
        """Not a passing grade -- a pin on a **known** remaining gap.

        The derivation clears the *inflation* (half the linker width) but not
        the source atom's **own** radius, and `IMP.bff` keeps the attachment
        atom in the obstacle set where FPS drops it
        (`av_routines.cpp:50`). At a fine grid the derived clearance therefore
        still leaves the source walled in by its own atom: measured on T4L
        chain A, four of the five sites 132/55/19/86 and 99 come back empty at
        `linker_width=4.5, grid=0.5`, and all four are non-empty under
        `r_source + width/2 + grid/2` (PRD-121, open item).

        Adding `r_source` would also grow volumes that already compute -- T4L
        132 at width 0.5 goes 66 805 -> 85 456 points -- so it is a change to a
        physics default, not a bug fix, and it is the owner's call. This test
        records today's answer so that changing it is deliberate.
        """
        av = _av(self.mdl, self.hier, linker_length=20.0, radii=(3.5, 0.0, 0.0),
                 linker_width=FPS_WIDTH, simulation_grid_resolution=0.5)
        av.resample()
        self.assertEqual(_n_points(av), 0)

    def test_both_doors_agree(self):
        """The decorator door and the fps.json door build the same volume.

        The derivation used to exist on one side only, which is the whole
        defect: two doors onto one volume, disagreeing about a parameter
        neither caller passed.
        """
        position = {
            "chain_identifier": "A", "residue_seq_number": 132,
            "atom_name": "CB", "linker_length": 20.0,
            "linker_width": FPS_WIDTH, "radius1": 3.5,
            "simulation_grid_resolution": 1.5,
        }
        through_json = IMP.bff.get_av_from_structure(
            PDB, json.dumps(position), 1.5)
        av = _av(self.mdl, self.hier, linker_length=20.0,
                 radii=(3.5, 0.0, 0.0), linker_width=FPS_WIDTH,
                 simulation_grid_resolution=1.5)
        av.resample()
        self.assertGreater(_n_points(av), 0)
        # Same clearance is the claim; the fps.json door additionally strips
        # the attachment residue's side chain, so the volumes are not identical
        # and only the parameter is compared.
        self.assertAlmostEqual(av.get_effective_allowed_sphere_radius(),
                               0.5 * FPS_WIDTH + 0.5 * 1.5, places=6)
        self.assertGreater(len(through_json.get_points()), 0)

    def test_empty_volume_reports_its_anchor_not_half_of_it(self):
        """An empty volume's mean is the source atom, not source/2.

        `sum` started at 1.0 while the source contributed 1.0 more, so the
        numerator held the source once and the denominator counted it twice.
        """
        av = _av(self.mdl, self.hier, linker_length=20.0, radii=(3.5, 0.0, 0.0),
                 linker_width=FPS_WIDTH, allowed_sphere_radius=1.0,
                 simulation_grid_resolution=0.5)
        av.resample()
        self.assertEqual(_n_points(av), 0)
        np.testing.assert_allclose(av.get_mean_position(),
                                   av.get_source_coordinates(), atol=1e-6)
        # and with the source excluded there is nothing to average at all
        self.assertTrue(np.all(np.isnan(np.asarray(
            av.get_mean_position(False)))))

    def test_mean_position_is_not_biased_toward_the_origin(self):
        """Including the source must not count it twice in the denominator.

        Checked against the arithmetic rather than a pin: the mean including
        the source is the weighted cloud mean and the source, at their true
        weights.
        """
        av = _av(self.mdl, self.hier, linker_length=20.0, radii=(3.5, 0.0, 0.0),
                 linker_width=0.5, simulation_grid_resolution=0.5)
        av.resample()
        pts = np.asarray(av.get_map().get_xyz_density(), dtype=float)
        self.assertGreater(len(pts), 0)
        xyz, w = pts[:, :3], pts[:, 3]
        cloud = (xyz * w[:, None]).sum(axis=0) / w.sum()
        source = np.asarray(av.get_source_coordinates(), dtype=float)
        expected = (cloud * w.sum() + source) / (w.sum() + 1.0)
        np.testing.assert_allclose(av.get_mean_position(), expected, atol=1e-4)
        np.testing.assert_allclose(av.get_mean_position(False), cloud,
                                   atol=1e-4)

    def test_av_export_writes_a_volume_at_the_fps_linker_width(self):
        """The user-facing symptom: an empty .xyz written with exit status 0.

        `imp_bff av-export --linker-width 4.5` used to print a mean position
        of exactly half the CB coordinate, write a file whose first line was
        `0`, and exit 0.
        """
        program = BIN / "imp_bff"
        if not program.exists():
            self.skipTest("bin/imp_bff not present")
        out = self.get_tmp_file_name("av_w45.xyz")
        subprocess.check_call(
            [sys.executable, str(program), "av-export", "-p", PDB, "-c", "A",
             "-r", "132", "-a", "CB", "--linker-width", str(FPS_WIDTH),
             "-o", out])
        with open(out) as fh:
            self.assertGreater(int(fh.readline().strip()), 0)


if __name__ == "__main__":
    IMP.test.main()
