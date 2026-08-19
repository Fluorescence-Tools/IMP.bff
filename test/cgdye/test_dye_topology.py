#!/usr/bin/env python3
"""Tests for Phase 2: dye topology, LJ scoring, exclusion lists."""

import math
import os
import sys
import tempfile
from pathlib import Path
from IMP.bff.tools import get_template_dir, get_structure_dir

import pytest


from IMP.bff.cgdye.topology import (
    build_angles,
    build_dihedrals,
    build_dye_topology,
    build_graph,
    find_cycles,
)
from IMP.bff.io.cif import as_forcefield_system, read_dye_forcefield_cif, write_dye_forcefield_cif
# CHARMM36_LJ, lj_cross and lj_params are scoring's, and were reached through
# `cgdye.topology` only because it imports them for its own use. Importing them
# from their owner is what lets that pass-through go.
from IMP.bff.scoring import CHARMM36_LJ, lj_cross, lj_params
from IMP.bff.scoring import (
    build_lj_type_table,
    compute_lj_pair_sites,
    lj_cross_params,
    lj_score,
    site_element_map,
)


class TestDyeTopologyBuilder:
    """Test 2.1: Dye topology generation."""

    def _make_atoms(self):
        return {
            1: {"serial": 1, "atom_name": "S1", "element": "S", "x": 0, "y": 0, "z": 0},
            2: {
                "serial": 2,
                "atom_name": "C1",
                "element": "C",
                "x": 1.8,
                "y": 0,
                "z": 0,
            },
            3: {
                "serial": 3,
                "atom_name": "C2",
                "element": "C",
                "x": 2.8,
                "y": 0,
                "z": 0,
            },
            4: {
                "serial": 4,
                "atom_name": "N1",
                "element": "N",
                "x": 3.8,
                "y": 0,
                "z": 0,
            },
        }

    def _make_bonds(self):
        return {(1, 2), (2, 3), (3, 4)}

    def _make_template(self):
        return {
            "name": "test",
            "features": {},
            "impropers": [],
        }

    def test_graph_construction(self):
        g = build_graph(self._make_bonds())
        assert 1 in g[2]
        assert 3 in g[2]
        assert 2 in g[3]

    def test_angles_from_chain(self):
        g = build_graph(self._make_bonds())
        angles = build_angles(g)
        assert (1, 2, 3) in angles
        assert (2, 3, 4) in angles

    def test_dihedrals_from_chain(self):
        g = build_graph(self._make_bonds())
        dihedrals = build_dihedrals(g)
        assert (1, 2, 3, 4) in dihedrals

    def test_build_dye_topology_bonds(self):
        atoms = self._make_atoms()
        bonds = self._make_bonds()
        template = self._make_template()
        topo = build_dye_topology(atoms, bonds, template)
        assert len(topo["bonds"]) == 3
        for a, b, d in topo["bonds"]:
            assert d > 0

    def test_build_dye_topology_angles(self):
        atoms = self._make_atoms()
        bonds = self._make_bonds()
        template = self._make_template()
        topo = build_dye_topology(atoms, bonds, template)
        assert len(topo["angles"]) == 2

    def test_build_dye_topology_dihedrals(self):
        atoms = self._make_atoms()
        bonds = self._make_bonds()
        template = self._make_template()
        topo = build_dye_topology(atoms, bonds, template)
        assert len(topo["dihedrals"]) == 1


class TestLJParameters:
    """Test 2.2: CHARMM36 LJ parameter table and cross terms."""

    def test_all_elements_present(self):
        for elem in ("C", "N", "O", "S", "H"):
            assert elem in CHARMM36_LJ
            assert "rmin_half" in CHARMM36_LJ[elem]
            assert "epsilon" in CHARMM36_LJ[elem]

    def test_rmin_positive(self):
        for elem, p in CHARMM36_LJ.items():
            assert p["rmin_half"] > 0, f"{elem} rmin_half must be positive"

    def test_epsilon_negative(self):
        for elem, p in CHARMM36_LJ.items():
            assert p["epsilon"] < 0, f"{elem} epsilon must be negative (well depth)"

    def test_lj_cross_combining_rules(self):
        rmin, eps = lj_cross("C", "N")
        assert (
            abs(rmin - (CHARMM36_LJ["C"]["rmin_half"] + CHARMM36_LJ["N"]["rmin_half"]))
            < 1e-10
        )
        expected_eps = math.sqrt(
            CHARMM36_LJ["C"]["epsilon"] * CHARMM36_LJ["N"]["epsilon"]
        )
        assert abs(eps - expected_eps) < 1e-10

    def test_lj_cross_symmetric(self):
        r1, e1 = lj_cross("C", "S")
        r2, e2 = lj_cross("S", "C")
        assert abs(r1 - r2) < 1e-10
        assert abs(e1 - e2) < 1e-10

    def test_lj_cross_unknown_element_falls_back_to_carbon(self):
        p = lj_params("X")
        assert p["rmin_half"] == CHARMM36_LJ["C"]["rmin_half"]


class TestLJTypeCifRoundtrip:
    """Test 2.2: LJ types read/write in system CIF."""

    def test_lj_types_roundtrip(self):
        system = {
            "name": "test_lj",
            "components": {},
            "sites": [],
            "groups": {},
            "rb_groups": {},
            "md_fixed_groups": {},
            "fixed_groups": [],
            "bond_types": {},
            "angle_types": {},
            "torsion_types": {},
            "improper_types": {},
            "bonds": [],
            "angles": [],
            "dihedrals": [],
            "impropers": [],
            "lj_types": {
                "LJ_C": {"element": "C", "rmin_half": 2.024, "epsilon": -0.064},
                "LJ_N": {"element": "N", "rmin_half": 1.893, "epsilon": -0.154},
            },
        }
        with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as f:
            path = f.name
        try:
            write_dye_forcefield_cif(path, system)
            loaded = read_dye_forcefield_cif(path)
            assert "LJ_C" in loaded.lj_types
            assert loaded.lj_types["LJ_C"].element == "C"
            assert abs(loaded.lj_types["LJ_C"].rmin_half - 2.024) < 1e-3
            assert "LJ_N" in loaded.lj_types
        finally:
            os.unlink(path)

    def test_build_lj_type_table(self):
        table = build_lj_type_table({"C", "S", "N"})
        assert len(table) == 3
        assert table["LJ_C"]["element"] == "C"
        assert table["LJ_S"]["element"] == "S"


class TestExclusionList:
    """Test 2.6: 1-2/1-3/1-4 exclusion list generation."""

    def test_bonded_pairs_excluded(self):
        system = {
            "bonds": [("A/1", "A/2", 1.5, "B1"), ("A/2", "A/3", 1.5, "B1")],
            "angles": [("A/1", "A/2", "A/3", 2.0, "A1")],
            "dihedrals": [],
            "impropers": [],
        }
        excluded = {frozenset(p) for p in as_forcefield_system(system).exclusions()}
        assert frozenset({"A/1", "A/2"}) in excluded
        assert frozenset({"A/2", "A/3"}) in excluded

    def test_angle_ends_excluded(self):
        system = {
            "bonds": [("A/1", "A/2", 1.5, "B1"), ("A/2", "A/3", 1.5, "B1")],
            "angles": [("A/1", "A/2", "A/3", 2.0, "A1")],
            "dihedrals": [],
            "impropers": [],
        }
        excluded = {frozenset(p) for p in as_forcefield_system(system).exclusions()}
        assert frozenset({"A/1", "A/3"}) in excluded

    def test_dihedral_ends_excluded(self):
        system = {
            "bonds": [],
            "angles": [],
            "dihedrals": [("A/1", "A/2", "A/3", "A/4", "T_PI")],
            "impropers": [],
        }
        excluded = {frozenset(p) for p in as_forcefield_system(system).exclusions()}
        assert frozenset({"A/1", "A/4"}) in excluded

    def test_nonbonded_pair_not_excluded(self):
        system = {
            "bonds": [("A/1", "A/2", 1.5, "B1")],
            "angles": [],
            "dihedrals": [],
            "impropers": [],
        }
        excluded = {frozenset(p) for p in as_forcefield_system(system).exclusions()}
        assert frozenset({"A/1", "A/3"}) not in excluded


class TestLJScoring:
    """Test 2.5: LJ pair score computation."""

    def test_lj_score_zero_at_rmin(self):
        e = lj_score(r=2.0, rmin=2.0, epsilon=-1.0)
        assert abs(e) < 1e-10

    def test_lj_score_zero_beyond_rmin(self):
        e = lj_score(r=3.0, rmin=2.0, epsilon=-1.0)
        assert e == 0.0

    def test_lj_score_repulsive_close(self):
        e = lj_score(r=1.0, rmin=2.0, epsilon=1.0)
        assert e > 0

    def test_lj_score_zero_beyond_rmin(self):
        e = lj_score(r=3.0, rmin=2.0, epsilon=1.0)
        assert e == 0.0

    def test_lj_cross_params_callable(self):
        rmin, eps = lj_cross_params("C", "N")
        assert rmin > 0
        assert eps > 0

    def test_compute_lj_pair_sites(self):
        system = {
            "sites": [
                {"id": "dye/C1", "atom_name": "C1", "component": "dye"},
                {"id": "dye/N1", "atom_name": "N1", "component": "dye"},
            ],
            "bonds": [],
            "angles": [],
            "dihedrals": [],
            "impropers": [],
        }
        pairs = compute_lj_pair_sites(system)
        assert len(pairs) == 1
        sa, sb, rmin, eps = pairs[0]
        assert rmin > 0
        assert eps > 0

    def test_site_element_map(self):
        system = {
            "sites": [
                {"id": "dye/C1", "atom_name": "C1"},
                {"id": "dye/S1", "atom_name": "S1"},
                {"id": "dye/N1", "atom_name": "N1"},
            ]
        }
        emap = site_element_map(system)
        assert emap["dye/C1"] == "C"
        assert emap["dye/S1"] == "S"
        assert emap["dye/N1"] == "N"



def _assert_close(a, b, places=9):
    assert abs(a - b) < 10 ** (-places), (a, b)


class TestTorsionConvention:
    """torsion_types are CHARMM V = k(1+cos(n phi - delta)); IMP.core.Cosine has the opposite sign."""

    def _score_at(self, ttype, phi_deg):
        import math
        import IMP
        import IMP.algebra
        import IMP.core
        from IMP.bff.scoring import torsion_cosine
        m = IMP.Model()
        def P(x):
            p = IMP.Particle(m)
            IMP.core.XYZ.setup_particle(p, IMP.algebra.Vector3D(*x))
            return p
        phi = math.radians(phi_deg)
        ps = [P((-1.5, 1, 0)), P((0, 0, 0)), P((1.5, 0, 0)), P((3.0, math.cos(phi), math.sin(phi)))]
        r = IMP.core.DihedralRestraint(m, torsion_cosine(ttype), *ps)
        return r.evaluate(False)

    def test_pi_torsion_is_planar_and_linker_staggered(self):
        import math
        t_pi = {"periodicity": 2, "phase_rad": math.pi, "k": 12.0}
        t_link = {"periodicity": 3, "phase_rad": 0.0, "k": 1.5}
        _assert_close(self._score_at(t_pi, 0.0), 0.0, places=9)
        _assert_close(self._score_at(t_pi, 180.0), 0.0, places=9)
        _assert_close(self._score_at(t_pi, 90.0), 24.0, places=9)
        _assert_close(self._score_at(t_link, 60.0), 0.0, places=9)
        _assert_close(self._score_at(t_link, 180.0), 0.0, places=9)
        _assert_close(self._score_at(t_link, 0.0), 3.0, places=9)

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
