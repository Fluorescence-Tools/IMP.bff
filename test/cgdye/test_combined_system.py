#!/usr/bin/env python3

import os
import sys
from pathlib import Path
from IMP.bff.tools.paths import get_template_dir, get_structure_dir

import IMP
import IMP.algebra
import IMP.atom
import IMP.core


from IMP.bff.cgdye.io.cif import read_dye_forcefield_cif, write_dye_forcefield_cif
from IMP.bff.cgdye.sim.dye_restraints import build_dye_restraints
from IMP.bff.cgdye.topology.combined import build_dye_protein_system


def test_build_combined_system_basic():
    system = build_dye_protein_system(
        str(get_structure_dir("cx4.mol2")),
        str(get_structure_dir("atto655.mol2")),
        "CX4",
        "atto655",
        protein_template=str(get_template_dir("cx4.template.cif")),
        dye_template=str(get_template_dir("atto655.template.cif")),
    )
    assert system["name"] == "CX4_atto655"
    assert "CX4" in system["components"]
    assert "atto655" in system["components"]
    assert len(system["sites"]) > 0
    assert len(system["bonds"]) > 0
    assert len(system["angles"]) > 0
    assert len(system["dihedrals"]) > 0
    assert "lj_types" in system
    assert "LJ_C" in system["lj_types"]


def test_combined_system_cif_roundtrip(tmp_path):
    system = build_dye_protein_system(
        str(get_structure_dir("cx4.mol2")),
        str(get_structure_dir("atto655.mol2")),
        "CX4",
        "atto655",
        protein_template=str(get_template_dir("cx4.template.cif")),
        dye_template=str(get_template_dir("atto655.template.cif")),
    )
    out = tmp_path / "combined.system.cif"
    write_dye_forcefield_cif(str(out), system)
    loaded = read_dye_forcefield_cif(str(out))
    assert loaded["name"] == "CX4_atto655"
    assert len(loaded["sites"]) == len(system["sites"])
    assert "LJ_C" in loaded.get("lj_types", {})


def test_build_dye_restraints_smoke():
    system = build_dye_protein_system(
        str(get_structure_dir("cx4.mol2")),
        str(get_structure_dir("atto655.mol2")),
        "CX4",
        "atto655",
        protein_template=str(get_template_dir("cx4.template.cif")),
        dye_template=str(get_template_dir("atto655.template.cif")),
    )

    model = IMP.Model()
    site_particles = {}
    for i, s in enumerate(system["sites"]):
        p = IMP.Particle(model)
        IMP.core.XYZ.setup_particle(
            p,
            IMP.algebra.Vector3D(float(i), 0.0, 0.0),
        )
        site_particles[s["id"]] = p

    restraints = build_dye_restraints(model, system, site_particles)
    assert len(restraints) > 0


# IMP runs every .py under test/ as a standalone script, and a file of bare
# pytest functions would import cleanly and exit 0 -- reporting success without
# running a single assertion. Hand the file to pytest explicitly so a failure
# here is a failure in ctest.
if __name__ == "__main__":
    import sys
    try:
        import pytest
    except ImportError:
        print("pytest not installed; skipping", __file__)
        sys.exit(0)
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
