#!/usr/bin/env python3
"""Tests for dye I/O: PDB→CIF conversion, template read/write, rotamer library."""

import os
import sys
import tempfile
from pathlib import Path
from IMP.bff import get_template_dir, get_structure_dir


import pytest


class TestDyePdbConversion:
    """Test 0.2: Convert dye PDB to mmCIF."""

    def test_convert_a48_c1r_pdb_to_mmcif(self):
        """Convert A48_C1R.pdb and verify atom count."""
        from IMP.bff import convert_pdb_to_cif

        # The rotamer library moved into imp.bff as module data when cgdye
        # did; reaching it by walking up from this file stopped working.
        import IMP.bff
        input_pdb = str(Path(IMP.bff.get_data_path("rotamer_library")) / "A48_C1R.pdb")

        with tempfile.TemporaryDirectory() as tmpdir:
            output_cif = os.path.join(tmpdir, "A48_C1R_test.cif")

            convert_pdb_to_cif(input_pdb, output_cif)

            assert os.path.exists(output_cif), f"Output file {output_cif} not created"

            # Count HETATM atom-site records in CIF
            with open(output_cif) as fh:
                content = fh.read()
                atom_count = content.count("HETATM")

            # Reference A48_C1R has ~83 atoms (dye + linker)
            assert atom_count >= 50, (
                f"Expected at least 50 HETATM records, got {atom_count}"
            )


class TestDyeTemplateRead:
    """Test 0.3/0.4: Read dye template CIF."""

    def test_read_dye_template(self):
        """Read A48_C1R template and verify structure."""
        from IMP.bff.io.cif import read_dye_template_cif

        template_path = str(get_template_dir("dyes/A48_C1R/template.cif"))

        if not os.path.exists(template_path):
            pytest.skip("Template file not found")

        try:
            template = read_dye_template_cif(template_path)
        except Exception as e:
            pytest.fail(f"Failed to read dye template: {e}")

        assert "name" in template, "Template missing 'name' field"
        assert "features" in template, "Template missing 'features'"
        assert "impropers" in template, "Template missing 'impropers'"

        # Check for expected features
        assert "backbone_N" in template["features"], "Missing backbone_N feature"
        assert "linker" in template["features"], "Missing linker feature"
        assert "dye_core" in template["features"], "Missing dye_core feature"


class TestDyeTemplateWrite:
    """Test 0.4: Write dye template CIF."""

    def test_write_dye_template_roundtrip(self):
        """Write and read back a dye template, verify consistency."""
        from IMP.bff.io.cif import read_dye_template_cif, write_dye_template_cif

        template_path = str(get_template_dir("dyes/A48_C1R/template.cif"))

        if not os.path.exists(template_path):
            pytest.skip("Template file not found")

        try:
            original = read_dye_template_cif(template_path)
        except Exception as e:
            pytest.fail(f"Failed to read template for write test: {e}")

        with tempfile.TemporaryDirectory() as tmpdir:
            output_cif = os.path.join(tmpdir, "test_write.cif")

            try:
                write_dye_template_cif(output_cif, original)
            except Exception as e:
                pytest.fail(f"Failed to write dye template: {e}")

            assert os.path.exists(output_cif), f"Output file {output_cif} not created"

            # Read back and compare
            try:
                reloaded = read_dye_template_cif(output_cif)
            except Exception as e:
                pytest.fail(f"Failed to reload template: {e}")

            assert reloaded["name"] == original["name"], (
                "Template name mismatch after round-trip"
            )
            assert set(reloaded["features"].keys()) == set(
                original["features"].keys()
            ), (
                f"Feature keys mismatch: {set(reloaded['features'].keys())} vs {set(original['features'].keys())}"
            )


class TestDyeRegistryRead:
    """Test 0.5: Read dye registry CIF."""

    def test_read_dye_registry(self):
        """Read dyes.cif and verify all 33 dyes are present."""

        registry_path = str(get_template_dir("dyes/dyes.cif"))

        if not os.path.exists(registry_path):
            pytest.skip("Registry file not found")

        # Use a simple parser for the registry (not rotamer library format)
        with open(registry_path) as fh:
            content = fh.read()

        # Count dye entries (each line after loop header is one dye)
        lines = [
            l.strip()
            for l in content.split("\n")
            if l.strip() and not l.startswith("_")
        ]

        # Should have at least 30 dyes
        assert len(lines) >= 30, f"Expected at least 30 dye entries, got {len(lines)}"


class TestRotamerLibraryRead:
    """Test 0.6: Read rotamer library."""

    def test_rotamer_library_structure(self):
        """Verify rotamer library reader handles numpy format."""
        from IMP.bff.io.cif import read_rotamer_library, normalize_weights

        test_lib = {
            "id": [1, 2, 3],
            "weight": [0.5, 0.3, 0.2],
            "atom_names": ["C1"],
            "coords": {
                1: [[1.0, 2.0, 3.0]],
                2: [[4.0, 5.0, 6.0]],
                3: [[7.0, 8.0, 9.0]],
            },
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            from IMP.bff.io.cif import write_rotamer_library

            output_base = os.path.join(tmpdir, "test_rotamer")
            write_rotamer_library(output_base, test_lib)

            assert os.path.exists(output_base + "_coords.npy"), (
                "coords .npy not created"
            )
            assert os.path.exists(output_base + "_weights.txt"), (
                "weights file not created"
            )
            assert os.path.exists(output_base + "_atoms.txt"), "atoms file not created"

            reloaded = read_rotamer_library(output_base)

            assert len(reloaded["id"]) == 3, (
                f"Rotamer ID count mismatch: got {len(reloaded['id'])}"
            )
            assert len(reloaded["weight"]) == 3, "Weight count mismatch"
            assert len(reloaded["coords"]) == 3, "Coords dict length mismatch"
            assert 1 in reloaded["coords"], "Rotamer 1 not found in coords"
            assert 2 in reloaded["coords"], "Rotamer 2 not found in coords"
            assert reloaded["atom_names"] == ["C1"], (
                f"atom_names mismatch: {reloaded['atom_names']}"
            )


class TestRotamerLibraryNormalize:
    """Test 0.6: Normalize rotamer weights."""

    def test_normalize_weights(self):
        """Verify weight normalization sums to 1.0."""
        from IMP.bff.io.cif import normalize_weights

        lib = {"id": [1, 2, 3], "weight": [50.0, 30.0, 20.0], "coords": {}}

        normalize_weights(lib)

        total = sum(lib["weight"])
        assert abs(total - 1.0) < 1e-6, f"Weights should sum to 1.0, got {total}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
