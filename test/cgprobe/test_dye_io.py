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

        # The rotamer library moved into imp.bff as module data when cgprobe
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
        from IMP.bff import read_probe_template_cif

        template_path = str(get_template_dir("dyes/A48_C1R/template.cif"))

        if not os.path.exists(template_path):
            pytest.skip("Template file not found")

        try:
            template = read_probe_template_cif(template_path)
        except Exception as e:
            pytest.fail(f"Failed to read dye template: {e}")

        assert template.name, "Template missing 'name' field"
        assert len(template.features) > 0, "Template missing 'features'"
        assert len(template.impropers) > 0, "Template missing 'impropers'"

        # Check for expected features
        assert "backbone_N" in template.features, "Missing backbone_N feature"
        assert "linker" in template.features, "Missing linker feature"
        assert "dye_core" in template.features, "Missing dye_core feature"


class TestProbeTemplateRead:
    """The shipped probe templates parse into the fields a caller reads."""

    def test_shipped_probe_template_carries_its_metadata(self):
        from IMP.bff import read_probe_template_cif

        template_path = str(get_template_dir("dyes/A48_C1R/template.cif"))
        if not os.path.exists(template_path):
            pytest.skip("Template file not found")

        template = read_probe_template_cif(template_path)
        assert template.name
        assert template.features, "no features parsed"
        # A feature is a degree of freedom or a named region; the reader gives
        # an unlabelled one "dof", so an empty type means the column was
        # misread rather than absent.
        for feature in template.features.values():
            assert feature.feature_type
        # the probe metadata is what separates this reader from the plain one
        assert template.center_atom or template.dipole_atom_1, (
            "probe metadata not read"
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
        """The numpy/text reader, against files numpy and the shell wrote."""
        import numpy as np
        from IMP.bff import read_rotamer_library

        coords = np.array([[[1.0, 2.0, 3.0]],
                           [[4.0, 5.0, 6.0]],
                           [[7.0, 8.0, 9.0]]])

        with tempfile.TemporaryDirectory() as tmpdir:
            base = os.path.join(tmpdir, "test_rotamer")
            np.save(base + "_coords.npy", coords)
            with open(base + "_weights.txt", "w") as f:
                f.write("0.5\n0.3\n0.2\n")
            with open(base + "_atoms.txt", "w") as f:
                f.write("C1\n")

            reloaded = read_rotamer_library(base)

            assert len(reloaded.weights) == 3, "Weight count mismatch"
            assert reloaded.n_rotamers == 3, "Rotamer count mismatch"
            assert list(reloaded.atom_names) == ["C1"], (
                f"atom_names mismatch: {list(reloaded.atom_names)}"
            )
            # `coords` is `(n_rotamers, n_atoms, 3)`; rotamer 2 is the second
            # slab, one atom wide here
            assert list(reloaded.coords[1][0]) == [4.0, 5.0, 6.0], (
                "Rotamer 2 not found in coords"
            )


class TestRotamerLibraryNormalize:
    """Test 0.6: Normalize rotamer weights."""

    def test_normalize_weights(self):
        """Verify weight normalization sums to 1.0."""
        from IMP.bff import normalize_weights
        from IMP.bff import RotamerLibrary

        lib = RotamerLibrary()
        lib.weights = [50.0, 30.0, 20.0]
        lib.n_rotamers = 3

        normalized = normalize_weights(lib)

        total = sum(normalized.weights)
        assert abs(total - 1.0) < 1e-6, f"Weights should sum to 1.0, got {total}"
        # value semantics: the original is untouched
        assert sum(lib.weights) == 100.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
