#!/usr/bin/env python3
"""Tests for Phase 1: backbone frame, attachment site, dye placement."""

import math
import os
import sys
import tempfile
from pathlib import Path
from IMP.bff import get_template_dir, get_structure_dir


import numpy as np
import pytest


class TestBackboneFrame:
    """Test 1.1: Backbone frame computation."""

    def test_frame_from_coords_identity(self):
        from IMP.bff import backbone_frame

        import IMP.algebra

        ca = IMP.algebra.Vector3D(10, 20, 30)
        n = IMP.algebra.Vector3D(11, 20, 30)
        c = IMP.algebra.Vector3D(10, 21, 30)
        frame = backbone_frame(ca, n, c)
        trans = frame          # `backbone_frame` answers the transformation itself
        assert abs(trans.get_translation()[0] - 10) < 1e-10
        assert abs(trans.get_translation()[1] - 20) < 1e-10
        assert abs(trans.get_translation()[2] - 30) < 1e-10

    def test_frame_axes_orthonormal(self):
        from IMP.bff import backbone_frame

        import IMP.algebra

        ca = IMP.algebra.Vector3D(5, 10, 15)
        n = IMP.algebra.Vector3D(5.5, 10.8, 15.3)
        c = IMP.algebra.Vector3D(4.8, 10.5, 15.7)
        frame = backbone_frame(ca, n, c)
        trans = frame          # `backbone_frame` answers the transformation itself
        rot = trans.get_rotation()

        x = rot.get_rotated(IMP.algebra.Vector3D(1, 0, 0))
        y = rot.get_rotated(IMP.algebra.Vector3D(0, 1, 0))
        z = rot.get_rotated(IMP.algebra.Vector3D(0, 0, 1))

        def mag(v):
            return math.sqrt(sum(v[i] ** 2 for i in range(3)))

        def dot(a, b):
            return sum(a[i] * b[i] for i in range(3))

        assert abs(mag(x) - 1.0) < 1e-10, f"|x| = {mag(x)}"
        assert abs(mag(y) - 1.0) < 1e-10, f"|y| = {mag(y)}"
        assert abs(mag(z) - 1.0) < 1e-10, f"|z| = {mag(z)}"
        assert abs(dot(x, y)) < 1e-10, f"x.y = {dot(x, y)}"
        assert abs(dot(x, z)) < 1e-10, f"x.z = {dot(x, z)}"
        assert abs(dot(y, z)) < 1e-10, f"y.z = {dot(y, z)}"

    def test_frame_matches_reference(self):
        """Verify that backbone frame matches reference transform."""
        from IMP.bff import backbone_frame

        import IMP.algebra

        ca_np = np.array([5.0, 10.0, 15.0])
        n_np = np.array([5.5, 10.8, 15.3])
        c_np = np.array([4.8, 10.5, 15.7])

        x_np = n_np - ca_np
        x_np /= np.linalg.norm(x_np)
        yt_np = c_np - ca_np
        yt_np /= np.linalg.norm(yt_np)
        z_np = np.cross(x_np, yt_np)
        z_np /= np.linalg.norm(z_np)
        y_np = np.cross(z_np, x_np)
        rot_np = np.vstack([x_np, y_np, z_np])

        frame = backbone_frame(
            IMP.algebra.Vector3D(*ca_np),
            IMP.algebra.Vector3D(*n_np),
            IMP.algebra.Vector3D(*c_np),
        )
        trans = frame          # `backbone_frame` answers the transformation itself
        rot = trans.get_rotation()

        for i in range(3):
            e = IMP.algebra.Vector3D(int(i == 0), int(i == 1), int(i == 2))
            col = rot.get_rotated(e)
            diff = np.sqrt(sum((col[j] - rot_np[i][j]) ** 2 for j in range(3)))
            assert diff < 1e-10, f"Column {i} diff = {diff}"

    def test_frame_x_along_ca_n(self):
        from IMP.bff import backbone_frame

        import IMP.algebra

        ca = IMP.algebra.Vector3D(0, 0, 0)
        n = IMP.algebra.Vector3D(2, 0, 0)
        c = IMP.algebra.Vector3D(0, 3, 0)
        frame = backbone_frame(ca, n, c)
        trans = frame          # `backbone_frame` answers the transformation itself
        rot = trans.get_rotation()
        x = rot.get_rotated(IMP.algebra.Vector3D(1, 0, 0))
        assert abs(x[0] - 1.0) < 1e-10
        assert abs(x[1]) < 1e-10
        assert abs(x[2]) < 1e-10


class TestAttachmentResolver:
    """Test 1.2: Resolve backbone atoms from hierarchy."""

    def test_resolve_site_from_pdb(self):
        import IMP
        import IMP.atom

        pdb_content = (
            "ATOM      1  N   ALA A   1       1.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA A   1       2.000   0.000   0.000  1.00  0.00           C\n"
            "ATOM      3  C   ALA A   1       2.000   1.000   0.000  1.00  0.00           C\n"
            "ATOM      4  O   ALA A   1       2.000   2.000   0.000  1.00  0.00           O\n"
            "END\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
            f.write(pdb_content)
            pdb_path = f.name

        try:
            m = IMP.Model()
            hier = IMP.atom.read_pdb(pdb_path, m, IMP.atom.AllPDBSelector())

            from IMP.bff import resolve_probe_site

            site = resolve_probe_site(hier, "A", 1)   # CA, N, C in that order
            assert len(site) == 3
            assert [IMP.atom.Atom(p).get_atom_type().get_string() for p in site] \
                == ["CA", "N", "C"]
        finally:
            os.unlink(pdb_path)

    def test_strip_sidechain_at_site(self):
        import IMP
        import IMP.atom

        pdb_content = (
            "ATOM      1  N   ALA A   1       1.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA A   1       2.000   0.000   0.000  1.00  0.00           C\n"
            "ATOM      3  C   ALA A   1       2.000   1.000   0.000  1.00  0.00           C\n"
            "ATOM      4  O   ALA A   1       2.000   2.000   0.000  1.00  0.00           O\n"
            "ATOM      5  CB  ALA A   1       3.000   0.000   0.000  1.00  0.00           C\n"
            "END\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
            f.write(pdb_content)
            pdb_path = f.name

        try:
            m = IMP.Model()
            hier = IMP.atom.read_pdb(pdb_path, m, IMP.atom.AllPDBSelector())

            from IMP.bff import strip_sidechain_at_site

            removed = strip_sidechain_at_site(hier, "A", 1)
            assert removed == 1

            names = {
                IMP.atom.Atom(a).get_atom_type().get_string()
                for a in IMP.atom.get_by_type(hier, IMP.atom.ATOM_TYPE)
            }
            assert "CB" not in names
            assert {"N", "CA", "C", "O"}.issubset(names)
        finally:
            os.unlink(pdb_path)


class TestDyePlacement:
    """Test 1.3: Transform dye coordinates to protein frame."""

    def test_place_dye_shifts_origin(self):
        import IMP
        import IMP.algebra
        import IMP.atom
        import IMP.core

        m = IMP.Model()
        root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(m, "root"))

        p = IMP.Particle(m)
        IMP.core.XYZ.setup_particle(p, IMP.algebra.Vector3D(1, 0, 0))
        at = IMP.atom.Atom.setup_particle(p, IMP.atom.AtomType("C"))
        IMP.atom.Hierarchy.setup_particle(p)
        root.add_child(IMP.atom.Hierarchy(p))

        ca = IMP.algebra.Vector3D(10, 20, 30)
        n = IMP.algebra.Vector3D(11, 20, 30)
        c = IMP.algebra.Vector3D(10, 21, 30)

        from IMP.bff import backbone_frame

        frame = backbone_frame(ca, n, c)
        from IMP.bff import place_probe

        place_probe(IMP.atom.Hierarchy(p), frame)

        pos = IMP.core.XYZ(p).get_coordinates()
        assert abs(pos[0] - 11) < 1e-6, f"Expected x~11, got {pos[0]}"
        assert abs(pos[1] - 20) < 1e-6, f"Expected y~20, got {pos[1]}"
        assert abs(pos[2] - 30) < 1e-6, f"Expected z~30, got {pos[2]}"

    def test_place_label_from_coords(self):
        import IMP
        import IMP.algebra
        import IMP.atom
        import IMP.core

        m = IMP.Model()
        p = IMP.Particle(m)
        IMP.core.XYZ.setup_particle(p, IMP.algebra.Vector3D(1, 0, 0))
        at = IMP.atom.Atom.setup_particle(p, IMP.atom.AtomType("C"))
        hier_p = IMP.atom.Hierarchy.setup_particle(p)

        from IMP.bff import place_probe_from_coords

        place_probe_from_coords(
            hier_p,
            IMP.algebra.Vector3D(10, 20, 30),
            IMP.algebra.Vector3D(11, 20, 30),
            IMP.algebra.Vector3D(10, 21, 30),
        )

        pos = IMP.core.XYZ(p).get_coordinates()
        assert abs(pos[0] - 11) < 1e-6
        assert abs(pos[1] - 20) < 1e-6
        assert abs(pos[2] - 30) < 1e-6


class TestMultiLabelAttachment:
    """Test 1.4: Attach multiple dyes."""

    def test_attach_two_dyes(self):
        import IMP
        import IMP.algebra
        import IMP.atom
        import IMP.core

        pdb_content = (
            "ATOM      1  N   ALA A   1       1.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA A   1       2.000   0.000   0.000  1.00  0.00           C\n"
            "ATOM      3  C   ALA A   1       2.000   1.000   0.000  1.00  0.00           C\n"
            "ATOM      4  N   ALA A   2      11.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      5  CA  ALA A   2      12.000   0.000   0.000  1.00  0.00           C\n"
            "ATOM      6  C   ALA A   2      12.000   1.000   0.000  1.00  0.00           C\n"
            "END\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
            f.write(pdb_content)
            pdb_path = f.name

        try:
            m = IMP.Model()
            prot = IMP.atom.read_pdb(pdb_path, m, IMP.atom.AllPDBSelector())

            def make_dye(model):
                p = IMP.Particle(model)
                IMP.core.XYZ.setup_particle(p, IMP.algebra.Vector3D(1, 0, 0))
                IMP.atom.Atom.setup_particle(p, IMP.atom.AtomType("C"))
                return IMP.atom.Hierarchy.setup_particle(p)

            dye1 = make_dye(m)
            dye2 = make_dye(m)

            from IMP.bff import ProbeAttachment, attach_probes

            results = attach_probes(prot, [ProbeAttachment(dye1, "A", 1), ProbeAttachment(dye2, "A", 2)])
            assert len(results) == 2
            assert results[0].residue == 1
            assert results[1].residue == 2

            pos1 = IMP.core.XYZ(dye1).get_coordinates()
            pos2 = IMP.core.XYZ(dye2).get_coordinates()
            assert abs(pos1[0] - 1) < 1e-6
            assert abs(pos2[0] - 11) < 1e-6
        finally:
            os.unlink(pdb_path)


class TestRotamerAttachment:
    """Test Trick 4: Attachment using IMP.rotamer."""

    def test_get_anchor_cb_position(self):
        import IMP
        import IMP.atom
        from IMP.bff import get_anchor_cb_position
        
        pdb_content = (
            "ATOM      1  N   ALA A   1       1.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA A   1       2.000   0.000   0.000  1.00  0.00           C\n"
            "ATOM      3  C   ALA A   1       2.000   1.000   0.000  1.00  0.00           C\n"
            "ATOM      4  CB  ALA A   1       3.000   0.000   0.000  1.00  0.00           C\n"
            "END\n"
        )
        import tempfile, os
        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
            f.write(pdb_content)
            pdb_path = f.name
            
        try:
            m = IMP.Model()
            prot = IMP.atom.read_pdb(pdb_path, m, IMP.atom.AllPDBSelector())
            
            # Since we don't have a rotamer library loaded in this simple test, 
            # we expect it to gracefully fall back to the crystallographic CB
            cb_pos = get_anchor_cb_position(prot, "A", 1)
            
            assert cb_pos is not None
            assert abs(cb_pos[0] - 3.0) < 1e-6
            assert abs(cb_pos[1] - 0.0) < 1e-6
            assert abs(cb_pos[2] - 0.0) < 1e-6
        finally:
            os.unlink(pdb_path)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
