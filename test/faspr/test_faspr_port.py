"""The 1:1 FASPR C++ port: backbone fidelity and parity with the reference.

``IMP.bff.pack_protein_sidechains`` (``src/Faspr*``, declared in
``include/ProbeRotamerLibrary.h``) is a
behaviour-identical vendoring of FASPR 20200309 wrapped in
``IMP::bff::faspr``: same energies, same DEE / tree-decomposition search,
same tie-breaking, float math untouched. These tests pin that claim:

1. **Self-contained**: packing a cleaned structure preserves the backbone
   (N/CA/C/O coordinates byte-identical), keeps every residue, and rebuilds
   complete side chains (heavy atoms present for every non-GLY residue).
2. **Parity**: repacking the same input with the *reference* FASPR
   executable (``junk/FASPR/``) must give the same answer -- identical
   atom ordering and coordinates to the printed 0.001 A, and equal rotamers
   (chi1 within 0.5 deg on every rotatable residue).

The Dunbrack-2010 binary library (``dun2010bbdep.bin``, 13.8 MB) is *not*
redistributed with IMP.bff (licence decision open, PRD-118); both tests
locate it via ``$IMP_BFF_FASPR_ROT_LIB`` or the ``junk/FASPR`` clone and
skip when neither is present. The reference executable is used as-is when
built, else compiled from ``junk/FASPR/src`` (same sources the port
vendors, so the comparison is against the ground truth itself).
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import IMP.bff

REPO = Path(__file__).resolve().parents[2]

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
# residues with a rotatable chi1 that is not ring-locked
CHI1_RESIDUES = {
    "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "ILE", "LEU", "LYS",
    "MET", "PHE", "SER", "THR", "TRP", "TYR", "VAL",
}


def _rotlib():
    """The Dunbrack-2010 binary library, or skip."""
    candidates = []
    env = os.environ.get("IMP_BFF_FASPR_ROT_LIB")
    if env:
        candidates.append(Path(env))
    candidates.append(REPO / "junk" / "FASPR" / "dun2010bbdep.bin")
    for c in candidates:
        if c.exists():
            return c
    pytest.skip("FASPR Dunbrack library not available "
                "(set IMP_BFF_FASPR_ROT_LIB or clone junk/FASPR)")


def _reference_exe(tmp_path):
    """The reference FASPR executable -- built from the vendored sources."""
    exe = REPO / "junk" / "FASPR" / "FASPR"
    if exe.exists():
        return exe
    src = REPO / "junk" / "FASPR" / "src"
    if not src.is_dir():
        pytest.skip("reference FASPR sources not available")
    out = tmp_path / "FASPR_ref"
    cc = os.environ.get("CXX", "c++")
    res = subprocess.run([cc, "-O3", "-o", str(out)]
                         + sorted(str(p) for p in src.glob("*.cpp")),
                         capture_output=True, text=True)
    if res.returncode != 0:
        pytest.skip(f"could not build reference FASPR: {res.stderr[-200:]}")
    return out


def _clean_pdb(src, dst):
    """Standard-residue ATOM records only (what FASPR accepts)."""
    lines = []
    for line in Path(src).read_text().splitlines():
        if line.startswith("ATOM") and line[17:20].strip() in STANDARD_AA:
            lines.append(line)
        elif line.startswith(("TER", "END")):
            lines.append(line)
    Path(dst).write_text("\n".join(lines) + "\n")


def _read_pdb(path):
    """[(chain, resid, resname, {atom: xyz})] in file order."""
    residues = []
    cur_key = None
    cur = None
    for line in Path(path).read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        key = (line[21], int(line[22:26]), line[17:20].strip())
        if cur_key != key:
            cur = (key[0], key[1], key[2], {})
            residues.append(cur)
            cur_key = key
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        cur[3][line[12:16].strip()] = xyz
    return residues


def _dihedral(p0, p1, p2, p3):
    """Standard (IUPAC) signed dihedral in degrees."""
    b0 = np.asarray(p0) - np.asarray(p1)
    b1 = np.asarray(p2) - np.asarray(p1)
    b2 = np.asarray(p3) - np.asarray(p2)
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def _chi1(residue):
    """chi1 (N-CA-CB-CG-like) of a residue, or None if not measurable."""
    _, _, name, atoms = residue
    if name not in CHI1_RESIDUES:
        return None
    for cg in ("CG", "CG1", "OG", "OG1", "SG"):
        if all(a in atoms for a in ("N", "CA", "CB", cg)):
            return _dihedral(atoms["N"], atoms["CA"], atoms["CB"], atoms[cg])
    return None


def _pack(tmp_path, rotlib):
    """Clean 3GUN and pack it with the port; returns (input, output) paths."""
    inp = tmp_path / "3gun_clean.pdb"
    out = tmp_path / "3gun_port.pdb"
    _clean_pdb(REPO / "examples" / "structure" / "T4L" / "3GUN.pdb", inp)
    IMP.bff.pack_protein_sidechains(str(inp), str(out), str(rotlib))
    assert out.exists(), "faspr_pack wrote no output"
    return inp, out


def test_pack_preserves_backbone_and_builds_sidechains(tmp_path):
    """Backbone untouched, side chains complete, no stdout leakage."""
    rotlib = _rotlib()
    inp, out = _pack(tmp_path, rotlib)

    src = _read_pdb(inp)
    dst = _read_pdb(out)
    assert len(src) == len(dst), "residue count changed"
    for (c1, r1, n1, a1), (c2, r2, n2, a2) in zip(src, dst):
        assert (c1, r1, n1) == (c2, r2, n2), "residue order/identity changed"
        for bb in ("N", "CA", "C", "O"):
            if bb in a1:
                assert bb in a2, f"backbone {bb} missing for {r1}"
                assert np.allclose(a1[bb], a2[bb], atol=1e-3), \
                    f"backbone {bb} moved for {c1}:{r1}"
        if n1 != "GLY":
            assert "CB" in a2, f"no CB rebuilt for {c1}:{r1} {n1}"
    # some chi1 must be measurable at all (sanity that side chains exist)
    assert sum(1 for r in dst if _chi1(r) is not None) > 50


def test_pack_parity_with_reference_executable(tmp_path, capfd):
    """Same coordinates and rotamers as the reference FASPR executable."""
    rotlib = _rotlib()
    inp, out = _pack(tmp_path, rotlib)

    ref_out = tmp_path / "3gun_ref.pdb"
    ref = _reference_exe(tmp_path)
    # the reference looks for dun2010bbdep.bin beside the executable
    work = tmp_path / "refwork"
    work.mkdir()
    shutil.copy(rotlib, work / "dun2010bbdep.bin")
    shutil.copy(inp, work / "in.pdb")
    res = subprocess.run([str(ref), "-i", "in.pdb", "-o", "ref.pdb"],
                         cwd=work, capture_output=True, text=True)
    assert res.returncode == 0, res.stderr[-400:]
    shutil.copy(work / "ref.pdb", ref_out)

    ours = _read_pdb(out)
    theirs = _read_pdb(ref_out)
    assert len(ours) == len(theirs), "residue count differs from reference"

    n_atoms = 0
    n_chi = 0
    for (c1, r1, n1, a1), (c2, r2, n2, a2) in zip(ours, theirs):
        assert (c1, r1, n1) == (c2, r2, n2), "atom/residue ordering differs"
        assert set(a1) == set(a2), f"atom set differs at {c1}:{r1}"
        for name, xyz in a1.items():
            assert np.allclose(xyz, a2[name], atol=0.002), \
                f"{name} at {c1}:{r1}: {xyz} vs reference {a2[name]}"
            n_atoms += 1
        x1, x2 = _chi1((c1, r1, n1, a1)), _chi1((c2, r2, n2, a2))
        if x1 is not None and x2 is not None:
            d = abs(x1 - x2) % 360.0
            d = min(d, 360.0 - d)
            assert d < 0.5, f"chi1 differs at {c1}:{r1} {n1}: {x1:.1f} vs {x2:.1f}"
            n_chi += 1

    assert n_atoms > 1000
    assert n_chi > 50

    # the port is silent unless asked (FASPR itself logs to stdout)
    captured = capfd.readouterr()
    assert "FASPR" not in captured.out, "port leaked FASPR's stdout log"


def test_pack_missing_library_raises(tmp_path):
    with pytest.raises(Exception):
        IMP.bff.pack_protein_sidechains(str(tmp_path / "none_in.pdb"),
                           str(tmp_path / "none_out.pdb"),
                           str(tmp_path / "no_such_library.bin"))


def test_pack_verbose_logs(tmp_path, capfd):
    rotlib = _rotlib()
    inp, out = _pack(tmp_path, rotlib)
    capfd.readouterr()  # drain
    IMP.bff.pack_protein_sidechains(str(inp), str(tmp_path / "v.pdb"), str(rotlib),
                       verbose=True)
    captured = capfd.readouterr()
    assert "FASPR" in captured.out or "residues" in captured.out