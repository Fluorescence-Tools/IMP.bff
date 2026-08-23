"""ZMatrix: the general internal-coordinate kernel behind PRD-118.

``IMP.bff.ZMatrix`` (``include/ZMatrix.h``) is the IMP-native kernel the
FASPR port and the ``.drot`` dye-rotamer store share: the standard (IUPAC)
dihedral, FASPR's ``Internal2Cartesian`` placement in double precision, bond
perception, and a molecule's spanning-tree Z-matrix with exact encode/decode.
It is general over the molecule -- proven here on a dye+linker template (83
atoms) and a protein residue, not just side chains.

What is pinned:

1. **The convention** -- ``dihedral_deg`` equals an independent numpy
   implementation of the praxeolitic IUPAC form. This test exists because the
   Python prototype shipped standard+180 for a day: its dihedral round-tripped
   perfectly against its own placement primitive while mirroring every
   rotamer off the world's convention (caught only by an A/B against FASPR
   and mdtraj; okf/log.md 2026-08-22). A convention that round-trips against
   itself proves nothing.
2. **Exactness** -- place a random atom, measure back: length, angle and
   dihedral agree to ~1e-13; encode(decode(chi)) == chi for arbitrary chi.
3. **FASPR agreement** -- the chi the kernel reads off a FASPR-repacked
   structure equals the mdtraj-verified reference value from the A/B.
4. **Generality** -- the dye template's perceived bonds equal its MOL2 bond
   table (87), and the template round-trips to itself.
"""

from pathlib import Path

import numpy as np
import pytest

import IMP.algebra
import IMP.bff

REPO = Path(__file__).resolve().parents[2]

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


def _np_dih(p0, p1, p2, p3):
    """Independent praxeolitic (IUPAC) dihedral, degrees."""
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def _V(x):
    return IMP.algebra.Vector3D(*[float(t) for t in x])


def _angmod(x):
    return np.abs((x + 180.0) % 360.0 - 180.0)


def test_dihedral_deg_is_iupac():
    rng = np.random.default_rng(11)
    for _ in range(500):
        pts = rng.normal(size=(4, 3))
        got = IMP.bff.dihedral_deg(*[_V(p) for p in pts])
        assert _angmod(got - _np_dih(*pts)) < 1e-8


def test_internal2cartesian_exact_and_iupac():
    """Place, measure back: the exact inverse, in the standard convention."""
    rng = np.random.default_rng(12)
    for _ in range(500):
        a, b, c = rng.normal(size=(3, 3))
        r = rng.uniform(0.8, 2.0)
        th = rng.uniform(60.0, 140.0)
        ph = rng.uniform(-180.0, 180.0)
        d = IMP.bff.internal2cartesian(_V(a), _V(b), _V(c), r, th, ph)
        dv = np.array([d[0], d[1], d[2]])
        assert abs(np.linalg.norm(dv - c) - r) < 1e-10
        cb, cd = b - c, dv - c
        ang = np.degrees(np.arccos(np.clip(
            cb @ cd / (np.linalg.norm(cb) * np.linalg.norm(cd)), -1, 1)))
        assert abs(ang - th) < 1e-8
        assert _angmod(_np_dih(a, b, c, dv) - ph) < 1e-8


def test_bond_angle_deg():
    rng = np.random.default_rng(13)
    for _ in range(200):
        a, b, c = rng.normal(size=(3, 3))
        u, w = a - b, c - b
        ref = np.degrees(np.arccos(np.clip(
            u @ w / (np.linalg.norm(u) * np.linalg.norm(w)), -1, 1)))
        assert abs(IMP.bff.bond_angle_deg(_V(a), _V(b), _V(c)) - ref) < 1e-8


def _dye_template():
    """The shipped Alexa488 C1R template: names, elements, coordinates."""
    path = Path(IMP.bff.get_data_path("rotamer_library")) / "A48_C1R.pdb"
    atoms = [l for l in path.read_text().splitlines() if l.startswith("ATOM")]
    elems = [l[76:78].strip() or "C" for l in atoms]
    xyz = [_V([float(l[30:38]), float(l[38:46]), float(l[46:54])])
           for l in atoms]
    return elems, xyz


def test_zmatrix_on_dye_template():
    """A dye+linker is a first-class molecule: 87 perceived bonds (the MOL2
    count), exact encode/decode round-trips."""
    elems, xyz = _dye_template()
    zm = IMP.bff.ZMatrix()
    zm.set_template(xyz, elems)

    n = len(xyz)
    assert n == 83
    assert len(zm.get_bonds()) == 87, "perceived bonds != MOL2 bond table"
    n_base = len(zm.get_base())
    n_rows = len(zm.get_rows()) // 4
    assert n_base + n_rows == n
    assert n_rows > 40  # a dye+linker has dozens of internal-coord rows

    # template round-trips to itself (torsion-exact)
    rec = zm.decode(list(zm.get_template_chi()))
    dev = max(np.linalg.norm(
        np.array([rec[i][k] for k in range(3)])
        - np.array([xyz[i][k] for k in range(3)])) for i in range(n))
    assert dev < 1e-8

    # arbitrary chi: encode(decode(chi)) == chi
    rng = np.random.default_rng(14)
    chi = list(rng.uniform(-180.0, 180.0, n_rows))
    back = zm.encode(zm.decode(chi))
    assert np.max(_angmod(np.array(back) - np.array(chi))) < 1e-8


def test_zmatrix_on_protein_residue():
    """A protein residue from a shipped example -- same kernel, no dye."""
    res = {}
    for l in (REPO / "examples" / "structure" / "T4L" / "3GUN.pdb") \
            .read_text().splitlines():
        if (l.startswith("ATOM") and l[21] == "A" and int(l[22:26]) == 36
                and l[12:16].strip() in ("N", "CA", "CB", "OG", "C", "O")):
            res[l[12:16].strip()] = (
                l[76:78].strip() or "C",
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    assert "OG" in res, "SER 36 not found"
    names = list(res)
    xyz = [_V(res[k][1]) for k in names]
    zm = IMP.bff.ZMatrix()
    zm.set_template(xyz, [res[k][0] for k in names])
    rec = zm.decode(zm.encode(xyz))
    dev = max(np.linalg.norm(
        np.array([rec[i][k] for k in range(3)])
        - np.array(res[names[i]][1])) for i in range(len(names)))
    assert dev < 1e-8


def test_zmatrix_rejects_disconnected():
    zm = IMP.bff.ZMatrix()
    with pytest.raises(Exception):
        zm.set_template([_V([0, 0, 0]), _V([50, 0, 0])], ["C", "C"])


def test_kernel_agrees_with_faspr_repack():
    """The chi the kernel reads off a FASPR-packed residue equals the
    mdtraj-verified value from the prototype's A/B (LYS A:19 = -65.9)."""
    rotlib = REPO / "junk" / "FASPR" / "dun2010bbdep.bin"
    if not rotlib.exists():
        pytest.skip("FASPR Dunbrack library not available")
    import tempfile
    lines = [l for l in (REPO / "examples" / "structure" / "T4L" / "3GUN.pdb")
             .read_text().splitlines()
             if (l.startswith("ATOM") and l[17:20].strip() in STANDARD_AA)
             or l.startswith(("TER", "END"))]
    with tempfile.TemporaryDirectory() as td:
        inp, out = Path(td) / "in.pdb", Path(td) / "out.pdb"
        inp.write_text("\n".join(lines) + "\n")
        IMP.bff.faspr_pack(str(inp), str(out), str(rotlib))
        atoms = {}
        for l in out.read_text().splitlines():
            if l.startswith("ATOM") and l[21] == "A" and int(l[22:26]) == 19:
                atoms[l[12:16].strip()] = np.array(
                    [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    chi1 = IMP.bff.dihedral_deg(*[_V(atoms[a]) for a in ("N", "CA", "CB", "CG")])
    assert _angmod(chi1 - (-65.9)) < 0.5