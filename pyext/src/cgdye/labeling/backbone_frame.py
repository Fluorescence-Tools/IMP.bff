"""Backbone reference frame computation for dye attachment.

Implements the standard reference frame:
  origin = CA position
  x-axis = (N - CA) normalized
  y-axis = cross(z, x) where z = cross(x, C-CA)
  z-axis = perpendicular to peptide plane
"""

import IMP
import IMP.atom
import IMP.algebra


def _find_atom(hierarchy, chain_id, resnum, atom_name):
    target_type = _atom_type_from_name(atom_name)
    for a in IMP.atom.get_by_type(hierarchy, IMP.atom.ATOM_TYPE):
        at = IMP.atom.Atom(a)
        if (
            target_type != IMP.atom.AtomType("UNK")
            and at.get_atom_type() == target_type
        ):
            pass
        elif _atom_name(a).upper() == atom_name.upper():
            pass
        else:
            continue
        res_p = a.get_parent()
        if not IMP.atom.Residue.get_is_setup(res_p):
            continue
        res = IMP.atom.Residue(res_p)
        if res.get_index() != resnum:
            continue
        chain_p = res_p.get_parent()
        if not IMP.atom.Chain.get_is_setup(chain_p):
            continue
        chain = IMP.atom.Chain(chain_p)
        if chain.get_id() != chain_id:
            continue
        return a
    return None


def _atom_type_from_name(name):
    n = name.upper()
    mapping = {
        "N": IMP.atom.AtomType("N"),
        "CA": IMP.atom.AtomType("CA"),
        "C": IMP.atom.AtomType("C"),
        "O": IMP.atom.AtomType("O"),
        "CB": IMP.atom.AtomType("CB"),
    }
    return mapping.get(n, IMP.atom.AtomType("UNK"))


def _atom_name(particle):
    at = IMP.atom.Atom(particle)
    tname = at.get_atom_type().get_string()
    return tname


def _xyz(particle):
    return IMP.core.XYZ(particle).get_coordinates()


def backbone_frame(hierarchy, chain_id, resnum):
    """Compute the backbone reference frame at a residue.

    Returns IMP.algebra.ReferenceFrame3D with:
      origin = CA
      x-axis = N - CA direction
      y-axis = in peptide plane, perpendicular to x
      z-axis = perpendicular to peptide plane
    """
    ca_p = _find_atom(hierarchy, chain_id, resnum, "CA")
    n_p = _find_atom(hierarchy, chain_id, resnum, "N")
    c_p = _find_atom(hierarchy, chain_id, resnum, "C")

    if ca_p is None:
        raise ValueError(f"CA atom not found for chain {chain_id} residue {resnum}")
    if n_p is None:
        raise ValueError(f"N atom not found for chain {chain_id} residue {resnum}")
    if c_p is None:
        raise ValueError(f"C atom not found for chain {chain_id} residue {resnum}")

    return backbone_frame_from_coords(_xyz(ca_p), _xyz(n_p), _xyz(c_p))


def _cross(a, b):
    return IMP.algebra.Vector3D(IMP.algebra.get_vector_product(a, b))


def backbone_frame_from_coords(ca, n, c):
    """Compute backbone frame from coordinate vectors.

    Args:
        ca: IMP.algebra.Vector3D — CA position (origin)
        n: IMP.algebra.Vector3D — N position
        c: IMP.algebra.Vector3D — C position

    Returns:
        IMP.algebra.ReferenceFrame3D
    """
    x = n - ca
    x = x / x.get_magnitude()

    yt = c - ca
    yt = yt / yt.get_magnitude()

    z = _cross(x, yt)
    z = z / z.get_magnitude()

    y = _cross(z, x)

    rot = IMP.algebra.get_rotation_from_x_y_axes(x, y)
    trans = IMP.algebra.Transformation3D(rot, ca)
    return IMP.algebra.ReferenceFrame3D(trans)


def frame_to_transformation(frame):
    """Extract the Transformation3D from a ReferenceFrame3D."""
    return frame.get_transformation_to()


def backbone_transformation(hierarchy, chain_id, resnum):
    """Compute the Transformation3D for a backbone residue."""
    frame = backbone_frame(hierarchy, chain_id, resnum)
    return frame_to_transformation(frame)


def backbone_transformation_from_coords(ca, n, c):
    """Compute the Transformation3D from backbone coordinates."""
    frame = backbone_frame_from_coords(ca, n, c)
    return frame_to_transformation(frame)
