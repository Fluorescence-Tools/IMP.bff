"""Photo-induced electron transfer (PET) quenching of a tethered dye.

The reference chemistry and the per-residue interaction table that describe how
much each amino-acid type quenches a xanthene dye in contact, and where on the
residue the quencher sits.

Moved here from QuEst (``quest/core/dye_diffusion.py``) by PRD-109: this is
general fluorescence machinery, not application logic, and ``IMP.bff`` had no
PET quenching of any kind before.

Four quantities describe a residue type:

``kQ``
    Quenching rate in 1/ns applied while the dye is in contact.
``quench_radius``
    Contact radius in Angstrom, measured from the dye **centre** to the
    residue's quenching centre. ``None`` means "inherit the model-wide
    critical distance".
``quench_atoms``
    Atom names whose centroid defines the quenching centre.
``slow_factor``
    Diffusion scaling in [0, 1] modelling unspecific stickiness near the
    residue.

Rates from residues whose contact spheres overlap **add up**; stickiness
factors **multiply**.
"""

from __future__ import annotations

from collections import OrderedDict

import numpy as np

import IMP.bff


__all__ = [
    "quenching_rate_per_frame",
    "STANDARD_AMINO_ACID_RESIDUES",
    "QUENCHER_ATOMS",
    "PET_QUENCHING_REFERENCE",
    "DEFAULT_DYE_RADIUS",
    "amino_acid_quenching_defaults",
    "normalize_amino_acid_quenching",
    "quencher_atom_indices",
    "quencher_centers",
]


STANDARD_AMINO_ACID_RESIDUES = (
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
)

#: Atoms carrying the redox-active moiety that mediates PET, per residue type.
#: The quenching centre is the centroid of whichever of these are present, which
#: puts the quencher on the reactive group -- indole, phenol, imidazole ring,
#: thioether or thiol sulfur -- rather than on CB. For tryptophan that moves the
#: quencher ~3.3 A away from CB, which matters at contact distances of a few
#: Angstrom. Residues with no known PET-active group fall back to a side-chain
#: stub, so a user-supplied rate still has a well-defined centre.
QUENCHER_ATOMS = OrderedDict((
    ("ALA", ("CB",)),
    ("ARG", ("CZ", "NE", "NH1", "NH2")),
    ("ASN", ("CG", "OD1", "ND2")),
    ("ASP", ("CG", "OD1", "OD2")),
    ("CYS", ("SG",)),
    ("GLN", ("CD", "OE1", "NE2")),
    ("GLU", ("CD", "OE1", "OE2")),
    ("GLY", ("CA",)),
    ("HIS", ("CG", "ND1", "CD2", "CE1", "NE2")),
    ("ILE", ("CB",)),
    ("LEU", ("CB",)),
    ("LYS", ("NZ",)),
    ("MET", ("SD",)),
    ("PHE", ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")),
    ("PRO", ("N", "CB", "CG", "CD")),
    ("SER", ("OG",)),
    ("THR", ("OG1",)),
    ("TRP", ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2", "NE1", "CG", "CD1")),
    ("TYR", ("CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH")),
    ("VAL", ("CB",)),
))

#: Default dye radius (Angstrom), used to turn surface contact distances into
#: the centre-to-centre radii the trajectory works with.
DEFAULT_DYE_RADIUS = 3.5

#: Reference PET parameters for a xanthene dye (Alexa488-like) in contact with
#: each residue type. Only residues with a redox-active side chain quench;
#: everything else stays at zero and can be enabled by the user.
#:
#: ``contact_distance`` is measured from the *dye surface* to the residue's
#: quenching centre -- roughly van der Waals contact (~3.5 A) plus the offset
#: from the moiety centroid to its outer atoms (~1.4 A for the aromatic rings,
#: zero for the single-atom sulfur centres). Expressing it relative to the dye
#: surface keeps the table transferable between dyes of different size;
#: :func:`amino_acid_quenching_defaults` adds the dye radius to obtain the
#: centre-to-centre ``quench_radius``.
#:
#: **Starting values to be calibrated against measured lifetimes, not
#: constants** -- Peulen et al., J. Phys. Chem. B 2017, 121, 8211.
PET_QUENCHING_REFERENCE = OrderedDict((
    ("TRP", {"kQ": 3.5, "contact_distance": 5.0}),
    ("TYR", {"kQ": 2.0, "contact_distance": 5.0}),
    ("MET", {"kQ": 1.67, "contact_distance": 3.5}),
    ("HIS", {"kQ": 1.0, "contact_distance": 4.7}),
    ("CYS", {"kQ": 0.8, "contact_distance": 3.5}),
    ("PRO", {"kQ": 2.0, "contact_distance": 4.0}),
))


_DEFAULT_TABLE = OrderedDict(
    (
        res,
        {
            "slow_factor": 1.0,
            "kQ": 0.0,
            # ``None`` means "inherit the model-wide critical distance".
            "quench_radius": None,
            "quench_atoms": list(QUENCHER_ATOMS[res]),
        },
    )
    for res in STANDARD_AMINO_ACID_RESIDUES
)


def _residue_name(value) -> str:
    if isinstance(value, bytes):
        return value.decode("ascii", errors="ignore").strip().upper()
    return str(value).strip().upper()


def _clamp_slow_factor(value) -> float:
    return min(1.0, max(0.0, float(value)))


def _quench_radius(value):
    """A positive quench radius, or ``None`` to inherit the global one."""
    if value is None:
        return None
    try:
        radius = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(radius) or radius <= 0.0:
        return None
    return radius


def _quench_atoms(value, default):
    """A de-duplicated list of atom names, falling back to *default*."""
    if value is None:
        return list(default)
    if isinstance(value, (str, bytes)):
        value = [value]
    atoms = []
    for atom in value:
        name = _residue_name(atom)
        if name and name not in atoms:
            atoms.append(name)
    return atoms or list(default)


def normalize_amino_acid_quenching(table=None):
    """The full per-residue interaction table, with defaults filled in.

    Accepts a partial table keyed by residue name, and the legacy shape where
    a bare number meant the slow factor. Unknown residue names are kept, so a
    non-standard residue can be given a rate.
    """
    normalized = OrderedDict(
        (res, dict(params, quench_atoms=list(params["quench_atoms"])))
        for res, params in _DEFAULT_TABLE.items()
    )
    if not table:
        return normalized
    for residue, params in table.items():
        name = _residue_name(residue)
        defaults = normalized.get(name) or {
            "slow_factor": 1.0,
            "kQ": 0.0,
            "quench_radius": None,
            "quench_atoms": list(QUENCHER_ATOMS.get(name, ("CB",))),
        }
        if isinstance(params, dict):
            slow_factor = params.get("slow_factor", defaults["slow_factor"])
            kQ = params.get("kQ", defaults["kQ"])
            quench_radius = params.get("quench_radius", defaults["quench_radius"])
            quench_atoms = params.get("quench_atoms", defaults["quench_atoms"])
        else:
            # Legacy format: a bare number meant the slow factor.
            slow_factor = params
            kQ = defaults["kQ"]
            quench_radius = defaults["quench_radius"]
            quench_atoms = defaults["quench_atoms"]
        normalized[name] = {
            "slow_factor": _clamp_slow_factor(slow_factor),
            "kQ": max(0.0, float(kQ)),
            "quench_radius": _quench_radius(quench_radius),
            "quench_atoms": _quench_atoms(
                quench_atoms, QUENCHER_ATOMS.get(name, ("CB",))
            ),
        }
    return normalized


def amino_acid_quenching_defaults(
    kQ_scale: float = 1.0,
    slow_factor: float = 1.0,
    dye_radius: float = DEFAULT_DYE_RADIUS,
):
    """A full interaction table built from :data:`PET_QUENCHING_REFERENCE`.

    :param kQ_scale:
        Dye-specific multiplier applied to every reference ``kQ``. Dyes that are
        harder to reduce or oxidise use a value below one.
    :param slow_factor:
        Diffusion scaling applied near every residue (unspecific stickiness).
    :param dye_radius:
        Radius of the dye sphere in Angstrom. The trajectory tracks the dye
        *centre*, so the reference surface contact distances are offset by this
        radius to give centre-to-centre quench radii.
    """
    radius = max(0.0, float(dye_radius))
    table = normalize_amino_acid_quenching()
    for residue, params in table.items():
        params["slow_factor"] = _clamp_slow_factor(slow_factor)
        reference = PET_QUENCHING_REFERENCE.get(residue)
        if reference is None:
            continue
        params["kQ"] = max(0.0, float(reference["kQ"]) * float(kQ_scale))
        params["quench_radius"] = _quench_radius(
            radius + float(reference["contact_distance"])
        )
    return table


def quencher_atom_indices(atoms, selection):
    """Atom indices per residue type, for a ``{residue: [atom names]}`` selection.

    *atoms* is a structured array with ``res_name`` and ``atom_name`` fields --
    the shape QuEst's structure reader and :mod:`IMP.bff.fret.av` both produce.
    Returns an ``OrderedDict`` keyed the same way as *selection*, with a
    ``uint32`` index array per residue type (possibly empty).
    """
    res_name = atoms["res_name"]
    atom_name = atoms["atom_name"]
    indices = OrderedDict()
    for residue in selection:
        found = [
            np.where((res_name == residue) & (atom_name == name))[0]
            for name in selection[residue]
        ]
        if found:
            indices[residue] = np.array(np.hstack(found), dtype=np.uint32)
        else:
            indices[residue] = np.array([], dtype=np.uint32)
    return indices


def quencher_centers(atoms, selection):
    """Quenching-centre coordinates per residue type.

    One row per selected atom, grouped by residue type -- the caller decides
    whether to average them into a per-residue centroid or to treat each atom as
    its own centre. Returns an ``OrderedDict`` of ``(n, 3)`` arrays.
    """
    indices = quencher_atom_indices(atoms, selection)
    coord = atoms["coord"]
    return OrderedDict(
        (residue, coord[indices[residue]]) for residue in selection
    )


def _rate_per_frame(collided, k_quench):
    """Total quenching rate per frame. **C++.**

    Sums the rate constants of the quenchers in contact at each frame -- rates
    add, because the channels are parallel.
    """
    collided = np.ascontiguousarray(collided)
    n_frames = int(collided.shape[0])
    return np.asarray(
        IMP.bff.quenching_rate_per_frame(
            [int(v) for v in collided.astype(np.int32).ravel()],
            n_frames,
            np.ascontiguousarray(k_quench, dtype=np.float64).ravel()),
        dtype=np.float64)


def quenching_rate_per_frame(collided, k_quench) -> np.ndarray:
    """Sum the rates of the quenching atoms the dye touched, frame by frame.

    The **per-atom** route to a quenching trace, as against sampling a stamped
    rate grid along the trajectory (:meth:`IMP.bff.DyeDiffusionSimulation.k_quench`).
    Use it when the contact flags are what you have -- from a distance
    calculation against explicit atoms rather than from a voxel map.

    Frames are independent, so each row reduces on its own thread and the
    ``(n_frames, n_atoms)`` product is never materialised -- for a long
    trajectory against a whole protein's quenching atoms that would be the
    largest array in the calculation.

    Moved here from ChiSurf (``chisurf/core/structure/av/dynamic.py``) by
    PRD-109, where it was a plain Python double loop whose docstring already
    claimed the threading it did not have.

    :param collided: ``(n_frames, n_atoms)`` flags, non-zero where the dye was
        within the critical distance of that atom in that frame.
    :param k_quench: ``(n_atoms,)`` quenching rate per atom, in 1/ns.
    :returns: ``(n_frames,)`` total quenching rate.
    """
    collided = np.ascontiguousarray(collided)
    k_quench = np.ascontiguousarray(k_quench, dtype=np.float64)
    if collided.ndim != 2:
        raise ValueError("`collided` must be (n_frames, n_atoms).")
    if collided.shape[1] != k_quench.shape[0]:
        raise ValueError(
            f"`collided` has {collided.shape[1]} atoms but `k_quench` has "
            f"{k_quench.shape[0]}."
        )
    return _rate_per_frame(collided.astype(np.uint8), k_quench)
