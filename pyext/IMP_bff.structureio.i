/*
 * Structures in and out: PDB, MOL2, mmCIF, DCD, and the tables beside them.
 *
 * The **RMF door stays Python**, and lazily so. Writing RMF needs `IMP.rmf`,
 * which is not one of this module's `required_modules`; the decision not to
 * widen that graph was already taken for the trajectory frame writers (see
 * `include/IMP/bff/DyeSampling.h`), and two library functions are not a reason
 * to reverse it. They are built on first use through the same `_LAZY` table
 * `AVNetworkRestraintWrapper` uses, so `import IMP.bff` still does not pull RMF.
 *
 * Everything else is C++, and the shims below reshape a flat buffer into the
 * `(n, 3)` or `(n_frames, n_atoms, 3)` the callers index.
 */

IMP_SWIG_VALUE(IMP::bff, DCDHeader, DCDHeaders);
IMP_SWIG_VALUE(IMP::bff, AtomBond, AtomBonds);
IMP_SWIG_VALUE(IMP::bff, CrossLink, CrossLinks);
IMP_SWIG_VALUE(IMP::bff, AtomReference, AtomReferences);

%rename(_read_dcd) IMP::bff::read_dcd;
%rename(_read_trajectory) IMP::bff::read_trajectory;
%rename(_read_dcd_header) IMP::bff::read_dcd_header;
%rename(_parse_conect_bonds) IMP::bff::parse_conect_bonds;
%rename(_infer_bonds) IMP::bff::infer_bonds;
%rename(_write_mol2) IMP::bff::write_mol2;
%rename(_read_score_series) IMP::bff::read_score_series;
%rename(_write_pdb) IMP::bff::write_pdb;
%rename(_apply_transform) IMP::bff::apply_transform;
%rename(_load_structure) IMP::bff::load_structure;
%rename(_structure_coordinates) IMP::bff::structure_coordinates;
%rename(_compute_rmsd) IMP::bff::compute_rmsd;
%rename(_read_xlink_table) IMP::bff::read_xlink_table;

%include "IMP/bff/TrajectoryIO.h"
%include "IMP/bff/StructureIO.h"

%template(AtomBondList) std::vector<IMP::bff::AtomBond>;
%template(CrossLinkList) std::vector<IMP::bff::CrossLink>;
%template(AtomReferenceList) std::vector<IMP::bff::AtomReference>;

%pythoncode %{
# The Python this replaces defined a `DCDFormatError(ValueError)` for the
# reader to raise. It is gone: `IMP_THROW` raises `IMP::ValueException`, which
# is itself a `ValueError`, and there is no way to raise a module-defined Python
# class from C++ without a custom exception typemap. Nothing but the reader's
# own test ever named it, and `except ValueError` catches what it caught.


def read_dcd_header(path):
    """Read a DCD header without loading coordinates.

    :returns: ``n_frames``, ``n_atoms``, ``first_step``, ``step_stride``,
        ``time_step``, ``has_unit_cell``, ``charmm_version``, ``endianness``
        and the byte ``offset`` at which frame data begins.
    """
    h = _IMP_bff._read_dcd_header(str(path))
    return {"n_frames": h.n_frames, "n_atoms": h.n_atoms,
            "first_step": h.first_step, "step_stride": h.step_stride,
            "time_step": h.time_step, "has_unit_cell": h.has_unit_cell,
            "charmm_version": h.charmm_version, "endianness": h.endianness,
            "offset": h.offset}


def read_dcd(path, max_frames=None):
    """``(n_frames, n_atoms, 3)`` coordinates from a DCD trajectory.

    In the file's own units -- Angstrom for the bundled libraries.
    """
    n_atoms = _IMP_bff._read_dcd_header(str(path)).n_atoms
    flat = _IMP_bff._read_dcd(str(path),
                              -1 if max_frames is None else int(max_frames))
    return flat.reshape(-1, n_atoms, 3)


def read_trajectory(path, n_atoms=None, max_frames=None):
    """Coordinates from a trajectory, whatever format it is in.

    **BinaryCIF is the format this package stores**; ``.dcd`` still reads,
    because the format is not gone from the world -- a user's own library may be
    one.

    :param n_atoms: required for BinaryCIF, which stores one row per
        (atom, frame) and cannot infer the split. Callers have it from the
        companion PDB; a DCD carries its own.
    """
    path = str(path)
    if path.lower().endswith(".dcd") and n_atoms is None:
        n_atoms = _IMP_bff._read_dcd_header(path).n_atoms
    flat = _IMP_bff._read_trajectory(
        path, -1 if n_atoms is None else int(n_atoms),
        -1 if max_frames is None else int(max_frames))
    return flat.reshape(-1, int(n_atoms), 3)


def parse_pdb_atoms(path):
    """``{serial: {...}}`` for the ATOM/HETATM records of a PDB.

    A view over :func:`read_pdb_records`, which is cached on (path, mtime,
    size). There were two readers of this one format until 2026-08-20, with two
    records behind them and only one of them cached.
    """
    return {a.serial: {"serial": a.serial, "name": a.atom_name,
                       "resname": a.res_name, "x": a.x, "y": a.y, "z": a.z,
                       "elem": a.element}
            for a in read_pdb_records(str(path))}


def parse_conect_bonds(path):
    """The bonds a PDB's ``CONECT`` records declare, as a set of pairs."""
    return {(b.a, b.b) for b in _IMP_bff._parse_conect_bonds(str(path))}


def _as_pdb_atoms(atoms):
    """A ``{serial: {...}}`` mapping, or an already-typed list, as records."""
    if not hasattr(atoms, "values"):
        return atoms
    out = PDBAtomRecords()
    for key in sorted(atoms):
        d = atoms[key]
        a = PDBAtomRecord()
        a.serial = int(d.get("serial", key))
        a.atom_name = str(d.get("name", ""))
        a.res_name = str(d.get("resname", ""))
        a.x, a.y, a.z = float(d["x"]), float(d["y"]), float(d["z"])
        a.element = str(d.get("elem", ""))
        out.append(a)
    return out


def infer_bonds(atoms, scale=1.22):
    """Bonds from geometry, by a covalent-radius cutoff.

    A guess, and that is why it is not the default: ``CONECT`` records are what
    a PDB says, and this is what a distance suggests.
    """
    return {(b.a, b.b)
            for b in _IMP_bff._infer_bonds(_as_pdb_atoms(atoms), float(scale))}


def write_mol2(path, atoms, bonds, mol_name="MOL"):
    """Write a Tripos MOL2 file."""
    typed = AtomBondList()
    for a, b in sorted(bonds):
        typed.append(AtomBond(int(a), int(b)))
    _IMP_bff._write_mol2(str(path), _as_pdb_atoms(atoms), typed, str(mol_name))


def read_score_series(stat_path):
    """``(frames, scores)`` parsed from a PMI stat file or a frame,score CSV."""
    flat = _IMP_bff._read_score_series(str(stat_path))
    pairs = flat.reshape(-1, 2)
    return list(pairs[:, 0]), list(pairs[:, 1])


def _transform_flat(transform):
    """A (4,4), (3,3) or (3,) transform as the flat vector the C++ takes."""
    if transform is None:
        return np.zeros(0)
    t = np.asarray(transform, dtype=np.float64)
    if t.shape not in ((4, 4), (3, 3), (3,)):
        raise ValueError("Unexpected transform shape %s" % (t.shape,))
    return np.ascontiguousarray(t).ravel()


def write_pdb(atoms, path, chain="A", res_name="ALA", transform=None,
              model_index=0):
    """Write ``(N, 3)`` or ``(N, 4)`` coordinates as a single-model PDB.

    A fourth column is an element hint and is ignored.
    """
    a = np.asarray(atoms, dtype=np.float64)
    if a.ndim != 2 or a.shape[1] < 3:
        raise ValueError("Expected (N,3) or (N,4) array, got %s" % (a.shape,))
    _IMP_bff._write_pdb(
        np.ascontiguousarray(a[:, :3]).ravel(), str(path), str(chain),
        str(res_name), _transform_flat(transform), int(model_index))


def load_structure(pdb_path):
    """Load a PDB and return ``(N, 3)`` coordinates."""
    return _IMP_bff._load_structure(str(pdb_path)).reshape(-1, 3)


def load_structure_with_particles(pdb_path):
    """Load a PDB and return ``(coords, leaves, model, hierarchy)``.

    The model is handed back because it owns the particles: drop it and the
    hierarchy is a decorator over nothing.
    """
    import IMP
    import IMP.atom
    model = IMP.Model()
    hier = read_pdb_hierarchy(str(pdb_path), model)
    coords = _IMP_bff._structure_coordinates(hier).reshape(-1, 3)
    return coords, IMP.atom.get_leaves(hier), model, hier


def compute_rmsd(coords_a, coords_b, selection_mask=None, superpose=False):
    """RMSD between two coordinate arrays, optionally after superposition.

    *selection_mask* selects the atoms that enter the RMSD **and** the
    superposition; the fitted rotation is applied to all of *coords_a*, which is
    what makes a selection a reference rather than a crop.
    """
    a = np.ascontiguousarray(np.asarray(coords_a, dtype=np.float64))
    b = np.ascontiguousarray(np.asarray(coords_b, dtype=np.float64))
    if a.shape != b.shape:
        raise ValueError("Shape mismatch: %s vs %s" % (a.shape, b.shape))
    mask = ([] if selection_mask is None else
            [int(bool(v)) for v in np.asarray(selection_mask).ravel()])
    return _IMP_bff._compute_rmsd(a.ravel(), b.ravel(), mask, bool(superpose))


def read_xlink_table(fn):
    """Read a crosslink table. ``{index: {protein_1, residue_1, ...}}``."""
    return {i: {"protein_1": x.protein_1, "protein_2": x.protein_2,
                "residue_1": x.residue_1, "residue_2": x.residue_2}
            for i, x in enumerate(_IMP_bff._read_xlink_table(str(fn)))}


def read_angle_file(hier, flex_dict):
    """Resolve the flexible residues and bonds a flexfit file names.

    :returns: ``(flexible_residues, bonds)`` as IMP particles and bonds.
    """
    import IMP
    import IMP.atom
    model = hier.get_model()

    residues = AtomReferenceList()
    for fr in flex_dict["Flexible residues"]:
        residues.append(AtomReference(
            str(fr["chain_identifier"]), int(fr["residue_seq_number"])))
    flexible = [model.get_particle(i)
                for i in select_flexible_residues(hier, residues)]

    pairs = AtomReferenceList()
    for bnd in flex_dict["Bonds"]:
        for end in (bnd[0], bnd[1]):
            pairs.append(AtomReference(
                str(end["chain_identifier"]), int(end["residue_seq_number"]),
                str(end["atom_name"])))
    bonds = [IMP.atom.Bond(model, i) for i in create_named_bonds(hier, pairs)]
    return flexible, bonds


def _build_write_rmf():
    """`write_rmf`, defined on first use.

    Lazy because `import RMF` at module scope would make `import IMP.bff`
    require a module this one does not depend on.
    """
    import json
    import IMP
    import IMP.atom
    import IMP.core
    import IMP.rmf
    import RMF

    def write_rmf(atoms, path, model_name="structure", transform=None,
                  metadata=None, radius=1.5):
        """Write ``(N, 3)`` coordinates to an RMF file via ``IMP.rmf``.

        Each coordinate becomes an ``XYZR`` ball under a single hierarchy named
        *model_name*; *metadata* is stored as the RMF file description.
        """
        coords = np.ascontiguousarray(
            np.asarray(atoms, dtype=np.float64)[:, :3])
        t = _transform_flat(transform)
        if t.size:
            coords = _IMP_bff._apply_transform(coords.ravel(), t).reshape(-1, 3)

        model = IMP.Model()
        root = IMP.atom.Hierarchy.setup_particle(
            IMP.Particle(model, str(model_name)))
        for i in range(coords.shape[0]):
            p = IMP.Particle(model, "p%d" % i)
            IMP.core.XYZR.setup_particle(
                p, IMP.algebra.Sphere3D(
                    IMP.algebra.Vector3D(*coords[i]), float(radius)))
            IMP.atom.Mass.setup_particle(p, 1.0)
            root.add_child(IMP.atom.Hierarchy.setup_particle(p))

        fh = RMF.create_rmf_file(str(path))
        try:
            if metadata:
                fh.set_description(json.dumps(metadata, default=str))
            IMP.rmf.add_hierarchy(fh, root)
            IMP.rmf.save_frame(fh, "frame_0")
        finally:
            del fh

    return write_rmf


def _build_write_rotamer_library_rmf():
    """`write_rotamer_library_rmf`, defined on first use. See `_build_write_rmf`."""
    import IMP
    import IMP.atom
    import IMP.core
    import IMP.rmf
    import RMF

    def write_rotamer_library_rmf(path, library):
        """Write a rotamer library to an RMF file, one frame per rotamer.

        Weights go on the root node as a per-frame value in the ``score``
        category; the transition matrix, when there is one, is a static value in
        the ``kinetic`` category, flattened because RMF's ``IntsTag`` is 1-D.
        """
        path = str(path)
        if not path.endswith(".rmf3"):
            path += ".rmf3"

        model = IMP.Model()
        root = IMP.atom.Hierarchy.setup_particle(
            IMP.Particle(model, "rotamers"))
        chain = IMP.atom.Chain.setup_particle(IMP.Particle(model, "A"), "A")
        root.add_child(chain)
        res = IMP.atom.Residue.setup_particle(
            IMP.Particle(model, "DYE"), IMP.atom.ResidueType("DYE"), 1)
        chain.add_child(res)

        particles = []
        for name in library["atom_names"]:
            p = IMP.Particle(model, name)
            # Mass rather than Atom, so the particle name survives into the RMF.
            IMP.atom.Mass.setup_particle(p, 1.0)
            IMP.core.XYZR.setup_particle(p)
            IMP.core.XYZR(p).set_radius(1.0)
            res.add_child(IMP.atom.Hierarchy.setup_particle(p))
            particles.append(p)

        fh = RMF.create_rmf_file(path)
        IMP.rmf.add_hierarchies(fh, [root])
        weight_key = fh.get_key(fh.get_category("score"), "weight",
                                RMF.FloatTag())
        trans_key = fh.get_key(fh.get_category("kinetic"), "transition_matrix",
                               RMF.IntsTag())
        root_node = fh.get_root_node()
        if "transitions" in library and library["transitions"] is not None:
            root_node.set_static_value(
                trans_key, np.asarray(library["transitions"]).ravel().tolist())

        for i, rid in enumerate(sorted(library["coords"])):
            for p, c in zip(particles, library["coords"][rid]):
                IMP.core.XYZ(p).set_coordinates(
                    IMP.algebra.Vector3D(c[0], c[1], c[2]))
            IMP.rmf.save_frame(fh, str(rid))
            root_node.set_frame_value(weight_key, float(library["weight"][i]))
        del fh

    return write_rotamer_library_rmf


def _build_read_rotamer_library_rmf():
    """`read_rotamer_library_rmf`, defined on first use. See `_build_write_rmf`."""
    import os
    import IMP
    import IMP.atom
    import IMP.core
    import IMP.rmf
    import RMF

    def read_rotamer_library_rmf(path):
        """Read a rotamer library from an RMF file.

        :returns: ``{"id", "weight", "atom_names", "coords", "transitions"}``.
        """
        path = str(path)
        if not os.path.exists(path):
            if os.path.exists(path + ".rmf3"):
                path += ".rmf3"
            else:
                raise FileNotFoundError("RMF library not found: %s" % path)

        model = IMP.Model()
        fh = RMF.open_rmf_file_read_only(path)
        root = IMP.rmf.create_hierarchies(fh, model)[0]
        atoms = IMP.atom.get_leaves(root)
        weight_key = fh.get_key(fh.get_category("score"), "weight",
                                RMF.FloatTag())
        trans_key = fh.get_key(fh.get_category("kinetic"), "transition_matrix",
                               RMF.IntsTag())
        root_node = fh.get_root_node()

        weights, coords = [], {}
        for i in range(fh.get_number_of_frames()):
            IMP.rmf.load_frame(fh, RMF.FrameID(i))
            weights.append(root_node.get_frame_value(weight_key))
            frame = np.zeros((len(atoms), 3))
            for j, a in enumerate(atoms):
                c = IMP.core.XYZ(a).get_coordinates()
                frame[j] = [c[0], c[1], c[2]]
            coords[i + 1] = frame

        transitions = None
        if root_node.get_has_value(trans_key):
            flat = root_node.get_static_value(trans_key)
            n = int(round(len(flat) ** 0.5))
            transitions = np.asarray(flat).reshape(n, n).tolist()

        return {"id": list(range(1, len(weights) + 1)), "weight": weights,
                "atom_names": [a.get_name() for a in atoms],
                "coords": coords, "transitions": transitions}

    return read_rotamer_library_rmf


_LAZY["write_rmf"] = _build_write_rmf
_LAZY["write_rotamer_library_rmf"] = _build_write_rotamer_library_rmf
_LAZY["read_rotamer_library_rmf"] = _build_read_rotamer_library_rmf
%}
