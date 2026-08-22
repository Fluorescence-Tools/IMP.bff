/*
 * Structures in and out: PDB, MOL2, mmCIF, DCD, and the tables beside them.
 *
 * The **RMF door stays Python**, and lazily so. Writing RMF needs `IMP.rmf`,
 * which is not one of this module's `required_modules`; the decision not to
 * widen that graph was already taken for the trajectory frame writers (see
 * `include/IMP/bff/DyeSampling.h`), and three library functions are not a
 * reason to reverse it. They are built on first use through the same `_LAZY`
 * table `AVNetworkRestraintWrapper` uses, so `import IMP.bff` still does not
 * pull RMF.
 *
 * `read_angle_file` likewise stays Python: it turns C++ `ParticleIndex` lists
 * into IMP `Particle`/`Bond` objects, which is IMP API glue, not a kernel.
 *
 * Everything else is C++ directly. The flat coordinate buffers come back as
 * numpy views; the `(n, 3)` / `(n_frames, n_atoms, 3)` reshape is the caller's,
 * which is where that shape belongs.
 */

IMP_SWIG_VALUE(IMP::bff, DCDHeader, DCDHeaders);
IMP_SWIG_VALUE(IMP::bff, AtomBond, AtomBonds);
IMP_SWIG_VALUE(IMP::bff, CrossLink, CrossLinks);
IMP_SWIG_VALUE(IMP::bff, AtomReference, AtomReferences);

%include "IMP/bff/TrajectoryIO.h"
%include "IMP/bff/StructureIO.h"

%template(AtomBondList) std::vector<IMP::bff::AtomBond>;
%template(CrossLinkList) std::vector<IMP::bff::CrossLink>;
%template(AtomReferenceList) std::vector<IMP::bff::AtomReference>;

%pythoncode %{
# The Python this replaces defined a `DCDFormatError(ValueError)` for the
# reader to raise. It is gone: `IMP_THROW` raises `IMP::ValueException`, which
# is itself a `ValueError`. Nothing but the reader's own test ever named it,
# and `except ValueError` catches what it caught.


def load_structure_with_particles(pdb_path):
    """Load a PDB and return ``(coords, leaves, model, hierarchy)``.

    The model is handed back because it owns the particles: drop it and the
    hierarchy is a decorator over nothing. The coords are the flat view,
    reshaped to ``(-1, 3)`` here because that is the one shape the hierarchy
    fixes (three per leaf).
    """
    import IMP
    import IMP.atom
    model = IMP.Model()
    hier = read_pdb_hierarchy(str(pdb_path), model)
    coords = structure_coordinates(hier).reshape(-1, 3)
    return coords, IMP.atom.get_leaves(hier), model, hier


def read_angle_file(hier, flex_dict):
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
    for bend in flex_dict["Bonds"]:
        for end in (bend[0], bend[1]):
            pairs.append(AtomReference(
                str(end["chain_identifier"]), int(end["residue_seq_number"]),
                str(end["atom_name"])))
    bonds = [IMP.atom.Bond(model, i) for i in create_named_bonds(hier, pairs)]
    return flexible, bonds


def _build_write_rmf():
    """`write_rmf`, defined on first use. Lazy because `import RMF` at module
    scope would make `import IMP.bff` require a module this one does not depend
    on."""
    import json
    import IMP
    import IMP.atom
    import IMP.core
    import IMP.rmf
    import RMF

    def write_rmf(atoms, path, model_name="structure", transform=None,
                  metadata=None, radius=1.5):
        """Write `(N, 3)` coordinates to an RMF file via `IMP.rmf`."""
        coords = np.ascontiguousarray(
            np.asarray(atoms, dtype=np.float64)[:, :3])
        t = np.ascontiguousarray(
            np.asarray(transform, dtype=np.float64)).ravel() if transform is not None else np.zeros(0)
        if t.size:
            coords = _IMP_bff.apply_transform(coords.ravel(), t).reshape(-1, 3)

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
        """Write a rotamer library to an RMF file, one frame per rotamer."""
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
        """Read a rotamer library from an RMF file."""
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