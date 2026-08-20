/*
 * Building an accessible volume: the two front doors, and the PDB read.
 *
 * `representation/av.py` was 1,233 lines around the two doors and the reads
 * they need -- a PDB parser with its own element-column rules, a van der Waals
 * table, an attachment-atom lookup, the FPS strip and its cache, and a
 * threading lock guarding SWIG object construction. None of it is Python-shaped
 * once the object it builds is a C++ value: the lock exists only because IMP's
 * decorators are constructed *through* SWIG, and on this side of the boundary
 * there is nothing to serialise.
 */

IMP_SWIG_VALUE(IMP::bff, PDBAtomRecord, PDBAtomRecords);
IMP_SWIG_VALUE(IMP::bff, StripSelection, StripSelections);

%apply(double* IN_ARRAY2, int DIM1, int DIM2) {
    (double* atoms_xyzr, int n_atoms, int n_cols)
};

%rename(_compute_av) IMP::bff::compute_av;
%rename(_compute_av_from_structure) IMP::bff::compute_av_from_structure;

%include "IMP/bff/StripMask.h"
%include "IMP/bff/AVBuilder.h"

%template(PDBAtomRecordVector) std::vector<IMP::bff::PDBAtomRecord>;
%template(MapIntDouble) std::map<int, double>;

%pythoncode %{
def compute_av(atoms_xyz, atoms_vdw, source_xyz, linker_length=20.0,
               linker_width=0.5, dye_radii=(3.5, 0.0, 0.0),
               grid_resolution=1.5,
               allowed_sphere_radius=DEFAULT_ALLOWED_SPHERE_RADIUS,
               search_stencil=0):
    """An accessible volume from raw atomic coordinates.

    The **array** front door: no PDB file and no fps position definition, which
    is what a restraint has. :func:`compute_av_from_structure` is the other.

    :param atoms_xyz: ``(N, 3)`` obstacle coordinates.
    :param atoms_vdw: ``(N,)`` van der Waals radii, A.
    :param source_xyz: ``(3,)`` where the linker is tied.
    :param dye_radii: ``(r1, r2, r3)``; ``(3.5, 0, 0)`` is the AV1 model.
    :param search_stencil: ``0`` keeps the AV's own default of **74**, the
        LabelLib reference metric. ``26`` is ~1.9x faster for ~20 % less volume.
    """
    xyz = np.ascontiguousarray(np.asarray(atoms_xyz, dtype=np.float64)).reshape(-1, 3)
    vdw = np.ascontiguousarray(np.asarray(atoms_vdw, dtype=np.float64)).ravel()
    r = (list(dye_radii) + [0.0, 0.0, 0.0])[:3]
    return _IMP_bff._compute_av(
        np.ascontiguousarray(np.column_stack([xyz, vdw])),
        np.asarray(source_xyz, dtype=np.float64).ravel(),
        float(linker_length), float(linker_width),
        float(r[0]), float(r[1]), float(r[2]), float(grid_resolution),
        float(allowed_sphere_radius), int(search_stencil))


def compute_av_from_structure(pdb_path, position, disc_step=None):
    """An accessible volume from a structure and an fps position definition.

    The **structure** front door. ``position`` is one entry of an fps.json
    ``Positions`` section; the FPS strip is applied before anything is
    measured, and the attachment atom is resolved by
    ``(chain_identifier, residue_seq_number, atom_name)``.

    A declared ``simulation_grid_resolution`` that disagrees with ``disc_step``
    raises. That field is *written into* the particle from ``disc_step``, so a
    caller who states the resolution in the position and leaves ``disc_step``
    at its default silently gets 1.5 A -- which is not hypothetical: both
    quenching identifiability benchmarks (PRD-110, PRD-111) passed it this way,
    so their ``--resolution`` flag was inert and every run was made at 1.5 A
    while the validation pages said 2.5 A. A resolution being wrong is
    invisible in the result, which is what makes it worth refusing loudly.
    """
    p = dict(position or {})
    declared = p.get("simulation_grid_resolution")
    if declared is not None and disc_step is not None:
        if abs(float(declared) - float(disc_step)) > 1e-9:
            raise ValueError(
                f"simulation_grid_resolution={float(declared)} in the position "
                f"disagrees with disc_step={float(disc_step)}")
    elif declared is not None:
        raise ValueError(
            f"the position declares simulation_grid_resolution="
            f"{float(declared)} but no disc_step was given, so the AV would be "
            f"built at disc_step=1.5")
    return _IMP_bff._compute_av_from_structure(
        str(pdb_path),
        str(p.get("chain_identifier") or ""),
        int(p.get("residue_seq_number", 0)),
        str(p.get("atom_name", "CA")),
        float(p.get("linker_length", 20.0)),
        float(p.get("linker_width", 1.0)),
        float(p.get("radius1", 3.5)),
        float(p.get("radius2", 0.0)),
        float(p.get("radius3", 0.0)),
        float(disc_step if disc_step is not None else 1.5),
        str(p.get("strip_mask") or ""),
        float(p["allowed_sphere_radius"]) if "allowed_sphere_radius" in p else -1.0,
        float(p.get("contact_volume_thickness", 0.0)),
        float(p.get("contact_volume_trapped_fraction", -1.0)))


def compute_avs_for_structure(positions, pdb_path, disc_step=None):
    """Every position of an fps.json ``Positions`` section, keyed by name.

    ``pdb_path`` may be one path or a list; a position's ``body_id`` indexes
    into it, which is how a docking run gives each rigid body its own
    structure. A position whose attachment atom is not in its structure comes
    back as an **empty** volume rather than one computed at a guessed
    coordinate -- a dye attached to an unrelated atom yields a plausible and
    entirely wrong cloud.
    """
    if isinstance(pdb_path, (list, tuple)):
        paths = list(pdb_path)
    elif isinstance(pdb_path, str) and "," in pdb_path:
        paths = [p.strip() for p in pdb_path.split(",")]
    else:
        paths = [pdb_path]

    out = {}
    for name, p in dict(positions).items():
        body = int(p.get("body_id", 0))
        path = str(paths[body] if body < len(paths) else paths[0])
        step = float(disc_step or p.get("simulation_grid_resolution", 1.5))
        found = find_attachment_point(
            path, str(p.get("chain_identifier") or ""),
            int(p.get("residue_seq_number", 0)), str(p.get("atom_name", "CA")))
        if found.size == 0:
            out[name] = AccessibleVolume(grid_step=step, position_name=name)
            continue
        av = compute_av_from_structure(path, p, disc_step=step)
        av.position_name = name
        av.params = p
        out[name] = av
    return out
%}
