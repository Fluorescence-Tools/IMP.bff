/*
 * Structures in and out: PDB, MOL2, mmCIF, DCD, and the tables beside them.
 *
 * There is no Python left here, and no lazy door. The four RMF functions --
 * `write_rmf`, `write_rotamer_library_rmf`, `read_rotamer_library_rmf` and
 * `protein_frames_from_rmf` -- were `%pythoncode` builders registered in a
 * `_LAZY` table, so that `import IMP.bff` would not require `IMP.rmf`. `rmf`
 * is one of this module's modules now (`dependencies.py`), and they are C++ in
 * `RmfIO.h`; the library the two rotamer ones carry is the same
 * #IMP::bff::RotamerLibrary every other reader returns, where it used to be a
 * dict from here and a value from everywhere else.
 *
 * `get_anchor_cb_position` went the same way, to `ProbeAttachment.h` beside the
 * rest of the labelling layer -- `rotamer` is a module of this one too. Its
 * companion `place_probe_from_rotamer_cb` is not ported: it computed the
 * backbone-dependent C-beta, logged it, and then called
 * `place_probe_from_coords(label, ca, n, c)`, which is what a caller gets by
 * calling that function. The C-beta was never used.
 *
 * `read_angle_file` and `load_structure_with_particles` were Python on the
 * argument that turning `ParticleIndex` lists into `Particle`/`Bond` objects
 * is IMP API glue -- but C++ decorates a particle as readily, and the flexfit
 * block is JSON, which C++ already parses everywhere else in this module.
 * Both are in `StructureIO.h`.
 *
 * The flat coordinate buffers come back as numpy views; the `(n, 3)` /
 * `(n_frames, n_atoms, 3)` reshape is the caller's, which is where that shape
 * belongs.
 *
 * The Python this replaces defined a `DCDFormatError(ValueError)` for the
 * reader to raise. It is gone: `IMP_THROW` raises `IMP::ValueException`, which
 * is itself a `ValueError`. Nothing but the reader's own test ever named it,
 * and `except ValueError` catches what it caught.
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

// The spheres `write_rmf` writes, straight from numpy.
%apply(double* IN_ARRAY2, int DIM1, int DIM2) {
    (double* coords, int n_atoms, int n_dim)
};

%include "IMP/bff/RmfIO.h"
