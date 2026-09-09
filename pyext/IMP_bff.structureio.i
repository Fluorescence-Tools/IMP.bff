/*
 * Structures in and out: PDB, MOL2, mmCIF, DCD, RMF, and the tables beside
 * them.
 *
 * Flat coordinate buffers come back as numpy views; the `(n, 3)` /
 * `(n_frames, n_atoms, 3)` reshape is the caller's, which is where that shape
 * belongs. A malformed trajectory raises `IMP::ValueException`, which is a
 * Python `ValueError`, so there is no format-specific exception type to catch.
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

/* The trajectory writer takes keyword arguments -- a caller names the frame
   or its metadata and not the other. */
%feature("kwargs") IMP::bff::RmfStructureWriter::RmfStructureWriter;
%feature("kwargs") IMP::bff::RmfStructureWriter::append;

%include "IMP/bff/RmfIO.h"

%attribute(IMP::bff::RmfStructureWriter, int, n_atoms, get_n_atoms);
%attribute(IMP::bff::RmfStructureWriter, int, n_frames, get_number_of_frames);

/* `with RmfStructureWriter(...) as w:` closes the file on the way out, which
   is what every caller of the writer this replaces was already doing. */
%extend IMP::bff::RmfStructureWriter {
    %pythoncode {
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            self.close()
            return False
    }
}
