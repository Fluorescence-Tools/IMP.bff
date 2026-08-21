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

%template(PDBAtomRecordVector) std::vector<IMP::bff::PDBAtomRecord>;
%template(MapIntDouble) std::map<int, double>;
%template(MapStringAccessibleVolume) std::map<std::string, IMP::bff::AccessibleVolume>;

%include "IMP/bff/StripMask.h"
%include "IMP/bff/AVBuilder.h"
