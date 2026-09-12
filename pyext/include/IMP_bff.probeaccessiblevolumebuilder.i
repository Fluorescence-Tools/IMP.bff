/*
 * Building an accessible volume: the two front doors, and the PDB read.
 *
 * The PDB parser with its element-column rules, the van der Waals table, the
 * attachment-atom lookup and the FPS strip with its cache are all C++, beside
 * the doors that use them. No lock guards the build: one is needed only when
 * IMP decorators are constructed *through* SWIG, and on this side of the
 * boundary there is nothing to serialise.
 */

IMP_SWIG_VALUE(IMP::bff, PDBAtomRecord, PDBAtomRecords);
IMP_SWIG_VALUE(IMP::bff, StripReport, StripReports);

%apply(double* IN_ARRAY2, int DIM1, int DIM2) {
    (double* atoms_xyzr, int n_atoms, int n_cols)
};

%template(PDBAtomRecordVector) std::vector<IMP::bff::PDBAtomRecord>;
%template(MapIntDouble) std::map<int, double>;

%include "IMP/bff/StripMask.h"
%include "IMP/bff/ProbeAccessibleVolumeBuilder.h"
