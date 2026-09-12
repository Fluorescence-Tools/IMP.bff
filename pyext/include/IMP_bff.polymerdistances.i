// PolymerDistances scientific graph producer.
%apply(double* IN_ARRAY1, int DIM1) {(double* in_axis, int n_axis)};
%shared_ptr(IMP::bff::PolymerDistances);
%include "IMP/bff/PolymerDistances.h"
