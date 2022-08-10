//%apply(double* IN_ARRAY1, int DIM1) {(double* curve1, int n_curve1), (double* curve2, int n_curve2)}
%attribute_py(IMP::bff::DecayCurve, double, acquisition_time, get_acquisition_time, set_acquisition_time);
%attribute_py(IMP::bff::DecayCurve, double, shift, get_shift, set_shift);
%attribute_np(IMP::bff::DecayCurve, std::vector<double>, x, get_x, set_x);
%attribute_np(IMP::bff::DecayCurve, std::vector<double>, y, get_y, set_y);
%attribute_np(IMP::bff::DecayCurve, std::vector<double>, ey, get_ey, set_ey);

%include "IMP/bff/DecayCurve.h"

