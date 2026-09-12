/*
 * A TCSPC decay as a node: the multi-exponential model curve a
 * time-correlated instrument produces, so a lifetime fit joins a parse fit
 * on the graph instead of returning to Python once per iteration.
 *
 * The kernels are tttrlib's, taken header-only from the vendored copy of
 * DecayConvolution.h; this class is the graph around them, exactly as
 * GraphExpression is the graph around tttrlib's expression engine.
 */
%apply(double* IN_ARRAY1, int DIM1) {(double* in_response, int n_response)};
%apply(double* IN_ARRAY1, int DIM1) {(double* in_data_y, int n_data_y)};
%apply(double* IN_ARRAY1, int DIM1) {(double* in_data_ey, int n_data_ey)};
%apply(double* IN_ARRAY1, int DIM1) {(double* in_table, int n_table)};
%shared_ptr(IMP::bff::TCSPCDecay);
%include "IMP/bff/TCSPCDecay.h"

// Complete multi-exponential TCSPC topology family over canonical owner
// ports, with native objectives and BIC scoring.
%shared_ptr(IMP::bff::TCSPCLifetimeSearchSpace);
%include "IMP/bff/TCSPCModelSearch.h"
