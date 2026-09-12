/*
 * The processes that deactivate or depolarise a probe.
 *
 * The kappa-squared kernels are the C++ surface: flat in, flat (or
 * concatenated) out, with the reshape and the tuple-splitting the caller's.
 * The exchange-kinetics half is `FRETExchange.h`: the master equation of an
 * exchanging labelled population, and the three averaging limits beside it.
 */

IMP_SWIG_VALUE(IMP::bff, FRETRegimes, FRETRegimesList);

// The three arrays a kappa^2 sweep or sample produces, as one value with
// numpy attributes -- it was a return plus two `std::vector<double>&`
// out-parameters, and the class a caller had to build one of is not this
// module's to give (see `IMP_bff.types.i`).
IMP_SWIG_VALUE(IMP::bff, Kappa2Distribution, Kappa2Distributions);

%include "IMP/bff/FRETOrientationFactor.h"

%attribute_np(IMP::bff::Kappa2Distribution, std::vector<double>, values,
              get_values);
%attribute_np(IMP::bff::Kappa2Distribution, std::vector<double>, scale,
              get_scale);
%attribute_np(IMP::bff::Kappa2Distribution, std::vector<double>, hist,
              get_hist);
%include "IMP/bff/FRETExchange.h"