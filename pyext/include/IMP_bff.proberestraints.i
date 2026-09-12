/*
 * Chi-squared scoring of labelling data, with and without volumes.
 *
 * Both restraints were Python classes holding a list and a dict and calling
 * `chi2_score` -- which was itself Python, four lines of scalar arithmetic in a
 * module of distance conventions. `ProbeAccessibleVolume` is a C++ value, so the
 * network restraint holds the volumes themselves rather than proxies for them,
 * and the distance it scores never crosses the boundary.
 */

IMP_SWIG_VALUE(IMP::bff, AVMeasurement, AVMeasurements);
IMP_SWIG_VALUE(IMP::bff, ProbeSite, ProbeSites);
IMP_SWIG_VALUE(IMP::bff, SimpleProbeNetworkRestraint, SimpleProbeNetworkRestraints);
IMP_SWIG_VALUE(IMP::bff, DirectProbeRestraint, DirectProbeRestraints);

// `(N, 3)` attachment coordinates, straight from numpy.
%apply(double* IN_ARRAY2, int DIM1, int DIM2) {
    (double* xyz, int n_atoms, int n_dim)
};
// The overloads without a dimension count cannot bind the array typemap
// (IN_ARRAY2 needs the two int dims to follow) -- refuse them at wrap time.
%ignore IMP::bff::DirectProbeRestraint::DirectProbeRestraint(double* xyz, int n_atoms);
%ignore IMP::bff::DirectProbeRestraint::DirectProbeRestraint(double* xyz);

%include "IMP/bff/ProbeRestraints.h"

%template(MapStringProbeAccessibleVolume) std::map<std::string, IMP::bff::ProbeAccessibleVolume>;
%template(AVMeasurementVector) std::vector<IMP::bff::AVMeasurement>;
%template(ProbeSiteVector) std::vector<IMP::bff::ProbeSite>;
