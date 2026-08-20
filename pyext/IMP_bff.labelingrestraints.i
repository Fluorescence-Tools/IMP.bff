/*
 * Chi-squared scoring of labelling data, with and without volumes.
 *
 * Both restraints were Python classes holding a list and a dict and calling
 * `chi2_score` -- which was itself Python, four lines of scalar arithmetic in a
 * module of distance conventions. `AccessibleVolume` is a C++ value, so the
 * network restraint holds the volumes themselves rather than proxies for them,
 * and the distance it scores never crosses the boundary.
 */

IMP_SWIG_VALUE(IMP::bff, AVMeasurement, AVMeasurements);
IMP_SWIG_VALUE(IMP::bff, LabelingSite, LabelingSites);
IMP_SWIG_VALUE(IMP::bff, SimpleAVNetworkRestraint, SimpleAVNetworkRestraints);
IMP_SWIG_VALUE(IMP::bff, DirectLabelingRestraint, DirectLabelingRestraints);

// `(N, 3)` attachment coordinates, straight from numpy.
%apply(double* IN_ARRAY2, int DIM1, int DIM2) {
    (double* xyz, int n_atoms, int n_dim)
};

%include "IMP/bff/LabelingRestraints.h"

%template(MapStringAccessibleVolume) std::map<std::string, IMP::bff::AccessibleVolume>;
%template(AVMeasurementVector) std::vector<IMP::bff::AVMeasurement>;
%template(LabelingSiteVector) std::vector<IMP::bff::LabelingSite>;
