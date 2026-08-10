IMP_SWIG_DECORATOR(IMP::bff, AV, AVs);
IMP_SWIG_OBJECT_SERIALIZE(IMP::bff, AVNetworkRestraint, AVNetworkRestraints);

%template(MapStringAVPairDistanceMeasurement) std::map<std::string, IMP::bff::AVPairDistanceMeasurement>;
%attribute_py(IMP::bff::AV, IMP::bff::PathMap, map, get_map);

IMP_SWIG_VALUE(IMP::bff, AVPairDistanceMeasurement, AVPairDistanceMeasurements)
IMP_SWIG_VALUE_SERIALIZE_IMPL(IMP::bff, AVPairDistanceMeasurement)

%include "IMP/bff/AV.h"
%include "IMP/bff/AVNetworkRestraint.h"
