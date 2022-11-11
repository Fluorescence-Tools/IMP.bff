IMP_SWIG_DECORATOR(IMP::bff, AV, AVs);

%template(MapStringAVPairDistanceMeasurement) std::map<std::string, IMP::bff::AVPairDistanceMeasurement>;
%attribute_py(IMP::bff::AV, IMP::bff::PathMap, map, get_map);

%include "IMP/bff/AV.h"
%include "IMP/bff/AVNetworkRestraint.h"
