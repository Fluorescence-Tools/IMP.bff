IMP_SWIG_DECORATOR(IMP::bff, AV, AVs);

%template(MapStringAVPairDistanceMeasurement) std::map<std::string, IMP::bff::AVPairDistanceMeasurement>;

%include "IMP/bff/AV.h"
%include "IMP/bff/AVNetworkRestraint.h"
