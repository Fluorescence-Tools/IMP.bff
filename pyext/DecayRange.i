%attribute(IMP::bff::DecayRange, int, start, get_start, set_start);
%attribute(IMP::bff::DecayRange, int, stop, get_stop, set_stop);

IMP_SWIG_VALUE_SERIALIZE_IMPL(IMP::bff, DecayRange)

%include "IMP/bff/DecayRange.h"
