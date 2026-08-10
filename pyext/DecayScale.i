%attribute(IMP::bff::DecayScale, double, number_of_photons, get_number_of_photons);
%attribute(IMP::bff::DecayScale, double, constant_background, get_constant_background, set_constant_background);

IMP_SWIG_VALUE_SERIALIZE_IMPL(IMP::bff, DecayScale)

%include "IMP/bff/DecayScale.h"
