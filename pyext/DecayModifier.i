%attribute(IMP::bff::DecayModifier, bool, active, is_active, set_active);
%attribute_py(IMP::bff::DecayModifier, IMP::bff::DecayCurve*, data, get_data, set_data);
%class_callable(IMP::bff::DecayModifier, add);

%include "IMP/bff/DecayModifier.h"
