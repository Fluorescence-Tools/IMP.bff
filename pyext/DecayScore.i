%attributestring(IMP::bff::DecayScore, std::string, score_type, get_score_type, set_score_type);
%attribute(IMP::bff::DecayScore, double, score, score);
%attribute_py(IMP::bff::DecayScore, IMP::bff::DecayCurve, model, get_model, set_model);
%attribute_py(IMP::bff::DecayScore, IMP::bff::DecayCurve, data, get_data, set_data);
%attribute_np(IMP::bff::DecayScore, std::vector<double>, weighted_residuals, get_weighted_residuals);

%include "IMP/bff/DecayScore.h"
