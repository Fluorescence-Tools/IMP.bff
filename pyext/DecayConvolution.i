%attribute(IMP::bff::DecayConvolution, IMP::bff::DecayCurve&, corrected_irf, get_corrected_irf);
%attribute(IMP::bff::DecayConvolution, int, convolution_method, get_convolution_method, set_convolution_method);
%attribute(IMP::bff::DecayConvolution, double, excitation_period, get_excitation_period, set_excitation_period);
%attribute(IMP::bff::DecayConvolution, double, irf_shift_channels, get_irf_shift_channels, set_irf_shift_channels);
%attribute(IMP::bff::DecayConvolution, double, irf_background_counts, get_irf_background_counts, set_irf_background_counts);
%attribute_py(IMP::bff::DecayConvolution, IMP::bff::DecayCurve, irf, get_irf, set_irf);

%include "IMP/bff/DecayConvolution.h"
