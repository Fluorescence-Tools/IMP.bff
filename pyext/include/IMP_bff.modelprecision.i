/*
 * What repeated docking says about a model's precision: the per-position
 * RMSF summaries as a C++ value (`ModelPrecision.h`, core). Wrapped where the
 * mean-position restraint used to carry it.
 */

IMP_SWIG_VALUE(IMP::bff, PositionUncertainty, PositionUncertainties);

%include "IMP/bff/ModelPrecision.h"

// The scalar RMSF summaries as read-only attributes; the array fields stay
// `get_rmsf()`/`get_mean_coords()` numpy views, reshaped by the caller.
%attribute(IMP::bff::PositionUncertainty, double, rmsf_mean, get_rmsf_mean);
%attribute(IMP::bff::PositionUncertainty, double, rmsf_max, get_rmsf_max);
%attribute(IMP::bff::PositionUncertainty, double, mobile_rmsf_mean,
           get_mobile_rmsf_mean);
