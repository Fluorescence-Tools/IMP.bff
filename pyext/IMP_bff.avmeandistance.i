/*
 * The mean-position FRET restraint, the precision of a docked model, and the
 * set an fps.json network builds.
 *
 * `AVMeanDistanceRestraint` had two Python copies -- one in
 * `restraints/network.py` without derivatives and one in `restraints/docking.py`
 * with them -- and the wrapper used the one without, so an fps.json-driven
 * docking could only be sampled, never minimised. There is one now, and it
 * computes the gradient; the branch costs nothing when no accumulator is passed.
 *
 * `AVNetworkRestraintWrapper` was a `%pythoncode` class subclassing
 * `IMP.pmi.restraints.RestraintBase`, built lazily so that `import IMP.bff`
 * would not require IMP.pmi. It is gone, and with it the last lazy door in
 * this file: what it did that is not PMI bookkeeping -- build the network
 * restraint, or one mean-position restraint per distance over rigid-body
 * volumes, and give the volumes a radius and a mass -- is
 * `probe_network_restraint_set` in `AVMeanDistanceRestraint.h`. A caller who
 * wants a `RestraintBase` writes one around that set, in the program that
 * already imports PMI; `bin/imp_bff` does not need one, because it adds the
 * set to the model directly.
 */

IMP_SWIG_OBJECT(IMP::bff, AVMeanDistanceRestraint, AVMeanDistanceRestraints);
IMP_SWIG_VALUE(IMP::bff, PositionUncertainty, PositionUncertainties);

// The constructor takes its particles as `IMP::ParticleIndexAdaptor` (a
// Particle, an AV decorator or an index all convert -- the same adaptor
// IMP::core::AngleRestraint uses) and its keyword arguments come from
// SWIG, so the shadow that did both by hand is gone.
%feature("kwargs") IMP::bff::AVMeanDistanceRestraint::AVMeanDistanceRestraint;

%include "IMP/bff/AVMeanDistanceRestraint.h"
%include "IMP/bff/ModelPrecision.h"

// The scalar RMSF summaries as read-only attributes; the array fields stay
// `get_rmsf()`/`get_mean_coords()` numpy views, reshaped by the caller.
%attribute(IMP::bff::PositionUncertainty, double, rmsf_mean, get_rmsf_mean);
%attribute(IMP::bff::PositionUncertainty, double, rmsf_max, get_rmsf_max);
%attribute(IMP::bff::PositionUncertainty, double, mobile_rmsf_mean,
           get_mobile_rmsf_mean);
