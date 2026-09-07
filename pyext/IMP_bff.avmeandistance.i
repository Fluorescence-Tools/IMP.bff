/*
 * The mean-position FRET restraint and the set an fps.json network builds
 * (the precision of a docked model is modelprecision.i, core).
 *
 * `AVMeanDistanceRestraint` computes its gradient, so an fps.json-driven
 * docking can be minimised and not only sampled; the branch costs nothing when
 * no accumulator is passed.
 *
 * There is no PMI restraint here. `IMP.pmi.restraints.RestraintBase` is a
 * Python class in a module this one does not depend on, and what a wrapper
 * around it would do apart from PMI's bookkeeping -- build the network
 * restraint, or one mean-position restraint per distance over rigid-body
 * volumes, and give the volumes a radius and a mass -- is
 * `probe_network_restraint_set` in `AVMeanDistanceRestraint.h`. A caller who
 * wants a `RestraintBase` writes one around that set in the program that
 * already imports PMI; `bin/imp_bff` adds the set to the model directly.
 */

IMP_SWIG_OBJECT(IMP::bff, AVMeanDistanceRestraint, AVMeanDistanceRestraints);
IMP_SWIG_OBJECT(IMP::bff, AVFlatBottomRestraint, AVFlatBottomRestraints);
IMP_SWIG_OBJECT(IMP::bff, AVRebuildOptimizerState, AVRebuildOptimizerStates);

// `ProbeParticle` and `MDRestraintSystem` are plain records, not IMP values:
// one carries a `IMP::Pointer` and neither has an ordering, so they are
// wrapped as structs and their vector members as std::vector templates.
%template(ProbeParticles) std::vector<IMP::bff::ProbeParticle>;

// `get_restraints` hands back an owned IMP object, so it gets IMP's ownership
// typemap rather than the default pointer wrapper.
%newobject IMP::bff::MDRestraintSystem::get_restraints;

// The constructor takes its particles as `IMP::ParticleIndexAdaptor` (a
// Particle, an AV decorator or an index all convert -- the same adaptor
// IMP::core::AngleRestraint uses) and its keyword arguments come from
// SWIG, so the shadow that did both by hand is gone.
%feature("kwargs") IMP::bff::AVMeanDistanceRestraint::AVMeanDistanceRestraint;

%include "IMP/bff/AVMeanDistanceRestraint.h"
