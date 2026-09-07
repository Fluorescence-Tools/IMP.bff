/*
 * Stage-2 scoring.
 *
 * The inner kernels are in RotamerEnergy.h, wrapped directly in swig.i-in;
 * the orchestration -- the CHARMM36 table, the Lorentz-Berthelot rules, the
 * Boltzmann weight, the AABB pre-filter, the masks and selectors, the
 * end-to-end rotamer score, the mean-field updates and the typed-system
 * walkers -- is Scoring.h. This file is the typemaps and the value types.
 */

%include "IMP/bff/Scoring.h"

// RotamerScoreResult is a value type with vector members.
IMP_SWIG_VALUE(IMP::bff, RotamerScoreResult, RotamerScoreResults);

%attribute_np(IMP::bff::RotamerScoreResult, std::vector<double>, weights, get_weights);
%attribute_np(IMP::bff::RotamerScoreResult, std::vector<double>, energies, get_energies);
%attribute_py(IMP::bff::RotamerScoreResult, double, partition, get_partition);

%template(LJSitePairList) std::vector<IMP::bff::LJSitePair>;
