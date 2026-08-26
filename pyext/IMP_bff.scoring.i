/*
 * Stage-2 scoring, from C++ to the Python surface it had as a module.
 *
 * The inner kernels are in RotamerEnergy.h, wrapped directly in swig.i-in.
 * The orchestration -- the CHARMM36 table, the Lorentz-Berthelot rules, the
 * Boltzmann weight, the AABB pre-filter, the masks and selectors, the
 * end-to-end rotamer score, the single-dye mean-field update and the
 * typed-system walkers -- is C++ now, in Scoring.h. Three things remain here:
 *
 * - `rotamer_mean_field_weights_multi_dye` is C++ too now: it was numpy
 *   iteration over the C++ pair-energy matrix, which is the same loop the
 *   single-dye update already ran in C++, written twice.
 * `internal_topology_system` moved to `TopologyBuild.h` with the rest of the
 * system building: it took `parse_dye_mol2`'s dicts, and the typed
 * `Mol2Component` the C++ reader already returns says the same thing without
 * a dict in the middle -- including the element per site, which the dict
 * version dropped on the floor.
 * `torsion_cosine` and `build_dye_restraints` were here too, on the argument
 * that creating IMP.core objects is API glue -- but C++ calls that API as
 * readily as Python does, and the CHARMM-to-Cosine phase shift is physics
 * that had no business living in an interface file. Both are in Scoring.h.
 */

%include "IMP/bff/Scoring.h"

// RotamerScoreResult is a value type with vector members.
IMP_SWIG_VALUE(IMP::bff, RotamerScoreResult, RotamerScoreResults);

%attribute_np(IMP::bff::RotamerScoreResult, std::vector<double>, weights, get_weights);
%attribute_np(IMP::bff::RotamerScoreResult, std::vector<double>, energies, get_energies);
%attribute_py(IMP::bff::RotamerScoreResult, double, partition, get_partition);

%template(LJSitePairList) std::vector<IMP::bff::LJSitePair>;
