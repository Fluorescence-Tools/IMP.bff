/*
 * A labelled site's donor decay: the field picture and the particle picture,
 * both as C++ objects.
 *
 * The orchestration was Python -- every field is `ng^3` doubles, and the
 * Python held all of them, handed them down to a kernel and took them back on
 * every call. All of that is C++ now: the constructors take the volume and the
 * obstacles as values, `update_*`/`simulate_*`/`run` methods do the work, and
 * every field is read back through a `get_*()` method returning a managed
 * numpy view. The scalar statistics are native attributes via `%attribute`.
 *
 * `donor_decay()` returns a #GridDiffusionResult (time axis, fluorescence,
 * final density) whose `time` is built from the *resolved* step -- no
 * lead-step buffer for Python to split.
 */

IMP_SWIG_VALUE(IMP::bff, ObstacleAtoms, ObstacleAtomsList);
IMP_SWIG_VALUE(IMP::bff, DynamicAccessibleVolume, DynamicAccessibleVolumes);
IMP_SWIG_VALUE(IMP::bff, QuenchedDonorDecay, QuenchedDonorDecays);

// Array-backed fields publish as 1-D managed numpy views.
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {
    (double** out_view, int* n_out_view)
};
%apply(int** ARGOUTVIEWM_ARRAY1, int* DIM1) {
    (int** out_view_i, int* n_out_view_i)
};

%apply(double* IN_ARRAY2, int DIM1, int DIM2) {
    (double* coords, int n_atoms, int n_dim)
};

// Scalar statistics as native attributes (no Python in the wrapper).
%attribute(IMP::bff::DynamicAccessibleVolume, double, tau0, get_tau0);
%attribute(IMP::bff::DynamicAccessibleVolume, double, dg, get_dg);
%attribute(IMP::bff::DynamicAccessibleVolume, int, ng, get_ng);
%attribute(IMP::bff::QuenchedDonorDecay, double, tau0, get_tau0);
%attribute(IMP::bff::QuenchedDonorDecay, double, dg, get_dg);
%attribute(IMP::bff::QuenchedDonorDecay, int, n_photons, get_n_photons);
%attribute(IMP::bff::QuenchedDonorDecay, double, t_step, get_t_step);
%attribute(IMP::bff::QuenchedDonorDecay, bool, has_walk, get_has_walk);
%attribute(IMP::bff::QuenchedDonorDecay, int, n_frames, get_n_frames);
%attribute(IMP::bff::QuenchedDonorDecay, double, mean_k_quench, get_mean_k_quench);
%attribute(IMP::bff::QuenchedDonorDecay, double, collision_fraction, get_collision_fraction);
%attribute(IMP::bff::QuenchedDonorDecay, double, quantum_yield, get_quantum_yield);
%attribute(IMP::bff::QuenchedDonorDecay, double, fluorescence_lifetime, get_fluorescence_lifetime);

%include "IMP/bff/PhotophysicsQuenchingModel.h"