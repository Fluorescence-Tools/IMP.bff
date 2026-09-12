/* The photon trace returns delay and emitted-flag views, and the decay curve
   accumulates in place into the caller's histogram. The %apply has to sit
   here, before the header that declares them. */
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_delays, int* n_delays)};
%apply(unsigned char** ARGOUTVIEWM_ARRAY1, int* DIM1) {(unsigned char** out_emitted, int* n_emitted)};
%apply(double* INPLACE_ARRAY1, int DIM1) {(double* decay, int n_decay)};
%include "IMP/bff/PhotophysicsPhotonSimulation.h"
