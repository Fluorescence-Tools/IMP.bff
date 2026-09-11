/*
 * Greedy Olga: which FRET pair to measure next.
 *
 * The kernel takes the two matrices as numpy arrays and publishes a
 * `(pairs, decay)` tuple of managed numpy views. The shapes are part of its
 * contract, stated once in GreedyOlga.h: `(n_frames, n_pairs)` efficiencies in,
 * square pairwise RMSDs in, two flat views of the selection out.
 *
 * `select_informative_sites` adds the homo-oligomer door: an `(n_pairs, 2)`
 * int array of the two site indices each pair connects, and out come the
 * chosen *sites* -- greedied over by the pair set they imply, which is the
 * unit a statistical labelling mix actually buys.
 */

%apply(double* IN_ARRAY2, int DIM1, int DIM2) {
    (double* effs, int n_frames, int n_pairs),
    (double* rmsds, int n_rmsd_rows, int n_rmsd_cols)
};
%apply(int* IN_ARRAY2, int DIM1, int DIM2) {
    (int* pair_sites, int n_site_rows, int n_site_cols)
};

%apply(int** ARGOUTVIEWM_ARRAY1, int* DIM1) {
    (int** out_pairs, int* n_out_pairs),
    (int** out_sites, int* n_out_sites)
};
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {
    (double** out_decay, int* n_out_decay)
};

%include "IMP/bff/GreedyOlga.h"
