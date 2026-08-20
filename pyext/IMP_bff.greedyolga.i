/*
 * Greedy Olga: which FRET pair to measure next.
 *
 * The kernel takes the two matrices as numpy arrays and returns a small value
 * carrying the selection and its precision decay. What the `.i` adds is the
 * shape contract callers had from the Python this replaced: two-dimensional
 * arrays in, a `(pairs, decay)` tuple of numpy arrays out.
 */

%apply(double* IN_ARRAY2, int DIM1, int DIM2) {
    (double* effs, int n_frames, int n_pairs),
    (double* rmsds, int n_rmsd_rows, int n_rmsd_cols)
};

// The C++ name is taken by the Python surface below, which is what callers use.
%rename(_select_informative_pairs) IMP::bff::select_informative_pairs;

%include "IMP/bff/GreedyOlga.h"

%pythoncode %{
def select_informative_pairs(effs, rmsds, err, max_pairs, unique_only=True,
                             diag_weight=0.99):
    """Olga-style greedy informative pair selection.

    Repeatedly adds the candidate pair that leaves the smallest expected mean
    RMSD between the true structure and the one the measurements would pick
    out, then reports the precision decay over that selection.

    Parameters
    ----------
    effs : numpy.ndarray
        FRET efficiencies, shape ``(n_frames, n_pairs)``. Expected finite:
        Olga's GUI does its NaN filtering before calling the selector.
    rmsds : numpy.ndarray
        Pairwise RMSD between frames, shape ``(n_frames, n_frames)``.
    err : float
        Expected absolute error in FRET efficiency.
    max_pairs : int
        How many pairs to select; capped at ``n_pairs``.
    unique_only : bool
        A candidate may be selected at most once.
    diag_weight : float
        Olga's diagonal correction on the denominator.

    Returns
    -------
    tuple of numpy.ndarray
        The selected pair indices in selection order, and the expected mean
        RMSD left after each of them.
    """
    out = _IMP_bff._select_informative_pairs(
        np.ascontiguousarray(effs, dtype=np.float64),
        np.ascontiguousarray(rmsds, dtype=np.float64),
        float(err), int(max_pairs), bool(unique_only), float(diag_weight))
    return (np.array(out.pairs, dtype=np.int64),
            np.array(out.decay, dtype=np.float64))
%}
