/*
 * The excitation and emission crosstalk matrices, moved here from ChiSurf's
 * `core/fluorescence/crosstalk.py` (owner, 2026-09-04: the definition lives
 * in bff, not in chisurf). A labelled matrix is the one thing both crosstalk
 * matrices share and the one thing two implementations disagree on silently:
 * what a row and a column mean, how a light-path payload is ordered into
 * values against labels, and what forward and inverse mixing do.
 *
 * The kernels' array arguments are `const std::vector<double>&`, so they ride
 * the memcpy input typemap of IMP_bff.types.i -- no `%apply` of their own --
 * and their outputs are named exactly `(out_view, n_out_view)`, which the
 * managed ARGOUTVIEWM_ARRAY1 pair in types.i already claims: the caller gets
 * an ndarray that owns its buffer.
 */
%include "IMP/bff/PhotophysicsCrosstalkMatrix.h"
