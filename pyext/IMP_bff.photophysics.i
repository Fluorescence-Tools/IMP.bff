/*
 * The processes that deactivate or depolarise a dye.
 *
 * `photophysics.py` was 934 lines around `OrientationFactor.h`, and almost
 * every function in it was a `.ravel()`, a call and a `.reshape()` -- the
 * shapes the kappa-squared kernels' callers want, stated in a module of their
 * own rather than on the kernels.
 *
 * The kernels are the C++ surface now: flat in, flat (or concatenated) out;
 * the reshape and the tuple-splitting are the caller's. The exchange-kinetics
 * half is `ExchangeFRET.h`: the master equation of an exchanging labelled
 * population, and the three averaging limits beside it.
 */

IMP_SWIG_VALUE(IMP::bff, FRETRegimes, FRETRegimesList);

%include "IMP/bff/OrientationFactor.h"
%include "IMP/bff/ExchangeFRET.h"