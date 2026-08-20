/*
 * The shapes a label distribution takes, as numpy arrays.
 *
 * The kernels return `std::vector<double>`, which SWIG hands back as a proxy
 * sequence; every consumer in this package wants an ndarray. That conversion
 * was a Python module of one-line wrappers -- `np.asarray(IMP.bff.f(...))`,
 * fifteen times. It is stated once here instead, on the functions themselves.
 */

%rename(_gaussian_chain) IMP::bff::gaussian_chain;
%rename(_worm_like_chain) IMP::bff::worm_like_chain;
%rename(_worm_like_chain_linker) IMP::bff::worm_like_chain_linker;
%rename(_poisson_0toN) IMP::bff::poisson_0toN;
%rename(_normal_distribution) IMP::bff::normal_distribution;
%rename(_generalized_normal_distribution) IMP::bff::generalized_normal_distribution;
%rename(_distance_between_gaussian) IMP::bff::distance_between_gaussian;
%rename(_polynomial_transfer) IMP::bff::polynomial_transfer;
%rename(_polynomial_transfer_vector) IMP::bff::polynomial_transfer_vector;
%rename(_gaussian_rmp_to_rda_mean) IMP::bff::gaussian_rmp_to_rda_mean;
%rename(_polynomial_transfer_ascending) IMP::bff::polynomial_transfer_ascending;

%include "IMP/bff/Distributions.h"
%include "IMP/bff/PolymerChain.h"

%pythoncode %{
def _dist_axis(x):
    return np.ascontiguousarray(np.asarray(x, dtype=np.float64)).ravel()


def gaussian_chain(distances, segment_length, number_of_segments):
    """Radial distribution of an ideal chain."""
    return np.asarray(_IMP_bff._gaussian_chain(
        _dist_axis(distances), float(segment_length), int(number_of_segments)),
        dtype=np.float64)


def worm_like_chain(distances, kappa, chain_length=0.0, normalize=True,
                    distance=True):
    """Radial distribution of a worm-like chain (Becker, Rosa & Everaers 2010).

    Values at or beyond the contour length stay zero -- a chain cannot be longer
    than itself, and the closed form diverges there.
    """
    return np.asarray(_IMP_bff._worm_like_chain(
        _dist_axis(distances), float(kappa), float(chain_length),
        bool(normalize), bool(distance)), dtype=np.float64)


def worm_like_chain_linker(distances, kappa, chain_length=0.0, sigma=6.0,
                           normalize=True):
    """A worm-like chain broadened by the dye linkers at each end.

    The chain distribution is between the *attachment points*; what a FRET
    experiment measures is between the *dyes*.
    """
    return np.asarray(_IMP_bff._worm_like_chain_linker(
        _dist_axis(distances), float(kappa), float(chain_length), float(sigma),
        bool(normalize)), dtype=np.float64)


def poisson_0toN(lam, N):
    """Poisson probabilities for ``k = 0 .. N-1``."""
    return np.asarray(_IMP_bff._poisson_0toN(float(lam), int(N)),
                      dtype=np.float64)


def normal_distribution(x, loc=0.0, scale=1.0, norm=True):
    """Normal density on ``x``; ``norm`` divides by the sum."""
    return np.asarray(_IMP_bff._normal_distribution(
        _dist_axis(x), float(loc), float(scale), bool(norm)), dtype=np.float64)


def generalized_normal_distribution(x, loc=0.0, scale=1.0, shape=0.0,
                                    norm=True):
    """Normal density with a skew, applied by transforming the axis."""
    return np.asarray(_IMP_bff._generalized_normal_distribution(
        _dist_axis(x), float(loc), float(scale), float(shape), bool(norm)),
        dtype=np.float64)


def distance_between_gaussian(x, sigma1, sigma2, distance):
    """The distance distribution between two isotropic Gaussians."""
    return np.asarray(_IMP_bff._distance_between_gaussian(
        _dist_axis(x), float(sigma1), float(sigma2), float(distance)),
        dtype=np.float64)


def polynomial_transfer(rmp, coeffs):
    """Evaluate a transfer polynomial, **highest power first**.

    That is the order ``numpy.polyfit`` returns and therefore what a fitted
    calibration carries. :func:`polynomial_transfer_ascending` takes the other
    one, and the difference is silent: at ``rmp = 45`` with ``[0, 1, 0.02]``
    they give 45.02 and 85.5.
    """
    c = _dist_axis(coeffs)
    if np.ndim(rmp) == 0:
        return _IMP_bff._polynomial_transfer(float(rmp), c)
    x = np.asarray(rmp, dtype=np.float64)
    y = np.asarray(_IMP_bff._polynomial_transfer_vector(
        np.ascontiguousarray(x.ravel()), c), dtype=np.float64)
    return y.reshape(x.shape)


def polynomial_transfer_ascending(rmp, coeffs):
    """Evaluate a transfer polynomial, **lowest power first**."""
    return polynomial_transfer(rmp, _dist_axis(coeffs)[::-1].copy())


def gaussian_rmp_to_rda_mean(rmp, sigma):
    """``Rmp + sigma^2 / Rmp`` -- the Gaussian mean-position correction.

    ``sigma`` is the per-component width of the *separation vector*. Zero where
    ``rmp`` is zero or negative: the expansion is in ``sigma/Rmp`` and says
    nothing at coincident mean positions.
    """
    if np.ndim(rmp) == 0:
        return _IMP_bff._gaussian_rmp_to_rda_mean(float(rmp), float(sigma))
    x = np.asarray(rmp, dtype=np.float64)
    return np.array([_IMP_bff._gaussian_rmp_to_rda_mean(float(v), float(sigma))
                     for v in x.ravel()]).reshape(x.shape)
%}
