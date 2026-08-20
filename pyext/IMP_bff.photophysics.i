/*
 * The processes that deactivate or depolarise a dye.
 *
 * `photophysics.py` was 934 lines around `OrientationFactor.h`, and almost
 * every function in it was a `.ravel()`, a call and a `.reshape()` -- the
 * shapes the kappa-squared kernels' callers want, stated in a module of their
 * own rather than on the kernels.
 *
 * The exchange-kinetics half is `ExchangeFRET.h`: the master equation of an
 * exchanging labelled population, and the three averaging limits beside it.
 */

IMP_SWIG_VALUE(IMP::bff, FRETRegimes, FRETRegimesList);

%rename(_s2_delta_from_anisotropy) IMP::bff::s2_delta_from_anisotropy;
%rename(_fret_efficiency_regimes) IMP::bff::fret_efficiency_regimes;
%rename(_fret_efficiency_exact_kinetic) IMP::bff::fret_efficiency_exact_kinetic;
%rename(_fret_efficiency_exact_kinetic_pair) IMP::bff::fret_efficiency_exact_kinetic_pair;

%include "IMP/bff/OrientationFactor.h"
%include "IMP/bff/ExchangeFRET.h"

%extend IMP::bff::FRETRegimes {
    %pythoncode %{
        def __getitem__(self, key):
            """The dictionary surface the Python result had."""
            return {"static": self.static_efficiency, "dynamic": self.dynamic,
                    "dynamic_plus": self.dynamic_plus,
                    "kappa2_avg": self.kappa2_avg}[key]
    %}
}

%pythoncode %{
def kappa2_from_dipoles(mu_donor, mu_acceptor, r_vectors):
    """``kappa^2`` for every donor/acceptor dipole pair.

    A zero-length separation contributes zero rather than dividing by zero:
    coincident states are reachable and are not an error.

    :param mu_donor: ``(n_d, 3)`` transition dipoles.
    :param mu_acceptor: ``(n_a, 3)``.
    :param r_vectors: ``(n_d, n_a, 3)`` donor-to-acceptor separations.
    """
    d = np.ascontiguousarray(np.asarray(mu_donor, dtype=np.float64))
    a = np.ascontiguousarray(np.asarray(mu_acceptor, dtype=np.float64))
    r = np.ascontiguousarray(np.asarray(r_vectors, dtype=np.float64))
    out = kappa2_dipole_matrix(d.ravel(), a.ravel(), r.ravel())
    return np.asarray(out).reshape(d.shape[0], a.shape[0])


def s2_delta_from_anisotropy(s2_donor, s2_acceptor, r_inf_AD, r_0=0.38):
    """``(s2_delta, delta)`` -- the order parameter of the inter-dye angle, and
    that angle in radians. Sindbert et al., JACS 133:2463 (2011), eq. 10."""
    out = _IMP_bff._s2_delta_from_anisotropy(
        float(s2_donor), float(s2_acceptor), float(r_inf_AD), float(r_0))
    return float(out[0]), float(out[1])


def kappa(donor_dipole, acceptor_dipole):
    """``(distance, kappa)`` for two dipoles given as two points each."""
    d = np.asarray(donor_dipole, dtype=np.float64)
    a = np.asarray(acceptor_dipole, dtype=np.float64)
    out = dipole_kappa_distance(d[0], d[1], a[0], a[1])
    return float(out[0]), float(out[1])


def kappa_distance(d1, d2, a1, a2):
    """``(distance, kappa)`` between the centres of two dipoles."""
    out = dipole_kappa_distance(
        np.asarray(d1, dtype=np.float64), np.asarray(d2, dtype=np.float64),
        np.asarray(a1, dtype=np.float64), np.asarray(a2, dtype=np.float64))
    return float(out[0]), float(out[1])


def kappasq(delta, sD2, sA2, beta1, beta2):
    """kappa^2 from a set of order parameters and angles.

    Sindbert et al., JACS 133:2463 (2011), eq. 9.
    """
    return wobbling_kappa2(float(delta), float(sD2), float(sA2), float(beta1),
                           float(beta2))


def kappa2_distribution_diffusion_with_traps(sD2, sA2, fret_efficiency,
                                             n_samples=10000, n_bins=31,
                                             k2_min=0.0, k2_max=4.0, seed=-1):
    """``p(kappa^2)`` for **diffusion with traps**, given a measured efficiency.

    Each dye is either freely diffusing -- rotating fast enough to average its
    orientation, contributing 2/3 -- or trapped, wobbling in a cone set by its
    order parameter. ``sD2`` and ``sA2`` are the trapped fractions, so the pair
    splits into four sub-populations.

    The name matters: this was ``kappa2_distribution_dynamic``, and "dynamic" in
    kappa^2 usage means the fast-rotation limit, which is only the free/free
    term here.

    :returns: ``(bin_edges, counts, samples)``.
    """
    out = np.asarray(sample_kappa2_diffusion_with_traps(
        float(sD2), float(sA2), float(fret_efficiency), int(n_samples),
        int(n_bins), float(k2_min), float(k2_max), int(seed)))
    n_edges = int(n_bins)
    return out[:n_edges], out[n_edges:2 * n_edges - 1], out[2 * n_edges - 1:]


def kappasq_all_delta(delta, sD2, sA2, step=0.25, n_bins=31, k2_min=0.0,
                      k2_max=4.0):
    """``(scale, hist, kappa2)`` for a wobbling-in-a-cone model at a fixed delta.

    ``beta1`` is swept over ``(0, pi/2)`` and ``phi`` over ``(0, 2 pi)``.
    """
    scale = VectorDouble()
    hist = VectorDouble()
    flat = wobbling_kappa2_distribution_delta(
        float(delta), float(sD2), float(sA2), float(step), int(n_bins),
        float(k2_min), float(k2_max), scale, hist)
    k2 = np.asarray(flat, dtype=np.float64)
    n_beta = max(1, int(np.floor((np.pi / 2.0 - 0.001) /
                                 (step * np.pi / 180.0))) + 1)
    if n_beta and k2.size % n_beta == 0:
        k2 = k2.reshape(n_beta, -1)
    return (np.asarray(scale, dtype=np.float64),
            np.asarray(hist, dtype=np.float64), k2)


def kappa2_distribution_wobbling_in_cone(sD2, sA2, n_bins=81, k2_min=0.0,
                                         k2_max=4.0, n_samples=10000, seed=0):
    """``(scale, hist, kappa2)`` for a wobbling-in-a-cone model, sampled."""
    scale = VectorDouble()
    hist = VectorDouble()
    k2 = wobbling_kappa2_distribution(
        float(sD2), float(sA2), int(n_bins), float(k2_min), float(k2_max),
        int(n_samples), int(seed), scale, hist)
    return (np.asarray(scale, dtype=np.float64),
            np.asarray(hist, dtype=np.float64),
            np.asarray(k2, dtype=np.float64))


def kappa2_isotropic_distribution(k2, normalize=True):
    """``p(kappa^2)`` for isotropically oriented, *static* dipoles.

    The closed form: singular at ``kappa^2 = 1`` and zero above 4. Static
    because each molecule keeps its orientation for the whole excited-state
    lifetime -- the dynamic limit is the delta function at 2/3 instead.
    """
    r = np.asarray(isotropic_kappa2_density(
        np.ascontiguousarray(np.asarray(k2, dtype=np.float64)).ravel()))
    if normalize:
        r = r / max(1.0, r.sum())
    return r.reshape(np.shape(k2))


def kappa2_to_distance_ratio(k2_amp, k2_val, n_bins=32):
    """``(r_ratio, weights, k2_mean)`` -- the distance-ratio distribution.

    A FRET measurement does not see ``kappa^2``; it sees an *apparent* distance,
    related by ``R_app/R_DA = (<kappa^2>/kappa^2)**(1/6)``. This is that change
    of variable, Jacobian included.
    """
    out = np.asarray(kappa2_distance_ratio_transform(
        np.ascontiguousarray(np.asarray(k2_amp, dtype=np.float64)).ravel(),
        np.ascontiguousarray(np.asarray(k2_val, dtype=np.float64)).ravel(),
        int(n_bins)))
    nb = int(n_bins)
    return out[:nb], out[nb:2 * nb], float(out[2 * nb])


def convolve_distance_with_k2_ratio(r_da, amp_r_da, r_ratio, weights_ratio,
                                    n_bins=256):
    """``R_app = R_DA * (R_app/R_DA)``, as a multiplicative convolution.

    One path, and no flag to pick between two. There used to be a ``use_fast``
    branch calling ``_fast_convolve_loop``, which **was never defined in this
    package** -- it came across in the kappa-squared migration as a call to
    something that stayed behind, so ``use_fast=True`` raised ``NameError``.
    """
    r_da = np.asarray(r_da, dtype=np.float64).ravel()
    amp_r_da = np.asarray(amp_r_da, dtype=np.float64).ravel()
    r_ratio = np.asarray(r_ratio, dtype=np.float64).ravel()
    weights_ratio = np.asarray(weights_ratio, dtype=np.float64).ravel()
    if r_da.size == 0 or r_ratio.size == 0:
        return np.empty(0), np.empty(0)

    lo = float(r_da.min() * r_ratio.min())
    hi = float(r_da.max() * r_ratio.max())
    if not hi > lo:
        # Every product identical -- a delta. numpy widens a zero-width range to
        # (a - 0.5, a + 0.5) rather than returning nothing, and the contract
        # here is the one numpy set.
        lo, hi = lo - 0.5, hi + 0.5
    hist = np.asarray(outer_product_histogram(
        np.ascontiguousarray(r_da), np.ascontiguousarray(amp_r_da),
        np.ascontiguousarray(r_ratio), np.ascontiguousarray(weights_ratio),
        int(n_bins), lo, hi), dtype=np.float64)

    edges = np.linspace(lo, hi, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    peak = hist.max() if hist.size else 0.0
    mask = hist > 1e-10 * peak
    return centers[mask], hist[mask]


def fret_efficiency_regimes(dist_matrix, kappa2_matrix, weights_d, weights_a,
                            R0=52.0):
    """The three averaging limits of a FRET pair, side by side."""
    return _IMP_bff._fret_efficiency_regimes(
        np.ascontiguousarray(np.asarray(dist_matrix, dtype=np.float64)).ravel(),
        np.ascontiguousarray(np.asarray(kappa2_matrix, dtype=np.float64)).ravel(),
        np.asarray(weights_d, dtype=np.float64).ravel(),
        np.asarray(weights_a, dtype=np.float64).ravel(), float(R0))


def fret_efficiency_exact_kinetic(p_matrix, fret_rates, tau0, dt=1.0,
                                  weights=None):
    """The exact FRET efficiency of an exchanging ensemble.

    Solves the master equation rather than assuming a limit: the integrated
    populations solve ``(k_rad I + K_fret - M)^T G = w`` with
    ``M = (P - I)/dt``, and the efficiency is ``sum_i k_fret_i G_i``.

    **The transpose is load-bearing.** ``p_matrix[i][j]`` is the probability of
    ``i -> j``, so populations evolve as a row vector. Solving ``A G = w``
    instead is right only for a symmetric ``M``, and for a real transition
    matrix it gives the wrong fast-exchange limit.
    """
    p = np.ascontiguousarray(np.asarray(p_matrix, dtype=np.float64))
    return _IMP_bff._fret_efficiency_exact_kinetic(
        p.ravel(), np.asarray(fret_rates, dtype=np.float64).ravel(),
        float(tau0), float(dt),
        np.zeros(0) if weights is None
        else np.asarray(weights, dtype=np.float64).ravel())


def fret_efficiency_exact_kinetic_pair(dist_matrix, kappa2_matrix, p_d, p_a,
                                       weights_d, weights_a, R0=52.0, tau0=4.0,
                                       dt=0.1):
    """The exact efficiency for two exchanging ensembles.

    The joint kinetics is the Kronecker product of the two transition matrices,
    which is ``(nd*na)^2`` -- expensive for a large library, and the reason this
    is not the default.
    """
    flat = lambda a: np.ascontiguousarray(np.asarray(a, dtype=np.float64)).ravel()
    return _IMP_bff._fret_efficiency_exact_kinetic_pair(
        flat(dist_matrix), flat(kappa2_matrix), flat(p_d), flat(p_a),
        flat(weights_d), flat(weights_a), float(R0), float(tau0), float(dt))
%}
