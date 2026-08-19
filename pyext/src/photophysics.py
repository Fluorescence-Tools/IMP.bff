"""``IMP.bff.photophysics`` -- the processes that deactivate or depolarise a dye.

Interaction terms with an arity, in the shape a force field has them: a
functional form plus parameters looked up by type.

* **1-body** -- radiative decay, internal conversion, rotational relaxation.
  These depend on the dye and the solvent, not on anything else in the system.
* **2-body** -- PET (a dye and a quencher), FRET (a dye and a dye). Their
  parameters are *pair* properties: ``kQ`` depends on the redox potentials of
  both partners, and R0 on both spectra plus the medium.
* **N-body** -- homo-FRET and multi-chromophore transfer, where the rate is not
  a sum of pairs.

Rate constants **add**, because the channels are parallel -- which is why they
are objects that compose rather than numbers computed inside an observable.

Everything here is **representation-agnostic**: a term consumes
:class:`~IMP.bff.representation.States` -- positions, weights, orientations -- so
one implementation serves an accessible volume, a rotamer library, a
coarse-grained model and an MD trajectory alike.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping, Optional, Sequence
import abc
import typing

import numpy as np

import IMP.bff

__all__ = [
    'FRETTerm',
    'InteractionTerm',
    'PETTerm',
    'RadiativeTerm',
    'kappa2_from_dipoles',
    'kappa2_isotropic',
    'orientation',
    'total_rate',
]

# --------------------------------------------------------------------------
# kappa2
# --------------------------------------------------------------------------
"""Orientation factor κ² of donor/acceptor transition dipoles.

Lives in ``fret`` (label-pair physics, no coordinate model attached); cgdye's
rotamer scoring and the ensemble pair kernels import it. Resolves the PRD-47
κ² placement on the bff side: κ² lives in imp.bff.
"""

def kappa2_from_dipoles(mu_donor: np.ndarray, mu_acceptor: np.ndarray,
                        r_vectors: np.ndarray) -> np.ndarray:
    """``kappa^2`` for every donor/acceptor dipole pair. **C++.**

    A zero-length separation contributes zero rather than a division by zero:
    coincident states are reachable and are not an error.

    :param mu_donor: ``(n_d, 3)`` transition dipoles.
    :param mu_acceptor: ``(n_a, 3)``.
    :param r_vectors: ``(n_d, n_a, 3)`` donor-to-acceptor separations.
    :returns: the ``(n_d, n_a)`` matrix.
    """
    d = np.ascontiguousarray(np.asarray(mu_donor, dtype=np.float64))
    a = np.ascontiguousarray(np.asarray(mu_acceptor, dtype=np.float64))
    r = np.ascontiguousarray(np.asarray(r_vectors, dtype=np.float64))
    out = IMP.bff.kappa2_dipole_matrix(d.ravel(), a.ravel(), r.ravel())
    return np.asarray(out).reshape(d.shape[0], a.shape[0])

def kappa2_isotropic() -> float:
    """⟨κ²⟩ for isotropically and rapidly reorienting dipoles (2/3)."""
    return 2.0 / 3.0


# --------------------------------------------------------------------------
# orientation
# --------------------------------------------------------------------------
"""Orientation: kappa^2 distributions, wobbling in a cone, and order parameters.

The forward model for the orientation factor. Given how freely each dye rotates
within the donor's excited-state lifetime -- expressed as the order parameters
``sD2`` and ``sA2``, which is what a time-resolved anisotropy measures -- this
gives the *distribution* of kappa^2 rather than the isotropic 2/3, and turns it
into the distance-ratio distribution a FRET measurement actually sees.

This module had **no consumers at all** before PRD-113 stage 4. It was never
exported through ``api.py``, so 933 lines of the anisotropy modelling this
package exists to do were unreachable, and it imported ``numba`` directly rather
than through ``IMP.bff._jit`` -- so on an installation without numba it would
not even import. Both were fixed then; the four jitted kernels became C++ in
stage 4d and now live in ``OrientationFactor.h``.

.. warning::
   **The port found that** :func:`kappasq_all` **was returning half the right
   answer.** It drew dipole directions with ``np.random.random(3)`` and
   normalised them, which does not sample the sphere: it fills only the
   positive octant of the unit cube, and non-uniformly at that. In the rigid
   isotropic limit ``sD2 = sA2 = 1``, where the average orientation factor must
   be :math:`2/3`, it gave **0.333**. The C++ draws three standard normals and
   normalises -- uniform on the sphere -- and gives **0.6663** at 4e5 samples.

   Nothing downstream carried the error, because the module had no callers.
   Anyone who took numbers from this function by hand before 2026-08-18 should
   take them again.

The isotropic average is the only analytic check this module has, so it is
worth stating why it exists: :func:`kappasq_all` at ``sD2 = sA2 = 0`` returns
exactly :math:`2/3` for every sample (freely rotating dyes have no orientation
preference at all), and at ``sD2 = sA2 = 1`` -- rigid dyes, random mutual
orientation -- the *mean* is :math:`2/3` while individual values span
:math:`[0, 4]`. A sampler that fails either is not sampling orientations.

The other half of kappa^2 -- the geometry, ``kappa2_from_dipoles``, which turns
two transition dipoles and a separation vector into a number -- is in
:mod:`IMP.bff.photophysics`. Distribution and geometry are different
questions and now live side by side rather than in unrelated packages.
"""

def kappasq_dwt(sD2, sA2, fret_efficiency, n_samples=10000, n_bins=31,
                k2_min=0.0, k2_max=4.0, seed=-1):
    """p(kappa^2) for a dynamic pair, conditioned on a measured efficiency.

    **C++** (:func:`IMP.bff.dynamic_kappa2_distribution`). The Python drew
    10 000 orientation pairs in an explicit loop and called ``kappasq`` three
    times per sample: 128.6 ms against 3.02 ms for 20 000 samples, 43x.

    :returns: ``(bin_edges, counts, samples)`` -- the three arrays the loop
        built, split out of the one buffer the kernel publishes.
    """
    out = np.asarray(IMP.bff.dynamic_kappa2_distribution(
        float(sD2), float(sA2), float(fret_efficiency), int(n_samples),
        int(n_bins), float(k2_min), float(k2_max), int(seed)))
    n_edges = int(n_bins)
    return out[:n_edges], out[n_edges:2 * n_edges - 1], out[2 * n_edges - 1:]

# ``kappasq_all_delta_new`` was **deleted**, not ported.
#
# It was an alternative implementation of ``kappasq_all_delta`` -- a nested
# Python loop over a beta1/beta2 grid -- with no consumer anywhere: not
# exported through ``api.py``, not imported by any module here, not covered by
# a test, and not reachable from chisurf, quest, imp-tricks or ucfret. It
# arrived with the PRD-107 migration and never had one. ChiSurf has its own
# copy of the same function, which is where the name is still live.
#
# The ported one is ``kappasq_all_delta``, over
# :func:`IMP.bff.wobbling_kappa2_distribution_delta`.

def kappasq_all_delta(
        delta: float,
        sD2: float,
        sA2: float,
        step: float = 0.25,
        n_bins: int = 31,
        k2_min: float = 0.0,
        k2_max: float = 4.0
) -> typing.Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Computes a orientation factor distribution for a wobbling in a cone model
    using parameters that can be estimated by experimental anisotropies.

    The function used second rank order parameter of the donor and acceptor
    and the angle delta between the symmetry axes of the dyes as input. These
    parameters can be estimated by the residual anisotropy the the dyes. The
    second rank order parameter of the donor and acceptor are estimated by the
    dye's residual anisotropies. The angle between the symmetry axes is estimated
    by the residual anisotropy of the FRET sensitized emission (see:
    `chisurf.fluorescence.anisotropy.kappa2.s2delta`).

    This function computes a orientation factor distribution, :math:`p(/kappa^2)`,
    for a wobbling in a cone model (WIC) for second rank structure factors of
    the donor and acceptor, and an angle :math:`delta`. The angle
    :math:`delta` is the angle between the symmetry axes of the dyes and can be
    estimated using experimental residual anisotropies [1]_.

    Parameters
    ----------
    delta : float
        The angle delta (in rad) for which the WIC orientation factor
        distribution is calculated.
    sD2 : float
        Second rank order parameter S2 of the donor dye. This can correspond
        to the fraction of trapped donor dye.
    sA2 : float
        Second rank order parameter S2 of the acceptor dye. This can correspond
        to the fraction of trapped acceptor dye.
    step : float
        The step size in degrees that is used to sample the
        angles.
    n_bins : int
        The number of bins in the kappa2 distribution that is
        generated.
    k2_max : float
        Upper kappa2 bound in the generated histogram
    k2_min : float
        Lower kappa2 bound in the generate histogram

    Returns
    -------
    k2scale : numpy-array
        A linear scale in the range of [0, 4] with *n_bins* elements
    k2hist : numpy-array
        The histogram of kappa2 values
    k2 : numpy-array
        A numpy-array containing all computed kappa2 values. The histogram
        corresponds to a histogram over all returned kappa2 values.

    Examples
    --------
    >>> from scikit_fluorescence.modeling.kappa2 import kappasq_all_delta
    >>> k2s, k2h, k2v = kappasq_all_delta(
    ...     delta=0.2,
    ...     sD2=0.15,
    ...     sA2=0.25,
    ...     step=2.0,
    ...     n_bins=31
    ... )
    >>> np.allclose(k2h, np.array([   0.        ,    0.        ,    0.        ,    0.        ,
    ...    3205.72877776, 1001.19048825,  611.44917432,  252.97166906,
    ...       0.        ,    0.        ,    0.        ,    0.        ,
    ...       0.        ,    0.        ,    0.        ,    0.        ,
    ...       0.        ,    0.        ,    0.        ,    0.        ,
    ...       0.        ,    0.        ,    0.        ,    0.        ,
    ...       0.        ,    0.        ,    0.        ,    0.        ,
    ...       0.        ,    0.        ]), rtol=0.3, atol=2.0)
    True

    Notes
    -----
    The angle :math:`/beta_1` is varied in the range (0,pi/2)
    The angle :math:`/phi` is varied in the range (0, 2 pi)

    References
    ----------
    .. [1] Simon Sindbert, Stanislav Kalinin, Hien Nguyen, Andrea Kienzler,
    Lilia Clima, Willi Bannwarth, Bettina Appel, Sabine Mueller, Claus A. M.
    Seidel, "Accurate Distance Determination of Nucleic Acids via Foerster
    Resonance Energy Transfer: Implications of Dye Linker Length and Rigidity"
    vol. 133, pp. 2463-2480, J. Am. Chem. Soc., 2011

    """
    scale = IMP.bff.VectorDouble()
    hist = IMP.bff.VectorDouble()
    flat = IMP.bff.wobbling_kappa2_distribution_delta(
        float(delta), float(sD2), float(sA2), float(step),
        int(n_bins), float(k2_min), float(k2_max), scale, hist)
    k2 = np.asarray(flat, dtype=np.float64)
    n_beta = max(1, int(np.floor((np.pi / 2.0 - 0.001) / (step * np.pi / 180.0))) + 1)
    if n_beta and k2.size % n_beta == 0:
        k2 = k2.reshape(n_beta, -1)
    return (np.asarray(scale, dtype=np.float64),
            np.asarray(hist, dtype=np.float64),
            k2)


def kappasq_all(
        sD2: float,
        sA2: float,
        n_bins: int = 81,
        k2_min: float = 0.0,
        k2_max: float = 4.0,
        n_samples: int = 10000,
        seed: int = 0
) -> typing.Tuple[np.array, np.array, np.array]:
    """Computes a orientation factor distribution for a wobbling in a cone model
    using specific second rank structure factors of the donor and acceptor.

    This function computes a orientation factor distribution, :math:`p(/kappa^2)`,
    for a wobbling in a cone model (WIC) for second rank structure factors of
    the donor and acceptor estimated using experimental residual anisotropies [1]_.

    Parameters
    ----------
    sD2 : float
        Second rank order parameter S2 of the donor dye
    sA2 : float
        Second rank order parameter S2 of the acceptor dye
    n_bins : int
        The number of bins in the kappa2 histogram that is
        generated.
    k2_max : float
        Upper kappa2 bound in the generated histogram
    k2_min : float
        Lower kappa2 bound in the generate histogram
    n_samples : int
        The number random vector pairs that are drawn (default: 10000)

    Returns
    -------
    k2scale : numpy-array
        A linear scale in the range of [0, 4] with *n_bins* elements
    k2hist : numpy-array
        The histogram of kappa2 values
    k2 : numpy-array
        A numpy-array containing all computed kappa2 values. The histogram
        corresponds to a histogram over all returned kappa2 values.

    Examples
    --------
    >>> from scikit_fluorescence.modeling.kappa2 import kappasq_all_delta
    >>> k2_scale, k2_hist, k2 = kappasq_all(
    ...     sD2=0.3,
    ...     sA2=0.5,
    ...     n_bins=31,
    ...     n_samples=100000
    ... )
    >>> k2_scale
    array([0.        , 0.13333333, 0.26666667, 0.4       , 0.53333333,
           0.66666667, 0.8       , 0.93333333, 1.06666667, 1.2       ,
           1.33333333, 1.46666667, 1.6       , 1.73333333, 1.86666667,
           2.        , 2.13333333, 2.26666667, 2.4       , 2.53333333,
           2.66666667, 2.8       , 2.93333333, 3.06666667, 3.2       ,
           3.33333333, 3.46666667, 3.6       , 3.73333333, 3.86666667,
           4.        ])
    >>> reference = np.array([0.0000e+00, 0.0000e+00, 0.0000e+00, 3.1920e+04, 4.3248e+04,
    ...    1.4842e+04, 5.8930e+03, 2.5190e+03, 1.0840e+03, 3.9700e+02,
    ...    9.4000e+01, 3.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
    ...    0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
    ...    0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
    ...    0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00])
    >>> np.allclose(reference, k2_hist, rtol=0.3, atol=2.0)
    True

    References
    ----------
    .. [1] Simon Sindbert, Stanislav Kalinin, Hien Nguyen, Andrea Kienzler,
    Lilia Clima, Willi Bannwarth, Bettina Appel, Sabine Mueller, Claus A. M.
    Seidel, "Accurate Distance Determination of Nucleic Acids via Foerster
    Resonance Energy Transfer: Implications of Dye Linker Length and Rigidity"
    vol. 133, pp. 2463-2480, J. Am. Chem. Soc., 2011

    """
    scale = IMP.bff.VectorDouble()
    hist = IMP.bff.VectorDouble()
    k2 = IMP.bff.wobbling_kappa2_distribution(
        float(sD2), float(sA2), int(n_bins), float(k2_min), float(k2_max),
        int(n_samples), int(seed), scale, hist)
    return (np.asarray(scale, dtype=np.float64),
            np.asarray(hist, dtype=np.float64),
            np.asarray(k2, dtype=np.float64))


def kappa_distance(
        d1: np.array,
        d2: np.array,
        a1: np.array,
        a2: np.array
) -> typing.Tuple[float, float]:
    r"""Calculates the distance between the center of two dipoles and the
    orientation-factor kappa of the dipoles

    Calculates for the vectors d1 and d2 pointing to the donors and the vectors
    a1 and a2 pointing to the ends of the acceptor dipole the orientation
    factor kappa.

    Parameters
    ----------
    d1 : numpy-array
        Vector pointing to the first point of the dipole D
    d2 : numpy-array
        Vector pointing to the second point of the dipole D
    a1 : numpy-array
        Vector pointing to the first point of the dipole A
    a2 : numpy-array
        Vector pointing to the second point of the dipole A

    Returns
    -------
    tuple
        distance between the center of the dipoles and the orientation factor
        for the two dipoles kappa

    Notes
    -----
    The four vectors defining the dipole of the donor :math:`\vec{r}_{D1}` and
    :math:`\vec{r}_{D2}` specified by the parameters `d1` and `d2` and
    :math:`\vec{r}_{A1}` and :math:`\vec{r}_{A2}` specified by the parameters
    `a1` and `a1` are used to compute orientation factor :math:`kappa^2`
    and the distance between the center of the two dipoles :math:`R_{DA}`.

    The distance :math:`R_{DA}` between the dipole centers and :math:`kappa`
    is calculated as follows:

    ..math::

        R_{D,21}=|\vec{r}_{D2} - \vec{r}_{D1}| \\
        R_{A,21}=|\vec{r}_{A2} - \vec{r}_{A1}| \\
        \hat{\mu}_{D}=1/R_{D,21} \cdot (\vec{r}_{D2}-\vec{r}_{D1}) \\
        \hat{\mu}_{A}=1/R_{A,21} \cdot (\vec{r}_{A2}-\vec{r}_{A1}) \\
        \vec{m}_{D}=\vec{r}_{D1}+1/2 \cdot \hat{\mu}_{D} \\
        \vec{m}_{A}=\vec{r}_{A1}+1/2 \cdot \hat{\mu}_{A} \\
        \vec{r}_{DA}=\vec{m}_{D}-\vec{m}_{A} \\
        R_{DA}=|\vec{m}_{D}-\vec{m}_{A}| \\
        \hat{\mu}_{DA}=\vec{r}_{DA} / R_{DA} \\
        \kappa=\langle\mu_A,\mu_D\rangle-3\cdot\langle\mu_D,\mu_{DA}\rangle \cdot \langle\mu_A,\mu_{DA}\rangle

    Examples
    --------
    >>> import scikit_fluorescence.modeling.kappa2
    >>> donor_dipole = np.array(
    ...      [
    ...          [0.0, 0.0, 0.0],
    ...          [1.0, 0.0, 0.0]
    ...      ], dtype=np.float64
    ... )
    >>> acceptor_dipole = np.array(
    ...     [
    ...         [0.0, 0.5, 0.0],
    ...         [0.0, 0.5, 1.0]
    ...     ], dtype=np.float64
    ... )
    >>> scikit_fluorescence.modeling.kappa2.kappa(
    ...     donor_dipole,
    ...     acceptor_dipole
    ... )
    (0.8660254037844386, 1.0000000000000002)

    r"""
    d, k = IMP.bff.dipole_kappa_distance(
        np.asarray(d1, dtype=np.float64).ravel(), np.asarray(d2, dtype=np.float64).ravel(),
        np.asarray(a1, dtype=np.float64).ravel(), np.asarray(a2, dtype=np.float64).ravel())
    return float(d), float(k)


def kappa(
        donor_dipole: np.ndarray,
        acceptor_dipole: np.ndarray
) -> typing.Tuple[float, float]:
    """Calculates the orientation-factor kappa

    :param donor_dipole: 2x3 vector of the donor-dipole
    :param acceptor_dipole: 2x3 vector of the acceptor-dipole
    :return: distance, kappa

    Example
    -------

    >>> import numpy as np
    >>> donor_dipole = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
    >>> acceptor_dipole = np.array([[0.0, 0.5, 0.0], [0.0, 0.5, 1.0]], dtype=np.float64)
    >>> kappa(donor_dipole, acceptor_dipole)
    (0.8660254037844386, 1.0000000000000002)
    """
    return kappa_distance(
        donor_dipole[0], donor_dipole[1],
        acceptor_dipole[0], acceptor_dipole[1]
    )


def s2delta(
        s2_donor: float,
        s2_acceptor: float,
        r_inf_AD: float,
        r_0: float = 0.38
) -> typing.Tuple[float, float]:
    r"""Calculate s2delta from the residual anisotropies of the donor and acceptor

    Parameters
    ----------
    r_0 : float
        Fundamental anisotropy, the anisotropy of the dyes at time zero (
        default value 0.4)
    s2_donor : float
        The second rank oder parameter of the donor dye. The second rank oder
        parameter can be computed using the dye's residual anisotropy (see
        Notes below)
    s2_acceptor : float
        The second rank oder parameter of the direct excited acceptor dye.
    r_inf_AD : float
        The residual anisotropy on the acceptor excited by the donor dye.

    Returns
    -------
    s2delta : float
         A second rank order parameter of the angle [1]_ eq. 10
    delta : float
        The angle between the two symmetry axes of the dipols in units of rad.

    Examples
    --------
    >>> from scikit_fluorescence.modeling.kappa2 import s2delta
    >>> r0 = 0.38
    >>> s2donor = 0.2
    >>> s2acceptor = 0.3
    >>> r_inf_AD = 0.01
    >>> s2delta(
    ...     r_0=r0,
    ...     s2_donor=s2donor,
    ...     s2_acceptor=s2acceptor,
    ...     r_inf_AD=r_inf_AD
    ... )
    (0.4385964912280701, 0.6583029208008411)

    Notes
    -----
    The parameters `s2_donor` and `s2_acceptor`, which correspond to :math:`S^{(2)}_D`
    and :math:`S^{(2)}_A` are calculated using the dye's residual anisotropy [1]_

    ..math::

        S^{(2)}_D = - /sqrt{/frac{r_{D,inf}}{r_{0,D}}} \\
        S^{(2)}_D = /sqrt{/frac{r_{A,inf}}{r_{0,A}}} \\

    References
    ----------

    .. [1] Simon Sindbert, Stanislav Kalinin, Hien Nguyen, Andrea Kienzler,
    Lilia Clima, Willi Bannwarth, Bettina Appel, Sabine Mueller, Claus A. M.
    Seidel, "Accurate Distance Determination of Nucleic Acids via Foerster
    Resonance Energy Transfer: Implications of Dye Linker Length and Rigidity"
    vol. 133, pp. 2463-2480, J. Am. Chem. Soc., 2011

    """
    s2_delta = r_inf_AD/(r_0 * s2_donor * s2_acceptor)
    delta = np.arccos(np.sqrt((2.0 * s2_delta + 1.0) / 3.0))
    return s2_delta, delta


def calculate_kappa_distance(
        xyz: np.array,
        aid1: int,
        aid2: int,
        aia1: int,
        aia2: int
) -> typing.Tuple[np.ndarray, np.ndarray]:
    """Calculates the orientation factor kappa2 and the distance of a
    trajectory given the atom-indices of the donor and the acceptor.

    :param xyz: numpy-array (frame, atom, xyz)
    :param aid1: int, atom-index of d-dipole 1
    :param aid2: int, atom-index of d-dipole 2
    :param aia1: int, atom-index of a-dipole 1
    :param aia2: int, atom-index of a-dipole 2

    :return: distances, kappa2
    """
    n_frames = xyz.shape[0]
    ks = np.empty(n_frames, dtype=np.float32)
    ds = np.empty(n_frames, dtype=np.float32)

    for i_frame in range(n_frames):
        try:
            d, k = kappa_distance(
                xyz[i_frame, aid1], xyz[i_frame, aid2],
                xyz[i_frame, aia1], xyz[i_frame, aia2]
            )
            ks[i_frame] = k
            ds[i_frame] = d
        except:
            print("Frame ", i_frame, "skipped, calculation error")

    return ds, ks


def kappasq(
        delta: float,
        sD2: float,
        sA2: float,
        beta1: float,
        beta2: float
) -> float:
    """Calculates kappa2 given a set of oder parameters and angles

    Parameters
    ----------
    delta : float
        The angle between the symmetry axis of rotation of the dyes in units
        of rad.
    sD2 : float
        The second rank oder parameter of the donor
    sA2 : float
        The second rank oder parameter of the acceptor
    beta1 : float
        The angle between the symmetry axes of the rotation of the dye and
        the distance vector RDA between the two dipoles
    beta2
        The angle between the symmetry axes of the rotation of the dye and
        the distance vector RDA between the two dipoles

    Returns
    -------
    kappa2 : float
        The orientation factor that corresponds to the provided angles.

    Notes
    -----

    This function corresponds to eq. 9 in [1]_.

    References
    ----------
    .. [1] Simon Sindbert, Stanislav Kalinin, Hien Nguyen, Andrea Kienzler,
    Lilia Clima, Willi Bannwarth, Bettina Appel, Sabine Mueller, Claus A. M.
    Seidel, "Accurate Distance Determination of Nucleic Acids via Foerster
    Resonance Energy Transfer: Implications of Dye Linker Length and Rigidity"
    vol. 133, pp. 2463-2480, J. Am. Chem. Soc., 2011

    """
    return IMP.bff.wobbling_kappa2(delta, sD2, sA2, beta1, beta2)


def p_isotropic_orientation_factor(k2: np.ndarray,
                                   normalize: bool = True) -> np.ndarray:
    """``p(kappa^2)`` for isotropically oriented, *static* dipoles. **C++.**

    The closed form: singular at ``kappa^2 = 1`` and zero above 4. Static
    because each molecule keeps its orientation for the whole excited-state
    lifetime -- the dynamic limit is the delta function at 2/3 instead.

    :param normalize: divide by the sum, which is what every caller wants and
        what makes the two branches comparable.
    """
    r = np.asarray(IMP.bff.isotropic_kappa2_density(
        np.ascontiguousarray(np.asarray(k2, dtype=np.float64)).ravel()))
    if normalize:
        r = r / max(1.0, r.sum())
    return r.reshape(np.shape(k2))

def kappa2_to_distance_ratio(k2_amp, k2_val, n_bins: int = 32) -> tuple:
    """Turn a kappa^2 distribution into the distance-ratio distribution. **C++.**

    A FRET measurement does not see ``kappa^2``; it sees an *apparent* distance,
    related by ``R_app/R_DA = (<kappa^2>/kappa^2)**(1/6)``. This is that change
    of variable, Jacobian included -- ``|dk2/dr| = 6 <k2> / r**7``, without
    which the distribution is wrong wherever the mapping is non-linear, which
    is everywhere.

    :returns: ``(r_ratio, weights, k2_mean)``.
    """
    out = np.asarray(IMP.bff.kappa2_distance_ratio_transform(
        np.ascontiguousarray(np.asarray(k2_amp, dtype=np.float64)).ravel(),
        np.ascontiguousarray(np.asarray(k2_val, dtype=np.float64)).ravel(),
        int(n_bins)))
    nb = int(n_bins)
    return out[:nb], out[nb:2 * nb], float(out[2 * nb])

def convolve_distance_with_k2_ratio(
    r_da: np.ndarray,
    amp_r_da: np.ndarray,
    r_ratio: np.ndarray,
    weights_ratio: np.ndarray,
    n_bins: int = 256,
) -> tuple:
    """Convolve distance distribution with κ² ratio distribution.
    
    Computes R_app = R_DA × (R_app/R_DA) via multiplicative convolution.
    
    Parameters
    ----------
    r_da : array
        R_DA distance values
    amp_r_da : array
        Amplitudes/weights for R_DA distribution
    r_ratio : array
        R_app/R_DA ratio values
    weights_ratio : array
        Weights for ratio distribution
    n_bins : int
        Number of bins for output histogram
    Returns
    -------
    tuple of (r_app, amp_app)
        r_app : array
            Apparent distance values
        amp_app : array
            Amplitudes for apparent distance distribution
    """
    r_da = np.asarray(r_da, dtype=np.float64).ravel()
    amp_r_da = np.asarray(amp_r_da, dtype=np.float64).ravel()
    r_ratio = np.asarray(r_ratio, dtype=np.float64).ravel()
    weights_ratio = np.asarray(weights_ratio, dtype=np.float64).ravel()
    if r_da.size == 0 or r_ratio.size == 0:
        return np.empty(0), np.empty(0)

    # One path, not two, and no flag to pick between them. There used to be a
    # `use_fast` branch calling
    # `_fast_convolve_loop`, which **was never defined in this package** -- it
    # came across in the kappa-squared migration as a call to something that
    # stayed behind, so `use_fast=True` raised NameError for every input large
    # enough to take it. The slow branch built the full outer product instead:
    # n x m twice over, once for the products and once for the weights, purely
    # to hand them to np.histogram.
    #
    # The kernel does the binning without forming either array, so there is
    # nothing left for a flag to choose between.
    lo = float(r_da.min() * r_ratio.min())
    hi = float(r_da.max() * r_ratio.max())
    if not hi > lo:
        # Every product identical -- a delta. numpy widens a zero-width range to
        # (a - 0.5, a + 0.5) rather than returning nothing, and the contract
        # here is the one numpy set, so it is matched rather than reinvented.
        lo, hi = lo - 0.5, hi + 0.5
    hist = np.asarray(IMP.bff.outer_product_histogram(
        np.ascontiguousarray(r_da), np.ascontiguousarray(amp_r_da),
        np.ascontiguousarray(r_ratio), np.ascontiguousarray(weights_ratio),
        int(n_bins), lo, hi), dtype=np.float64)

    edges = np.linspace(lo, hi, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    peak = hist.max() if hist.size else 0.0
    mask = hist > 1e-10 * peak
    return centers[mask], hist[mask]


# --------------------------------------------------------------------------
# terms
# --------------------------------------------------------------------------
"""Interaction terms: the processes that deactivate a dye, and their rate constants.

Shaped the way a force field's terms are: a functional form, an **arity**, and
parameters looked up by type. Each term answers one question -- *given the
states of the participants, what is the rate constant per state* -- and the
answers **add**, because the channels are parallel. That additivity is the whole
reason these are objects rather than numbers computed inside an observable: the
field solver already relies on summing a quenching map and a FRET map, and
burying either inside its own function is what left ``fret_rate_trace`` and
``fret_rate_map`` as two incompatible calls.

Terms are **representation-agnostic**. They consume
:class:`~IMP.bff.representation.States` -- positions, weights, orientations -- so
one implementation serves an accessible volume, a rotamer library, a
coarse-grained model and an MD trajectory alike. A term that needs orientations
says so, and a representation that cannot supply them (an AV) is told, rather
than silently averaged.

Arity is what distinguishes state-dependent from emergent:

* **1-body** -- :class:`RadiativeTerm`. Depends on the dye and the solvent only.
* **2-body** -- :class:`PETTerm` (a dye against quencher atoms),
  :class:`FRETTerm` (a dye against a dye).
* **N-body** -- homo-FRET and multi-chromophore transfer, where the rate is not
  a sum over pairs. Not implemented; the interface is shaped to take it.

These currently delegate to the kernels in :mod:`IMP.bff.quenching`, and the
tests assert they reproduce them exactly. The kernels move underneath the terms
in a later stage; wrapping first and moving second is what keeps the change
reviewable.
"""

class InteractionTerm(abc.ABC):
    """One deactivation channel.

    :cvar arity: how many participants the term is defined over. 1 for a
        process of the dye alone, 2 for a pair, more for a genuinely many-body
        process.
    """

    arity: int = 1
    #: True when the term needs transition dipoles, so a caller can tell whether
    #: a representation without them (an accessible volume) forces an isotropic
    #: assumption instead of resolving the orientation.
    needs_orientations: bool = False

    @abc.abstractmethod
    def rate_constants(self, *participants, **kwargs) -> np.ndarray:
        """Rate constants in 1/ns, one per state of the first participant."""

    def __repr__(self) -> str:
        return f"{type(self).__name__}(arity={self.arity})"


@dataclass(repr=False)
class RadiativeTerm(InteractionTerm):
    """1-body: the dye's own decay, ``1/tau0``.

    The floor every other channel adds to. Depends on the dye and its medium,
    not on where it is or what is near it -- which is exactly what "1-body"
    means here.
    """

    lifetime: float
    arity: int = field(default=1, init=False)

    def rate_constants(self, states, **_) -> np.ndarray:
        if self.lifetime <= 0.0:
            raise ValueError(f"lifetime must be > 0, not {self.lifetime}")
        return np.full(int(states.n_points), 1.0 / float(self.lifetime),
                       dtype=np.float64)


@dataclass(repr=False)
class PETTerm(InteractionTerm):
    """2-body: photoinduced electron transfer between a dye and quencher atoms.

    ``k(r) = sum_a kQ_a * exp(-(|r - r_a| - r_dye) / rC_a)``

    The rate is a **pair** property -- ``kQ`` depends on the redox potentials of
    both partners -- so the parameters come in keyed by ``(dye, comp_id)`` as
    :class:`~IMP.bff.label.PETParameters` rather than as a table of the
    quencher alone.

    :param parameters: ``{comp_id: PETParameters}`` for **one** dye.
    :param dye_radius: subtracted from the centre-to-centre distance, because
        the tabulated contact distances are measured from the dye *surface*.
    """

    parameters: Mapping[str, "PETParameters"]
    dye_radius: float = 3.5
    arity: int = field(default=2, init=False)

    def rate_constants(self, states, atoms, **_) -> np.ndarray:
        """:param atoms: the structured atom array ``IMP.bff.quenching`` uses."""
        import IMP.bff.quenching.model as maps

        table = {
            comp_id: {
                atom: (p.rate_constant,
                       p.attenuation_length if p.attenuation_length else 1.0)
                for atom in _quencher_atoms(comp_id)
            }
            for comp_id, p in self.parameters.items()
        }
        kQ, rC = maps.atomic_quenching_parameters(atoms, table)
        xyz = np.ascontiguousarray(np.asarray(atoms["coord"], dtype=np.float64))
        positions = np.ascontiguousarray(
            np.asarray(states.positions, dtype=np.float64))
        return _pet_rates(positions, xyz, kQ, rC, float(self.dye_radius))


@dataclass(repr=False)
class FRETTerm(InteractionTerm):
    """2-body: Förster transfer between two dyes.

    ``k(r) = (1/tau0) * (R0/r)^6 * kappa^2 / (2/3)``

    **R0 is derived, not supplied** -- from the two dyes' spectra, the donor's
    quantum yield, the medium's refractive index and kappa^2. Passing it in was
    how ``forster_radius=52.0`` came to be a default in a dozen signatures.

    :param donor, acceptor: :class:`~IMP.bff.dye.Dye` species.
    :param refractive_index: of the medium between them.
    :param kappa2: orientation factor. ``None`` resolves it from the
        participants' orientations when both have them, and falls back to the
        isotropic 2/3 otherwise -- reporting which, through
        :attr:`used_isotropic_kappa2`.
    """

    donor: "Dye"
    acceptor: "Dye"
    refractive_index: float = 1.4
    kappa2: Optional[float] = None
    arity: int = field(default=2, init=False)
    needs_orientations: bool = field(default=True, init=False)

    def __post_init__(self):
        self.used_isotropic_kappa2 = self.kappa2 is None

    @property
    def forster_radius(self) -> float:
        """R0 in Angstrom, derived from the pair and the medium."""
        from IMP.bff.dye import forster_radius
        k2 = 2.0 / 3.0 if self.kappa2 is None else float(self.kappa2)
        return 10.0 * forster_radius(
            self.donor, self.acceptor, k2, self.refractive_index)

    def rate_constants(self, donor_states, acceptor_states, r_min: float = 7.0,
                       **_) -> np.ndarray:
        from IMP.bff.quenching.fret_trace import fret_rate_trace

        tau0 = self.donor.lifetime
        if tau0 is None or tau0 <= 0.0:
            raise ValueError(
                f"{self.donor.name} has no lifetime, so a FRET rate cannot be "
                "expressed as 1/tau0 * (R0/r)^6")
        return fret_rate_trace(
            np.ascontiguousarray(np.asarray(donor_states.positions, dtype=np.float64)),
            np.ascontiguousarray(np.asarray(acceptor_states.positions, dtype=np.float64)),
            R0=self.forster_radius, tau0=float(tau0), r_min=r_min,
            kappa2=self.kappa2,
        )


def total_rate(terms: Sequence[InteractionTerm], *participants, **kwargs) -> np.ndarray:
    """Sum the channels.

    Parallel deactivation channels add, which is the property that makes a
    quenching map and a FRET map summable on the same grid. Doing it here rather
    than inside each observable is what keeps that true when a channel is added.
    """
    if not terms:
        raise ValueError("no interaction terms to sum")
    out = None
    for term in terms:
        r = np.asarray(term.rate_constants(*participants[:term.arity], **kwargs),
                       dtype=np.float64)
        out = r.copy() if out is None else out + r
    return out


def _quencher_atoms(comp_id: str):
    from IMP.bff.quenching.pet import QUENCHER_ATOMS
    return QUENCHER_ATOMS.get(comp_id, ())


def _pet_rates(positions, atoms_xyz, kQ, rC, dye_radius) -> np.ndarray:
    """Per-state PET rate, the same law ``quenching_rate_map`` stamps on a grid.

    Written over *states* rather than voxels so it serves a rotamer library and
    an MD trajectory as well as an accessible volume.
    """
    active = (kQ > 0.0) & (rC > 0.0)
    if not np.any(active):
        return np.zeros(positions.shape[0], dtype=np.float64)
    xyz, k, r = atoms_xyz[active], kQ[active], rC[active]
    out = np.zeros(positions.shape[0], dtype=np.float64)
    for i in range(positions.shape[0]):
        d = np.sqrt(np.sum((xyz - positions[i]) ** 2, axis=1)) - dye_radius
        out[i] = float(np.sum(k * np.exp(-d / r)))
    return out
