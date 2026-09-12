/**
 * \file IMP/bff/PhotophysicsPolarisation.h
 * \brief What a polariser sees: mixing an isotropic signal and an anisotropy
 *        into parallel and perpendicular channels, and taking them apart again.
 *
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_PHOTOPHYSICSPOLARISATION_H
#define IMPBFF_PHOTOPHYSICSPOLARISATION_H

#include <IMP/bff/IMPCompatibility.h>
#include <cmath>
#include <cstddef>

IMPBFF_BEGIN_NAMESPACE

//! \name Polarised detection
//! @{

//! Which polarised channel a signal is detected in.
enum PolarisedChannel {
  POL_VV = 0,   //!< vertical excitation, vertical detection ("parallel")
  POL_VH,       //!< vertical excitation, horizontal detection ("perpendicular")
  POL_MAGIC     //!< the magic angle: anisotropy-free by construction
};

/**
 * \brief The weight the anisotropy enters a channel with.
 *
 * For ideal polarisers the parallel channel carries `1 + 2 r(t)` and the
 * perpendicular `1 - r(t)`, both times a third of the total intensity. Real
 * optics mix the two, and the mixing is described by `l1` and `l2` -- the
 * fractions of the wrong polarisation each detection arm admits -- which turn
 * the factors into `2 - 3 l1` and `-1 + 3 l2`. At `l1 = l2 = 0` they are the
 * textbook 2 and -1.
 *
 * The magic angle carries zero by construction: it is the polariser setting
 * at which the anisotropy cancels, so that the decay of the total intensity can
 * be measured without knowing the rotation at all.
 *
 * These are optics, not FRET. The same two numbers describe a
 * polarisation-resolved FCS curve, an anisotropy image and a time-resolved
 * anisotropy decay.
 */
template <typename T = double>
inline T anisotropy_weight(PolarisedChannel channel, const T& l1, const T& l2) {
  switch (channel) {
    case POL_VV: return T(2.0) - T(3.0) * l1;
    case POL_VH: return T(-1.0) + T(3.0) * l2;
    default: return T(0.0);
  }
}

/**
 * \brief Mix an isotropic signal and its anisotropy-weighted partner into one
 *        polarised channel.
 *
 * `out[i] = (iso[i] + weight * r0 * aniso[i]) / g_channel`, where `g_channel`
 * is the `g` factor for the perpendicular arm and one for the parallel one.
 *
 * **The `g` factor is a property of the instrument, not of the sample**: the
 * ratio of the two detection arms' efficiencies. It is measured once, on a
 * sample with no anisotropy or by tail matching, and it multiplies rather than
 * shifts -- which is why getting it wrong looks like a rotational correlation
 * time that is wrong by a constant factor rather than like a bad fit.
 */
template <typename T = double>
inline void mix_polarised(const T* iso, const T* aniso, std::size_t n,
                          const T& r0, const T& weight, const T& g_channel, T* out) {
  for (std::size_t i = 0; i < n; ++i)
    out[i] = (iso[i] + weight * r0 * aniso[i]) / g_channel;
}

/**
 * \brief Remove an additive contribution from a channel before inverting.
 *
 * **The inverse is a ratio of differences, so anything additive biases it.**
 * Uncorrelated background -- dark counts, room light, afterpulsing -- and
 * scattered excitation light both add to `VV` and `VH` and neither belongs to
 * the fluorescence. Leaving them in pulls the recovered anisotropy toward the
 * anisotropy of the contaminant, which is zero for background and nearly one
 * for scatter, so the bias does not even have a consistent sign.
 *
 * Scatter deserves its own subtraction rather than being lumped into the
 * background: it is the excitation pulse's own polarisation, so it is close to
 * fully polarised and enters the two channels with the ideal weights, not the
 * sample's. Treating it as isotropic is the mistake that makes a rotational
 * correlation time come out short.
 */
template <typename T = double>
inline void subtract_additive(const T* channel, const T* additive, std::size_t n, T* out) {
  for (std::size_t i = 0; i < n; ++i) out[i] = channel[i] - additive[i];
}

/**
 * \brief The isotropic (anisotropy-free) signal from two corrected channels.
 *
 * With `VV = S + w_vv R` and `g VH = S + w_vh R`, eliminating `R` gives
 *
 *     S = (w_vv * g * VH - w_vh * VV) / (w_vv - w_vh)
 *
 * For ideal polarisers `w_vv = 2` and `w_vh = -1`, and this is the familiar
 * magic-angle combination `(VV + 2 g VH) / 3`. **For real optics it is not**:
 * using the ideal weights on data taken with `l1` and `l2` leaves a residual
 * anisotropy in what is supposed to be the anisotropy-free signal.
 *
 * The channels must already have background and scatter removed.
 */
template <typename T = double>
inline void isotropic_from_polarised(const T* vv, const T* vh, std::size_t n, const T& g,
                                     const T& w_vv, const T& w_vh, T* out) {
  const T den = w_vv - w_vh;
  for (std::size_t i = 0; i < n; ++i) out[i] = (w_vv * g * vh[i] - w_vh * vv[i]) / den;
}

//! The ideal-polariser case of \ref isotropic_from_polarised: `(VV + 2 g VH) / 3`.
template <typename T = double>
inline void magic_angle(const T* vv, const T* vh, std::size_t n, const T& g, T* out) {
  isotropic_from_polarised(vv, vh, n, g, T(2.0), T(-1.0), out);
}

/**
 * \brief The anisotropy from two corrected channels, for real optics.
 *
 *     r = (VV - g VH) / (w_vv * g * VH - w_vh * VV)
 *
 * -- the difference over the isotropic signal. At `w_vv = 2`, `w_vh = -1` this
 * is the textbook `(VV - g VH) / (VV + 2 g VH)`.
 *
 * Note what it does NOT need: the total intensity, the quantum yield, the
 * concentration, or the excitation power. They divide out, which is why
 * anisotropy is measurable on a sample whose brightness is unknown -- and why
 * the corrections it DOES need are the ones worth insisting on. There are
 * three: `g`, because it multiplies one channel and not the other; `l1` and
 * `l2` through the weights, because ideal weights on real optics bias `r`
 * toward zero; and the additive terms, which must be removed first with
 * \ref subtract_additive.
 *
 * **This is tttrlib's `fit23` expression**, reached from the other direction.
 * `DecayFit23.cpp`'s `decay23_r_ad` writes the denominator as
 * `Fp (1 - 3 l2) + (2 - 3 l1) g Fs`; expanding `w_vv` and `w_vh` above gives
 * `(1 - 3 l2) VV + (2 - 3 l1) g VH`, the same thing. The two agree to 1.2e-16
 * in this project's test, which is worth more than either derivation alone:
 * they are independent implementations and they do not share an error.
 *
 * `fit23` corrects its background as `(S - gamma B) / (1 - gamma)` with
 * `gamma` the background FRACTION rather than by subtracting counts. For the
 * anisotropy the difference is nothing: the `1 / (1 - gamma)` multiplies both
 * channels and cancels in the ratio. It matters for `Fp` and `Fs` as
 * intensities, which is why that form is the one to use when the isotropic
 * signal is wanted on an absolute scale.
 *
 * \param floor denominators at or below this give zero rather than a division
 */
template <typename T = double>
inline void anisotropy_from_polarised(const T* vv, const T* vh, std::size_t n,
                                      const T& g, const T& w_vv, const T& w_vh,
                                      T* out, double floor = 1e-300) {
  for (std::size_t i = 0; i < n; ++i) {
    //  (w_vv - w_vh) divides out of numerator and denominator alike, so this
    //  is r = R / S written without forming either
    const T scaled_iso = w_vv * g * vh[i] - w_vh * vv[i];
    out[i] = (scaled_iso > T(floor)) ? (vv[i] - g * vh[i]) / scaled_iso : T(0.0);
  }
}

//! The ideal-polariser case of \ref anisotropy_from_polarised.
template <typename T = double>
inline void anisotropy_from_polarised(const T* vv, const T* vh, std::size_t n,
                                      const T& g, T* out, double floor = 1e-300) {
  anisotropy_from_polarised(vv, vh, n, g, T(2.0), T(-1.0), out, floor);
}

//! @}

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_PHOTOPHYSICSPOLARISATION_H
