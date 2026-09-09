/**
 * \file IMP/bff/Polarisation.h
 * \brief What a polariser sees: mixing an isotropic signal and an anisotropy
 *        into parallel and perpendicular channels, and taking them apart again.
 *
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_POLARISATION_H
#define IMPBFF_POLARISATION_H

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
 * \brief The magic-angle (total) intensity from the two polarised channels:
 *        `(I_vv + 2 g I_vh) / 3`.
 *
 * Anisotropy-free whatever the rotation is doing, which is the point: it is
 * the combination in which the polarisation cancels exactly, so a lifetime can
 * be fitted from polarised data without modelling the rotation.
 */
template <typename T = double>
inline void magic_angle(const T* vv, const T* vh, std::size_t n, const T& g, T* out) {
  for (std::size_t i = 0; i < n; ++i) out[i] = (vv[i] + T(2.0) * g * vh[i]) / T(3.0);
}

/**
 * \brief The anisotropy from the two polarised channels:
 *        `r = (I_vv - g I_vh) / (I_vv + 2 g I_vh)`.
 *
 * The difference over the total. Note what it does NOT need: the total
 * intensity, the quantum yield, the concentration, or the excitation power.
 * They divide out, which is why anisotropy is measurable on a sample whose
 * brightness is unknown -- and why a wrong `g` is the one instrument error it
 * cannot survive.
 */
template <typename T = double>
inline void anisotropy_from_polarised(const T* vv, const T* vh, std::size_t n,
                                      const T& g, T* out, double floor = 1e-300) {
  for (std::size_t i = 0; i < n; ++i) {
    const T tot = vv[i] + T(2.0) * g * vh[i];
    out[i] = (tot > T(floor)) ? (vv[i] - g * vh[i]) / tot : T(0.0);
  }
}

//! @}

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_POLARISATION_H
