/**
 * \file IMP/bff/internal/FftLength.h
 * \brief Choosing a transform length that the FFT is actually fast at.
 *
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_FFTLENGTH_H
#define IMPBFF_FFTLENGTH_H

#include <IMP/bff/IMPCompatibility.h>
#include <cstddef>

IMPBFF_BEGIN_NAMESPACE

//! \name FFT lengths
//! @{

/**
 * \brief The smallest length `>= n_min` whose only prime factors are 2, 3 and
 *        5 (a "5-smooth" or regular number).
 *
 * **Why it is worth a function.** A mixed-radix FFT is fast at lengths that
 * factor into small primes and falls back to Bluestein's algorithm otherwise,
 * which is several times slower. The cost is not marginal: measured with
 * pocketfft on an M-series Mac, four rows of length **1563 = 3 x 521** took
 * **0.1245 ms** and four rows of length **1600** took **0.0174 ms** -- seven
 * times faster for a transform two per cent longer. A signal whose natural
 * length happens to carry a large prime factor pays that every evaluation.
 *
 * So: pad. For a linear convolution the padding is required anyway, and for a
 * periodic one -- where the natural length is the period and looks
 * mandatory -- a circular shift or a wrapped convolution can be computed as a
 * linear one at a padded length and folded back, which is exactly the trick
 * this makes cheap.
 *
 * 7 is deliberately excluded even though many libraries handle it: the gain
 * over the next 5-smooth length is small and the guarantee is simpler to state.
 */
inline std::size_t good_fft_length(std::size_t n_min) {
  if (n_min <= 1) return 1;
  for (std::size_t n = n_min;; ++n) {
    std::size_t m = n;
    for (std::size_t p : {std::size_t(2), std::size_t(3), std::size_t(5)})
      while (m % p == 0) m /= p;
    if (m == 1) return n;
  }
}

//! @}

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_FFTLENGTH_H
