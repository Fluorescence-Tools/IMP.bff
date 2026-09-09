/**
 * \file IMP/bff/Rebin.h
 * \brief Rebinning as the partition it usually is, with the general matrix as
 *        the fallback.
 *
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_REBIN_H
#define IMPBFF_REBIN_H

#include <cmath>
#include <cstddef>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! \name Rebinning
//! @{

/**
 * \brief A linear map from a fine axis to a coarse one, applied as segment
 *        sums when it is a partition and as a matrix product when it is not.
 *
 * **Why this exists.** A rebinning is almost always "each coarse bin is the
 * sum of a contiguous run of fine bins": every weight one, the runs disjoint
 * and covering. Written as a dense `(n_coarse, n_fine)` matrix product that is
 * `n_coarse * n_fine` multiply-adds; written as what it is, it is `n_fine`
 * additions. In the model this came from, 226 by 1563, that is 353 thousand
 * products against 1563 sums -- and it was **77 % of an evaluation** before
 * anyone looked, because a matrix product is what a rebinning looks like when
 * it is written down.
 *
 * **The structure is discovered, not assumed.** `Rebin` inspects the matrix
 * once: if every nonzero is exactly one, each row's nonzeros are contiguous,
 * and the rows are disjoint, it stores runs; otherwise it keeps the matrix.
 * So a smoothing rebin, an overlapping one, or one with fractional edge
 * weights still gives the right answer, only at the matrix's price. A
 * component that were merely fast would be a trap the first time someone
 * passed it a real interpolation.
 */
class Rebin {
 public:
  //! \param R row-major `(n_coarse, n_fine)`.
  Rebin(const double* R, std::size_t n_coarse, std::size_t n_fine)
      : n_coarse_(n_coarse), n_fine_(n_fine) {
    beg_.assign(n_coarse, 0);
    len_.assign(n_coarse, 0);
    std::size_t covered = 0;
    for (std::size_t b = 0; b < n_coarse && partition_; ++b) {
      const double* row = &R[b * n_fine];
      std::size_t first = n_fine, last = 0, cnt = 0;
      for (std::size_t i = 0; i < n_fine; ++i) {
        if (row[i] == 0.0) continue;
        if (row[i] != 1.0) { partition_ = false; break; }
        if (cnt == 0) first = i;
        last = i;
        ++cnt;
      }
      if (!partition_) break;
      if (cnt && last - first + 1 != cnt) { partition_ = false; break; }
      beg_[b] = cnt ? first : 0;
      len_[b] = cnt;
      covered += cnt;
    }
    if (partition_ && covered > n_fine) partition_ = false;    // rows overlap
    if (!partition_) dense_.assign(R, R + n_coarse * n_fine);
  }

  //! True when the map was recognised as a partition and is applied as sums.
  bool is_partition() const { return partition_; }
  std::size_t n_coarse() const { return n_coarse_; }
  std::size_t n_fine() const { return n_fine_; }

  //! `out = R * in`. `in` has `n_fine` entries, `out` has `n_coarse`.
  template <typename T>
  void apply(const T* in, T* out) const {
    if (partition_) {
      for (std::size_t b = 0; b < n_coarse_; ++b) {
        T s = T(0.0);
        const T* p = in + beg_[b];
        for (std::size_t i = 0; i < len_[b]; ++i) s += p[i];
        out[b] = s;
      }
    } else {
      for (std::size_t b = 0; b < n_coarse_; ++b) {
        T s = T(0.0);
        const double* row = &dense_[b * n_fine_];
        for (std::size_t i = 0; i < n_fine_; ++i)
          if (row[i] != 0.0) s += T(row[i]) * in[i];
        out[b] = s;
      }
    }
  }

 private:
  std::size_t n_coarse_, n_fine_;
  bool partition_ = true;
  std::vector<std::size_t> beg_, len_;
  std::vector<double> dense_;
};

//! @}

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_REBIN_H
