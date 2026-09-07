/**
 * \file SpecialFunctions.cpp
 * \brief Special functions the distance distributions need.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/SpecialFunctions.h>
#include <IMP/bff/internal/OutputView.h>

#include <cmath>

IMPBFF_BEGIN_NAMESPACE

double i0(double x) {
    // Abramowitz & Stegun 9.8.1 / 9.8.2, branching at |x| = 3.75. The
    // coefficients are transcribed digit for digit from the form ChiSurf uses,
    // because parity with it has to be exact rather than close.
    //
    // NOTE the first one: A&S and Numerical Recipes both print 3.5156229, and
    // ChiSurf has 3.5156299 -- two digits transposed, at some point long ago.
    // It is reproduced here ON PURPOSE. The contract of this port is "the same
    // answer, faster", and every published worm-like-chain fit in the stack was
    // made with 3.5156299; correcting it here would change those numbers
    // silently, which is the one outcome a port must not have. It costs about
    // 1e-6 relative near |x| = 3.75 -- slightly worse than the ~1e-7 the
    // approximation is otherwise good for. If it is ever fixed it should be
    // fixed in one place, deliberately, with the refits that implies.
    const double ax = std::fabs(x);
    if (ax < 3.75) {
        const double y = (ax / 3.75) * (ax / 3.75);
        return 1.0 + y * (3.5156299 + y * (3.0899424 + y * (1.2067492
                 + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
    }
    const double y = 3.75 / ax;
    return (std::exp(ax) / std::sqrt(ax))
           * (0.39894228 + y * (0.1328592e-1 + y * (0.225319e-2
              + y * (-0.157565e-2 + y * (0.916281e-2 + y * (-0.2057706e-1
              + y * (0.2635537e-1 + y * (-0.1647633e-1
              + y * 0.392377e-2))))))));
}

void i0_array(const std::vector<double>& x, double** out_view, int* n_out_view) {
    std::vector<double> out(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) out[i] = i0(x[i]);
    internal::copy_to_view(out, out_view, n_out_view);
}

IMPBFF_END_NAMESPACE
