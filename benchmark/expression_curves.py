"""bff.GraphExpression against numpy, at the curve lengths that actually occur.

The engine exists to replace a Python ``eval`` in ChiSurf's parse models,
evaluated once per fit iteration over a curve of 100-4096 points. That is a
*small*-n problem, and small n is where a vector engine's fixed costs show:
what wins at a million rows can lose at five hundred.

The equations below are the shapes the catalogue actually contains -- a sum
of exponentials, a mixed transcendental, a near-trivial line, and the FCS
form. The FCS row used to be the loss case, at 0.26x numpy here: its
`**(-1)` parses as a *negated* constant, which the constant-power fold did
not recognise, so every point paid a std::pow() call where numpy rewrites
the same thing into a reciprocal. It is now folded at compile time, and the
row is 3.9x faster in absolute terms.

Scalars are passed as scalars, not as ``np.full(n, v)``. Materialising them
as columns makes the engine memcpy a whole array per operand, is not what a
real model does, and flattered numpy by ~30% in an earlier round of this.

    PYTHONPATH=~/dev/imp/cmake-build-arm64/lib python benchmark/expression_curves.py
"""

import os
import time

import numpy as np

from IMP import bff

REPEATS = 11
LENGTHS = [128, 512, 1024, 2048, 4096]

# (label, equation, numpy twin). `x` is the curve axis; everything else is a
# scalar parameter, exactly as a fit supplies them.
CASES = [
    ("two exponentials", "exp(-x/1.5)+0.2*exp(-x/4.0)",
     lambda x: np.exp(-x / 1.5) + 0.2 * np.exp(-x / 4.0)),
    ("mixed",            "sqrt(abs(x))+1/(x*x+1)",
     lambda x: np.sqrt(np.abs(x)) + 1 / (x * x + 1)),
    ("a line",           "0.3+2.0*x",
     lambda x: 0.3 + 2.0 * x),
    ("FCS",              "0.3+1/1.7*(1+x/1.2)**(-1)/sqrt(1+1/2.1**2*x/1.2)",
     lambda x: 0.3 + 1 / 1.7 * (1 + x / 1.2) ** (-1)
     / np.sqrt(1 + 1 / 2.1 ** 2 * x / 1.2)),
]


def best(fn):
    times = []
    for _ in range(REPEATS):
        t0 = time.perf_counter()
        fn()
        times.append(time.perf_counter() - t0)
    return min(times) * 1e6  # microseconds; these curves are short


def main():
    print(f"load average {os.getloadavg()[0]:.1f}   min of {REPEATS}   "
          f"times in us\n")
    header = f"{'equation':<20}" + "".join(f"{n:>16}" for n in LENGTHS)
    print(header)
    print("-" * len(header))

    for label, text, twin in CASES:
        cells = []
        for n in LENGTHS:
            x = np.linspace(0.01, 10.0, n)
            block = np.ascontiguousarray(x.reshape(1, n))
            ex = bff.GraphExpression("model")
            ex.set_expression(text)

            got = ex.compute_columns(["x"], block)
            np.testing.assert_allclose(got, twin(x), rtol=1e-12)

            t_ours = best(lambda: ex.compute_columns(["x"], block))
            t_numpy = best(lambda: twin(x))
            cells.append(f"{t_ours:>7.1f}{t_numpy / t_ours:>8.2f}x")
        print(f"{label:<20}" + "".join(cells))

    print("\ncolumns are: our time (us), then our speed as a multiple of "
          "numpy's\nequations:")
    for label, text, _ in CASES:
        print(f"  {label:<20}{text}")


if __name__ == "__main__":
    main()
