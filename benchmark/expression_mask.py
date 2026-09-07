"""How much a gate costs when its answer is bytes instead of doubles.

A selection asks a yes/no question of every row, and the answer to it is
one bit. Carrying that back as an array of doubles spends eight bytes a
row to say one bit, and at the row counts a burst table reaches -- 1e5 to
1e7 -- that write is a real fraction of the work. This measures the
byte-mask path against the double path, against pandas' ``DataFrame.eval``
(the thing it replaces), and against numpy.

Run with the bff build on PYTHONPATH::

    PYTHONPATH=~/dev/imp/cmake-build-arm64/lib python benchmark/expression_mask.py

Timings are min-of-N, because the interesting quantity is the floor and on
a loaded machine the mean measures the other processes instead. The load
average is printed with the results for exactly that reason -- above about
2 the numbers are not worth quoting.
"""

import os
import time

import numpy as np

from IMP import bff

REPEATS = 7
ROW_COUNTS = [100_000, 1_000_000, 5_000_000]

# Gates of the shape a burst selection actually uses: one bound, a window,
# several terms, and one with real arithmetic inside it. Each carries its
# numpy twin rather than a syntax rewrite, so what numpy is timed on is
# visible rather than inferred.
GATES = [
    ("one bound",    "x > 0.3",
     lambda c: c["x"] > 0.3),
    ("a window",     "x > 0.3 and x < 0.7",
     lambda c: (c["x"] > 0.3) & (c["x"] < 0.7)),
    ("three terms",  "x > 0.3 and y < 500 and z > 1",
     lambda c: (c["x"] > 0.3) & (c["y"] < 500) & (c["z"] > 1)),
    ("arithmetic",   "(x-0.5)*(x-0.5) + z*z < 0.25",
     lambda c: (c["x"] - 0.5) * (c["x"] - 0.5) + c["z"] * c["z"] < 0.25),
]


def best(fn):
    """Min-of-REPEATS wall time, in milliseconds."""
    times = []
    for _ in range(REPEATS):
        t0 = time.perf_counter()
        fn()
        times.append(time.perf_counter() - t0)
    return min(times) * 1e3


def columns(n, seed=7):
    rng = np.random.default_rng(seed)
    return {"x": rng.uniform(0, 1, n),
            "y": rng.uniform(0, 1000, n),
            "z": rng.normal(0, 1, n)}


def main():
    try:
        import pandas as pd
    except ImportError:
        pd = None

    print(f"load average {os.getloadavg()[0]:.1f}   min of {REPEATS}   "
          f"times in ms\n")
    header = (f"{'gate':<14}{'rows':>10}{'mask':>9}{'double':>9}{'numpy':>9}"
              f"{'pandas':>9}{'vs dbl':>8}{'vs np':>8}{'vs pd':>8}")
    print(header)
    print("-" * len(header))

    for n in ROW_COUNTS:
        cols = columns(n)
        names = list(cols)
        block = np.ascontiguousarray(np.vstack([cols[k] for k in names]))
        frame = pd.DataFrame(cols) if pd is not None else None

        for label, text, twin in GATES:
            ex = bff.Expression("gate")
            ex.set_expression(text)

            # Agreement first: a benchmark of a wrong answer is worthless.
            assert np.array_equal(ex.compute_mask(names, block).astype(bool),
                                  twin(cols)), text

            t_mask = best(lambda: ex.compute_mask(names, block))
            t_double = best(lambda: ex.compute_columns(names, block))
            t_numpy = best(lambda: twin(cols))
            t_pandas = (best(lambda: frame.eval(text))
                        if frame is not None else float("nan"))

            print(f"{label:<14}{n:>10,}{t_mask:>9.2f}{t_double:>9.2f}"
                  f"{t_numpy:>9.2f}{t_pandas:>9.2f}"
                  f"{t_double / t_mask:>7.2f}x{t_numpy / t_mask:>7.2f}x"
                  f"{t_pandas / t_mask:>7.2f}x")
        print()

    print("gates:")
    for label, text, _ in GATES:
        print(f"  {label:<14}{text}")


if __name__ == "__main__":
    main()
