"""Wall time of the accessible-volume solver, per volume.

The gate for the PathMap work (step 2 of the independent-core plan): that step
replaces `PathMap`'s `IMP::em::SampledDensityMap` base with a bff-owned grid,
which rewrites the AV hot path. The standing rule is that performance only
improves, and until now there was nothing to measure it against --
`benchmark_av_screening.py` times `ProbeNetworkRestraint`, not the solver.

Two doors are timed because they cost differently and can regress
independently:

* **structure** -- `get_av_from_structure`, which parses the PDB, applies
  the strip mask and resamples. This is what a caller with a file pays.
* **array** -- `get_av`, handed the obstacle array directly. This isolates
  the path search from the parsing, so a regression can be attributed.

Grid resolution dominates: cost goes as the voxel count, so a 0.5 A step is
roughly 27x the work of 1.5 A. The sweep exists so a change that is neutral on
a coarse grid and quadratic on a fine one cannot hide.

    python benchmark/benchmark_av.py                     # default sweep
    python benchmark/benchmark_av.py --repeats 7 --json baseline.json
    python benchmark/benchmark_av.py --compare baseline.json

Baselines are recorded in `okf/validation/` with the commit they were taken
at; a bare number with no commit is not a baseline.
"""
import argparse
import json
import os
import statistics
import time

import IMP.bff as bff

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# (label, pdb, chain, resseq, atom, linker_length, linker_width, radius)
SITES = [
    ("T4L-132", "examples/structure/T4L/3GUN.pdb", "", 132, "CB", 22.0, 3.5, 3.5),
    ("T4L-150", "examples/structure/T4L/3GUN.pdb", "", 150, "CB", 22.0, 2.5, 3.5),
    ("GBP-400", "examples/structure/GBP/mGBP2A.pdb", "A", 400, "CB", 22.0, 4.5, 3.5),
]
STEPS = [1.5, 1.0, 0.5]


def timeit(fn, repeats):
    """Median of `repeats` runs, after one untimed warm-up.

    Median rather than mean: the first timed run on a cold page cache is a
    different measurement from the rest, and one outlier should not move the
    number a regression gate reads.
    """
    fn()
    ts = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        ts.append((time.perf_counter() - t0) * 1e3)
    return {
        "median_ms": statistics.median(ts),
        "min_ms": min(ts),
        "n": len(ts),
    }


def obstacles_for(pdb):
    """The (N, 4) x/y/z/radius array the array door takes."""
    recs = bff.read_pdb_records(pdb)
    import numpy as np
    a = np.empty((len(recs), 4), dtype=np.float64)
    for i, r in enumerate(recs):
        a[i] = (r.x, r.y, r.z, r.vdw_radius)
    return a


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--steps", type=float, nargs="*", default=STEPS)
    ap.add_argument("--json", help="write results here")
    ap.add_argument("--compare", help="read a baseline and print deltas")
    ap.add_argument("--threshold", type=float, default=1.08,
                    help="ratio above which the gate fails (default 1.08)")
    ap.add_argument("--noise", action="store_true",
                    help="run the sweep twice and report the self-ratio; the "
                         "gate is only trustworthy if this is under --threshold")
    args = ap.parse_args()

    if args.noise:
        # A performance gate is only as good as the machine it runs on. This
        # measures the floor directly: the same binary against itself. If the
        # self-ratio exceeds the threshold, the box is too busy and any verdict
        # from --compare is noise. Observed on a loaded machine (load 14):
        # 1.19x. Take baselines when this reads close to 1.00x.
        a = run_sweep(args)
        b = run_sweep(args)
        worst = max(b[k][d]["min_ms"] / a[k][d]["min_ms"]
                    for k in a for d in ("structure", "array")
                    if a[k][d]["min_ms"])
        print("\nnoise floor (same binary, twice): %.2fx  -> %s"
              % (worst, "TOO NOISY TO GATE" if worst > args.threshold else "quiet enough"))
        return 1 if worst > args.threshold else 0

    out = run_sweep(args)
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=2, sort_keys=True)
        print("\nwrote %s" % args.json)
    return compare(out, args)


def run_sweep(args):
    out = {}
    for label, rel, chain, resseq, atom, ll, lw, r1 in SITES:
        pdb = os.path.join(ROOT, rel)
        if not os.path.exists(pdb):
            print("  skip %s (no %s)" % (label, rel))
            continue
        arr = obstacles_for(pdb)
        # the attachment coordinate, so the array door measures the same volume
        try:
            src = bff.get_attachment_point(pdb, chain, resseq, atom)
        except Exception as exc:
            print("  skip %s (%s)" % (label, exc))
            continue
        for step in args.steps:
            k = "%s@%.1f" % (label, step)
            s = timeit(lambda: bff.get_av_from_structure(
                pdb, chain, resseq, atom, ll, lw, r1, 0.0, 0.0, step), args.repeats)
            a = timeit(lambda: bff.get_av(
                arr, list(src), ll, lw, r1, 0.0, 0.0, step), args.repeats)
            out[k] = {"structure": s, "array": a, "n_atoms": int(arr.shape[0])}
            print("%-18s structure %8.2f ms   array %8.2f ms   (%d atoms)"
                  % (k, s["median_ms"], a["median_ms"], arr.shape[0]))
    return out


def compare(out, args):
    if args.compare:
        base = json.load(open(args.compare))
        # Compared on the MINIMUM, not the median. The minimum is the least
        # noisy estimate of how fast the code can run: it is the sample least
        # contaminated by scheduler preemption and page faults, which only ever
        # add time. Measured here, a self-comparison at --repeats 3 swung 1.07x
        # on medians -- a gate that fires on its own noise is worse than none.
        print("\n%-20s %10s %10s %8s" % ("case/door", "base ms", "now ms", "ratio"))
        worst, worst_k = 0.0, ""
        for k in sorted(set(base) & set(out)):
            for door in ("structure", "array"):
                b = base[k][door]["min_ms"]
                n = out[k][door]["min_ms"]
                ratio = n / b if b else float("nan")
                if ratio > worst:
                    worst, worst_k = ratio, k + "/" + door
                print("%-20s %10.2f %10.2f %8.2fx  %s"
                      % (k + "/" + door[:4], b, n, ratio,
                         "SLOWER" if ratio > args.threshold else ""))
        bad = worst > args.threshold
        print("\nworst %.2fx at %s (threshold %.2fx)  -> %s"
              % (worst, worst_k, args.threshold, "REGRESSION" if bad else "ok"))
        return 1 if bad else 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
