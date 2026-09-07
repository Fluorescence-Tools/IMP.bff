"""A/B two builds on a machine that is never idle.

The usual way to gate a performance change -- measure before, change, measure
after -- assumes the two measurements are comparable. On a shared box they are
not: this one runs several agent sessions at once, and `benchmark_av.py
--noise` has read 1.12x comparing a build with *itself*. A change measured that
way is indistinguishable from the load that happened to be running.

Paired measurement removes the assumption instead of waiting for quiet. Both
builds are measured in the same round, alternately, so a burst of load lands on
both; the statistic is the per-round *ratio*, and the report is the median of
those ratios. Absolute times still wander -- the ratio does not.

Each build is a self-contained directory holding `IMP/`, `_IMP_bff.so` and
`libimp_bff.0.dylib` with its install names rewritten to point inside itself
(see `make_variant` in the plan notes), so neither run touches the shared build
tree and no lock is needed to measure.

    python benchmark/paired_ab.py --a /tmp/v_before --b /tmp/v_after --rounds 7

Reports, per case, median(t_b / t_a). Above 1 means B is slower.
"""
import argparse
import json
import os
import statistics
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def measure(variant, python, repeats):
    """Run the AV sweep inside `variant` and return its JSON."""
    out = os.path.join(tempfile.mkdtemp(), "r.json")
    env = dict(os.environ)
    env["PYTHONPATH"] = variant
    # The variant's own libimp_bff is found through the rewritten install name,
    # so DYLD only has to supply the IMP libraries both builds share.
    # stderr is captured rather than discarded: when the sweep fails the
    # reason is the only useful thing the harness has, and a driver that hides
    # it turns a one-line error into a debugging session.
    r = subprocess.run(
        [python, os.path.join(HERE, "benchmark_av.py"),
         "--repeats", str(repeats), "--json", out],
        cwd=ROOT, env=env, stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        tail = "\n".join(r.stderr.strip().splitlines()[-8:])
        raise SystemExit("sweep failed in %s (exit %d):\n%s"
                         % (variant, r.returncode, tail))
    with open(out) as fh:
        return json.load(fh)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--a", required=True, help="baseline build directory")
    ap.add_argument("--b", required=True, help="candidate build directory")
    # 11, measured rather than guessed: at load 13-18 a seven-round run
    # reported GBP-400@1.0/array as a 1.081x regression, slower in every one of
    # its rounds, and eleven rounds put the same case at 0.988. Seven is not
    # enough samples for the median to settle on a busy machine, and a gate
    # that flips verdict with sample count teaches people to re-roll until it
    # passes -- which is worse than no gate.
    ap.add_argument("--rounds", type=int, default=11)
    ap.add_argument("--repeats", type=int, default=3,
                    help="inner repeats per case per round")
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--threshold", type=float, default=1.05)
    args = ap.parse_args()

    ratios = {}
    for r in range(args.rounds):
        # A then B, adjacent in time. Alternating the order every other round
        # cancels any systematic advantage from going first (a warm cache, a
        # scheduler that has just freed a core).
        first, second = (args.a, args.b) if r % 2 == 0 else (args.b, args.a)
        m1 = measure(first, args.python, args.repeats)
        m2 = measure(second, args.python, args.repeats)
        ma, mb = (m1, m2) if r % 2 == 0 else (m2, m1)
        for k in sorted(set(ma) & set(mb)):
            for door in ("structure", "array"):
                ta, tb = ma[k][door]["min_ms"], mb[k][door]["min_ms"]
                if ta:
                    ratios.setdefault(k + "/" + door, []).append(tb / ta)
        print("round %d/%d done" % (r + 1, args.rounds), flush=True)

    # A case counts as a regression only when its median exceeds the threshold
    # AND its *best* round was still slower than the baseline. A real slowdown
    # moves the whole distribution; load moves the median around while leaving
    # some rounds faster. Without the second half of that test a self-
    # comparison at load 17 reports a 1.07x "regression" on one case out of
    # eighteen, which is how a useless gate teaches people to ignore it.
    print("\n%-24s %8s %8s %8s" % ("case/door", "median", "min", "max"))
    worst, worst_k, flagged = 0.0, "", []
    for k in sorted(ratios):
        v = ratios[k]
        med = statistics.median(v)
        bad = med > args.threshold and min(v) > 1.0
        if bad:
            flagged.append((k, med))
        if med > worst:
            worst, worst_k = med, k
        print("%-24s %8.3f %8.3f %8.3f %s"
              % (k, med, min(v), max(v), "SLOWER" if bad else ""))
    print("\nworst median %.3f at %s (threshold %.2f)" % (worst, worst_k, args.threshold))
    if flagged:
        print("REGRESSION in %d case(s): %s"
              % (len(flagged), ", ".join("%s %.3fx" % kv for kv in flagged)))
        return 1
    print("ok -- no case is slower in every round")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
