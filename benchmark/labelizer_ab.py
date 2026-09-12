"""Regenerate the numbers in ``okf/validation/labelizer_ab.md``.

The native Labelizer replaces two external programs (DSSP, MSMS) and one banned
library (LabelLib). Whether that changes a number is a measurement, and this
prints it: per-parameter agreement against the reference implementation's own
published output for 1DDB, the depth distribution behind the solvent-exposure
tolerance, and the cost and consequence of the two-tier dye placement.

    python benchmark/labelizer_ab.py
    python benchmark/labelizer_ab.py --points 200 --probe 1.5
    python benchmark/labelizer_ab.py --refine 10

`test/label/test_labelizer_ab.py` pins the conclusions; this shows the
distributions they were drawn from, which is what a tolerance should be set
against. Recorded 2026-08-24 (arm64): `cr` and `ss` exact on 195/195, `se`
74.9 % bin-exact with +0.015 A bias and r = 0.976.
"""
import argparse
import csv
import os
import time

import numpy as np

import IMP.bff as bff

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "test", "input", "labelizer")

#: The reference's two-letter tag -> the dictionary `score_type` it becomes.
TAGS = {
    "cr": "cysteine_resemblance",
    "cs": "conservation",
    "se": "solvent_exposure",
    "ss": "secondary_structure",
}


def _path(name):
    return os.path.join(DATA, name)


def _reference(tag):
    with open(_path("1DDB-39_%s.csv" % tag)) as fh:
        return {row[0]: float(row[1]) for row in list(csv.reader(fh))[1:]}


def _agreement(args):
    """Per-parameter agreement against the reference's published scores."""
    options = bff.LabelizerOptions()
    options.probe_radius = args.probe
    options.n_sphere_points = args.points

    started = time.time()
    scores = bff.labelizer_score_structure(_path("1DDB-39.pdb"), bff.labelizer_model_paper(),
                                    options, _path("1DDB-39_cs.pdb"))
    elapsed = time.time() - started

    got = {}
    for row in scores:
        key = bff.labelizer_residue_key(row.asym_id, row.seq_id)
        got.setdefault(row.score_type, {})[key] = (
            row.value if row.status == "scored" else None)

    print("scored 195 residues x 6 parameters in %.2f s "
          "(probe %.1f A, %d sphere points)"
          % (elapsed, args.probe, args.points))
    print()
    print("%-4s %5s %12s %10s %10s" % ("tag", "n", "bin-exact", "mean|d|",
                                       "max|d|"))
    for tag in sorted(TAGS):
        reference = _reference(tag)
        ours, refs = [], []
        for key, value in reference.items():
            mine = got[TAGS[tag]].get(key)
            if mine is None:
                continue
            ours.append(mine)
            refs.append(value)
        delta = np.abs(np.asarray(ours) - np.asarray(refs))
        exact = int((delta < 1e-9).sum())
        print("%-4s %5d %7d/%-4d %10.5f %10.5f"
              % (tag, len(delta), exact, len(delta), delta.mean(), delta.max()))
    print()
    print("`cs` disagreeing by a constant is expected and is not a defect here:")
    print("  the shipped example writes its conservation output back over its")
    print("  own input, and the lookup has a two-cycle, so feeding it the")
    print("  shipped file yields one more iteration of that cycle -- the two")
    print("  values swapped, hence a constant |d| of %.5f. The lookup itself"
          % abs(2.3553071957924936 - 1.6322095472510827))
    print("  is exact to sixteen digits; see okf/validation/labelizer_ab.md.")
    return got


def _depth(args):
    """The distribution the solvent-exposure tolerance is set from.

    The reference depth is only recoverable to half a bin, because the
    published score is binned; the comparison is therefore against bin centres
    and part of the scatter below is that quantisation, not disagreement.
    """
    table = bff.labelizer_load_table("N_SE11_MEAN_SURFACE_DIST")
    centre = {round(table.values[i], 9): table.bins[i]
              for i in range(len(table.bins))}
    reference = _reference("se")

    structure = bff.labelizer_read_structure(_path("1DDB-39.pdb"))
    started = time.time()
    depth = np.asarray(bff.labelizer_residue_depth(structure, args.probe, args.points))
    elapsed = time.time() - started

    ours, refs = [], []
    for i, residue in enumerate(structure.residues):
        key = bff.labelizer_residue_key(residue.chain, residue.seq_id)
        want = centre.get(round(reference[key], 9))
        if want is None:
            continue
        ours.append(min(depth[i], 4.0))
        refs.append(want)
    ours, refs = np.asarray(ours), np.asarray(refs)
    residual = ours - refs
    width = float(np.diff(sorted(table.bins)).mean())

    print()
    print("residue depth against MSMS (%.2f s, bin width %.3f A)"
          % (elapsed, width))
    print("  bias           %+.3f A mean, %+.3f A median"
          % (residual.mean(), np.median(residual)))
    print("  scatter        %.3f A sd" % residual.std())
    print("  quantisation   %.3f A sd (bin width / sqrt 12) -- part of the above"
          % (width / np.sqrt(12.0)))
    print("  correlation    r = %.4f, slope %.4f"
          % (np.corrcoef(ours, refs)[0, 1], np.polyfit(refs, ours, 1)[0]))
    print("  within 1 bin   %.1f %%" % (100.0 * (np.abs(residual) < width).mean()))


def _pairs(args, got):
    """What the cheap dye model costs, and what it gets wrong."""
    combined = {k: v for k, v in got["combined"].items() if v is not None}

    options = bff.LabelizerFRETOptions()
    options.n_refine = 0
    started = time.time()
    screened = list(bff.labelizer_fret_pair_scores(_path("1DDB-39.pdb"), combined, options))
    screen_time = time.time() - started

    above = sum(1 for v in combined.values()
                if v >= options.label_score_threshold)
    print()
    print("pair screen: %d pairs from %d labelable sites in %.3f s (alpha cone)"
          % (len(screened), above, screen_time))
    if not args.refine:
        return

    options.n_refine = args.refine
    started = time.time()
    refined = list(bff.labelizer_fret_pair_scores(_path("1DDB-39.pdb"), combined, options))
    refine_time = time.time() - started
    print("  + rebuilding the top %d with real accessible volumes: %.2f s"
          % (args.refine, refine_time - screen_time))

    before = {(p.seq_id_1, p.seq_id_2): p for p in screened}
    changed = [p for p in refined
               if p.probe_model == bff.PROBE_MODEL_ACCESSIBLE_VOLUME]
    print()
    print("  %-14s %-22s %-22s %s" % ("pair", "cone", "AV", "dd"))
    for p in changed:
        was = before[(p.seq_id_1, p.seq_id_2)]
        print("  %-14s %8.4f @ %5.1f A %8.4f @ %5.1f A %+6.1f A"
              % ("%s%d-%s%d" % (p.asym_id_1, p.seq_id_1,
                                p.asym_id_2, p.seq_id_2),
                 was.value, was.distance, p.value, p.distance,
                 p.distance - was.distance))

    top_before = [(q.seq_id_1, q.seq_id_2) for q in screened[:args.refine]]
    top_after = [(q.seq_id_1, q.seq_id_2) for q in refined[:args.refine]]
    kept = len(set(top_before) & set(top_after))
    print()
    print("  of the top %d by cone score, %d are still in the top %d after "
          "refinement" % (args.refine, kept, args.refine))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--probe", type=float, default=1.4,
                    help="rolling-probe radius, Angstrom (default 1.4)")
    ap.add_argument("--points", type=int, default=590,
                    help="unit-sphere samples per atom (default 590)")
    ap.add_argument("--refine", type=int, default=5,
                    help="how many top pairs to rebuild with accessible "
                         "volumes; 0 skips the pair section (default 5)")
    args = ap.parse_args(argv)

    got = _agreement(args)
    _depth(args)
    _pairs(args, got)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
