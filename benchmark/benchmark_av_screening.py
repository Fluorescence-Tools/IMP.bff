"""Per-frame wall time of ``ProbeNetworkRestraint`` on the T4L docking trajectory.

This is the PRD-105 recording benchmark (acceptance criterion 8): it scores
the ``chi2_C1_33p`` distance set of ``examples/structure/T4L`` against the
first ``--frames`` frames of ``t4l_docking.rmf3`` and prints mean / median
wall time per frame together with the restraint's diagnostics counters,
once per requested mode.

    python benchmark_av_screening.py                # all modes, 100 frames
    python benchmark_av_screening.py --modes legacy default --frames 50

The baseline recorded on 2026-08-17 (arm64 build, before PRD-105) was
38.5 ms/frame mean for the legacy path.
"""
import argparse
import json
import statistics
import time

import RMF
import IMP
import IMP.rmf
import IMP.atom
import IMP.bff

MODES = {
    # name: (space_fixed, shared_map, distance)
    "default": (True, True, "quad"),
    "shared-mc": (True, True, "mc"),
    "lattice-private": (True, False, "quad"),
    "lattice-private-mc": (True, False, "mc"),
    "legacy-quad": (False, False, "quad"),
    "legacy": (False, False, "mc"),
}


def build(mode, hier, fps_json_path, score_set, n_samples, quad_k):
    kwargs = {}
    if mode is not None:
        space_fixed, shared_map, distance = MODES[mode]
        kwargs = dict(space_fixed=space_fixed, shared_map=shared_map,
                      distance=distance, quad_k=quad_k)
    return IMP.bff.ProbeNetworkRestraint(
        hier, fps_json_path, "ProbeNetworkRestraint",
        score_set, n_samples, **kwargs)


def run(mode, frames, n_samples, quad_k, repeats):
    m = IMP.Model()
    rmf_fn = IMP.bff.get_example_path("structure/T4L/t4l_docking.rmf3")
    f = RMF.open_rmf_file_read_only(rmf_fn)
    hier = IMP.rmf.create_hierarchies(f, m)[0]
    IMP.rmf.load_frame(f, RMF.FrameID(0))
    fps_json_path = IMP.bff.get_example_path("structure/T4L/fret.fps.json")
    r = build(mode, hier, fps_json_path, "chi2_C1_33p", n_samples, quad_k)
    r.unprotected_evaluate(None)  # cold

    times = []
    scores = []
    all_frames = list(f.get_root_frames())[:frames]
    for frame in all_frames:
        IMP.rmf.load_frame(f, frame)
        t0 = time.perf_counter()
        v = r.unprotected_evaluate(None)
        times.append((time.perf_counter() - t0) * 1e3)
        scores.append(v)

    # sampling-repeat proxy: the same frame evaluated `repeats` times
    rep_times = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        r.unprotected_evaluate(None)
        rep_times.append((time.perf_counter() - t0) * 1e3)

    out = {
        "mode": mode or "(constructor defaults)",
        "frames": len(times),
        "ms_per_frame_mean": statistics.mean(times),
        "ms_per_frame_median": statistics.median(times),
        "ms_per_repeat_mean": statistics.mean(rep_times) if rep_times else None,
        "score_first": scores[0],
        "score_last": scores[-1],
    }
    if hasattr(r, "get_diagnostics_json"):
        out["diagnostics"] = json.loads(r.get_diagnostics_json())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--modes", nargs="*", default=list(MODES))
    ap.add_argument("--frames", type=int, default=100)
    ap.add_argument("--n-samples", type=int, default=50000)
    ap.add_argument("--quad-k", type=int, default=100)
    ap.add_argument("--repeats", type=int, default=20)
    ap.add_argument("--baseline", action="store_true",
                    help="use the constructor defaults only (pre-PRD-105 builds)")
    args = ap.parse_args()
    modes = [None] if args.baseline else args.modes
    for mode in modes:
        res = run(mode, args.frames, args.n_samples, args.quad_k, args.repeats)
        print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
