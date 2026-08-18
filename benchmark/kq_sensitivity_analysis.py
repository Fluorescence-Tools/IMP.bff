"""kQ sensitivity analysis: what is the quenching parameter vector measurable from?

    <arm64>/bin/python benchmark/kq_sensitivity_analysis.py [--resolution 2.5]

PRD-110 stopped because theta is not identifiable from one decay: two of five
directions carry the information, condition number 1.4e10
(``okf/validation/quenching_identifiability.md``). The question that replaces it
is not how to fit faster but **what has to be measured**.

The cheapest candidate answer, and the one that needs no new physics and no new
instrument, is *more sites*. theta is global -- it is dye chemistry and dye
mobility, the same molecule at every position -- while the geometry is per site.
So N decays share one theta, and the joint Fisher information is the **sum** of
the per-site ones. If the blind directions rotate with the geometry, the sum is
better conditioned than any term in it and the answer is "label more sites". If
every site is blind in the *same* direction, no number of sites helps, the
degeneracy is structural rather than geometric, and the model has to be reduced
to what a decay can carry.

There is a mechanism to hope for, and it is worth naming before measuring so the
result can disappoint it: **diffusion only matters when the dye has to travel.**
At a site with a tryptophan in the dye's immediate contact shell the population
is quenched where it already sits, and the mobility never enters. At a site
whose only quencher is 10 A away, reaching it *is* the rate. Those two should
weight `free_diffusion` differently -- if they do, multi-site is the answer.

Sites are chosen on T4 lysozyme to span that axis, not for their pedigree.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

import IMP.bff

_spec = importlib.util.spec_from_file_location(
    "quenching_identifiability", Path(__file__).with_name("quenching_identifiability.py"))
_ident = importlib.util.module_from_spec(_spec)
sys.modules["quenching_identifiability"] = _ident
_spec.loader.exec_module(_ident)

PARAMETERS = _ident.PARAMETERS
Site = _ident.Site
STRUCTURE = _ident.STRUCTURE
N_PEAK_COUNTS = _ident.N_PEAK_COUNTS

#: Six chain-A sites spanning "quencher in contact" to "no quencher at all".
#: The last is a control: a site the decay cannot speak about is worth including
#: precisely because a joint fit is only as good as what each term contributes.
SITES = (
    dict(chain="A", residue=132, atom="CB", note="PRD-110 reference"),
    dict(chain="A", residue=124, atom="CB", note="TRP126 at 3.5 A -- contact"),
    dict(chain="A", residue=93, atom="CB", note="CYS97 6.1, TRP158 7.1 -- crowded"),
    dict(chain="A", residue=19, atom="CB", note="TYR18/TYR25 at ~6 A -- many, medium"),
    dict(chain="A", residue=65, atom="CB", note="HIS31 at 9.7 A -- sparse, distant"),
    dict(chain="A", residue=53, atom="CB", note="CYS54 at 6.4 A, otherwise bare"),
)


def geometry(site_obj, quencher_xyz) -> dict:
    """The covariate the hypothesis is about: how far must the dye travel?"""
    dg, x0 = site_obj.dg, site_obj.x0
    idx = np.argwhere(site_obj.bounds > 0)
    centre = (np.array(site_obj.bounds.shape) - 1) // 2
    points = x0[None, :] + (idx - centre[None, :]) * dg
    d = np.linalg.norm(points[:, None, :] - quencher_xyz[None, :, :], axis=2)
    nearest = d.min(axis=1)
    return {
        "accessible_voxels": int(idx.shape[0]),
        "nearest_quencher_from_attachment": float(
            np.linalg.norm(quencher_xyz - x0, axis=1).min()),
        "min_dye_quencher_distance": float(nearest.min()),
        "median_dye_quencher_distance": float(np.median(nearest)),
        "fraction_within_6A": float((nearest < 6.0).mean()),
    }


def quencher_positions(atoms) -> np.ndarray:
    """Coordinates of the redox moieties -- what the PET rate is measured from."""
    keep = np.zeros(len(atoms), dtype=bool)
    for i, row in enumerate(atoms):
        names = IMP.bff.QUENCHER_ATOMS.get(str(row["res_name"]), ())
        keep[i] = str(row["atom_name"]) in names
    return np.ascontiguousarray(atoms["coord"][keep], dtype=np.float64)


def site_fisher(site_obj, theta):
    """Noise-normalised, relative-scaled Jacobian and Fisher matrix for one site."""
    reference = site_obj.decay(theta)
    counts = reference * N_PEAK_COUNTS
    sigma = np.sqrt(np.maximum(counts, 1.0)) / N_PEAK_COUNTS
    J_rel, eigenvalues, eigenvectors = _ident.analyse(site_obj, theta, sigma)
    return reference, J_rel, eigenvalues, eigenvectors


def describe(vector, names, threshold=0.15) -> str:
    terms = " ".join(
        f"{w:+.2f}*{n}" for n, w in
        sorted(zip(names, vector), key=lambda kv: -abs(kv[1]))
        if abs(w) > threshold)
    return terms or "(no dominant term)"


def spectrum(fisher, names, label):
    values, vectors = np.linalg.eigh(fisher)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    condition = values[0] / max(values[-1], 1e-300)
    print(f"\n  {label}")
    print(f"    eigenvalues: {np.array2string(values, precision=3)}")
    print(f"    condition number: {condition:.4g}")
    for rank in range(len(names)):
        print(f"      lambda={values[rank]:11.4g}  {describe(vectors[:, rank], names)}")
    try:
        errors = np.sqrt(np.diag(np.linalg.inv(fisher)))
    except np.linalg.LinAlgError:
        errors = None
        print("    Fisher matrix is singular.")
    if errors is not None:
        print("    marginal uncertainty of each parameter (relative):")
        for name, error in zip(names, errors):
            verdict = "identifiable" if error < 0.5 else "POORLY DETERMINED"
            print(f"      {name:<20} +/- {error * 100:9.1f} %   {verdict}")
    return values, vectors, condition, errors


def principal_angles(a, b) -> np.ndarray:
    """Angles between two subspaces, in degrees. 0 means they share a direction.

    **The rank matters.** Two subspaces of dimension `r` in `R^n` are forced to
    intersect once `2r > n`, so with five parameters a rank-3 comparison reports
    a 0-degree first angle for every pair regardless of the geometry -- an
    arithmetic fact wearing the costume of a result. Keep `2r <= n`.
    """
    qa, _ = np.linalg.qr(a)
    qb, _ = np.linalg.qr(b)
    singular = np.linalg.svd(qa.T @ qb, compute_uv=False)
    return np.degrees(np.arccos(np.clip(singular, -1.0, 1.0)))


def joint_fit(site_objects, theta_true, observed, sigmas):
    """The confirmation the Fisher analysis is worth: one theta, six decays.

    PRD-110's gate was decided by a fit as well as by eigenvalues, and for the
    same reason -- a condition number is a linearisation about one point, and an
    optimiser walking a real surface is the thing that has to work.
    """
    from scipy.optimize import least_squares

    bounds = (np.array([p[2] for p in PARAMETERS]),
              np.array([p[3] for p in PARAMETERS]))
    start = np.clip(
        np.array([theta_true[i] * (1.6 if i % 2 == 0 else 0.6)
                  for i in range(len(theta_true))]),
        bounds[0], bounds[1])

    def residuals(theta):
        return np.concatenate([
            (obj.decay(theta) - data) / sig
            for obj, data, sig in zip(site_objects, observed, sigmas)])

    print("\n  joint least-squares through the solver "
          f"({len(site_objects)} sites, one theta):", flush=True)
    for obj in site_objects:
        obj.n_evaluations = 0
    t0 = time.perf_counter()
    result = least_squares(residuals, start, bounds=bounds,
                           diff_step=0.05, xtol=1e-6)
    elapsed = time.perf_counter() - t0
    solves = sum(obj.n_evaluations for obj in site_objects)

    print(f"    parameter evaluations: {solves // len(site_objects)}, "
          f"forward solves: {solves}, wall clock: {elapsed:.1f} s")
    recovery = {}
    worst = 0.0
    for i, (name, *_rest) in enumerate(PARAMETERS):
        true, found = theta_true[i], result.x[i]
        error = (found - true) / true * 100.0
        worst = max(worst, abs(error))
        recovery[name] = {"true": float(true), "found": float(found),
                          "percent": float(error), "start": float(start[i])}
        print(f"    {name:<20} true {true:8.3f}  start {start[i]:8.3f}  "
              f"found {found:8.3f}  ({error:+7.1f} %)")
    print(f"    worst recovery error: {worst:.1f} %")
    return {"forward_solves": int(solves), "seconds": elapsed,
            "cost": float(result.cost), "worst_percent": float(worst),
            "recovery": recovery}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolution", type=float, default=2.5,
                        help="AV grid resolution in Angstrom (default 2.5).")
    parser.add_argument("--blind-rank", type=int, default=2,
                        help="how many weakest directions form the blind subspace. "
                             "Must satisfy 2*rank <= n_parameters: two subspaces "
                             "of dimension r in R^n intersect whenever 2r > n, so a "
                             "larger rank reports a 0-degree angle by dimension "
                             "counting alone and measures nothing.")
    parser.add_argument("--fit", action="store_true",
                        help="also run the joint least-squares fit (slow).")
    parser.add_argument("--out", type=Path, default=None,
                        help="write the numbers as JSON.")
    args = parser.parse_args()

    pdb_path = IMP.bff.get_example_path(STRUCTURE)
    names = [p[0] for p in PARAMETERS]
    theta0 = np.array([p[1] for p in PARAMETERS], dtype=np.float64)

    print(f"structure {STRUCTURE}, resolution {args.resolution} A, "
          f"{len(SITES)} sites, theta = "
          + ", ".join(f"{n}={v:g}" for n, v in zip(names, theta0)), flush=True)

    records, fishers, blind = [], [], []
    site_objects, observed, sigmas = [], [], []
    rng = np.random.default_rng(0)
    for site in SITES:
        tag = f"{site['chain']}{site['residue']}.{site['atom']}"
        print(f"\n== {tag}  ({site['note']})", flush=True)
        t0 = time.perf_counter()
        site_obj = Site(pdb_path, args.resolution, site=site)
        geom = geometry(site_obj, quencher_positions(site_obj.atoms))
        print(f"  grid {site_obj.bounds.shape}, {geom['accessible_voxels']} voxels; "
              f"nearest quencher {geom['min_dye_quencher_distance']:.1f} A from the "
              f"volume, {geom['fraction_within_6A'] * 100:.1f} % of it within 6 A",
              flush=True)

        reference, J_rel, values, vectors = site_fisher(site_obj, theta0)
        sigma = np.sqrt(np.maximum(reference * N_PEAK_COUNTS, 1.0)) / N_PEAK_COUNTS
        site_objects.append(site_obj)
        sigmas.append(sigma)
        observed.append(reference + rng.normal(0.0, sigma))
        fisher = J_rel.T @ J_rel
        fishers.append(fisher)
        blind.append(vectors[:, -args.blind_rank:])

        condition = values[0] / max(values[-1], 1e-300)
        print(f"    lambda_max {values[0]:.4g}, condition {condition:.4g}; "
              f"best direction {describe(vectors[:, 0], names)}")
        # How much of this site's information is about the mobility at all?
        column = dict(zip(names, np.linalg.norm(J_rel, axis=0)))
        print("    per-parameter |dF/dlog p| (relative units): "
              + "  ".join(f"{n}={column[n]:.3g}" for n in names))
        records.append({
            "site": tag, "note": site["note"], "geometry": geom,
            "eigenvalues": values.tolist(), "condition_number": float(condition),
            "column_norms": {k: float(v) for k, v in column.items()},
            "best_direction": dict(zip(names, vectors[:, 0].tolist())),
            "final_fluorescence": float(reference[-1]),
            "seconds": time.perf_counter() - t0,
        })

    print("\n" + "=" * 72)
    joint = sum(fishers)
    _v, _e, joint_condition, joint_errors = spectrum(
        joint, names, f"joint Fisher information, {len(SITES)} sites")

    best_single = min(r["condition_number"] for r in records)
    print(f"\n  best single site: condition {best_single:.4g}")
    print(f"  all sites jointly: condition {joint_condition:.4g}"
          f"   ({best_single / joint_condition:.3g}x better)")

    print("\n  do the blind directions rotate between sites?"
          f"  (principal angles between the {args.blind_rank} weakest directions;"
          " 0 deg = the same blindness, 90 deg = complementary)")
    pairs = []
    for i in range(len(SITES)):
        for j in range(i + 1, len(SITES)):
            angles = principal_angles(blind[i], blind[j])
            pairs.append({"a": records[i]["site"], "b": records[j]["site"],
                          "angles_deg": angles.tolist()})
            print(f"    {records[i]['site']:<9} vs {records[j]['site']:<9} "
                  + "  ".join(f"{a:6.2f}" for a in angles))
    largest = max(max(p["angles_deg"]) for p in pairs)
    print(f"    largest principal angle over all pairs: {largest:.2f} deg")

    # The decisive scalar, and the one that cannot be flattered by a choice of
    # rank: if every site were blind in the *same* direction, the joint matrix
    # would be about as blind as the sum of the per-site blindnesses. If the
    # directions rotate, the sum is far better than its terms.
    per_site_minima = [float(np.linalg.eigvalsh(f)[0]) for f in fishers]
    shared_blindness = float(sum(per_site_minima))
    joint_minimum = float(np.linalg.eigvalsh(joint)[0])
    print("\n  is the blindness shared or rotated?")
    print(f"    sum of per-site smallest eigenvalues:  {shared_blindness:.4g}"
          "   (what N sites give if they are all blind the same way)")
    print(f"    smallest joint eigenvalue:             {joint_minimum:.4g}")
    print(f"    gain from rotation:                    "
          f"{joint_minimum / max(shared_blindness, 1e-300):.4g}x")

    identifiable = (joint_errors is not None
                    and bool(np.all(joint_errors < 0.5)))
    print("\n  VERDICT: theta is "
          + ("identifiable" if identifiable else "NOT identifiable")
          + f" from {len(SITES)} decays jointly.")

    payload = {
        "structure": STRUCTURE, "resolution": args.resolution,
        "parameters": names, "theta0": theta0.tolist(),
        "sites": records,
        "joint_eigenvalues": _v.tolist(),
        "joint_condition_number": float(joint_condition),
        "joint_relative_errors": None if joint_errors is None else joint_errors.tolist(),
        "best_single_condition_number": float(best_single),
        "blind_subspace_angles": pairs,
        "largest_principal_angle_deg": float(largest),
        "blind_rank": args.blind_rank,
        "per_site_smallest_eigenvalue": per_site_minima,
        "sum_of_per_site_smallest": shared_blindness,
        "joint_smallest_eigenvalue": joint_minimum,
        "rotation_gain": float(joint_minimum / max(shared_blindness, 1e-300)),
        "identifiable": identifiable,
    }
    if args.fit:
        payload["joint_fit"] = joint_fit(site_objects, theta0, observed, sigmas)

    if args.out:
        args.out.write_text(json.dumps(payload, indent=2) + "\n")
        print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
