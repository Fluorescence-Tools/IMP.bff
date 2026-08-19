"""AV ↔ rotamer-ensemble comparison on the bundled hGBP1 and T4L systems (PRD-108).

Writes the authoritative table ``okf/validation/av_vs_rotamer.md`` and the
pin file ``test/references/cgdye_av_vs_rotamer_pins.json`` when run from the
repository (``--okf-dir``/``--pins`` override the paths); otherwise prints
the tables.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, Optional, Sequence, Tuple
import json

import numpy as np

from IMP.bff.representation.av import AccessibleVolume, compute_av
from IMP.bff.representation.distance import av_pair_statistics, chi2_score, fret_pair_efficiencies, fret_pair_geometry
from IMP.bff.representation.rotamer import RotamerEnsemble

__all__ = [
    'av_for_position',
    'compare_av_and_rotamer_pairs',
    'compare_av_and_rotamer_positions',
    'ensemble_for_position',
    'ensemble_pair_statistics',
    'markdown_table',
    'summary_numbers',
]

# --------------------------------------------------------------------------
# compare
# --------------------------------------------------------------------------
"""AV ↔ rotamer-ensemble cross-validation (PRD-108 stage 2).

For each labelling position: the accessible volume the fps.json position
describes (``fret.av.compute_av``, the AV1/AV3 model with its linker
parameters) next to the screened rotamer ensemble of a named library at the
same site (:class:`RotamerEnsemble`); per pair: ⟨R_DA⟩, ⟨R_DA⟩_E, σ_R and R_mp
from both, and -- when the fps.json carries experimental distances -- the
χ² each model gives. The numbers are *recorded* (the table in
``okf/validation/av_vs_rotamer.md`` is authoritative); a test pins them
against drift and asserts loose sanity bounds only: an AV with a 20 Å linker
and a screened rotamer library are two different physical models of the same
label, and their difference is data, not a failure.

Filed under ``fret`` rather than ``av`` or ``representation.rotamer``, which is
where it reads more naturally: a tool that compares two things has to sit above
both. It uses the AV builder in :mod:`IMP.bff.representation.av` and the rotamer ensemble
in :mod:`IMP.bff.representation.rotamer`, and ``fret`` is the lowest package
above both. Putting it in ``av`` closed the loop
``av -> fret -> restraints -> av``.
"""

def _av_position(pos: Dict[str, Any]) -> Dict[str, Any]:
    """The fps position with its authored ``strip_mask`` dropped.

    The shipped T4L/TG2 masks (``... or resname HOH SOL WAT ...``) are outside
    the fps dialect the strip engine reads today (PRD-106 open item); the AV
    is computed with the default strip (attachment residue side chain minus
    the attachment atom).
    """
    out = dict(pos)
    out.pop("strip_mask", None)
    return out


def av_for_position(pdb: str, pos: Dict[str, Any], disc_step: Optional[float] = None) -> AccessibleVolume:
    """The accessible volume of one fps.json position on ``pdb``."""
    p = _av_position(pos)
    ds = float(disc_step or p.get("simulation_grid_resolution", 1.5))
    return compute_av(
        np.zeros((0, 4)), np.zeros(3),
        float(p.get("linker_length", 20.0)), float(p.get("linker_width", 1.0)),
        (float(p.get("radius1", 3.5)), float(p.get("radius2", 0.0)), float(p.get("radius3", 0.0))),
        disc_step=ds, pdb_path=pdb, source_info=p)


def ensemble_for_position(pdb: str, pos: Dict[str, Any], library: str, **kwargs) -> RotamerEnsemble:
    """The screened rotamer ensemble of ``library`` at the position's residue."""
    return RotamerEnsemble.from_site(
        pdb, pos.get("chain_identifier") or None, int(pos["residue_seq_number"]), library, **kwargs)


def ensemble_pair_statistics(e1: AccessibleVolume, e2: AccessibleVolume, forster_radius: float) -> Tuple[float, float, float, float, float]:
    """(R_mp, ⟨R_DA⟩, ⟨R_DA⟩_E, σ_R, ⟨κ²⟩) from the full pair matrix (no sampling).

    ⟨R_DA⟩_E is taken with κ² = 2/3 for every pair -- the AV convention -- so
    it is comparable with ``av_pair_statistics``; the ensembles' own ⟨κ²⟩ is
    reported separately.
    """
    g_iso = fret_pair_geometry(e1.points[:, :3], e1.points[:, 3], e2.points[:, :3], e2.points[:, 3])
    w = g_iso["weight"]; r = g_iso["R"]
    rmp = float(np.linalg.norm(e1.mean_position - e2.mean_position))
    rda = float(np.sum(r * w))
    sigma = float(np.sqrt(max(np.sum(r * r * w) - rda * rda, 0.0)))
    e_static = fret_pair_efficiencies(g_iso, forster_radius)["static"]
    if e_static <= 0:
        rda_e = rda
    elif e_static >= 1:
        rda_e = 0.0
    else:
        rda_e = float(forster_radius * (1.0 / e_static - 1.0) ** (1.0 / 6.0))
    mu1, mu2 = getattr(e1, "mu", None), getattr(e2, "mu", None)
    kappa2_avg = 2.0 / 3.0
    if mu1 is not None and mu2 is not None and len(mu1) == e1.n_points and len(mu2) == e2.n_points:
        kappa2_avg = fret_pair_geometry(e1.points[:, :3], e1.points[:, 3], e2.points[:, :3], e2.points[:, 3], mu1, mu2)["kappa2_avg"]
    return rmp, rda, rda_e, sigma, kappa2_avg


def compare_av_and_rotamer_positions(
    pdb: str,
    positions: Dict[str, Dict[str, Any]],
    libraries: Dict[str, str],
    *,
    disc_step: Optional[float] = None,
    n_samples: int = 50000,
    **ensemble_kwargs,
) -> Dict[str, Dict[str, Any]]:
    """AV and rotamer ensemble per position with per-position numbers.

    ``libraries`` maps position name → rotamer library name. Returns
    ``{name: {"av", "ensemble", "d_mean_position", "av_n_points",
    "n_rotamers", "av_extent", "rot_extent"}}`` (extents = rms distance of
    the weighted cloud from its mean, Å).
    """
    from IMP.bff.representation.rotamer import load_protein_frames

    frame = load_protein_frames(pdb)[0]
    p_names = [str(n).upper() for n in frame["atom_names"]]
    p_res = np.asarray(frame["residue_indices"])
    p_chain = np.asarray([str(c).upper() for c in frame.get("chain_ids", [""] * len(p_names))])
    p_heavy = np.array([not n.startswith("H") for n in p_names])
    p_coords = np.asarray(frame["coords"], dtype=np.float64)

    def _interpenetration_weight(ens: RotamerEnsemble, cutoff: float = 2.5, min_weight: float = 1e-3) -> Tuple[float, float]:
        """(weight on rotamers whose heavy atoms come within ``cutoff`` Å of a
        protein heavy atom outside the labelled residue, min such distance over
        rotamers carrying more than ``min_weight``)."""
        same = (p_res == ens.residue) & ((p_chain == "") | (p_chain == (ens.chain or "").upper()) | (ens.chain == ""))
        prot = p_coords[p_heavy & ~same]
        keep = np.array([(not str(n).upper().startswith("H")) and str(n).upper() not in {"N", "CA", "C", "O", "OXT"}
                         for n in ens.atom_names])
        w_clash = 0.0
        d_min = float("inf")
        for i in range(ens.n_rotamers):
            d = np.linalg.norm(ens.atoms[i][keep][:, None, :] - prot[None, :, :], axis=2).min()
            if ens.weights[i] > min_weight:
                d_min = min(d_min, float(d))
            if d < cutoff:
                w_clash += float(ens.weights[i])
        return w_clash, d_min

    out: Dict[str, Dict[str, Any]] = {}
    for name, pos in positions.items():
        av = av_for_position(pdb, pos, disc_step)
        ens = ensemble_for_position(pdb, pos, libraries[name], position_name=name, **ensemble_kwargs)
        w_clash, d_min = _interpenetration_weight(ens)
        def _extent(cloud: AccessibleVolume) -> float:
            if cloud.n_points == 0:
                return float("nan")
            w = cloud.points[:, 3] / cloud.points[:, 3].sum()
            d2 = np.sum((cloud.points[:, :3] - cloud.mean_position) ** 2, axis=1)
            return float(np.sqrt(np.sum(w * d2)))
        out[name] = {
            "av": av, "ensemble": ens,
            "d_mean_position": float(np.linalg.norm(av.mean_position - ens.mean_position)),
            "av_n_points": int(av.n_points), "n_rotamers": int(ens.n_rotamers),
            "av_extent": _extent(av), "rot_extent": _extent(ens),
            "partition": float(ens.partition),
            "interpenetration_weight_2.5A": w_clash,
            "min_heavy_atom_distance": d_min,
        }
    return out


#: FRETpredict's partition-function cutoff: below it every rotamer clashes and
#: the ensemble is uniform-fallback, i.e. meaningless (the site is buried for
#: this library).
Z_CUTOFF = 0.05


def compare_av_and_rotamer_pairs(
    per_position: Dict[str, Dict[str, Any]],
    pairs: Iterable[Tuple[str, str]],
    forster_radius: float,
    *,
    experimental: Optional[Dict[Tuple[str, str], Dict[str, float]]] = None,
    n_samples: int = 50000,
) -> list[Dict[str, Any]]:
    """Per pair: AV and rotamer R_mp/⟨R_DA⟩/⟨R_DA⟩_E/σ_R (+ χ² against experiment).

    ``rotamer_valid`` is False when either ensemble has Z < ``Z_CUTOFF``
    (its numbers are the uniform fallback and are excluded from Σχ²).
    """
    rows = []
    for a, b in pairs:
        av1, av2 = per_position[a]["av"], per_position[b]["av"]
        e1, e2 = per_position[a]["ensemble"], per_position[b]["ensemble"]
        rmp_av, rda_av, rdae_av, sig_av = av_pair_statistics(av1, av2, forster_radius, n_samples=n_samples)
        rmp_r, rda_r, rdae_r, sig_r, k2_r = ensemble_pair_statistics(e1, e2, forster_radius)
        row: Dict[str, Any] = {
            "pair": f"{a}-{b}", "p1": a, "p2": b,
            "Rmp_av": rmp_av, "Rmp_rot": rmp_r,
            "RDAMean_av": rda_av, "RDAMean_rot": rda_r,
            "RDAMeanE_av": rdae_av, "RDAMeanE_rot": rdae_r,
            "sigma_av": sig_av, "sigma_rot": sig_r,
            "kappa2_rot": k2_r,
            "rotamer_valid": bool(per_position[a]["partition"] >= Z_CUTOFF and per_position[b]["partition"] >= Z_CUTOFF),
        }
        if experimental and (a, b) in experimental:
            exp = experimental[(a, b)]
            dtype = exp.get("distance_type", "RDAMean")
            key = {"RDAMean": "RDAMean", "RDAMeanE": "RDAMeanE", "Rmp": "Rmp"}[dtype]
            row["exp_distance"] = float(exp["distance"])
            row["exp_type"] = dtype
            row["chi2_av"] = chi2_score(row[f"{key}_av"], exp["distance"], exp["error_neg"], exp["error_pos"])
            row["chi2_rot"] = chi2_score(row[f"{key}_rot"], exp["distance"], exp["error_neg"], exp["error_pos"])
        rows.append(row)
    return rows


def markdown_table(per_position: Dict[str, Dict[str, Any]], rows: Sequence[Dict[str, Any]], title: str) -> str:
    """The comparison as two markdown tables."""
    lines = [f"### {title}", "", "| position | AV points | rotamers | Z | AV rms extent (Å) | rot rms extent (Å) | Δ mean position (Å) | weight within 2.5 Å of protein | min heavy-atom distance, w > 1e-3 (Å) |", "|---|---|---|---|---|---|---|---|---|"]
    for name, r in per_position.items():
        flag = " ‡" if r["partition"] < Z_CUTOFF else ""
        lines.append(f"| {name}{flag} | {r['av_n_points']} | {r['n_rotamers']} | {r['partition']:.3f} | {r['av_extent']:.1f} | {r['rot_extent']:.1f} | {r['d_mean_position']:.1f} | {r['interpenetration_weight_2.5A']:.2f} | {r['min_heavy_atom_distance']:.2f} |")
    if any(r["partition"] < Z_CUTOFF for r in per_position.values()):
        lines.append(f"\n‡ Z < {Z_CUTOFF}: every rotamer of the library clashes at this site (buried for the rotamer model); the ensemble is FRETpredict's uniform fallback and its pairs are excluded from Σχ².")
    has_exp = any("chi2_av" in r for r in rows)
    lines += ["", "| pair | R_mp AV / rot | ⟨R_DA⟩ AV / rot | ⟨R_DA⟩_E AV / rot (κ²=2/3) | σ_R AV / rot | ⟨κ²⟩ rot" + (" | exp (type) | χ² AV / rot |" if has_exp else " |"),
              "|---|---|---|---|---|---" + ("|---|---|" if has_exp else "|")]
    for r in rows:
        line = (f"| {r['pair']}{'' if r.get('rotamer_valid', True) else ' ‡'} | {r['Rmp_av']:.1f} / {r['Rmp_rot']:.1f} | {r['RDAMean_av']:.1f} / {r['RDAMean_rot']:.1f} "
                f"| {r['RDAMeanE_av']:.1f} / {r['RDAMeanE_rot']:.1f} | {r['sigma_av']:.1f} / {r['sigma_rot']:.1f} | {r['kappa2_rot']:.2f}")
        if has_exp:
            line += (f" | {r.get('exp_distance', float('nan')):.1f} ({r.get('exp_type', '')}) | "
                     f"{r.get('chi2_av', float('nan')):.2f} / {r.get('chi2_rot', float('nan')):.2f} |")
        else:
            line += " |"
        lines.append(line)
    if has_exp:
        valid = [r for r in rows if r.get("rotamer_valid", True)]
        chi_av = sum(r.get("chi2_av", 0.0) for r in valid)
        chi_rot = sum(r.get("chi2_rot", 0.0) for r in valid)
        chi_av_all = sum(r.get("chi2_av", 0.0) for r in rows)
        lines += ["", f"Σχ² over the {len(valid)} pairs with valid ensembles: AV {chi_av:.2f}, rotamer {chi_rot:.2f} "
                  f"(AV over all {len(rows)} pairs: {chi_av_all:.2f})."]
    return "\n".join(lines) + "\n"


def summary_numbers(per_position: Dict[str, Dict[str, Any]], rows: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    """The pin-able numbers (no objects)."""
    return {
        "positions": {n: {k: v for k, v in r.items() if k not in ("av", "ensemble")} for n, r in per_position.items()},
        "pairs": [dict(r) for r in rows],
    }
