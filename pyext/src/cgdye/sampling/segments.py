"""Segmentation and FP detection ported from fpsimp."""

from __future__ import annotations

import numpy as np
from typing import Dict, List, Tuple, Optional
from pathlib import Path
import logging

from .fp_lib import load_fp_library

logger = logging.getLogger(__name__)

def _smith_waterman(query: str, templ: str,
                    match: float = 2.0, mismatch: float = -1.0,
                    gap_open: float = -5.0, gap_extend: float = -1.0):
    """Local alignment with affine gaps, returning the aligned blocks.

    Replaces Biopython's PairwiseAligner in local mode with the same scoring, so
    that IMP.bff needs nothing beyond what IMP itself brings. Gotoh's three-matrix
    formulation: `main` for a residue pair, `gap_q`/`gap_t` for a gap continuing
    in one sequence, which is what makes an affine penalty exact rather than a
    per-position approximation.

    Returns
    -------
    tuple
        ``(query_blocks, templ_blocks)``, each a list of ``(start, end)``
        half-open index pairs, in the spelling Biopython's ``aligned`` uses.
    """
    n, m = len(query), len(templ)
    if n == 0 or m == 0:
        return [], []

    neg = float("-inf")
    main = np.zeros((n + 1, m + 1), dtype=np.float64)
    gap_q = np.full((n + 1, m + 1), neg, dtype=np.float64)   # gap in the template
    gap_t = np.full((n + 1, m + 1), neg, dtype=np.float64)   # gap in the query
    # 0 = stop, 1 = diagonal, 2 = up (gap in template), 3 = left (gap in query)
    trace = np.zeros((n + 1, m + 1), dtype=np.int8)

    best, best_ij = 0.0, (0, 0)
    for i in range(1, n + 1):
        qi = query[i - 1]
        for j in range(1, m + 1):
            gap_q[i, j] = max(main[i - 1, j] + gap_open, gap_q[i - 1, j] + gap_extend)
            gap_t[i, j] = max(main[i, j - 1] + gap_open, gap_t[i, j - 1] + gap_extend)
            diag = main[i - 1, j - 1] + (match if qi == templ[j - 1] else mismatch)
            cell = max(0.0, diag, gap_q[i, j], gap_t[i, j])
            main[i, j] = cell
            if cell == 0.0:
                trace[i, j] = 0
            elif cell == diag:
                trace[i, j] = 1
            elif cell == gap_q[i, j]:
                trace[i, j] = 2
            else:
                trace[i, j] = 3
            if cell > best:
                best, best_ij = cell, (i, j)

    if best <= 0.0:
        return [], []

    # Walk back to the first zero, collecting runs of diagonal steps as blocks.
    i, j = best_ij
    q_blocks, t_blocks = [], []
    run_q_end, run_t_end = i, j
    while i > 0 and j > 0 and trace[i, j] != 0:
        step = trace[i, j]
        if step == 1:
            i, j = i - 1, j - 1
        else:
            if run_q_end != i or run_t_end != j:
                q_blocks.append((i, run_q_end))
                t_blocks.append((j, run_t_end))
            if step == 2:
                i -= 1
            else:
                j -= 1
            run_q_end, run_t_end = i, j
    if run_q_end != i or run_t_end != j:
        q_blocks.append((i, run_q_end))
        t_blocks.append((j, run_t_end))

    return q_blocks[::-1], t_blocks[::-1]


def _align_best_window(query: str, templ: str) -> tuple[int, int, float]:
    """Local alignment; return (q_start_1based, q_end_1based, identity)."""
    q_blocks, t_blocks = _smith_waterman(query, templ)
    if not q_blocks:
        return (0, 0, 0.0)

    q_start = int(q_blocks[0][0])
    q_end = int(q_blocks[-1][1])

    matches = 0
    aligned_len = 0
    for (qs, qe), (ts, te) in zip(q_blocks, t_blocks):
        qseg = query[qs:qe]
        tseg = templ[ts:te]
        aligned_len += len(qseg)
        for qc, tc in zip(qseg, tseg):
            if qc == tc:
                matches += 1
    identity = (matches / aligned_len) if aligned_len else 0.0
    return (int(q_start + 1), int(q_end), float(identity))

def find_fp_domains(seq: str, min_identity: float = 0.35) -> List[Tuple[str, int, int, float]]:
    """Detect FP domains in a fusion protein sequence."""
    lib = load_fp_library()
    hits: List[Tuple[str, int, int, float]] = []

    for name, data in lib.items():
        templ = data["sequence"]
        motifs = data.get("motifs", [])
        
        best_seed = (-1, -1.0)
        for m in motifs:
            i = seq.find(m)
            if i >= 0 and len(m) > best_seed[1]:
                best_seed = (i, len(m))

        if best_seed[0] >= 0:
            i = best_seed[0]
            start = max(0, i - 50)
            end = min(len(seq), i + 300)
            qs, qe, idy = _align_best_window(seq[start:end], templ)
            if idy >= min_identity and qs > 0:
                hits.append((name, int(start + qs), int(start + qe), float(idy)))
            continue

        qs, qe, idy = _align_best_window(seq, templ)
        if idy >= min_identity and qs > 0:
            hits.append((name, int(qs), int(qe), float(idy)))

    hits.sort(key=lambda x: (-x[3], x[2] - x[1], x[1]))
    covered = set()
    non_overlap: List[Tuple[str, int, int, float]] = []
    for nm, s, e, idy in hits:
        if any(i in covered for i in range(s, e + 1)):
            continue
        for i in range(s, e + 1):
            covered.add(i)
        non_overlap.append((nm, s, e, idy))
    return non_overlap

def parse_plddt_from_pdb(pdb_path: str | Path, chain_id: str = "A") -> Dict[int, float]:
    """Parse pLDDT from AlphaFold PDB B-factor column."""
    plddt_sum: Dict[int, float] = {}
    plddt_cnt: Dict[int, int] = {}
    
    with open(pdb_path, "r") as f:
        for line in f:
            if not line.startswith("ATOM"): continue
            ch = line[21].strip()
            if chain_id and ch and ch != chain_id: continue
            try:
                resseq = int(line[22:26])
                b = float(line[60:66])
                plddt_sum[resseq] = plddt_sum.get(resseq, 0.0) + b
                plddt_cnt[resseq] = plddt_cnt.get(resseq, 0) + 1
            except ValueError: continue

    if not plddt_sum:
        raise ValueError(f"No pLDDT data found for chain '{chain_id}' in {pdb_path}")
    return {r: plddt_sum[r] / plddt_cnt[r] for r in sorted(plddt_sum.keys())}

def segments_from_plddt(
    seq_len: int,
    plddt: Dict[int, float],
    fp_domains: List[Tuple[str, int, int, float]],
    rigid_threshold: float = 70.0,
    min_rb_len: int = 12
) -> List[Dict]:
    """Create rigid/linker segments from pLDDT and FP detections."""
    labels = ['R' if plddt.get(i+1, 0.0) >= rigid_threshold else 'L' for i in range(seq_len)]
    
    # FP domains are always rigid
    for _, s, e, _ in fp_domains:
        for i in range(s-1, e):
            labels[i] = 'R'
            
    # Small rigid islands become linkers
    i = 0
    while i < seq_len:
        if labels[i] == 'R':
            j = i
            while j < seq_len and labels[j] == 'R': j += 1
            if (j - i) < min_rb_len:
                # Only if not part of an FP? Let's keep it simple for now
                for k in range(i, j): labels[k] = 'L'
            i = j
        else: i += 1
        
    segs = []
    i = 0
    while i < seq_len:
        curr = labels[i]
        j = i
        while j < seq_len and labels[j] == curr: j += 1
        
        kind = "core" if curr == 'R' else "linker"
        name = kind
        
        # Check if it's an FP
        for nm, s, e, _ in fp_domains:
            if max(i+1, s) <= min(j, e):
                if (min(j, e) - max(i+1, s) + 1) > 0.5 * (j - i):
                    kind = "fp"
                    name = nm
                    break
                    
        segs.append({"kind": kind, "name": name, "start": i+1, "end": j})
        i = j
    return segs
