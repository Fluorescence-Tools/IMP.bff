"""Fluorescent Protein (FP) library accessor (``data/cgdye/fp_library.json``).

Copy of fpsim's FP registry (imp.bff has 4 FPs, fpsim has drifted); tracked in
chisurf PRD-96. Do not extend it here.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json
import logging
import math
import random

import numpy as np

import IMP.algebra

# --------------------------------------------------------------------------
# fp_lib
# --------------------------------------------------------------------------
"""Fluorescent Protein (FP) library accessor (``data/cgdye/fp_library.json``).

Copy of fpsim's FP registry (imp.bff has 4 FPs, fpsim has drifted); tracked in
chisurf PRD-96. Do not extend it here.
"""

def load_fp_library() -> Dict[str, Dict]:
    """Load the FP library from the bundled JSON file."""
    from IMP.bff.tools import _data_root
    path = _data_root() / "fp_library.json"
    with open(path, "r") as f:
        return json.load(f)

def get_fp_names() -> List[str]:
    """Return a list of all supported FP names."""
    return list(load_fp_library().keys())

def get_fp_data(name: str) -> Dict:
    """Return the data for a specific FP."""
    lib = load_fp_library()
    if name not in lib:
        raise ValueError(f"FP '{name}' not found in library. Available: {', '.join(lib.keys())}")
    return lib[name]


# --------------------------------------------------------------------------
# rrt
# --------------------------------------------------------------------------
"""RRT samplers: rigid transforms (Transformation3D) and periodic torsion vectors."""

@dataclass
class IMPNode:
    node_id: int
    transform: IMP.algebra.Transformation3D
    parent_id: int | None


class IMPRRTTree:
    def __init__(self, root_transform):
        self.nodes = {0: IMPNode(0, root_transform, None)}
        self.next_id = 1

    def add_node(self, transform, parent_id):
        nid = self.next_id
        self.nodes[nid] = IMPNode(nid, transform, parent_id)
        self.next_id += 1
        return nid

    def nearest(self, transform, rot_weight=0.25):
        best_id = None
        best_d = float("inf")
        for nid, n in self.nodes.items():
            d = transform_distance(n.transform, transform, rot_weight=rot_weight)
            if d < best_d:
                best_d = d
                best_id = nid
        return best_id

    def path_to_root(self, node_id):
        path = []
        cur = node_id
        while cur is not None:
            path.append(cur)
            cur = self.nodes[cur].parent_id
        return list(reversed(path))


def make_transform(tx, ty, tz, rx, ry, rz):
    rot = IMP.algebra.get_rotation_from_fixed_xyz(rx, ry, rz)
    tr = IMP.algebra.Vector3D(tx, ty, tz)
    return IMP.algebra.Transformation3D(rot, tr)


def transform_to_tuple(t):
    tr = t.get_translation()
    e = IMP.algebra.get_fixed_xyz_from_rotation(t.get_rotation())
    return (tr[0], tr[1], tr[2], e.get_x(), e.get_y(), e.get_z())


def transform_distance(a, b, rot_weight=0.25):
    at = a.get_translation()
    bt = b.get_translation()
    dt = math.sqrt((at[0] - bt[0]) ** 2 + (at[1] - bt[1]) ** 2 + (at[2] - bt[2]) ** 2)

    ae = IMP.algebra.get_fixed_xyz_from_rotation(a.get_rotation())
    be = IMP.algebra.get_fixed_xyz_from_rotation(b.get_rotation())
    dr = math.sqrt(
        (ae.get_x() - be.get_x()) ** 2
        + (ae.get_y() - be.get_y()) ** 2
        + (ae.get_z() - be.get_z()) ** 2
    )
    return dt + rot_weight * dr


def sample_transform(bounds, rng):
    vals = [rng.uniform(lo, hi) for lo, hi in bounds]
    return make_transform(*vals)


def steer_transform(from_t, to_t, step_size, rot_weight=0.25):
    d = transform_distance(from_t, to_t, rot_weight=rot_weight)
    if d == 0:
        return from_t
    if d <= step_size:
        return to_t
    alpha = step_size / d
    fa = transform_to_tuple(from_t)
    ta = transform_to_tuple(to_t)
    vals = [fa[i] + alpha * (ta[i] - fa[i]) for i in range(6)]
    return make_transform(*vals)


def is_collision_sphere(transform, obstacles, probe_radius):
    tr = transform.get_translation()
    x, y, z = tr[0], tr[1], tr[2]
    for ox, oy, oz, orad in obstacles:
        d2 = (x - ox) ** 2 + (y - oy) ** 2 + (z - oz) ** 2
        md = probe_radius + orad
        if d2 < md * md:
            return True
    return False


def rrt_grow_step(tree, target_t, step_size, collision_fn, rot_weight=0.25):
    near_id = tree.nearest(target_t, rot_weight=rot_weight)
    near_t = tree.nodes[near_id].transform
    new_t = steer_transform(near_t, target_t, step_size, rot_weight=rot_weight)
    if collision_fn(new_t):
        return None
    return tree.add_node(new_t, near_id)


def run_rigid_body_rrt(
    start_transform,
    bounds,
    n_iter,
    step_size,
    collision_fn,
    goal_transform=None,
    goal_bias=0.1,
    goal_tolerance=1.0,
    rot_weight=0.25,
    seed=0,
):
    """Grow an RRT over rigid-body transforms (``IMP.algebra.Transformation3D``).

    ``bounds`` are six (lo, hi) pairs for (tx, ty, tz, rx, ry, rz);
    ``collision_fn(transform) -> bool`` gates growth. Returns ``(tree, goal_node_id)``.
    """
    rng = random.Random(seed)
    tree = IMPRRTTree(start_transform)
    goal_node_id = None

    for _ in range(n_iter):
        if goal_transform is not None and rng.random() < goal_bias:
            target = goal_transform
        else:
            target = sample_transform(bounds, rng)

        new_id = rrt_grow_step(
            tree, target, step_size, collision_fn, rot_weight=rot_weight
        )
        if new_id is None:
            continue
        if goal_transform is not None:
            d = transform_distance(
                tree.nodes[new_id].transform, goal_transform, rot_weight=rot_weight
            )
            if d <= goal_tolerance:
                goal_node_id = new_id
                break

    return tree, goal_node_id


# ---------------------------------------------------------------------------
# RRT over a vector of periodic torsion angles (internal linker DOFs)
# ---------------------------------------------------------------------------

def _wrap_angle(a):
    """Wrap an angle into [-pi, pi)."""
    return (a + math.pi) % (2.0 * math.pi) - math.pi


def torsion_distance(a, b):
    """Euclidean distance between two torsion vectors on the torus."""
    return math.sqrt(sum(_wrap_angle(x - y) ** 2 for x, y in zip(a, b)))


def steer_torsions(from_cfg, to_cfg, step_size):
    """Move from ``from_cfg`` towards ``to_cfg`` by at most ``step_size`` (rad)."""
    d = torsion_distance(from_cfg, to_cfg)
    if d == 0.0:
        return list(from_cfg)
    alpha = min(1.0, step_size / d)
    return [_wrap_angle(f + alpha * _wrap_angle(t - f)) for f, t in zip(from_cfg, to_cfg)]


class TorsionRRTTree:
    """RRT tree whose nodes are torsion vectors (list of angles in rad)."""

    def __init__(self, root_cfg):
        self.configs = [list(root_cfg)]
        self.parents = [None]

    def add_node(self, cfg, parent_id):
        self.configs.append(list(cfg))
        self.parents.append(parent_id)
        return len(self.configs) - 1

    def nearest(self, cfg):
        best_id, best_d = 0, float("inf")
        for i, c in enumerate(self.configs):
            d = torsion_distance(c, cfg)
            if d < best_d:
                best_id, best_d = i, d
        return best_id

    def __len__(self):
        return len(self.configs)


def run_torsion_rrt(
    n_dofs,
    n_iter,
    step_size,
    collision_fn,
    *,
    start_cfg=None,
    goal_cfg=None,
    goal_bias=0.1,
    goal_tolerance=0.1,
    seed=0,
):
    """Grow an RRT over ``n_dofs`` periodic torsions.

    ``collision_fn(cfg) -> bool`` returns True when the configuration clashes
    (it is expected to apply ``cfg`` to the molecule as a side effect if the
    caller needs the coordinates). Returns ``(tree, goal_node_id)``.
    """
    rng = random.Random(seed)
    root = list(start_cfg) if start_cfg is not None else [0.0] * n_dofs
    tree = TorsionRRTTree(root)
    goal_node_id = None
    for _ in range(n_iter):
        if goal_cfg is not None and rng.random() < goal_bias:
            target = list(goal_cfg)
        else:
            target = [rng.uniform(-math.pi, math.pi) for _ in range(n_dofs)]
        near_id = tree.nearest(target)
        new_cfg = steer_torsions(tree.configs[near_id], target, step_size)
        if collision_fn(new_cfg):
            continue
        new_id = tree.add_node(new_cfg, near_id)
        if goal_cfg is not None and torsion_distance(new_cfg, goal_cfg) <= goal_tolerance:
            goal_node_id = new_id
            break
    return tree, goal_node_id


# --------------------------------------------------------------------------
# segments
# --------------------------------------------------------------------------
"""Segmentation and FP detection ported from fpsim (Biopython-free copy).

This is a copy of fpsim's segmenter kept for ``dye label-fusion``; the
duplication and its drift are tracked in chisurf PRD-96 (fpsim <-> imp.bff).
Do not extend it here -- fpsim is the home of the algorithm.
"""

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
