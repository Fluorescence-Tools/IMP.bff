#!/usr/bin/env python
"""LJ pair scoring with CHARMM36 parameters and exclusion list generation.

Tasks 2.5/2.6: per-element LJ scoring, 1-2/1-3/1-4 exclusions.

Algorithmic tricks ported from IMP RotamerCalculator.cpp:
  - Trick 1: Axis-aligned bounding-box (AABB) clash pre-filter, analogous to
    ResidueRotamer::create_bounding_boxes / intersect().  Skips expensive pairwise
    LJ evaluation for rotamer frames whose atom cloud cannot possibly clash with
    a reference (protein) coordinate set.
  - Trick 3: Vectorized batch LJ evaluation using NumPy broadcasting — replaces
    the Python-level for-loop over frames with a single (n_frames, n_pairs)
    tensor operation.
"""


import numpy as np

from IMP.bff.cgdye.topology.dye import CHARMM36_LJ, lj_cross, lj_energy


# ---------------------------------------------------------------------------
# Scalar helpers (kept for single-frame Metropolis acceptance in LinkerSampler)
# ---------------------------------------------------------------------------

def lj_score(r, rmin, epsilon):
    """Repulsive-only Lennard-Jones energy for one pair (scalar).

    ``epsilon * ((rmin/r)**12 - 2*(rmin/r)**6)`` for ``r < rmin``, 0 beyond
    (no attractive tail in the coarse-grained linker sampler); ``r`` is clamped
    to 0.01 A. Thin scalar front to :func:`IMP.bff.cgdye.topology.dye.lj_energy`.
    """
    if rmin == 0 or epsilon == 0:
        return 0.0
    return lj_energy(r, rmin, epsilon, repulsive_only=True)


# ---------------------------------------------------------------------------
# Exclusion / topology helpers
# ---------------------------------------------------------------------------

def compute_exclusions(system):
    """Compute 1-2, 1-3, and 1-4 excluded pairs from system topology.

    Returns a set of frozenset({site_id_a, site_id_b}) pairs.
    """
    excluded = set()

    bond_pairs = set()
    for a, b, _length, _tid in system.get("bonds", []):
        bond_pairs.add(frozenset({a, b}))
    excluded.update(bond_pairs)

    angle_ends = set()
    for a, b, c, _theta, _tid in system.get("angles", []):
        pair = frozenset({a, c})
        if pair not in excluded:
            angle_ends.add(pair)
    excluded.update(angle_ends)

    dihedral_ends = set()
    for a, b, c, d, _tid in system.get("dihedrals", []):
        pair = frozenset({a, d})
        if pair not in excluded:
            dihedral_ends.add(pair)
    excluded.update(dihedral_ends)

    for a, b, c, d, _tid in system.get("impropers", []):
        for x, y in [(a, c), (a, d), (b, d)]:
            excluded.add(frozenset({x, y}))

    return excluded


def dye_internal_system(atoms_dict, bonds):
    """Minimal force-field ``system`` for a lone dye from ``parse_mol2`` output.

    Sites are the MOL2 atoms in serial order (ids ``dye:<atom_name>``, the
    order ``LinkerSampler.get_coords`` uses); bonds, angles and dihedrals are
    derived from the MOL2 connectivity so that :func:`compute_exclusions`
    removes the 1-2, 1-3 and 1-4 pairs from the LJ pair list. Without them the
    linker sampler scored bonded neighbours with the repulsive LJ, which is
    unphysical and biases torsion/angle sampling. Site ids are
    ``dye:<serial>:<atom_name>`` (unique; see below).
    """
    from IMP.bff.cgdye.topology.dye import build_angles, build_dihedrals, build_graph

    ordered = sorted(atoms_dict.values(), key=lambda x: x["serial"])
    # Site ids carry the serial: MOL2 atom names are not unique across the
    # dye + linker residues (alexa488_r48 has 83 atoms, 68 distinct names), and
    # name-only ids made same-named atoms one site, i.e. self-pairs at r = 0
    # with 1e31 energies in every frame.
    sid = {a["serial"]: f"dye:{a['serial']}:{a['atom_name']}" for a in ordered}
    sites = [{"id": sid[a["serial"]], "atom_name": a["atom_name"]} for a in ordered]
    graph = build_graph(bonds)
    return {
        "sites": sites,
        "bonds": [(sid[a], sid[b], 0.0, None) for a, b in sorted(bonds)],
        "angles": [(sid[a], sid[b], sid[c], 0.0, None) for a, b, c in build_angles(graph)],
        "dihedrals": [(sid[a], sid[b], sid[c], sid[d], None) for a, b, c, d in build_dihedrals(graph)],
    }


def build_lj_type_table(elements):
    """Build _ff_lj_type entries for a set of elements.

    Returns dict {type_id: {element, rmin_half, epsilon}}.
    """

    lj_types = {}
    for elem in sorted(elements):
        params = CHARMM36_LJ.get(elem, CHARMM36_LJ["C"])
        lj_types[f"LJ_{elem}"] = {
            "element": elem,
            "rmin_half": params["rmin_half"],
            "epsilon": params["epsilon"],
        }
    return lj_types


def site_element_map(system):
    """Extract {site_id: element} from system sites.

    Derives element from atom_name prefix (first letter capitalized).
    """
    import re

    result = {}
    for s in system.get("sites", []):
        sid = s["id"]
        aname = s.get("atom_name", "")
        m = re.match(r"([A-Za-z])", aname)
        elem = m.group(1).upper() if m else "C"
        result[sid] = elem
    return result


def lj_cross_params(elem_i, elem_j):
    """Lorentz-Berthelot combining rules for two elements.

    Returns (rmin_ij, epsilon_ij).
    """
    return lj_cross(elem_i, elem_j)


def compute_lj_pair_sites(system, excluded=None):
    """Compute all non-excluded site pairs with LJ cross parameters.

    Returns list of (site_id_a, site_id_b, rmin_ij, epsilon_ij).
    """
    if excluded is None:
        excluded = compute_exclusions(system)
    elem_map = site_element_map(system)
    sites = system.get("sites", [])
    pairs = []
    for i in range(len(sites)):
        for j in range(i + 1, len(sites)):
            sa = sites[i]["id"]
            sb = sites[j]["id"]
            pair = frozenset({sa, sb})
            if pair in excluded:
                continue
            ea = elem_map.get(sa, "C")
            eb = elem_map.get(sb, "C")
            rmin, eps = lj_cross_params(ea, eb)
            pairs.append((sa, sb, rmin, eps))
    return pairs


# ---------------------------------------------------------------------------
# Trick 1 — Axis-Aligned Bounding Box (AABB) pre-filter
# Mirrors ResidueRotamer::create_bounding_boxes / intersect() from IMP.
# Boxes are padded by half the maximum LJ rmin (≈ 3.5 Å default) so that
# any overlapping box pair *could* have a clashing atom pair.
# ---------------------------------------------------------------------------

class BoundingBoxFilter:
    """AABB clash pre-filter for rotamer coordinate batches.

    Algorithmic source: IMP RotamerCalculator.cpp,
    ResidueRotamer::create_bounding_boxes() and intersect().

    Usage::

        bbf = BoundingBoxFilter(pad=3.5)
        # boxes for all rotamer frames: shape (n_frames, 6)  [xmin,ymin,zmin,xmax,ymax,zmax]
        rot_boxes = bbf.build(rotamer_coords)       # (n_frames, n_atoms, 3)
        prot_box  = bbf.build_single(protein_coords)  # (n_prot_atoms, 3)
        mask = bbf.intersects_reference(rot_boxes, prot_box)  # bool (n_frames,)
        surviving = rotamer_coords[mask]
    """

    def __init__(self, pad: float = 3.5):
        """
        Args:
            pad: Padding in Angstroms added to each face of the AABB.
                 IMP uses 3.5 Å (≈ C–C LJ rmin/2).
        """
        self.pad = float(pad)

    def build(self, coords: np.ndarray) -> np.ndarray:
        """Build AABBs for a batch of coordinate frames.

        Args:
            coords: shape (n_frames, n_atoms, 3)

        Returns:
            boxes: shape (n_frames, 6) — columns are
                   [xmin, ymin, zmin, xmax, ymax, zmax] with padding applied.
        """
        mins = coords.min(axis=1) - self.pad   # (n_frames, 3)
        maxs = coords.max(axis=1) + self.pad   # (n_frames, 3)
        return np.concatenate([mins, maxs], axis=1)  # (n_frames, 6)

    def build_single(self, coords: np.ndarray) -> np.ndarray:
        """Build a single AABB for a static coordinate set.

        Args:
            coords: shape (n_atoms, 3)

        Returns:
            box: shape (6,) — [xmin, ymin, zmin, xmax, ymax, zmax]
        """
        return np.concatenate([
            coords.min(axis=0) - self.pad,
            coords.max(axis=0) + self.pad,
        ])

    @staticmethod
    def intersects(box_a: np.ndarray, box_b: np.ndarray) -> bool:
        """Test AABB intersection between two single boxes (shape (6,)).

        Direct port of ResidueRotamer::intersect() from IMP.
        """
        return (
            box_a[3] >= box_b[0] and box_a[0] <= box_b[3] and
            box_a[4] >= box_b[1] and box_a[1] <= box_b[4] and
            box_a[5] >= box_b[2] and box_a[2] <= box_b[5]
        )

    def intersects_reference(
        self,
        rot_boxes: np.ndarray,
        ref_box: np.ndarray,
    ) -> np.ndarray:
        """Vectorized AABB intersection test of many rotamer frames vs one reference.

        Args:
            rot_boxes: shape (n_frames, 6)
            ref_box:   shape (6,)

        Returns:
            mask: bool array (n_frames,) — True if frame *could* clash.
        """
        # rot_boxes columns: [xmin0, ymin1, zmin2, xmax3, ymax4, zmax5]
        # ref_box:           [xmin0, ymin1, zmin2, xmax3, ymax4, zmax5]
        overlap_x = (rot_boxes[:, 3] >= ref_box[0]) & (rot_boxes[:, 0] <= ref_box[3])
        overlap_y = (rot_boxes[:, 4] >= ref_box[1]) & (rot_boxes[:, 1] <= ref_box[4])
        overlap_z = (rot_boxes[:, 5] >= ref_box[2]) & (rot_boxes[:, 2] <= ref_box[5])
        return overlap_x & overlap_y & overlap_z

    def filter_frames(
        self,
        coords: np.ndarray,
        reference_coords: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Filter rotamer frames that cannot clash with a reference coordinate set.

        Args:
            coords: shape (n_frames, n_atoms, 3) — rotamer frames
            reference_coords: shape (n_ref_atoms, 3) — static reference (e.g. protein)

        Returns:
            surviving_coords: subset of coords where AABB overlaps reference
            mask:             bool (n_frames,) survivor mask
        """
        rot_boxes = self.build(coords)
        ref_box = self.build_single(reference_coords)
        mask = self.intersects_reference(rot_boxes, ref_box)
        return coords[mask], mask


# ---------------------------------------------------------------------------
# Trick 3 — Vectorized batch LJ scoring
# Replaces the Python for-loop in the old evaluate_batch with NumPy
# (n_frames, n_pairs) broadcasting — O(1) Python loops.
# ---------------------------------------------------------------------------

class InternalEnergyEvaluator:
    """Evaluates the internal LJ energy of dye conformations.

    evaluate()       — scalar, used inside Metropolis loop (unchanged API).
    evaluate_batch() — vectorized NumPy over all frames simultaneously
                       (Trick 3 from IMP RotamerCalculator inner-loop pattern).
    evaluate_batch_filtered() — AABB pre-filter (Trick 1) then vectorized LJ.
    """

    def __init__(self, system, excluded=None):
        self.pairs = compute_lj_pair_sites(system, excluded=excluded)
        # Map site_id to index in the coordinate array
        self.id_to_idx = {s["id"]: i for i, s in enumerate(system["sites"])}

        # Pre-extract index arrays and parameter arrays for vectorized eval
        if self.pairs:
            self._idx_a = np.array([self.id_to_idx[p[0]] for p in self.pairs], dtype=np.intp)
            self._idx_b = np.array([self.id_to_idx[p[1]] for p in self.pairs], dtype=np.intp)
            self._rmin  = np.array([p[2] for p in self.pairs], dtype=np.float64)
            self._eps   = np.array([p[3] for p in self.pairs], dtype=np.float64)
        else:
            self._idx_a = np.array([], dtype=np.intp)
            self._idx_b = np.array([], dtype=np.intp)
            self._rmin  = np.array([], dtype=np.float64)
            self._eps   = np.array([], dtype=np.float64)

    # ------------------------------------------------------------------
    # Scalar evaluate — kept for Metropolis single-step acceptance
    # ------------------------------------------------------------------
    def evaluate(self, coords):
        """Evaluate LJ energy for a single frame.

        Args:
            coords: Array of shape (n_atoms, 3)
        """
        # One frame through the vectorised kernel: the same arithmetic as
        # evaluate_batch (regression-tested equal to 1e-10), no Python pair loop.
        return float(self.evaluate_batch(np.asarray(coords, dtype=np.float64)[None, :, :])[0])

    # ------------------------------------------------------------------
    # Trick 3 — vectorized batch LJ evaluation
    # ------------------------------------------------------------------
    def evaluate_batch(self, coords_batch: np.ndarray) -> np.ndarray:
        """Vectorized LJ evaluation for all frames simultaneously.

        Trick 3: instead of a Python for-loop over frames, we broadcast
        the full (n_frames, n_pairs) distance tensor in one NumPy call.

        Args:
            coords_batch: shape (n_frames, n_atoms, 3)

        Returns:
            energies: shape (n_frames,)
        """
        if len(self.pairs) == 0:
            return np.zeros(coords_batch.shape[0])

        # Gather positions for all pairs — shape (n_frames, n_pairs, 3)
        pos_a = coords_batch[:, self._idx_a, :]   # (F, P, 3)
        pos_b = coords_batch[:, self._idx_b, :]   # (F, P, 3)

        # Pairwise distances — (F, P)
        delta = pos_a - pos_b
        r = np.linalg.norm(delta, axis=-1)        # (F, P)

        # Repulsive-only LJ through the shared kernel (r clamped away from zero)
        lj_val = lj_energy(r, self._rmin[None, :], self._eps[None, :], repulsive_only=True)
        return lj_val.sum(axis=-1)                # (F,)

    # ------------------------------------------------------------------
    # Combined Trick 1 + 3: filter then score
    # ------------------------------------------------------------------
    def evaluate_batch_filtered(
        self,
        coords_batch: np.ndarray,
        reference_coords: np.ndarray,
        pad: float = 3.5,
    ) -> tuple[np.ndarray, np.ndarray]:
        """AABB pre-filter (Trick 1) then vectorized LJ (Trick 3).

        Frames whose bounding box does not overlap the reference (protein)
        bounding box are assigned energy=0 without any pairwise computation.

        Args:
            coords_batch:     shape (n_frames, n_atoms, 3) — dye rotamers
            reference_coords: shape (n_ref_atoms, 3)       — protein atoms
            pad:              AABB padding in Å (default 3.5, same as IMP)

        Returns:
            energies: shape (n_frames,) — 0 for non-overlapping frames
            clash_mask: bool (n_frames,) — True where AABB overlapped
        """
        bbf = BoundingBoxFilter(pad=pad)
        surviving, clash_mask = bbf.filter_frames(coords_batch, reference_coords)

        energies = np.zeros(coords_batch.shape[0])
        if surviving.shape[0] > 0:
            energies[clash_mask] = self.evaluate_batch(surviving)

        return energies, clash_mask
